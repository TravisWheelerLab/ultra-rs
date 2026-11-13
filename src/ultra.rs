//! The ULTRA driver (`ultra.hpp` / `ultra.cpp`): reads sequence windows,
//! runs the model on each, joins repeats found in overlapping windows, and
//! writes them out.
//!
//! # Threads
//!
//! The C++ runs `-t N` pthreads that each pull windows, analyse them and
//! hand their repeats over through unsynchronised shared state, so its output
//! with more than one thread varies from run to run. This port keeps the
//! single-threaded schedule for every thread count: windows are read in the
//! same batches (the whole window pool, refilled when the ready heap runs
//! dry), analysed by up to `N` worker threads, and the results are consumed
//! strictly in read-ID order by the calling thread, which does exactly what
//! the C++ does after analysing a window with one thread. Every quantity in
//! that bookkeeping that depends on the thread count is evaluated as for one
//! thread ([`LOGICAL_THREADS`]), so `-t N` writes the same output as `-t 1`.
//! With `-t 1` no thread is spawned.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::mpsc;
use std::sync::Mutex;

use crate::cli::{Settings, ULTRA_VERSION_STRING};
use crate::file_reader::FileReader;
use crate::mask::MRegion;
use crate::repeat::{get_next_repeat, joint_repeat_region, repeats_overlap, RepeatRegion};
use crate::repeat_splitter::{validate_splits, SplitWindow};
use crate::sequence_window::SequenceWindow;
use crate::std_sort::std_sort;
use crate::symbol::symbol_for_char;
use crate::umodel::UModel;
use crate::writers::{
    log2_pval_for_score, pval_for_score, BEDFileWriter, JSONFileWriter, OutSink, RepeatFileWriter,
    TabFileWriter,
};

/// The thread count used wherever the C++ bookkeeping depends on
/// `numberOfThreads` (`3 * numberOfThreads` in `OutputRepeats`, the
/// round-robin `primaryThread`). Results are always consumed as by one
/// thread; see the module documentation.
pub const LOGICAL_THREADS: i32 = 1;

/// Per-thread bookkeeping of the C++ `uthread`. There is one, for the single
/// logical thread that consumes window results.
#[derive(Default)]
pub struct UThread {
    pub id: i32,
    /// Repeats found in the current window and not yet moved to
    /// `out_repeats` (they stay here for `find_repeats_in_string`).
    pub repeats: Vec<Box<RepeatRegion>>,
    pub active_read_id: i32,
    pub smallest_read_id: i32,
}

/// The per-thread model and splitter of the C++ `uthread`: everything needed
/// to analyse one window. Analysis does not depend on which analyzer ran the
/// previous window, so any analyzer may take any window.
pub struct WindowAnalyzer {
    pub model: UModel,
    pub splitter: SplitWindow,
}

impl WindowAnalyzer {
    /// The model-and-split part of C++ `AnalyzeSequenceWindow`: runs the model
    /// over the window and returns its repeats, with sequence, traceback,
    /// scores and logo numbers stored when `store_trace_and_sequence` is set,
    /// and splits/consensi filled in for periods up to `max_split`.
    pub fn analyze_window(
        &mut self,
        sequence: &SequenceWindow,
        settings: &Settings,
        store_trace_and_sequence: bool,
    ) -> Vec<Box<RepeatRegion>> {
        let mut repeats = Vec::new();
        let sleng: i32 = (sequence.length as i32).wrapping_add(sequence.overlap as i32);

        let model = &mut self.model;
        model.matrix.restart_matrix();

        for i in 0..sleng {
            model.calculate_current_column(sequence, i, None);
        }
        model.calculate_current_column_without_emission(None);

        let prev = model.matrix.previous_column_index;
        model.matrix.calculate_traceback(prev);
        let mut i: i32 = 0;
        let mut r = get_next_repeat(sequence, &model.matrix, &mut i, settings.upstream_consensus);

        while let Some(mut rr) = r {
            if store_trace_and_sequence {
                rr.store_sequence(sequence);
                rr.store_traceback(&model.matrix);
                rr.store_scores(&model.matrix);

                rr.get_logo_numbers();
            }
            // std::max(int * unsigned long long, unsigned long long), narrowed to int.
            let split_window: i32 = std::cmp::max(
                (rr.repeat_period as i64 as u64).wrapping_mul(settings.split_depth),
                settings.min_split_window,
            ) as i32;
            if rr.repeat_period <= settings.max_split {
                // unsigned long > int: the int converts to unsigned long.
                if rr.repeat_length > (2i32.wrapping_mul(split_window)) as i64 as u64 {
                    let mut splits = self.splitter.splits_for_region(
                        &rr,
                        split_window,
                        settings.split_threshold,
                    );
                    let mut consensi = self.splitter.consensi_for_split(&rr, &splits, 0.65);
                    validate_splits(&mut consensi, &mut splits, 0.85);
                    rr.splits = Some(splits);
                    rr.consensi = Some(consensi);
                } else {
                    rr.splits = Some(Vec::new());
                    rr.consensi = Some(vec![rr.string_consensus.clone()]);
                }
            }

            repeats.push(rr);
            r = get_next_repeat(sequence, &model.matrix, &mut i, settings.upstream_consensus);
        }

        repeats
    }
}

/// C++ `CompareRepeatOrder`: descending by sequence ID, then by start, so
/// the earliest repeat ends up at the back of the vector.
pub fn compare_repeat_order(lhs: &RepeatRegion, rhs: &RepeatRegion) -> bool {
    if lhs.sequence_id != rhs.sequence_id {
        return lhs.sequence_id > rhs.sequence_id;
    }
    lhs.sequence_start > rhs.sequence_start
}

/// C++ `Ultra`.
pub struct Ultra {
    pub settings: Settings,

    /// `None` when `settings.run_without_reader` is set.
    pub reader: Option<FileReader>,
    pub writers: Vec<Box<dyn RepeatFileWriter>>,
    /// `writers[i]` writes to `outs[i]`.
    pub outs: Vec<OutSink>,

    /// The requested thread count (`settings.threads`): the number of window
    /// analyzers. Output bookkeeping uses [`LOGICAL_THREADS`].
    pub number_of_threads: i32,
    pub min_reader_size: i32,
    pub primary_thread: i32,

    /// The `.settings` file; `None` means stdout (or already closed).
    pub settings_out: Option<BufWriter<File>>,

    pub score_threshold: f32,

    pub output_repeat_sequence: bool,
    pub output_read_id: bool,
    pub multithreading: bool,
    pub shuffle_sequence: bool,

    pub store_sequence: bool,
    pub store_traceback: bool,
    pub store_scores: bool,
    pub store_profile_numbers: bool,

    pub can_output: bool,
    pub first_repeat: bool,
    pub repeat_buffer: u64,

    pub total_coverage: u64,
    /// Upstream quirk: never updated, so `last_rep_end` resets for every
    /// repeat with a read ID above 0.
    pub last_seq_id: u64,
    pub last_rep_end: u64,

    pub store_trace_and_sequence: bool,

    pub masks_for_seq: HashMap<u64, Vec<MRegion>>,
    pub out_repeats: Vec<Box<RepeatRegion>>,
    /// One per requested thread (C++ `models` plus each thread's splitter).
    pub analyzers: Vec<WindowAnalyzer>,

    pub count: i32,
    pub pass_id: i32,
    pub threads: Vec<UThread>,
    pub corrected_read_ids: i32,
}

/// `fopen(path, "w")`, or the C++ error message and `exit(-1)`.
fn open_output(path: &str, what: &str) -> File {
    match File::create(path) {
        Ok(f) => f,
        Err(_) => {
            let _ = io::stdout().flush();
            eprintln!("Unable to open {} {}", what, path);
            std::process::exit(255);
        }
    }
}

impl Ultra {
    /// C++ `Ultra(Settings *)`. Opens the output files (and the `.settings`
    /// file), creates the writers, the input reader (unless
    /// `run_without_reader`) and one model per thread. Exits with status 255
    /// if a file cannot be opened.
    pub fn new(s: &Settings) -> Ultra {
        let settings = s.clone();
        let mut writers: Vec<Box<dyn RepeatFileWriter>> = Vec::new();
        let mut outs: Vec<OutSink> = Vec::new();
        let mut settings_out: Option<BufWriter<File>> = None;

        if !settings.out_file.is_empty() {
            let mut c = 0;
            if settings.ultra_out {
                c += 1;
            }
            if settings.json_out {
                c += 1;
            }
            if settings.bed_out {
                c += 1;
            }
            if c > 1 {
                if settings.ultra_out {
                    let ultra_path = format!("{}.tsv", settings.out_file);
                    outs.push(OutSink::file(open_output(&ultra_path, "output file")));
                    writers.push(Box::new(TabFileWriter::new()));
                }

                if settings.json_out {
                    let json_path = format!("{}.json", settings.out_file);
                    outs.push(OutSink::file(open_output(&json_path, "output file")));
                    writers.push(Box::new(JSONFileWriter::new()));
                }

                if settings.bed_out {
                    let bed_path = format!("{}.bed", settings.out_file);
                    outs.push(OutSink::file(open_output(&bed_path, "output file")));
                    writers.push(Box::new(BEDFileWriter::new()));
                }
            } else {
                outs.push(OutSink::file(open_output(
                    &settings.out_file,
                    "output file",
                )));
                if settings.ultra_out {
                    writers.push(Box::new(TabFileWriter::new()));
                } else if settings.json_out {
                    writers.push(Box::new(JSONFileWriter::new()));
                } else if settings.bed_out {
                    writers.push(Box::new(BEDFileWriter::new()));
                }
            }

            let settings_file_path = format!("{}.settings", settings.out_file);
            if !settings.hide_settings {
                settings_out = Some(BufWriter::new(open_output(
                    &settings_file_path,
                    "settings output file",
                )));
            }
        } else {
            if settings.ultra_out {
                writers.push(Box::new(TabFileWriter::new()));
            } else if settings.json_out {
                writers.push(Box::new(JSONFileWriter::new()));
            } else if settings.bed_out {
                writers.push(Box::new(BEDFileWriter::new()));
            }
            outs.push(OutSink::stdout());
        }

        let number_of_threads = settings.threads;
        let score_threshold = settings.min_score;
        let output_read_id = settings.show_wid;
        let output_repeat_sequence = settings.show_seq;
        let shuffle_sequence = false;

        let mut reader = None;
        if !settings.run_without_reader {
            // int maxWindows, unsigned long long lengths.
            let mut r = FileReader::new(
                &settings.in_file,
                settings.windows as i32,
                settings.window_size as u64,
                settings.overlap as u64,
                settings.threads > 1,
            );
            r.fasta_reader.shuffle = shuffle_sequence;
            reader = Some(r);
        }

        // long long narrowed to int.
        let leng: i32 = settings
            .window_size
            .wrapping_add(settings.overlap.wrapping_add(2)) as i32;
        // Upstream sets this unconditionally, then again under conditions.
        let store_trace_and_sequence = true;

        let mut analyzers = Vec::new();
        for _ in 0..number_of_threads.max(0) {
            // max_period + 1 so that a max period of 10 finds period-10 repeats.
            let mut model = UModel::new(
                settings.max_period.wrapping_add(1) as i32,
                settings.max_insert as i32,
                settings.max_delete as i32,
                leng.wrapping_add(1),
            );

            // Assigned before CalculateScores, which turns it into its log2.
            model.period_decay = settings.period_decay;

            model.set_match_probabilities(settings.match_probability);
            model.set_atcg_probabilities(
                settings.a_freq,
                settings.t_freq,
                settings.c_freq,
                settings.g_freq,
            );

            model.tp_zero_to_match = settings.transition_nr;
            model.tp_match_to_zero = settings.transition_rn;

            model.tp_match_to_insertion = settings.transition_ri;
            model.tp_match_to_deletion = settings.transition_rd;

            model.tp_consecutive_insertion = settings.transition_ii;
            model.tp_consecutive_deletion = settings.transition_dd;

            model.calculate_scores();

            let mut splitter = SplitWindow::new();
            if !settings.no_split {
                splitter.allocate_split_window(4, settings.max_period.wrapping_add(1) as i32);
            } else {
                splitter.max_period = 0;
            }
            analyzers.push(WindowAnalyzer { model, splitter });
        }

        let mut threads = Vec::new();
        if number_of_threads > 0 {
            let mut t = UThread::default();
            t.repeats.reserve(64);
            threads.push(t);
        }

        let multithreading = number_of_threads != 1;
        if let Some(r) = reader.as_mut() {
            r.multithread = multithreading;
        }

        Ultra {
            settings,
            reader,
            writers,
            outs,
            number_of_threads,
            min_reader_size: 100,
            primary_thread: 0,
            settings_out,
            score_threshold,
            output_repeat_sequence,
            output_read_id,
            multithreading,
            shuffle_sequence,
            store_sequence: false,
            store_traceback: false,
            store_scores: false,
            store_profile_numbers: false,
            can_output: false,
            first_repeat: true,
            repeat_buffer: 2000,
            total_coverage: 0,
            last_seq_id: 0,
            last_rep_end: 0,
            store_trace_and_sequence,
            masks_for_seq: HashMap::new(),
            out_repeats: Vec::new(),
            analyzers,
            count: 0,
            pass_id: 0,
            threads,
            corrected_read_ids: -1,
        }
    }

    fn reader(&mut self) -> &mut FileReader {
        self.reader
            .as_mut()
            .expect("Ultra has no reader (run_without_reader)")
    }

    /// C++ `AnalyzeFile`: reads and analyses the whole input, writing the
    /// settings block and the writers' headers first. Call
    /// [`output_repeats`](Self::output_repeats)`(true)` afterwards to flush.
    pub fn analyze_file(&mut self) {
        if self.number_of_threads <= 0 {
            let _ = io::stdout().flush();
            eprintln!("Number of threads must be greater than 0.");
            std::process::exit(255);
        }
        let shuffle = self.shuffle_sequence;
        self.reader().fasta_reader.shuffle = shuffle;
        self.reader().fill_windows();

        self.output_ultra_settings();
        self.initialize_writer();

        if self.number_of_threads == 1 {
            self.analyze_file_with_thread();
        } else {
            self.analyze_file_with_workers();
        }
    }

    /// C++ `InitializeWriter`.
    pub fn initialize_writer(&mut self) {
        for i in 0..self.writers.len() {
            let _ = self.writers[i].initialize_writer(&self.settings, &mut self.outs[i]);
        }
    }

    /// C++ `GetSequenceWindow` (single-threaded): returns `seq` to the waiting
    /// pool, refills the pool when no ready window is left, and takes the
    /// ready window with the smallest read ID.
    pub fn get_sequence_window(
        &mut self,
        seq: Option<Box<SequenceWindow>>,
    ) -> Option<Box<SequenceWindow>> {
        let reader = self.reader();
        let mut should_read = false;

        if let Some(s) = seq {
            reader.add_waiting_window(s);
        }

        // Check to see if we need to read more sequence first
        if !reader.done_reading_file() && !reader.is_reading() && reader.ready_windows_size() == 0 {
            should_read = true;
            reader.set_is_reading(true);
        }

        if should_read {
            reader.fill_windows();
            reader.set_is_reading(false);
        }

        reader.get_ready_window()
    }

    /// C++ `SmallestReadID`.
    pub fn smallest_read_id(&self) -> i32 {
        let mut smallest: i32 = 100000000;
        for t in &self.threads {
            if t.smallest_read_id < smallest {
                smallest = t.smallest_read_id;
            }
        }
        smallest
    }

    /// C++ `AnalyzeFileWithThread` for thread 0 with one thread: the plain
    /// sequential loop.
    pub fn analyze_file_with_thread(&mut self) {
        let mut current_window = self.get_sequence_window(None);
        while current_window.is_some() || !self.reader().done_reading_file() {
            if let Some(w) = current_window.as_deref() {
                self.analyze_sequence_window(w, 0);
            }
            current_window = self.get_sequence_window(current_window);
        }
    }

    /// The `-t N` replacement for the C++ pthread loop. Reproduces
    /// [`analyze_file_with_thread`](Self::analyze_file_with_thread) exactly:
    /// the window pool is refilled at the same points, each batch of ready
    /// windows is analysed in parallel, and the results are applied in
    /// read-ID order, returning each window to the pool after its results,
    /// as the sequential loop does.
    fn analyze_file_with_workers(&mut self) {
        let mut analyzers = std::mem::take(&mut self.analyzers);
        let worker_settings = self.settings.clone();
        let store = self.store_trace_and_sequence;

        type Job = (usize, Box<SequenceWindow>);
        // `Err` is sent by a worker that is unwinding from a panic, so the
        // coordinator fails instead of waiting forever for its result.
        type Done = Result<(usize, Box<SequenceWindow>, Vec<Box<RepeatRegion>>), ()>;

        struct PanicNotice(mpsc::Sender<Done>);
        impl Drop for PanicNotice {
            fn drop(&mut self) {
                if std::thread::panicking() {
                    let _ = self.0.send(Err(()));
                }
            }
        }

        let (job_tx, job_rx) = mpsc::channel::<Job>();
        let job_rx = Mutex::new(job_rx);
        let (done_tx, done_rx) = mpsc::channel::<Done>();

        std::thread::scope(|scope| {
            for analyzer in analyzers.iter_mut() {
                let notice = PanicNotice(done_tx.clone());
                let job_rx = &job_rx;
                let settings = &worker_settings;
                scope.spawn(move || loop {
                    // The lock is held only while taking a job.
                    let job = job_rx.lock().expect("job queue poisoned").recv();
                    let Ok((index, window)) = job else { break };
                    let repeats = analyzer.analyze_window(&window, settings, store);
                    if notice.0.send(Ok((index, window, repeats))).is_err() {
                        break;
                    }
                });
            }
            drop(done_tx);

            loop {
                // As in GetSequenceWindow, refill only when no window is
                // ready. Within a batch the sequential loop always finds a
                // ready window, so it never refills there either.
                {
                    let reader = self.reader();
                    if !reader.done_reading_file()
                        && !reader.is_reading()
                        && reader.ready_windows_size() == 0
                    {
                        reader.set_is_reading(true);
                        reader.fill_windows();
                        reader.set_is_reading(false);
                    }
                }

                let mut batch = 0usize;
                while let Some(w) = self.reader().get_ready_window() {
                    job_tx.send((batch, w)).expect("worker threads stopped");
                    batch += 1;
                }
                if batch == 0 {
                    if self.reader().done_reading_file() {
                        break;
                    }
                    continue;
                }

                let mut pending: BTreeMap<usize, (Box<SequenceWindow>, Vec<Box<RepeatRegion>>)> =
                    BTreeMap::new();
                for next in 0..batch {
                    let (window, repeats) = loop {
                        if let Some(done) = pending.remove(&next) {
                            break done;
                        }
                        // Unwinding from here drops `job_tx`, which stops the
                        // remaining workers so the scope can join them.
                        let (i, w, r) = done_rx
                            .recv()
                            .expect("worker threads stopped")
                            .expect("a worker thread panicked");
                        pending.insert(i, (w, r));
                    };
                    self.begin_window(&window);
                    self.finish_window(0, repeats);
                    self.reader().add_waiting_window(window);
                }
            }
            drop(job_tx);
        });

        self.analyzers = analyzers;
    }

    /// C++ `Log2PvalForScore`.
    pub fn log2_pval_for_score(&self, score: f32, period: f32) -> f64 {
        log2_pval_for_score(&self.settings, score, period)
    }

    /// C++ `PvalForScore`.
    pub fn pval_for_score(&self, score: f32) -> f64 {
        pval_for_score(&self.settings, score)
    }

    /// C++ `FindRepeatsInString`: library entry point. Analyses `seq` as one
    /// window (it must fit the model, so set `window_size` to at least the
    /// longest string) and returns its repeats, which are not filtered,
    /// joined or written. Returns `None`, after a message on stderr, if the
    /// string is too long or repeats are still pending.
    pub fn find_repeats_in_string(&mut self, seq: &str) -> Option<Vec<Box<RepeatRegion>>> {
        let bytes = seq.as_bytes();
        // Make sure that the seq window can fit in the DP matrix
        let matrix_length = self.analyzers[0].model.matrix.length;
        if bytes.len() > matrix_length {
            eprintln!(
                "ULTRA model has maximum size {} but string has length {}",
                matrix_length,
                bytes.len()
            );
            return None;
        }

        // Make sure that uth isn't holding any repeats right now
        if !self.threads[0].repeats.is_empty() {
            eprintln!("ULTRA repeat array is not empty.");
            return None;
        }

        // Create sequence window from the string
        let mut seq_window = SequenceWindow::new(bytes.len(), 0);
        seq_window.length = bytes.len();
        seq_window.start = 0;
        seq_window.end = bytes.len() as u64;

        for (i, &c) in bytes.iter().enumerate() {
            seq_window.seq_mem[i] = symbol_for_char(c);
        }

        // Store and change primary thread so the repeats stay with the thread.
        let pthread = self.primary_thread;
        self.primary_thread = -10;

        self.analyze_sequence_window(&seq_window, 0);

        let repeats = std::mem::take(&mut self.threads[0].repeats);

        self.primary_thread = pthread;

        Some(repeats)
    }

    /// C++ `AnalyzeSequenceWindow` on thread 0 using analyzer `analyzer`.
    pub fn analyze_sequence_window(&mut self, sequence: &SequenceWindow, analyzer: usize) {
        self.begin_window(sequence);
        let repeats = self.analyzers[analyzer].analyze_window(
            sequence,
            &self.settings,
            self.store_trace_and_sequence,
        );
        self.finish_window(0, repeats);
    }

    /// The read-ID bookkeeping at the top of C++ `AnalyzeSequenceWindow`.
    fn begin_window(&mut self, sequence: &SequenceWindow) {
        let uth = &mut self.threads[0];
        // unsigned long narrowed to int.
        uth.active_read_id = sequence.read_id as i32;
        if uth.active_read_id < uth.smallest_read_id {
            uth.smallest_read_id = uth.active_read_id;
        }
    }

    /// The end of C++ `AnalyzeSequenceWindow`: hands the window's repeats to
    /// `out_repeats` and writes out what is complete once the buffer is
    /// over `repeat_buffer`.
    fn finish_window(&mut self, tid: usize, repeats: Vec<Box<RepeatRegion>>) {
        self.threads[tid].repeats.extend(repeats);

        if self.primary_thread == self.threads[tid].id {
            let uth = &mut self.threads[tid];
            self.out_repeats.append(&mut uth.repeats);
            uth.smallest_read_id = uth.active_read_id;

            if self.out_repeats.len() as u64 > self.repeat_buffer
                && !self.settings.disable_streaming_out
            {
                self.output_repeats(false);
            }

            // Always 0 with one logical thread; kept to mirror the C++ round robin.
            #[allow(clippy::modulo_one)]
            {
                self.primary_thread = (self.primary_thread + 1) % LOGICAL_THREADS;
            }
        }
    }

    /// C++ `OutputRepeats`: sorts `out_repeats`, then joins and writes every
    /// repeat from windows more than 3 reads back (all of them when
    /// `flush`, which also ends the writers). Flushes the outputs.
    pub fn output_repeats(&mut self, flush: bool) {
        if !flush && self.settings.disable_streaming_out {
            return;
        }
        let mut max_read_id: i32 = self.smallest_read_id() - (3 * LOGICAL_THREADS);
        if max_read_id < 0 {
            max_read_id = 0;
        }

        if flush {
            for t in self.threads.iter_mut() {
                if !t.repeats.is_empty() {
                    self.out_repeats.append(&mut t.repeats);
                }
            }
            max_read_id = 100000000;
        }

        self.sort_repeat_regions();

        // unsigned long >= int: the int converts to unsigned long.
        let max_read = max_read_id as i64 as u64;

        while let Some(back) = self.out_repeats.last() {
            if back.read_id >= max_read {
                break;
            }

            let mut r = self.out_repeats.pop().expect("non-empty");

            // Check if we need to correct overlap
            while let Some(back) = self.out_repeats.last() {
                if back.read_id >= max_read {
                    // Upstream returns here without ending the writers or
                    // flushing (unreachable when flushing).
                    self.out_repeats.push(r);
                    return;
                }

                if repeats_overlap(&r, back, true) {
                    let mut joined =
                        joint_repeat_region(&r, back).expect("overlapping repeats join");
                    if let Some(splits) = joined.splits.as_mut() {
                        let consensi = joined.consensi.as_mut().expect("consensi is null");
                        validate_splits(consensi, splits, 0.85);
                    }
                    r = joined;
                    self.out_repeats.pop();
                } else {
                    break;
                }
            }

            r.sort_consensi();
            // unsigned long / int, compared with unsigned long long.
            if r.region_score >= self.score_threshold
                && r.repeat_length / (r.repeat_period as i64 as u64) >= self.settings.min_units
            {
                self.output_repeat(&r, false);
            }
        }

        if flush {
            for i in 0..self.writers.len() {
                let _ = self.writers[i].end_writer(&self.settings, &mut self.outs[i]);
            }
        }

        for out in self.outs.iter_mut() {
            let _ = out.flush();
        }
    }

    /// C++ `OutputULTRASettings`: the settings block, to the `.settings`
    /// file (which is then closed) or stdout, unless `hide_settings`.
    pub fn output_ultra_settings(&mut self) {
        let text = format!(
            "{{\"Version\": \"{}\", \n\"Parameters\": {{\n{}}}}}\n",
            ULTRA_VERSION_STRING,
            self.settings.json_string()
        );
        match self.settings_out.take() {
            Some(mut f) => {
                if !self.settings.hide_settings {
                    let _ = f.write_all(text.as_bytes());
                }
                let _ = f.flush();
            }
            None => {
                if !self.settings.hide_settings {
                    let _ = io::stdout().write_all(text.as_bytes());
                }
            }
        }
    }

    /// C++ `OutputRepeat`: adds the repeat to the coverage count, writes it
    /// unless `suppress_out`, and records it for the mask.
    pub fn output_repeat(&mut self, r: &RepeatRegion, _is_sub_rep: bool) {
        // Coverage is counted per sequence. Upstream 1.2.2 compares the read
        // id with a counter it never updates, so it stops merging overlapping
        // repeats after the first window.
        if r.sequence_id != self.last_seq_id {
            self.last_seq_id = r.sequence_id;
            self.last_rep_end = 0;
        }

        let end = r.sequence_start.wrapping_add(r.repeat_length);
        let mut rep_start = r.sequence_start;
        if rep_start < self.last_rep_end {
            rep_start = self.last_rep_end;
        }

        if rep_start < end {
            self.total_coverage = self
                .total_coverage
                .wrapping_add(end.wrapping_sub(rep_start));
        }

        if end > self.last_rep_end {
            self.last_rep_end = end;
        }

        if !self.settings.suppress_out {
            for i in 0..self.writers.len() {
                let _ = self.writers[i].write_repeat(&self.settings, r, &mut self.outs[i]);
            }
        }

        if self.settings.produce_mask {
            self.store_mask_for_region(r);
        }
    }

    /// C++ `StoreMaskForRegion`.
    pub fn store_mask_for_region(&mut self, r: &RepeatRegion) {
        self.masks_for_seq
            .entry(r.sequence_id)
            .or_default()
            .push(MRegion {
                start: r.sequence_start,
                end: r.repeat_length.wrapping_add(r.sequence_start),
            });
    }

    /// C++ `SortRepeatRegions`: libstdc++ `std::sort` with
    /// [`compare_repeat_order`]; the order of tied repeats matters.
    pub fn sort_repeat_regions(&mut self) {
        std_sort(&mut self.out_repeats, |a, b| compare_repeat_order(a, b));
    }

    /// C++ `Coverage`.
    pub fn coverage(&self) -> u64 {
        self.total_coverage
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(seq_id: u64, start: u64, read_id: u64) -> Box<RepeatRegion> {
        let mut r = Box::new(RepeatRegion::new());
        r.sequence_id = seq_id;
        r.sequence_start = start;
        r.read_id = read_id;
        r
    }

    #[test]
    fn compare_repeat_order_is_descending() {
        let a = region(1, 10, 0);
        let b = region(1, 20, 0);
        let c = region(2, 0, 0);
        assert!(compare_repeat_order(&b, &a));
        assert!(!compare_repeat_order(&a, &b));
        assert!(compare_repeat_order(&c, &b));
        assert!(!compare_repeat_order(&a, &a));
    }

    fn library_settings() -> Settings {
        let mut s = Settings::new();
        s.run_without_reader = true;
        s.window_size = 1000;
        s.overlap = 0;
        // The upstream example's "-p 25 -i 3 -d 3" loses "-p" as argv[0].
        s.max_insert = 3;
        s.max_delete = 3;
        s.assign_settings();
        s
    }

    // Expected output from upstream examples/library_example.cpp.
    #[test]
    fn find_repeats_in_string_matches_library_example() {
        let s = library_settings();
        let mut ultra = Ultra::new(&s);
        let mut find = |s: &str| -> Vec<String> {
            ultra
                .find_repeats_in_string(s)
                .unwrap()
                .iter()
                .map(|r| {
                    format!(
                        "{} {} {} {}",
                        r.window_start,
                        r.repeat_length,
                        r.repeat_period,
                        String::from_utf8(r.get_consensus()).unwrap()
                    )
                })
                .collect()
        };
        assert_eq!(find("aaaaaaaaaaaaaaaaaaaaaaaaaa"), ["0 26 1 A"]);
        assert_eq!(
            find("aggtaaggtaaggtaaggtaaggtaagcggtataacatacagatctgactactactactactactactactactac"),
            ["0 27 5 AGGTA", "47 29 3 ACT"]
        );
        assert!(ultra.threads[0].repeats.is_empty());
    }

    #[test]
    fn find_repeats_in_string_rejects_long_strings() {
        let s = library_settings();
        let mut ultra = Ultra::new(&s);
        let long = "a".repeat(5000);
        assert!(ultra.find_repeats_in_string(&long).is_none());
    }
}
