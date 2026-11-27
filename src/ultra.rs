//! The ULTRA driver (`ultra.hpp` / `ultra.cpp`): reads sequence windows,
//! runs the model on each, joins repeats found in overlapping windows, and
//! writes them out.
//!
//! # Threads
//!
//! The C++ runs `-t N` pthreads that each pull windows, analyse them and
//! hand their repeats over through unsynchronised shared state, so its output
//! with more than one thread varies from run to run. This port gives the
//! single-threaded output for every thread count: windows are read in file
//! order (a few at a time, as pool windows come free, so workers do not wait
//! on a pool refill), analysed by up to `N` worker threads, and the results
//! are consumed strictly in read-ID order by the calling thread, which does
//! exactly what the C++ does after analysing a window with one thread. What a
//! window contains and its read id depend only on file order, not on when it
//! is read. Every quantity in
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
    log2_pval_for_score, pval_for_score, BEDFileWriter, BPAFFileWriter, JSONFileWriter, OutSink,
    RepeatFileWriter,
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
    pub active_read_id: u64,
    pub smallest_read_id: u64,
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
        let mut r = get_next_repeat(sequence, &model.matrix, &mut i, settings.upstream_compat);

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
            r = get_next_repeat(sequence, &model.matrix, &mut i, settings.upstream_compat);
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

/// The order [`Ultra::output_repeats`] takes repeats in, outside
/// upstream-compat mode: `Less` when `lhs` is taken first. It is the order of
/// [`compare_repeat_order`] (sequence ID, then start) with its ties broken:
/// by period, then longest first, then by read ID.
///
/// `compare_repeat_order` leaves repeats with one start tied, and the sort
/// that uses it is not stable, so which of them came first depended on what
/// else was being sorted: a run with `--disable_streaming_out` could join or
/// write them in another order than a streaming run. Two repeats from one
/// window differ in start or period, because each starts one period and one
/// column before its first repeat column. A joined repeat keeps the start,
/// period and read ID of its first piece, so a streaming run that puts it
/// back to wait sorts it by the same three.
///
/// Two repeats with the same start and period come from two windows that
/// both cover that start, so it lies in their overlap. Longest first makes
/// the longer one the first piece of their join, which keeps it whole (see
/// [`joint_repeat_region`]). They are next to each other in this order, so
/// they are joined to each other before either is joined to anything else.
pub fn repeat_output_order(lhs: &RepeatRegion, rhs: &RepeatRegion) -> std::cmp::Ordering {
    lhs.sequence_id
        .cmp(&rhs.sequence_id)
        .then(lhs.sequence_start.cmp(&rhs.sequence_start))
        .then(lhs.repeat_period.cmp(&rhs.repeat_period))
        .then(rhs.repeat_length.cmp(&lhs.repeat_length))
        .then(lhs.read_id.cmp(&rhs.read_id))
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

/// Where `path` points once links and `..` are followed, for a file that may not exist
/// yet (then its directory is resolved).
fn resolved(path: &str) -> std::path::PathBuf {
    let p = std::path::Path::new(path);
    if let Ok(c) = std::fs::canonicalize(p) {
        return c;
    }
    let dir = p.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(std::path::Path::new("."));
    match (std::fs::canonicalize(dir), p.file_name()) {
        (Ok(d), Some(f)) => d.join(f),
        _ => p.to_path_buf(),
    }
}

/// The files a run with settings `s` writes: (path, what it is).
fn output_paths(s: &Settings) -> Vec<(String, &'static str)> {
    let mut v = Vec::new();
    if !s.out_file.is_empty() {
        let kinds = [(s.ultra_out, "tsv"), (s.json_out, "json"), (s.bed_out, "bed"), (s.bpaf_out, "bpaf")];
        if kinds.iter().filter(|k| k.0).count() > 1 {
            v.extend(kinds.iter().filter(|k| k.0).map(|k| (format!("{}.{}", s.out_file, k.1), "output file")));
        } else {
            v.push((s.out_file.clone(), "output file"));
        }
        if !s.hide_settings {
            v.push((format!("{}.settings", s.out_file), "settings output file"));
        }
    }
    if !s.mask_file.is_empty() {
        v.push((s.mask_file.clone(), "mask file"));
    }
    v
}

/// Before any output is created (which empties an existing file): the input must open,
/// and no output may be the input or another output. Exits with status 255 otherwise.
fn check_paths(s: &Settings) {
    let fail = |msg: String| -> ! {
        let _ = io::stdout().flush();
        eprintln!("{}", msg);
        std::process::exit(255)
    };
    let outs: Vec<(std::path::PathBuf, String, &str)> =
        output_paths(s).into_iter().map(|(p, what)| (resolved(&p), p, what)).collect();
    if !s.run_without_reader {
        // a pipe is not opened here: that would take from what the reader is to read
        let ok = match std::fs::metadata(&s.in_file) {
            Ok(m) if m.is_file() => File::open(&s.in_file).is_ok(),
            other => other.is_ok(),
        };
        if !ok {
            fail(format!("Unable to open file: {}", s.in_file));
        }
        // these read the input again from the start, which a pipe cannot give
        let again = [(!s.mask_file.is_empty(), "--mask"), (s.estimate_fdr, "--fdr"), (s.tune, "--tune")];
        if let Some(opt) = again.iter().find(|a| a.0) {
            if !std::fs::metadata(&s.in_file).is_ok_and(|m| m.is_file()) {
                fail(format!("{} reads the input more than once and needs a regular file: {}", opt.1, s.in_file));
            }
        }
        let input = resolved(&s.in_file);
        if let Some(o) = outs.iter().find(|o| o.0 == input) {
            fail(format!("The {} {} is the input file", o.2, o.1));
        }
    }
    for (i, a) in outs.iter().enumerate() {
        if let Some(b) = outs[..i].iter().find(|b| b.0 == a.0) {
            fail(format!("The {} and the {} are the same file: {}", b.2, a.2, a.1));
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
        check_paths(&settings);
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
            if settings.bpaf_out {
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

                if settings.bpaf_out {
                    let bpaf_path = format!("{}.bpaf", settings.out_file);
                    outs.push(OutSink::file(open_output(&bpaf_path, "output file")));
                    writers.push(Box::new(BPAFFileWriter::new()));
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
                } else if settings.bpaf_out {
                    writers.push(Box::new(BPAFFileWriter::new()));
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
            model.upstream_indel_windows = settings.upstream_compat;

            model.tp_match_to_insertion = settings.transition_ri;
            model.tp_match_to_deletion = settings.transition_rd;

            model.tp_consecutive_insertion = settings.transition_ii;
            model.tp_consecutive_deletion = settings.transition_dd;

            model.calculate_scores();

            let mut splitter = SplitWindow::new();
            splitter.upstream_compat = settings.upstream_compat;
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
    /// Returns the first error from writing the output; the analysis stops
    /// there.
    pub fn analyze_file(&mut self) -> io::Result<()> {
        if self.number_of_threads <= 0 {
            let _ = io::stdout().flush();
            eprintln!("Number of threads must be greater than 0.");
            std::process::exit(255);
        }
        let shuffle = self.shuffle_sequence;
        self.reader().fasta_reader.shuffle = shuffle;
        if self.number_of_threads == 1 {
            self.reader().fill_windows();
        } else {
            // The workers start on the first windows while the rest are
            // read; see analyze_file_with_workers.
            let n = 2 * self.number_of_threads as usize;
            self.reader().fill_windows_up_to(n);
        }

        self.output_ultra_settings()?;
        self.initialize_writer()?;

        if self.number_of_threads == 1 {
            self.analyze_file_with_thread()
        } else {
            self.analyze_file_with_workers()
        }
    }

    /// C++ `InitializeWriter`.
    pub fn initialize_writer(&mut self) -> io::Result<()> {
        for i in 0..self.writers.len() {
            self.writers[i].initialize_writer(&self.settings, &mut self.outs[i])?;
        }
        Ok(())
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
    ///
    /// The C++ starts from 100,000,000 and keeps read ids in an int, so past
    /// that many windows nothing more was written; here there is no cap.
    pub fn smallest_read_id(&self) -> u64 {
        let mut smallest: u64 = u64::MAX;
        for t in &self.threads {
            if t.smallest_read_id < smallest {
                smallest = t.smallest_read_id;
            }
        }
        smallest
    }

    /// C++ `AnalyzeFileWithThread` for thread 0 with one thread: the plain
    /// sequential loop.
    pub fn analyze_file_with_thread(&mut self) -> io::Result<()> {
        let mut current_window = self.get_sequence_window(None);
        while current_window.is_some() || !self.reader().done_reading_file() {
            if let Some(w) = current_window.as_deref() {
                self.analyze_sequence_window(w, 0)?;
            }
            current_window = self.get_sequence_window(current_window);
        }
        Ok(())
    }

    /// The `-t N` replacement for the C++ pthread loop. Gives the output of
    /// [`analyze_file_with_thread`](Self::analyze_file_with_thread): windows
    /// are read in file order and their results are applied in read-ID
    /// order, each window returning to the pool after its results, as the
    /// sequential loop does. Reading is not tied to the sequential loop's
    /// refill points: free windows are refilled a few at a time as they come
    /// back, keeping two windows per worker queued, so the workers never wait
    /// for a whole pool to be applied. Window contents and read ids depend
    /// only on file order, so this does not change the output.
    fn analyze_file_with_workers(&mut self) -> io::Result<()> {
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

        let workers = analyzers.len().max(1);
        let queue_target = 2 * workers;

        let result = std::thread::scope(|scope| {
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

            // Read order of the next window to hand out and of the next one
            // whose results are applied.
            let mut next_sent = 0usize;
            let mut next_applied = 0usize;
            let mut pending: BTreeMap<usize, (Box<SequenceWindow>, Vec<Box<RepeatRegion>>)> =
                BTreeMap::new();
            loop {
                // Keep two windows per worker queued. Windows are read as
                // soon as the pool has free ones, so workers do not wait for
                // a whole pool to be applied and refilled. Reading order and
                // read ids do not depend on when this happens.
                while next_sent - next_applied < queue_target {
                    if self.reader().ready_windows_size() == 0 {
                        let reader = self.reader();
                        if reader.done_reading_file() {
                            break;
                        }
                        reader.set_is_reading(true);
                        reader.fill_windows_up_to(workers);
                        reader.set_is_reading(false);
                    }
                    let mut sent = false;
                    while let Some(w) = self.reader().get_ready_window() {
                        job_tx.send((next_sent, w)).expect("worker threads stopped");
                        next_sent += 1;
                        sent = true;
                    }
                    if !sent {
                        // The pool has no free window, or the file is done.
                        break;
                    }
                }

                if next_applied == next_sent {
                    if self.reader().done_reading_file() && self.reader().ready_windows_size() == 0 {
                        break;
                    }
                    continue;
                }

                // Apply results strictly in read order, as the sequential
                // loop does, returning each window to the pool afterwards.
                let (window, repeats) = loop {
                    if let Some(done) = pending.remove(&next_applied) {
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
                let written = self.finish_window(0, repeats);
                self.reader().add_waiting_window(window);
                next_applied += 1;
                // Returning drops `job_tx`, which stops the workers once
                // the queued windows are done, so the scope can join them.
                written?;
            }
            drop(job_tx);
            Ok(())
        });

        self.analyzers = analyzers;
        result
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
        // A window of n symbols fills matrix columns 1..=n and then column
        // n + 1 without emission, so n must be below the matrix length (the
        // C++ only rejects n above it, and n equal to it ran off the end).
        let matrix_length = self.analyzers[0].model.matrix.length;
        if bytes.len() >= matrix_length {
            eprintln!(
                "ULTRA model has maximum size {} but string has length {}",
                matrix_length.saturating_sub(1),
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

        // With primary_thread at -10 the repeats are not passed on for
        // output, so nothing is written here.
        let written = self.analyze_sequence_window(&seq_window, 0);

        let repeats = std::mem::take(&mut self.threads[0].repeats);

        self.primary_thread = pthread;

        if let Err(e) = written {
            eprintln!("Error writing output: {e}");
            return None;
        }
        Some(repeats)
    }

    /// C++ `AnalyzeSequenceWindow` on thread 0 using analyzer `analyzer`.
    /// Returns the error from writing out finished repeats, if that fails.
    pub fn analyze_sequence_window(&mut self, sequence: &SequenceWindow, analyzer: usize) -> io::Result<()> {
        self.begin_window(sequence);
        let repeats = self.analyzers[analyzer].analyze_window(
            sequence,
            &self.settings,
            self.store_trace_and_sequence,
        );
        self.finish_window(0, repeats)
    }

    /// The read-ID bookkeeping at the top of C++ `AnalyzeSequenceWindow`.
    fn begin_window(&mut self, sequence: &SequenceWindow) {
        let uth = &mut self.threads[0];
        uth.active_read_id = sequence.read_id;
        if uth.active_read_id < uth.smallest_read_id {
            uth.smallest_read_id = uth.active_read_id;
        }
    }

    /// The end of C++ `AnalyzeSequenceWindow`: hands the window's repeats to
    /// `out_repeats` and writes out what is complete once the buffer is
    /// over `repeat_buffer`.
    fn finish_window(&mut self, tid: usize, repeats: Vec<Box<RepeatRegion>>) -> io::Result<()> {
        self.threads[tid].repeats.extend(repeats);

        if self.primary_thread == self.threads[tid].id {
            let uth = &mut self.threads[tid];
            self.out_repeats.append(&mut uth.repeats);
            uth.smallest_read_id = uth.active_read_id;

            if self.out_repeats.len() as u64 > self.repeat_buffer
                && !self.settings.disable_streaming_out
            {
                self.output_repeats(false)?;
            }

            // Always 0 with one logical thread; kept to mirror the C++ round robin.
            #[allow(clippy::modulo_one)]
            {
                self.primary_thread = (self.primary_thread + 1) % LOGICAL_THREADS;
            }
        }
        Ok(())
    }

    /// C++ `OutputRepeats`: sorts `out_repeats`, then joins and writes every
    /// repeat from windows more than 3 reads back (all of them when
    /// `flush`, which also ends the writers). Flushes the outputs. Stops at
    /// the first write error and returns it.
    pub fn output_repeats(&mut self, flush: bool) -> io::Result<()> {
        if !flush && self.settings.disable_streaming_out {
            return Ok(());
        }
        // Repeats from windows at least this far on may still be joined
        // with the next windows' repeats, so they wait (unless `flush`).
        let max_read = self.smallest_read_id().saturating_sub(3 * LOGICAL_THREADS as u64);

        if flush {
            for t in self.threads.iter_mut() {
                if !t.repeats.is_empty() {
                    self.out_repeats.append(&mut t.repeats);
                }
            }
            // A flush writes every repeat, so the read-id tests below are
            // skipped. (The C++ sets max_read to 100,000,000 here, which
            // dropped every repeat from a later window and could leave the
            // JSON array unclosed. A u64::MAX marker would drop the repeats
            // of find_repeats_in_string, whose read id is UNSET_ID.)
        }

        self.sort_repeat_regions();

        while let Some(back) = self.out_repeats.last() {
            if !flush && back.read_id >= max_read {
                break;
            }

            let mut r = self.out_repeats.pop().expect("non-empty");

            // Check if we need to correct overlap
            while let Some(back) = self.out_repeats.last() {
                if !flush && back.read_id >= max_read {
                    // Upstream returns here without ending the writers or
                    // flushing (unreachable when flushing).
                    self.out_repeats.push(r);
                    return Ok(());
                }

                // The repeat to join r with: the next one if it overlaps r
                // with the same period. Otherwise (not as the C++, which
                // stops there) the first later one that still starts within
                // r and has r's period: a repeat of another period found
                // inside r, such as a microsatellite in a satellite's unit,
                // sorts between r and its continuation from the next window
                // and kept the two from being joined, so they were written
                // as two overlapping records.
                let mut take = None;
                if repeats_overlap(&r, back, true) {
                    take = Some(self.out_repeats.len() - 1);
                } else if !self.settings.upstream_compat {
                    let r_end = r.sequence_start.wrapping_add(r.repeat_length);
                    for k in (0..self.out_repeats.len() - 1).rev() {
                        let c = &self.out_repeats[k];
                        if c.sequence_id != r.sequence_id || c.sequence_start > r_end {
                            break;
                        }
                        if !flush && c.read_id >= max_read {
                            // Not final yet: wait, as for the next repeat.
                            self.out_repeats.push(r);
                            return Ok(());
                        }
                        if repeats_overlap(&r, c, true) {
                            take = Some(k);
                            break;
                        }
                    }
                }
                let Some(k) = take else { break };
                let mut joined = joint_repeat_region(&r, &self.out_repeats[k], self.settings.upstream_compat)
                    .expect("overlapping repeats join");
                if let Some(splits) = joined.splits.as_mut() {
                    let consensi = joined.consensi.as_mut().expect("consensi is null");
                    validate_splits(consensi, splits, 0.85);
                }
                r = joined;
                self.out_repeats.remove(k);
            }

            r.sort_consensi();
            // unsigned long / int, compared with unsigned long long.
            // --min_length: documented, but the C++ never applies it.
            let long_enough = self.settings.upstream_compat || r.repeat_length >= self.settings.min_length;
            if r.region_score >= self.score_threshold
                && long_enough
                && r.repeat_length / (r.repeat_period as i64 as u64) >= self.settings.min_units
            {
                self.output_repeat(&r, false)?;
            }
        }

        if flush {
            for i in 0..self.writers.len() {
                self.writers[i].end_writer(&self.settings, &mut self.outs[i])?;
            }
        }

        for out in self.outs.iter_mut() {
            out.flush()?;
        }
        Ok(())
    }

    /// C++ `OutputULTRASettings`: the settings block, to the `.settings`
    /// file (which is then closed) or stdout, unless `hide_settings`.
    pub fn output_ultra_settings(&mut self) -> io::Result<()> {
        let text = format!(
            "{{\"Version\": \"{}\", \n\"Parameters\": {{\n{}}}}}\n",
            ULTRA_VERSION_STRING,
            self.settings.json_string()
        );
        match self.settings_out.take() {
            Some(mut f) => {
                if !self.settings.hide_settings {
                    f.write_all(text.as_bytes())?;
                }
                f.flush()
            }
            None => {
                if !self.settings.hide_settings {
                    io::stdout().lock().write_all(text.as_bytes())?;
                }
                Ok(())
            }
        }
    }

    /// C++ `OutputRepeat`: adds the repeat to the coverage count, writes it
    /// unless `suppress_out`, and records it for the mask. Returns the first
    /// writer's error, if any; the repeat is then not recorded for the mask.
    pub fn output_repeat(&mut self, r: &RepeatRegion, _is_sub_rep: bool) -> io::Result<()> {
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
                self.writers[i].write_repeat(&self.settings, r, &mut self.outs[i])?;
            }
        }

        if self.settings.produce_mask {
            self.store_mask_for_region(r);
        }
        Ok(())
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

    /// C++ `SortRepeatRegions`: puts the repeat to take first at the back.
    /// In upstream-compat mode this is libstdc++ `std::sort` with
    /// [`compare_repeat_order`], where the order of tied repeats matters and
    /// depends on the arrangement before the sort. Otherwise the order is
    /// [`repeat_output_order`], reversed.
    pub fn sort_repeat_regions(&mut self) {
        if self.settings.upstream_compat {
            std_sort(&mut self.out_repeats, |a, b| compare_repeat_order(a, b));
        } else {
            self.out_repeats.sort_by(|a, b| repeat_output_order(b, a));
        }
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

    /// Repeats with one start are taken by period, then longest first, then
    /// by read ID, whatever order the vector held them in.
    #[test]
    fn repeat_output_order_breaks_ties() {
        use std::cmp::Ordering::{Equal, Less};
        let piece = |start: u64, period: i32, length: u64, read_id: u64| {
            let mut r = region(1, start, read_id);
            r.repeat_period = period;
            r.repeat_length = length;
            r
        };
        let short = piece(50, 3, 20, 4);
        let long = piece(50, 3, 90, 5);
        let other_period = piece(50, 2, 10, 5);
        let later = piece(51, 2, 500, 0);
        assert_eq!(repeat_output_order(&long, &short), Less);
        assert_eq!(repeat_output_order(&other_period, &long), Less);
        assert_eq!(repeat_output_order(&short, &later), Less);
        assert_eq!(repeat_output_order(&short, &piece(50, 3, 20, 5)), Less);
        assert_eq!(repeat_output_order(&short, &short), Equal);
        assert_eq!(repeat_output_order(&piece(50, 9, 1, 9), &region(2, 0, 0)), Less);

        let mut s = library_settings();
        s.suppress_out = true;
        // The vector after the sort: the repeat to take first is at the back.
        let sorted = |arrangement: [&Box<RepeatRegion>; 4]| {
            let mut ultra = Ultra::new(&s);
            ultra.out_repeats = arrangement.into_iter().cloned().collect();
            ultra.sort_repeat_regions();
            ultra.out_repeats.iter().map(|r| (r.sequence_start, r.repeat_period, r.repeat_length)).collect::<Vec<_>>()
        };
        let want = vec![(51, 2, 500), (50, 3, 20), (50, 3, 90), (50, 2, 10)];
        assert_eq!(sorted([&short, &long, &other_period, &later]), want);
        assert_eq!(sorted([&later, &other_period, &long, &short]), want);
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

    /// --min_length drops shorter repeats (the C++ never applied it).
    #[test]
    fn min_length_filters_repeats() {
        for (compat, want) in [(false, 50), (true, 58)] {
            let mut s = library_settings();
            s.suppress_out = true;
            s.min_length = 10;
            s.upstream_compat = compat;
            let mut ultra = Ultra::new(&s);
            for (k, len) in [8u64, 10, 40].into_iter().enumerate() {
                let mut r = region(1, 1000 * k as u64, 0);
                r.repeat_period = 2;
                r.repeat_length = len;
                r.region_score = 100.0;
                ultra.out_repeats.push(r);
            }
            ultra.output_repeats(true).unwrap();
            assert_eq!(ultra.coverage(), want, "compat {compat}");
        }
    }

    /// A string as long as the matrix does not fit (the column after the
    /// last symbol is needed too); one shorter does.
    #[test]
    fn find_repeats_in_string_needs_one_spare_column() {
        let s = library_settings();
        let mut ultra = Ultra::new(&s);
        let n = ultra.analyzers[0].model.matrix.length;
        assert!(ultra.find_repeats_in_string(&"a".repeat(n)).is_none());
        let found = ultra.find_repeats_in_string(&"a".repeat(n - 1)).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].repeat_length as usize, n - 1);
    }

    /// Repeats from windows numbered 100,000,000 and up are written. The C++
    /// capped read ids there, so the final flush dropped them.
    #[test]
    fn repeats_from_late_windows_are_written() {
        let mut s = library_settings();
        s.suppress_out = true;
        let mut ultra = Ultra::new(&s);
        for (k, read_id) in [99_999_999u64, 100_000_000, 3_000_000_000].into_iter().enumerate() {
            let mut r = region(1, 1000 * k as u64, read_id);
            r.repeat_period = 2;
            r.repeat_length = 10;
            r.region_score = 100.0;
            ultra.out_repeats.push(r);
        }
        ultra.threads[0].smallest_read_id = 3_000_000_000;
        // Streaming writes the first; the second waits because the repeat
        // after it is not final yet and might join it.
        ultra.output_repeats(false).unwrap();
        assert_eq!(ultra.coverage(), 10);
        // The flush writes the rest (the C++ wrote neither: their read ids
        // are not below its 100,000,000 cap).
        ultra.output_repeats(true).unwrap();
        assert_eq!(ultra.coverage(), 30);
        assert!(ultra.out_repeats.is_empty());
    }

    /// Counts the repeats and trailers written; fails every write after the
    /// first `ok` with an error of kind `fail`.
    struct Counter {
        repeats: std::sync::Arc<Mutex<(usize, usize)>>,
        ok: usize,
        fail: io::ErrorKind,
    }

    impl RepeatFileWriter for Counter {
        fn write_repeat(&mut self, _: &Settings, _: &RepeatRegion, _: &mut dyn Write) -> io::Result<()> {
            let mut n = self.repeats.lock().unwrap();
            if n.0 == self.ok {
                return Err(io::Error::new(self.fail, "test sink failed"));
            }
            n.0 += 1;
            Ok(())
        }
        fn end_writer(&mut self, _: &Settings, _: &mut dyn Write) -> io::Result<()> {
            self.repeats.lock().unwrap().1 += 1;
            Ok(())
        }
    }

    fn counted(ok: usize, fail: io::ErrorKind) -> (Ultra, std::sync::Arc<Mutex<(usize, usize)>>) {
        let mut ultra = Ultra::new(&library_settings());
        let n = std::sync::Arc::new(Mutex::new((0, 0)));
        ultra.writers = vec![Box::new(Counter { repeats: n.clone(), ok, fail })];
        ultra.outs.truncate(1);
        (ultra, n)
    }

    const TWO_REPEATS: &str =
        "aggtaaggtaaggtaaggtaaggtaagcggtataacatacagatctgactactactactactactactactactac";

    /// The repeats of find_repeats_in_string, whose read id is UNSET_ID, are
    /// written by a flush. A flush used to stop at read id u64::MAX, which is
    /// UNSET_ID, so it wrote none of them but still ended the writer.
    #[test]
    fn a_flush_writes_repeats_from_find_repeats_in_string() {
        let (mut ultra, n) = counted(usize::MAX, io::ErrorKind::Other);
        let found = ultra.find_repeats_in_string(TWO_REPEATS).unwrap();
        assert!(found.iter().all(|r| r.read_id == crate::sequence_window::UNSET_ID));
        ultra.out_repeats.extend(found);
        ultra.output_repeats(true).unwrap();
        assert_eq!(*n.lock().unwrap(), (2, 1));
        assert!(ultra.out_repeats.is_empty());
        assert_eq!(ultra.coverage(), 27 + 29);
    }

    /// A writer's error comes back to the caller, for a closed pipe as for
    /// any other error. It used to end the process (status 141 or 255).
    #[test]
    fn write_errors_are_returned() {
        for kind in [io::ErrorKind::BrokenPipe, io::ErrorKind::Other] {
            let (mut ultra, n) = counted(1, kind);
            let found = ultra.find_repeats_in_string(TWO_REPEATS).unwrap();
            ultra.out_repeats.extend(found);
            let e = ultra.output_repeats(true).unwrap_err();
            assert_eq!(e.kind(), kind);
            // One repeat was written; the writer was not ended.
            assert_eq!(*n.lock().unwrap(), (1, 0));
        }
    }
}
