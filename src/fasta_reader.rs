//! Streams a FASTA file into [`SequenceWindow`]s. Port of `FASTAReader.hpp` /
//! `FASTAReader.cpp`.
//!
//! Windows are owned by the reader's two pools while idle and handed out as
//! `Box<SequenceWindow>`. The pools have no internal locking; see
//! [`crate::file_reader`] for the locking contract.

use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::sequence_window::{shuffle_sequence_window, SequenceWindow, READ_ALL};
use crate::symbol::Symbol;

/// Reads a FASTA file window by window.
pub struct FASTAReader {
    /// The open input. `None` once the C++ would have closed the stream (or,
    /// for the random-windows constructor, never opened it).
    pub file: Option<BufReader<File>>,

    /// The current line, without its `\n` (a `\r` is kept, as with
    /// `std::getline`).
    pub line: Vec<u8>,
    /// Position in `line` to read from next, `READ_ALL` when the line is used
    /// up, or `-1` when it is a header that has not been processed yet.
    pub line_place: i64,
    pub total_seq_length: u64,

    /// Min-heap on `read_id`, laid out exactly as libstdc++'s
    /// `std::push_heap`/`std::pop_heap` with `CompareSequenceWindows`.
    pub ready_windows: Vec<Box<SequenceWindow>>,
    /// Idle windows, used as a stack (push and pop at the back).
    pub waiting_windows: Vec<Box<SequenceWindow>>,

    /// If set (`max_windows == -1`), `max_windows` is ignored and the whole
    /// file is read in the first `fill_windows` call.
    pub read_whole_file: bool,
    pub shuffle: bool,
    pub max_windows: i32,
    pub max_seq_length: u64,
    pub max_overlap_length: u64,

    pub a_pctg: f64,
    pub t_pctg: f64,
    pub c_pctg: f64,
    pub g_pctg: f64,

    /// `max_overlap_length + 1` symbols; the first `overlap_length` are the
    /// tail of the previous window.
    pub overlap_buffer: Vec<Symbol>,
    pub overlap_length: u64,

    pub sequence_name: String,
    pub sequence_id: u64,
    /// Read ids are assigned to every window read, including empty ones that
    /// are returned to the waiting pool, so they may not be contiguous.
    pub read_id: u64,
    pub symbols_read_in_seq: u64,

    pub done_reading_file: bool,
    pub is_reading: bool,

    /// Kept for parity with the C++; it only switched the internal pthread
    /// locks on and off, which this port does not have.
    pub multithread: bool,
}

/// `std::getline(file, line, '\n')`: returns false only when nothing could be
/// extracted. The delimiter is dropped; everything else (including `\r`) is
/// kept. A read error is treated like end of file, as the stream would go bad
/// and fail.
fn getline(file: &mut Option<BufReader<File>>, line: &mut Vec<u8>) -> bool {
    line.clear();
    let f = match file.as_mut() {
        Some(f) => f,
        None => return false,
    };
    match f.read_until(b'\n', line) {
        Ok(0) | Err(_) => false,
        Ok(_) => {
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            true
        }
    }
}

/// `std::push_heap` (libstdc++ `__push_heap`) with comparator
/// `lhs.read_id > rhs.read_id`, done with swaps instead of a hole.
fn sift_up(heap: &mut [Box<SequenceWindow>], mut hole: usize, top: usize) {
    while hole > top {
        let parent = (hole - 1) / 2;
        if heap[parent].read_id > heap[hole].read_id {
            heap.swap(parent, hole);
            hole = parent;
        } else {
            break;
        }
    }
}

fn push_heap(heap: &mut [Box<SequenceWindow>]) {
    let len = heap.len();
    if len > 1 {
        sift_up(heap, len - 1, 0);
    }
}

/// `std::pop_heap` (libstdc++ `__pop_heap` / `__adjust_heap`): moves the
/// front to the back and restores the heap on the rest.
fn pop_heap(heap: &mut [Box<SequenceWindow>]) {
    let len = heap.len();
    if len < 2 {
        return;
    }
    heap.swap(0, len - 1);
    let len = len - 1;
    let comp = |h: &[Box<SequenceWindow>], a: usize, b: usize| h[a].read_id > h[b].read_id;

    let mut hole = 0usize;
    let mut second_child = 0usize;
    while second_child < len.saturating_sub(1) / 2 {
        second_child = 2 * (second_child + 1);
        if comp(heap, second_child, second_child - 1) {
            second_child -= 1;
        }
        heap.swap(hole, second_child);
        hole = second_child;
    }
    if len & 1 == 0 && second_child == (len - 2) / 2 {
        second_child = 2 * (second_child + 1);
        heap.swap(hole, second_child - 1);
        hole = second_child - 1;
    }
    sift_up(heap, hole, 0);
}

impl FASTAReader {
    /// Handles a line starting with `>`: the rest of the line (including any
    /// trailing `\r`) becomes the sequence name. Returns true if it was a
    /// header.
    pub fn read_special_line(&mut self) -> bool {
        // C++ `line[0]` on an empty string reads the terminating NUL.
        if self.line.first() == Some(&b'>') {
            self.sequence_name = String::from_utf8_lossy(&self.line[1..]).into_owned();
            self.sequence_id += 1;
            self.symbols_read_in_seq = 0;
            self.overlap_length = 0;
            self.line_place = READ_ALL;
            return true;
        }

        false
    }

    /// Saves the last `ol` new symbols of `window` (fewer if the window is
    /// shorter) as the overlap for the next window.
    pub fn copy_overlap_buffer_from_window(&mut self, window: &SequenceWindow, ol: u64) -> bool {
        let mut ol = ol;
        if ol > window.length as u64 {
            ol = window.length as u64;
        }

        self.overlap_length = ol;

        let ol = ol as usize;
        let end = window.new_seq_offset + window.length;
        self.overlap_buffer[..ol].copy_from_slice(&window.seq_mem[end - ol..end]);

        true
    }

    /// Reads input into every waiting window (or, in read-whole-file mode,
    /// into as many new windows as the file needs) and moves them to the ready
    /// pool. Empty windows go back to the waiting pool.
    pub fn fill_windows(&mut self) -> bool {
        if self.done_reading_file {
            self.file = None;
            return true;
        }

        if self.read_whole_file {
            while !self.done_reading_file || self.line_place != READ_ALL {
                let mut win = if self.waiting_windows.is_empty() {
                    // As in the C++, these windows get no windowID.
                    Box::new(SequenceWindow::new(
                        self.max_seq_length as usize,
                        self.max_overlap_length as usize,
                    ))
                } else {
                    self.get_waiting_window().unwrap()
                };

                if !self.read_window(&mut win) {
                    return false;
                }

                self.add_ready_window(win);
            }
        } else {
            while !self.waiting_windows.is_empty()
                && (!self.done_reading_file || self.line_place != READ_ALL)
            {
                let mut win = self.get_waiting_window().unwrap();
                if !self.read_window(&mut win) {
                    return false;
                }

                if win.length > 0 {
                    self.add_ready_window(win);
                } else {
                    self.add_waiting_window(win);
                }
            }
        }

        true
    }

    /// Fills one window from the input, starting with the saved overlap.
    pub fn read_window(&mut self, window: &mut SequenceWindow) -> bool {
        let mut reset_symbol_count = false;

        window.prepare_window(
            &self.sequence_name,
            self.sequence_id,
            self.symbols_read_in_seq,
            self.overlap_length as usize,
        );

        if self.overlap_length > 0 {
            window.copy_overlap(&self.overlap_buffer);
        }

        loop {
            if self.line_place >= 0 && self.line_place < READ_ALL {
                self.line_place =
                    window.read_line(&self.line, self.line_place, &mut self.total_seq_length);

                // The window is full partway through the line.
                if self.line_place >= 0 && self.line_place < READ_ALL {
                    break;
                }
            }

            if self.line_place == -1 {
                self.read_special_line();

                if window.length > 0 {
                    // This window ends the previous sequence.
                    reset_symbol_count = true;
                    break;
                } else {
                    // Start the window over at the new sequence.
                    self.symbols_read_in_seq = 0;
                    window.prepare_window(
                        &self.sequence_name,
                        self.sequence_id,
                        self.symbols_read_in_seq,
                        self.overlap_length as usize,
                    );
                }
            }

            if self.line_place >= READ_ALL {
                if getline(&mut self.file, &mut self.line) {
                    self.line_place = 0;
                } else {
                    self.done_reading_file = true;
                    self.file = None;
                    break;
                }
            }
        }

        self.symbols_read_in_seq += window.length as u64;

        if reset_symbol_count {
            self.symbols_read_in_seq = 0;
        }

        if self.shuffle {
            shuffle_sequence_window(window);
        }

        // Quirk: no overlap is saved until at least max_overlap_length symbols
        // of the sequence have been read, and then only from this window.
        if self.symbols_read_in_seq >= self.max_overlap_length {
            self.copy_overlap_buffer_from_window(window, self.max_overlap_length);
        }

        window.read_id = self.read_id;
        self.read_id += 1;

        true
    }

    /// Takes the ready window with the smallest `read_id`.
    pub fn get_ready_window(&mut self) -> Option<Box<SequenceWindow>> {
        if self.ready_windows.is_empty() {
            return None;
        }
        pop_heap(&mut self.ready_windows);
        self.ready_windows.pop()
    }

    /// Takes the most recently returned waiting window.
    pub fn get_waiting_window(&mut self) -> Option<Box<SequenceWindow>> {
        self.waiting_windows.pop()
    }

    pub fn add_ready_window(&mut self, window: Box<SequenceWindow>) -> bool {
        self.ready_windows.push(window);
        push_heap(&mut self.ready_windows);
        true
    }

    pub fn add_waiting_window(&mut self, window: Box<SequenceWindow>) -> bool {
        self.waiting_windows.push(window);
        true
    }

    /// Opens `file_path` and creates `mx_windows` empty windows (none if
    /// `mx_windows == -1`, which selects read-whole-file mode). If the file
    /// cannot be opened, prints `Unable to open file: <path>` to stderr and
    /// exits with status 255 (C++ `exit(-1)`).
    pub fn new(
        file_path: &str,
        mx_windows: i32,
        mx_seq_length: u64,
        mx_overlap_length: u64,
    ) -> Self {
        let mut reader = FASTAReader {
            file: None,
            line: Vec::new(),
            line_place: 0,
            total_seq_length: 0,
            ready_windows: Vec::new(),
            waiting_windows: Vec::new(),
            read_whole_file: false,
            shuffle: false,
            max_windows: mx_windows,
            max_seq_length: mx_seq_length,
            max_overlap_length: mx_overlap_length,
            a_pctg: 0.3,
            t_pctg: 0.3,
            c_pctg: 0.2,
            g_pctg: 0.2,
            overlap_buffer: vec![0; mx_overlap_length as usize + 1],
            overlap_length: 0,
            sequence_name: String::new(),
            sequence_id: 0,
            read_id: 0,
            symbols_read_in_seq: 0,
            done_reading_file: false,
            is_reading: false,
            multithread: true,
        };

        if mx_windows == -1 {
            reader.read_whole_file = true;
        } else {
            let cap = (mx_windows.max(0) as usize) + 1;
            reader.waiting_windows.reserve(cap);
            reader.ready_windows.reserve(cap);
        }

        reader.create_windows(mx_windows, mx_seq_length, mx_overlap_length);

        match File::open(file_path) {
            Ok(f) => reader.file = Some(BufReader::new(f)),
            Err(_) => {
                eprint!("Unable to open file: {}\n", file_path);
                std::process::exit(-1);
            }
        }

        reader.line.reserve(4086);
        reader.line.clear();
        reader.line_place = 0;
        reader
    }

    /// The "random windows" constructor. Nothing in ULTRA calls it.
    ///
    /// The C++ seeds `rand()` from `time(NULL) + getpid()` and names the
    /// sequence after the seed, but never draws a random number (the only
    /// consumer, `CreateRandomWindow`, is commented out), so no generator is
    /// kept here. Its parameters shadow the members, so `maxWindows`,
    /// `maxSeqLength`, `maxOverlapLength`, `overlapBuffer` and (unless
    /// `max_windows == -1`) `readWholeFile` are left uninitialized and no file
    /// is opened; this port uses zeros/false for those members. The windows
    /// themselves are created with `max_seq_length` and an overlap of 0.
    pub fn new_random(
        _random_windows: u64,
        max_windows: i32,
        max_seq_length: u64,
        _max_overlap_length: u64,
    ) -> Self {
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0);
        // (unsigned int)(time(NULL) + getpid())
        let seed = now.wrapping_add(std::process::id() as i64) as u32;

        let mut reader = FASTAReader {
            file: None,
            line: Vec::new(),
            line_place: 0,
            total_seq_length: 0,
            ready_windows: Vec::new(),
            waiting_windows: Vec::new(),
            read_whole_file: false,
            shuffle: false,
            max_windows: 0,
            max_seq_length: 0,
            max_overlap_length: 0,
            a_pctg: 0.3,
            t_pctg: 0.3,
            c_pctg: 0.2,
            g_pctg: 0.2,
            overlap_buffer: vec![0; 1],
            overlap_length: 0,
            sequence_name: seed.to_string(),
            sequence_id: 0,
            read_id: 0,
            symbols_read_in_seq: 0,
            done_reading_file: false,
            is_reading: false,
            multithread: true,
        };

        if max_windows == -1 {
            reader.read_whole_file = true;
        } else {
            let cap = (max_windows.max(0) as usize) + 1;
            reader.waiting_windows.reserve(cap);
            reader.ready_windows.reserve(cap);
        }

        // The shadowed parameter was set to 0 before the windows were made.
        reader.create_windows(max_windows, max_seq_length, 0);
        reader
    }

    /// Creates `count` windows with ids `0..count` and pushes them onto the
    /// waiting pool in that order.
    fn create_windows(&mut self, count: i32, seq_length: u64, overlap_length: u64) {
        for i in 0..count.max(0) {
            let mut new_window = Box::new(SequenceWindow::new(
                seq_length as usize,
                overlap_length as usize,
            ));
            new_window.window_id = i as u64;
            self.add_waiting_window(new_window);
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::symbol::char_for_symbol;
    use std::io::Write;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// A temporary file removed on drop.
    pub(crate) struct TempFasta(pub PathBuf);

    impl TempFasta {
        pub(crate) fn new(tag: &str, contents: &[u8]) -> Self {
            static COUNTER: AtomicUsize = AtomicUsize::new(0);
            let n = COUNTER.fetch_add(1, Ordering::SeqCst);
            let path = std::env::temp_dir().join(format!(
                "ultra_reader_test_{}_{}_{}.fa",
                std::process::id(),
                tag,
                n
            ));
            let mut f = File::create(&path).unwrap();
            f.write_all(contents).unwrap();
            TempFasta(path)
        }

        pub(crate) fn path(&self) -> &str {
            self.0.to_str().unwrap()
        }
    }

    impl Drop for TempFasta {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    /// One line per window in the same format as the C++ harness.
    pub(crate) fn describe(w: &SequenceWindow) -> String {
        let mut s = format!(
            "r{} s{} w{} [{}] {}-{} len{} ov{} ",
            w.read_id, w.seq_id, w.window_id, w.sequence_name, w.start, w.end, w.length, w.overlap
        );
        for &c in &w.seq()[..w.overlap] {
            s.push(char_for_symbol(c) as char);
        }
        s.push('|');
        for &c in w.new_seq() {
            s.push(char_for_symbol(c) as char);
        }
        s
    }

    /// Drives the reader as single-threaded `Ultra::AnalyzeFile` /
    /// `GetSequenceWindow` do.
    fn dump(contents: &[u8], max_windows: i32, max_seq: u64, max_overlap: u64) -> Vec<String> {
        let tmp = TempFasta::new("fasta", contents);
        let mut r = FASTAReader::new(tmp.path(), max_windows, max_seq, max_overlap);
        r.multithread = false;
        r.fill_windows();
        let mut out = Vec::new();
        let get = |r: &mut FASTAReader, seq: Option<Box<SequenceWindow>>| {
            if let Some(s) = seq {
                r.add_waiting_window(s);
            }
            if !r.done_reading_file && !r.is_reading && r.ready_windows.is_empty() {
                r.is_reading = true;
                r.fill_windows();
                r.is_reading = false;
            }
            r.get_ready_window()
        };
        let mut cur = get(&mut r, None);
        while cur.is_some() || !r.done_reading_file {
            if let Some(w) = &cur {
                out.push(describe(w));
            }
            cur = get(&mut r, cur);
        }
        out.push(format!("total {}", r.total_seq_length));
        out
    }

    const EDGE1: &[u8] = b"acgtNNxx\n\n>seq one desc\nACGTACGTAC\nacgtacgtac\n;comment ACGT\nGG-TT*CC12\n\n>empty\n>third\r\nAC\r\nGT\r\n\r\nTTTTAAAACCCCGGGG\n>last\nACGTTGCA";
    const EDGE2: &[u8] = b">a\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>b\nAC\n>c\nACGTACGTACGTACGTACGTACGTACGTACGTACG\n\n\n";

    fn lines(s: &str) -> Vec<String> {
        s.lines().map(|l| l.to_string()).collect()
    }

    // Expected values below are dumps from the C++ FASTAReader driven the same
    // way.

    #[test]
    fn edge1_two_windows_matches_cpp() {
        let expected = "r0 s0 w1 [] 0-5 len5 ov0 |ACGTN
r1 s0 w0 [] 5-8 len3 ov3 GTN|NNN
r2 s1 w0 [seq one desc] 0-5 len5 ov0 |ACGTA
r3 s1 w1 [seq one desc] 5-10 len5 ov3 GTA|CGTAC
r4 s1 w1 [seq one desc] 10-15 len5 ov3 TAC|ACGTA
r5 s1 w0 [seq one desc] 15-20 len5 ov3 GTA|CGTAC
r6 s1 w0 [seq one desc] 20-24 len4 ov3 TAC|GGTT
r7 s1 w1 [seq one desc] 24-28 len4 ov3 GTT|CCNN
r8 s3 w1 [third\r] 0-5 len5 ov0 |ACGTT
r9 s3 w0 [third\r] 5-10 len5 ov3 GTT|TTTAA
r10 s3 w0 [third\r] 10-15 len5 ov3 TAA|AACCC
r11 s3 w1 [third\r] 15-20 len5 ov3 CCC|CGGGG
r12 s4 w1 [last] 0-5 len5 ov0 |ACGTT
r13 s4 w0 [last] 5-8 len3 ov3 GTT|GCA
total 64";
        assert_eq!(dump(EDGE1, 2, 5, 3), lines(expected));
    }

    #[test]
    fn overlap_longer_than_window_matches_cpp() {
        let expected = "r0 s0 w0 [] 0-3 len3 ov0 |ACG
r1 s0 w0 [] 3-6 len3 ov0 |TNN
r2 s0 w0 [] 6-8 len2 ov3 TNN|NN
r3 s1 w0 [seq one desc] 0-3 len3 ov0 |ACG
r4 s1 w0 [seq one desc] 3-6 len3 ov0 |TAC
r5 s1 w0 [seq one desc] 6-9 len3 ov3 TAC|GTA
r6 s1 w0 [seq one desc] 9-12 len3 ov3 GTA|CAC
r7 s1 w0 [seq one desc] 12-15 len3 ov3 CAC|GTA
r8 s1 w0 [seq one desc] 15-18 len3 ov3 GTA|CGT
r9 s1 w0 [seq one desc] 18-21 len3 ov3 CGT|ACG
r10 s1 w0 [seq one desc] 21-23 len2 ov3 ACG|GT
r11 s1 w0 [seq one desc] 23-25 len2 ov2 GT|TC
r12 s1 w0 [seq one desc] 25-28 len3 ov2 TC|CNN
r13 s3 w0 [third\r] 0-3 len3 ov0 |ACG
r14 s3 w0 [third\r] 3-6 len3 ov0 |TTT
r15 s3 w0 [third\r] 6-9 len3 ov3 TTT|TTA
r16 s3 w0 [third\r] 9-12 len3 ov3 TTA|AAA
r17 s3 w0 [third\r] 12-15 len3 ov3 AAA|CCC
r18 s3 w0 [third\r] 15-18 len3 ov3 CCC|CGG
r19 s3 w0 [third\r] 18-20 len2 ov3 CGG|GG
r20 s4 w0 [last] 0-3 len3 ov0 |ACG
r21 s4 w0 [last] 3-6 len3 ov0 |TTG
r22 s4 w0 [last] 6-8 len2 ov3 TTG|CA
total 64";
        assert_eq!(dump(EDGE1, 1, 3, 5), lines(expected));
    }

    #[test]
    fn zero_overlap_three_windows_matches_cpp() {
        let expected = "r0 s0 w2 [] 0-7 len7 ov0 |ACGTNNN
r1 s0 w1 [] 7-8 len1 ov0 |N
r2 s1 w0 [seq one desc] 0-7 len7 ov0 |ACGTACG
r3 s1 w0 [seq one desc] 7-14 len7 ov0 |TACACGT
r4 s1 w1 [seq one desc] 14-21 len7 ov0 |ACGTACG
r5 s1 w2 [seq one desc] 21-26 len5 ov0 |GTTCC
r6 s1 w2 [seq one desc] 26-28 len2 ov0 |NN
r7 s3 w1 [third\r] 0-7 len7 ov0 |ACGTTTT
r8 s3 w0 [third\r] 7-14 len7 ov0 |TAAAACC
r9 s3 w0 [third\r] 14-20 len6 ov0 |CCGGGG
r10 s4 w1 [last] 0-7 len7 ov0 |ACGTTGC
r11 s4 w2 [last] 7-8 len1 ov0 |A
total 64";
        assert_eq!(dump(EDGE1, 3, 7, 0), lines(expected));
    }

    #[test]
    fn read_whole_file_matches_cpp() {
        let expected = "r0 s1 w18446744073709551615 [a] 0-16 len16 ov0 |ACGTACGTACGTACGT
r1 s1 w18446744073709551615 [a] 16-32 len16 ov4 ACGT|ACGTACGTACGTACGT
r2 s1 w18446744073709551615 [a] 32-40 len8 ov4 ACGT|ACGTACGT
r3 s2 w18446744073709551615 [b] 0-2 len2 ov0 |AC
r4 s3 w18446744073709551615 [c] 0-16 len16 ov0 |ACGTACGTACGTACGT
r5 s3 w18446744073709551615 [c] 16-32 len16 ov4 ACGT|ACGTACGTACGTACGT
r6 s3 w18446744073709551615 [c] 32-35 len3 ov4 ACGT|ACG
total 77";
        assert_eq!(dump(EDGE2, -1, 16, 4), lines(expected));
    }

    #[test]
    fn example_1_matches_cpp() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/upstream/ULTRA/examples/example_1.fa"
        );
        let contents = match std::fs::read(path) {
            Ok(c) => c,
            Err(_) => return, // upstream checkout not present
        };
        let d = dump(&contents, 2, 300, 20);
        let summary: Vec<String> = d
            .iter()
            .map(|l| l.split(' ').take(7).collect::<Vec<_>>().join(" "))
            .collect();
        let expected = "r0 s1 w1 [Rand_seq_1] 0-300 len300 ov0
r1 s1 w0 [Rand_seq_1] 300-600 len300 ov20
r2 s1 w0 [Rand_seq_1] 600-900 len300 ov20
r3 s1 w1 [Rand_seq_1] 900-1200 len300 ov20
r4 s1 w1 [Rand_seq_1] 1200-1260 len60 ov20
total 1260";
        assert_eq!(summary, lines(expected));
        assert!(d[1].starts_with(
            "r1 s1 w0 [Rand_seq_1] 300-600 len300 ov20 CTCGCTGTTTTCGAAATTTC|CCTTTAAGCG"
        ));
    }

    #[test]
    fn empty_and_comment_only_files_produce_nothing() {
        assert_eq!(dump(b"", 2, 4, 2), vec!["total 0".to_string()]);
        assert_eq!(
            dump(b"\n\n;only comments\n", 2, 4, 2),
            vec!["total 0".to_string()]
        );
    }

    #[test]
    fn crlf_header_and_trailing_header() {
        let expected = "r0 s1 w1 [x\r] 0-4 len4 ov0 |ACGT
total 4";
        assert_eq!(dump(b">x\r\nACGT\r\n>y\r\n", 2, 4, 2), lines(expected));
    }

    #[test]
    fn ready_heap_pops_in_read_id_order() {
        let tmp = TempFasta::new("heap", b"");
        let mut r = FASTAReader::new(tmp.path(), 0, 4, 1);
        let ids = [5u64, 3, 9, 1, 7, 2, 8, 0, 6, 4];
        for &id in &ids {
            let mut w = Box::new(SequenceWindow::new(4, 1));
            w.read_id = id;
            r.add_ready_window(w);
        }
        // Interleave a push to exercise sift-up after pops.
        let mut got = vec![];
        for _ in 0..3 {
            got.push(r.get_ready_window().unwrap().read_id);
        }
        let mut w = Box::new(SequenceWindow::new(4, 1));
        w.read_id = 10;
        r.add_ready_window(w);
        while let Some(w) = r.get_ready_window() {
            got.push(w.read_id);
        }
        assert_eq!(got, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        assert!(r.get_ready_window().is_none());
    }

    #[test]
    fn waiting_pool_is_lifo_and_ids_follow_creation_order() {
        let tmp = TempFasta::new("lifo", b"");
        let mut r = FASTAReader::new(tmp.path(), 3, 4, 1);
        assert_eq!(r.get_waiting_window().unwrap().window_id, 2);
        assert_eq!(r.get_waiting_window().unwrap().window_id, 1);
        let w0 = r.get_waiting_window().unwrap();
        assert_eq!(w0.window_id, 0);
        assert!(r.get_waiting_window().is_none());
        r.add_waiting_window(w0);
        assert_eq!(r.get_waiting_window().unwrap().window_id, 0);
    }

    #[test]
    fn non_utf8_header_is_converted_lossily() {
        let d = dump(b">n\xffm\nAC\n", 1, 4, 1);
        assert_eq!(d[0], "r0 s1 w0 [n\u{fffd}m] 0-2 len2 ov0 |AC");
    }

    #[test]
    fn random_constructor_builds_windows_without_a_file() {
        let mut r = FASTAReader::new_random(5, 2, 8, 3);
        assert_eq!(r.waiting_windows.len(), 2);
        assert_eq!(r.waiting_windows[0].max_length, 8);
        assert_eq!(r.waiting_windows[0].max_overlap, 0);
        assert!(r.sequence_name.parse::<u32>().is_ok());
        assert!(r.file.is_none());
        r.fill_windows();
        assert!(r.done_reading_file);
        assert!(r.ready_windows.is_empty());
    }
}
