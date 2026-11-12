//! Format-dispatching front end to the input readers. Port of
//! `FileReader.hpp` / `FileReader.cpp`. Only FASTA input exists; the JSON
//! constructor prints an error and exits, as in the C++.
//!
//! # Locking contract
//!
//! In the C++, `FASTAReader` guards its ready pool with `readyLock`
//! (`GetReadyWindow`, `AddReadyWindow`) and its waiting pool with
//! `waitingLock` (`GetWaitingWindow`, `AddWaitingWindow`), both only when
//! `multithread` is set. `Ultra::GetSequenceWindow` calls `AddWaitingWindow`
//! outside its `outerLock`, then takes `outerLock` around `DoneReadingFile`,
//! `IsReading`, `ReadyWindowsSize`, `FillWindows`, `SetIsReading` and
//! `GetReadyWindow`. The only overlap the C++ allows is one thread returning a
//! window to the waiting pool while another is inside `FillWindows`.
//!
//! This port has no internal locks: every method takes `&mut self` (or
//! `&self`). A multithreaded caller wraps the whole [`FileReader`] in one
//! `Mutex` and holds it for each call, and for the whole
//! check-fill-get sequence of `GetSequenceWindow`. A thread returning a window
//! then waits for an in-progress `fill_windows` instead of slipping its window
//! into the pool mid-fill. That cannot deadlock (no other lock is held while
//! waiting) and does not change what any window contains: windows are filled
//! strictly in file order and `read_id`s are assigned in that order, so only
//! which physical window (`window_id`) carries a given read changes, which was
//! already timing-dependent in the C++. Nothing reads per-window state that
//! survives reuse (`symbol_counts` / `symbol_freqs` are unused downstream).
//!
//! The C++ `AnalyzeFileWithThread` loop also reads `DoneReadingFile()`
//! without any lock; in Rust that read goes through the mutex too.

use crate::fasta_reader::FASTAReader;
use crate::sequence_window::SequenceWindow;

/// C++ `file_type`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileType {
    Unknown,
    Fasta,
    Json,
    Bed,
}

pub struct FileReader {
    pub max_windows: i32,
    pub max_seq_length: u64,
    pub max_overlap_length: u64,

    /// Copied into the FASTA reader at construction. As in the C++, setting it
    /// later does not propagate (and neither flag does anything here).
    pub multithread: bool,
    pub format: FileType,

    pub fasta_reader: FASTAReader,
}

impl FileReader {
    /// Opens a FASTA file. `max_windows == -1` reads the whole file in the
    /// first `fill_windows` call. Exits with status 255 if the file cannot be
    /// opened (see [`FASTAReader::new`]).
    pub fn new(
        fasta_file_path: &str,
        m_windows: i32,
        m_seq_length: u64,
        m_overlap_length: u64,
        mthread: bool,
    ) -> Self {
        let mut fasta_reader =
            FASTAReader::new(fasta_file_path, m_windows, m_seq_length, m_overlap_length);
        fasta_reader.multithread = mthread;
        FileReader {
            max_windows: m_windows,
            max_seq_length: m_seq_length,
            max_overlap_length: m_overlap_length,
            multithread: mthread,
            format: FileType::Fasta,
            fasta_reader,
        }
    }

    /// The JSON constructor: prints `JSON input no longer allowed.` to stderr
    /// and exits with status 255. Nothing in ULTRA calls it.
    pub fn new_json(
        _json_file_path: &str,
        _m_seq_length: u64,
        _m_overlap_length: u64,
        _mthread: bool,
    ) -> Self {
        eprint!("JSON input no longer allowed.\n");
        std::process::exit(-1);
    }

    pub fn get_ready_window(&mut self) -> Option<Box<SequenceWindow>> {
        if self.format == FileType::Fasta {
            self.fasta_reader.get_ready_window()
        } else {
            eprint!("SequenceWindow Error, invalid format type\n");
            std::process::exit(-1);
        }
    }

    pub fn add_ready_window(&mut self, window: Box<SequenceWindow>) -> bool {
        if self.format == FileType::Fasta {
            return self.fasta_reader.add_ready_window(window);
        }
        false
    }

    pub fn get_waiting_window(&mut self) -> Option<Box<SequenceWindow>> {
        if self.format == FileType::Fasta {
            return self.fasta_reader.get_waiting_window();
        }
        None
    }

    pub fn add_waiting_window(&mut self, window: Box<SequenceWindow>) -> bool {
        if self.format == FileType::Fasta {
            return self.fasta_reader.add_waiting_window(window);
        }
        false
    }

    pub fn is_reading(&self) -> bool {
        if self.format == FileType::Json {
            return false;
        }
        self.fasta_reader.is_reading
    }

    pub fn done_reading_file(&self) -> bool {
        if self.format == FileType::Json {
            return true;
        }
        self.fasta_reader.done_reading_file
    }

    pub fn set_is_reading(&mut self, value: bool) {
        if self.format == FileType::Fasta {
            self.fasta_reader.is_reading = value;
        }
    }

    pub fn ready_windows_size(&self) -> u64 {
        if self.format == FileType::Fasta {
            return self.fasta_reader.ready_windows.len() as u64;
        }
        0
    }

    pub fn fill_windows(&mut self) -> bool {
        if self.format == FileType::Fasta {
            self.fasta_reader.fill_windows()
        } else {
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta_reader::tests::{describe, TempFasta};

    /// `Ultra::GetSequenceWindow`, single-threaded.
    fn get_sequence_window(
        reader: &mut FileReader,
        seq: Option<Box<SequenceWindow>>,
    ) -> Option<Box<SequenceWindow>> {
        let mut should_read = false;
        if let Some(s) = seq {
            reader.add_waiting_window(s);
        }
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

    fn dump(contents: &[u8], max_windows: i32, max_seq: u64, max_overlap: u64) -> Vec<String> {
        let tmp = TempFasta::new("file", contents);
        let mut reader = FileReader::new(tmp.path(), max_windows, max_seq, max_overlap, false);
        reader.fasta_reader.shuffle = false;
        reader.fill_windows();
        let mut out = Vec::new();
        let mut cur = get_sequence_window(&mut reader, None);
        while cur.is_some() || !reader.done_reading_file() {
            if let Some(w) = &cur {
                out.push(describe(w));
            }
            cur = get_sequence_window(&mut reader, cur);
        }
        out.push(format!("total {}", reader.fasta_reader.total_seq_length));
        out
    }

    // Expected values are dumps from the C++ FileReader driven the same way.

    #[test]
    fn single_window_recycled_matches_cpp() {
        let input = b">a\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>b\nAC\n>c\nACGTACGTACGTACGTACGTACGTACGTACGTACG\n\n\n";
        let expected = [
            "r0 s1 w0 [a] 0-16 len16 ov0 |ACGTACGTACGTACGT",
            "r1 s1 w0 [a] 16-32 len16 ov4 ACGT|ACGTACGTACGTACGT",
            "r2 s1 w0 [a] 32-40 len8 ov4 ACGT|ACGTACGT",
            "r3 s2 w0 [b] 0-2 len2 ov0 |AC",
            "r4 s3 w0 [c] 0-16 len16 ov0 |ACGTACGTACGTACGT",
            "r5 s3 w0 [c] 16-32 len16 ov4 ACGT|ACGTACGTACGTACGT",
            "r6 s3 w0 [c] 32-35 len3 ov4 ACGT|ACG",
            "total 77",
        ];
        assert_eq!(dump(input, 1, 16, 4), expected);
    }

    #[test]
    fn state_queries_follow_the_fasta_reader() {
        let tmp = TempFasta::new("state", b">s\nACGT\n");
        let mut reader = FileReader::new(tmp.path(), 2, 3, 1, true);
        assert_eq!(reader.format, FileType::Fasta);
        assert!(reader.fasta_reader.multithread);
        assert!(!reader.done_reading_file());
        assert!(!reader.is_reading());
        reader.set_is_reading(true);
        assert!(reader.is_reading());
        reader.set_is_reading(false);
        assert!(reader.fill_windows());
        assert!(reader.done_reading_file());
        assert_eq!(reader.ready_windows_size(), 2);
        let w = reader.get_ready_window().unwrap();
        assert_eq!(describe(&w), "r0 s1 w1 [s] 0-3 len3 ov0 |ACG");
        reader.add_waiting_window(w);
        let w = reader.get_ready_window().unwrap();
        assert_eq!(describe(&w), "r1 s1 w0 [s] 3-4 len1 ov1 G|T");
        assert!(reader.get_ready_window().is_none());
        assert_eq!(reader.get_waiting_window().unwrap().window_id, 1);
    }
}
