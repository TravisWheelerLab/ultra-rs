//! A window of encoded sequence. Port of `SequenceWindow.hpp` /
//! `SequenceWindow.cpp`.
//!
//! The C++ class keeps four raw pointers into one allocation (`seqMem`,
//! `overlapSeq`, `seq`, `newSeq`). The first three always point at index 0, so
//! here the buffer is a single `Vec` and `new_seq_offset` records where the
//! non-overlap part begins.

use crate::symbol::{symbol_for_char, Symbol, NUM_SYMBOLS};

/// Returned by [`SequenceWindow::read_line`] when the whole line was consumed.
pub const READ_ALL: i64 = 1 << 29;

/// Value the C++ constructor assigns to the id fields (`-1` stored in an
/// `unsigned long`).
pub const UNSET_ID: u64 = u64::MAX;

#[derive(Debug, Clone)]
pub struct SequenceWindow {
    /// A unique number assigned to this window.
    pub window_id: u64,
    /// A unique number representing the whole sequence.
    pub seq_id: u64,
    /// A unique number representing this particular read.
    pub read_id: u64,

    pub sequence_name: String,

    /// Where the window starts in the sequence.
    pub start: u64,
    /// Where the window ends in the sequence.
    pub end: u64,

    /// Backing storage, `max_length + max_overlap + 1` symbols. The overlap
    /// with the previous window occupies `[0, overlap)` and the new symbols
    /// occupy `[new_seq_offset, new_seq_offset + length)`.
    pub seq_mem: Vec<Symbol>,
    /// Index into `seq_mem` of the first non-overlap symbol (C++ `newSeq`).
    pub new_seq_offset: usize,

    /// How many new (non-overlap) symbols are in this window.
    pub length: usize,
    /// Number of symbols shared with the previous window.
    pub overlap: usize,

    /// Running symbol counts. As in the C++ these are never reset by
    /// [`prepare_window`](Self::prepare_window), only by construction.
    pub symbol_counts: [u64; NUM_SYMBOLS],
    pub symbol_freqs: [f64; NUM_SYMBOLS],

    /// The largest number of new symbols the window can store.
    pub max_length: usize,
    /// The largest overlap the window can store.
    pub max_overlap: usize,
}

impl SequenceWindow {
    pub fn new(max_length: usize, max_overlap: usize) -> Self {
        SequenceWindow {
            window_id: UNSET_ID,
            seq_id: UNSET_ID,
            read_id: UNSET_ID,
            sequence_name: String::new(),
            start: 0,
            end: 0,
            seq_mem: vec![0; max_length + max_overlap + 1],
            new_seq_offset: max_overlap,
            length: 0,
            overlap: 0,
            symbol_counts: [0; NUM_SYMBOLS],
            symbol_freqs: [0.0; NUM_SYMBOLS],
            max_length,
            max_overlap,
        }
    }

    /// The whole window, overlap first (C++ `seq` / `overlapSeq`). Valid data
    /// runs for `overlap + length` symbols.
    #[inline]
    pub fn seq(&self) -> &[Symbol] {
        &self.seq_mem
    }

    /// The non-overlap part of the window (C++ `newSeq`), `length` symbols.
    #[inline]
    pub fn new_seq(&self) -> &[Symbol] {
        &self.seq_mem[self.new_seq_offset..self.new_seq_offset + self.length]
    }

    #[inline]
    pub fn new_seq_mut(&mut self) -> &mut [Symbol] {
        let (start, end) = (self.new_seq_offset, self.new_seq_offset + self.length);
        &mut self.seq_mem[start..end]
    }

    pub fn prepare_window(&mut self, seq_name: &str, seq_id: u64, start: u64, overlap: usize) {
        self.overlap = overlap;

        self.sequence_name.clear();
        self.sequence_name.push_str(seq_name);

        self.start = start;
        self.end = start;
        self.seq_id = seq_id;

        self.length = 0;
        self.new_seq_offset = overlap;

        self.symbol_freqs = [0.0; NUM_SYMBOLS];
    }

    /// Reads symbols from `line` starting at byte `place`.
    ///
    /// Returns [`READ_ALL`] if the rest of the line was consumed (or the line
    /// is a `;` comment), `-1` if the line starts a new sequence (`>`), and
    /// otherwise the position in `line` at which reading stopped because the
    /// window is full.
    ///
    /// As in the C++, the number of bytes examined is capped by the space left
    /// in the window, and bytes that are not ASCII alphanumerics are skipped
    /// without being stored but still count toward that cap.
    pub fn read_line(&mut self, line: &[u8], place: i64, total_seq_length: &mut u64) -> i64 {
        let window_remaining = self.max_length - self.length;
        let mut place = place as usize;
        let characters_to_read = window_remaining.min(line.len().saturating_sub(place));

        match line.first() {
            Some(b';') => return READ_ALL,
            Some(b'>') => return -1,
            _ => {}
        }

        let mut characters_read: u64 = 0;
        for _ in 0..characters_to_read {
            let c = line[place];
            place += 1;
            if c.is_ascii_alphanumeric() {
                characters_read += 1;
                let s = symbol_for_char(c);
                self.seq_mem[self.new_seq_offset + self.length] = s;
                self.length += 1;
                self.symbol_counts[s as usize] += 1;
                *total_seq_length += 1;
            }
        }

        self.end += characters_read;

        if place >= line.len() {
            return READ_ALL;
        }

        place as i64
    }

    /// Copies `overlap` symbols from `b` into the front of the window.
    pub fn copy_overlap(&mut self, b: &[Symbol]) {
        let n = self.overlap;
        self.seq_mem[..n].copy_from_slice(&b[..n]);
        self.new_seq_offset = n;
    }

    pub fn calculate_symbol_frequencies(&mut self) {
        let total_symbols: u64 = self.symbol_counts.iter().sum();

        // Check to see if there are any symbols at all.
        if total_symbols == 0 {
            return;
        }

        for i in 0..NUM_SYMBOLS {
            self.symbol_freqs[i] = self.symbol_counts[i] as f64 / total_symbols as f64;
        }
    }
}

/// Small non-cryptographic generator (xoshiro256++) used for shuffling.
///
/// The C++ seeds a fresh `std::mt19937` from `std::random_device` on every
/// shuffle, so shuffled output was never reproducible between runs. Only the
/// statistical behaviour needs to match.
#[derive(Debug, Clone)]
pub struct ShuffleRng {
    s: [u64; 4],
}

impl ShuffleRng {
    /// Seeds the generator from operating-system entropy.
    pub fn from_entropy() -> Self {
        use std::collections::hash_map::RandomState;
        use std::hash::{BuildHasher, Hasher};

        // `RandomState` is keyed from OS randomness; mix in the clock and pid
        // so that repeated calls in one process still diverge.
        let mut hasher = RandomState::new().build_hasher();
        hasher.write_u32(std::process::id());
        if let Ok(d) = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
            hasher.write_u128(d.as_nanos());
        }
        Self::from_seed(hasher.finish())
    }

    /// Seeds the generator deterministically (SplitMix64 expansion).
    pub fn from_seed(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        ShuffleRng {
            s: [next(), next(), next(), next()],
        }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }
}

/// Shuffles the non-overlap symbols of `window` in place (Fisher-Yates), using
/// a freshly seeded generator as the C++ does.
pub fn shuffle_sequence_window(window: &mut SequenceWindow) {
    let mut rng = ShuffleRng::from_entropy();
    shuffle_sequence_window_with(window, &mut rng);
}

/// Shuffles the non-overlap symbols of `window` in place with the given
/// generator.
pub fn shuffle_sequence_window_with(window: &mut SequenceWindow, rng: &mut ShuffleRng) {
    let length = window.length as u64;
    let new_seq = window.new_seq_mut();
    for i in 0..length {
        let j = (rng.next_u64() % (length - i)) + i;
        new_seq.swap(i as usize, j as usize);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbol::{N_A, N_C, N_G, N_T, N_UNK};

    #[test]
    fn new_window_matches_cpp_constructor() {
        let w = SequenceWindow::new(10, 4);
        assert_eq!(w.seq_mem.len(), 15);
        assert_eq!(w.new_seq_offset, 4);
        assert_eq!(w.window_id, UNSET_ID);
        assert_eq!(w.seq_id, UNSET_ID);
        assert_eq!(w.read_id, UNSET_ID);
        assert_eq!((w.length, w.overlap, w.start, w.end), (0, 0, 0, 0));
    }

    #[test]
    fn read_line_consumes_whole_line() {
        let mut w = SequenceWindow::new(10, 0);
        w.prepare_window("s", 0, 0, 0);
        let mut total = 0;
        assert_eq!(w.read_line(b"ACGTn", 0, &mut total), READ_ALL);
        assert_eq!(w.new_seq(), &[N_A, N_C, N_G, N_T, N_UNK]);
        assert_eq!((w.length, w.end, total), (5, 5, 5));
        assert_eq!(w.symbol_counts[N_A as usize], 1);
        assert_eq!(w.symbol_counts[N_UNK as usize], 1);
    }

    #[test]
    fn read_line_stops_when_window_is_full() {
        let mut w = SequenceWindow::new(3, 0);
        w.prepare_window("s", 0, 100, 0);
        let mut total = 0;
        assert_eq!(w.read_line(b"ACGTAC", 0, &mut total), 3);
        assert_eq!((w.length, w.start, w.end), (3, 100, 103));

        w.prepare_window("s", 0, 103, 0);
        assert_eq!(w.read_line(b"ACGTAC", 3, &mut total), READ_ALL);
        assert_eq!(w.new_seq(), &[N_T, N_A, N_C]);
        assert_eq!(total, 6);
    }

    #[test]
    fn read_line_skips_non_alphanumerics_but_counts_them_against_the_cap() {
        let mut w = SequenceWindow::new(4, 0);
        w.prepare_window("s", 0, 0, 0);
        let mut total = 0;
        // Four bytes are examined; the two spaces are skipped.
        assert_eq!(w.read_line(b"A C GT", 0, &mut total), 4);
        assert_eq!(w.new_seq(), &[N_A, N_C]);
        assert_eq!((w.end, total), (2, 2));
    }

    #[test]
    fn read_line_flags_headers_comments_and_empty_lines() {
        let mut w = SequenceWindow::new(4, 0);
        let mut total = 0;
        assert_eq!(w.read_line(b">name", 0, &mut total), -1);
        assert_eq!(w.read_line(b";comment", 0, &mut total), READ_ALL);
        assert_eq!(w.read_line(b"", 0, &mut total), READ_ALL);
        assert_eq!((w.length, total), (0, 0));
    }

    #[test]
    fn prepare_window_keeps_symbol_counts() {
        let mut w = SequenceWindow::new(8, 2);
        w.prepare_window("a", 0, 0, 0);
        let mut total = 0;
        w.read_line(b"AAAA", 0, &mut total);
        w.prepare_window("b", 1, 0, 2);
        assert_eq!(w.symbol_counts[N_A as usize], 4);
        assert_eq!((w.length, w.new_seq_offset, w.overlap), (0, 2, 2));
        assert_eq!(w.sequence_name, "b");
    }

    #[test]
    fn copy_overlap_fills_the_front() {
        let mut w = SequenceWindow::new(4, 2);
        w.prepare_window("a", 0, 10, 2);
        w.copy_overlap(&[N_G, N_T, N_A]);
        let mut total = 0;
        w.read_line(b"CC", 0, &mut total);
        assert_eq!(&w.seq()[..4], &[N_G, N_T, N_C, N_C]);
        assert_eq!(w.new_seq(), &[N_C, N_C]);
    }

    #[test]
    fn symbol_frequencies_sum_to_one() {
        let mut w = SequenceWindow::new(8, 0);
        w.prepare_window("a", 0, 0, 0);
        let mut total = 0;
        w.read_line(b"AACG", 0, &mut total);
        w.calculate_symbol_frequencies();
        assert_eq!(w.symbol_freqs[N_A as usize], 0.5);
        assert_eq!(w.symbol_freqs[N_C as usize], 0.25);
        assert_eq!(w.symbol_freqs.iter().sum::<f64>(), 1.0);
    }

    #[test]
    fn shuffle_preserves_composition_and_overlap() {
        let mut w = SequenceWindow::new(64, 2);
        w.prepare_window("a", 0, 0, 2);
        w.copy_overlap(&[N_G, N_G]);
        let mut total = 0;
        w.read_line(&b"ACGT".repeat(16), 0, &mut total);
        let mut before = w.new_seq().to_vec();

        let mut rng = ShuffleRng::from_seed(7);
        shuffle_sequence_window_with(&mut w, &mut rng);

        assert_eq!(&w.seq()[..2], &[N_G, N_G]);
        assert_ne!(w.new_seq(), &before[..]);
        let mut after = w.new_seq().to_vec();
        before.sort_unstable();
        after.sort_unstable();
        assert_eq!(before, after);
    }
}
