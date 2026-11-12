//! Repeat regions extracted from a DP traceback. Port of `repeat.hpp` /
//! `repeat.cpp`.
//!
//! A [`RepeatRegion`] owns everything it holds: it copies the fields it needs
//! out of the [`SequenceWindow`] and the [`UMatrix`] and keeps no reference to
//! either, so it can outlive the window's reuse.
//!
//! Byte strings (`traceback`, `sequence`, `string_consensus`, the `consensi`
//! entries) are `Vec<u8>`: they hold ASCII and the code indexes them by byte.
//! Use the `*_str` accessors to print them.
//!
//! C++ pointer members that some code tests against null are `Option`s here;
//! `None` is the C++ null pointer.

use crate::sequence_window::SequenceWindow;
use crate::symbol::{char_for_symbol, symbol_for_char, Symbol, NUM_SYMBOLS};
use crate::umatrix::{UMatrix, CT_DELETION, CT_INSERTION, CT_MATCH};

pub const OC_NONE: i32 = 0;
pub const OC_TRUE: i32 = 1;

/// Byte `i` of `s`, or 0 when `i` is at or past the end. Mirrors
/// `std::string::operator[]`, which returns `'\0'` at `size()`; past that
/// the C++ is undefined and this keeps returning 0.
#[inline]
fn byte_at(s: &[u8], i: i64) -> u8 {
    if i < 0 {
        // Undefined in C++; treat as the terminator.
        return 0;
    }
    s.get(i as usize).copied().unwrap_or(0)
}

/// C++ `s.substr(start, len)` for non-negative `start <= s.len()`.
fn substr(s: &[u8], start: usize, len: usize) -> Vec<u8> {
    let end = start.saturating_add(len).min(s.len());
    s[start..end].to_vec()
}

fn ascii_str(s: &[u8]) -> &str {
    std::str::from_utf8(s).expect("repeat strings are ASCII")
}

/// One tandem repeat (C++ `RepeatRegion`).
#[derive(Debug, Clone, PartialEq)]
pub struct RepeatRegion {
    pub sequence_name: String,

    /// `'.'`, `'I'`, `'D'` per position; filled by
    /// [`store_traceback`](Self::store_traceback). Empty until then.
    pub traceback: Vec<u8>,
    /// Upper-case `ACGTN`; filled by [`store_sequence`](Self::store_sequence).
    pub sequence: Vec<u8>,

    pub sequence_id: u64,
    pub read_id: u64,

    pub win_overlap_size: u64,
    /// Left uninitialised by the C++ constructor (set by `GetNextRepeat`,
    /// `SubRepeat` and `joint_repeat_region`); 0 here.
    pub win_total_length: u64,

    pub combined_repeat: bool,

    pub false_start: i32,
    pub window_start: u64,
    pub sequence_start: u64,

    pub repeat_length: u64,

    pub repeat_period: i32,

    pub start_score: f32,
    pub end_score: f32,
    pub region_score: f32,
    /// Per-position score deltas, `repeat_length` entries
    /// ([`store_scores`](Self::store_scores) or a join). `None` = C++ NULL.
    pub scores: Option<Vec<f32>>,

    /// C++ `logoMemory` and `logo` together: a flat buffer of
    /// `repeat_period * (NUM_SYMBOLS + 1)` ints; C++ `logo[i][j]` is
    /// `logo[i * NUM_SYMBOLS + j]` (see [`logo_at`](Self::logo_at)). The tail
    /// past `repeat_period * NUM_SYMBOLS` is unused, as in the C++. `None` =
    /// NULL (joined regions, or no logo built).
    pub logo: Option<Vec<i32>>,
    /// `repeat_period` symbols. `None` = NULL.
    pub consensus: Option<Vec<Symbol>>,

    /// [`get_logo_numbers`](Self::get_logo_numbers) allocates
    /// `repeat_length + 1` entries (the last is never written; 0 here, garbage
    /// in C++); a join allocates `repeat_length`. `None` = not allocated (the
    /// C++ constructor leaves this pointer uninitialised).
    pub logo_numbers: Option<Vec<i32>>,

    /// Only allocated by the dead `SplitRepeats` path.
    pub look_back: Option<Vec<i32>>,
    pub look_forward: Option<Vec<i32>>,
    /// Never allocated anywhere upstream (C++ `char *`, signed 8-bit).
    pub forward_counts: Option<Vec<i8>>,
    pub back_counts: Option<Vec<i8>>,

    /// `None` = C++ nullptr, distinct from an empty vector.
    pub splits: Option<Vec<i32>>,
    /// `None` = C++ nullptr, distinct from an empty vector.
    pub consensi: Option<Vec<Vec<u8>>>,

    pub mismatches: i32,
    pub insertions: i32,
    pub deletions: i32,

    pub overlap_correction: i32,

    pub string_consensus: Vec<u8>,
}

impl Default for RepeatRegion {
    fn default() -> Self {
        Self::new()
    }
}

impl RepeatRegion {
    /// C++ `RepeatRegion()`.
    pub fn new() -> Self {
        RepeatRegion {
            sequence_name: String::new(),
            traceback: Vec::new(),
            sequence: Vec::new(),
            sequence_id: 0,
            read_id: 0,
            win_overlap_size: 0,
            win_total_length: 0,
            combined_repeat: false,
            false_start: 0,
            window_start: 0,
            sequence_start: 0,
            repeat_length: 0,
            repeat_period: 0,
            start_score: 0.0,
            end_score: 0.0,
            region_score: 0.0,
            scores: None,
            logo: None,
            consensus: None,
            logo_numbers: None,
            look_back: None,
            look_forward: None,
            forward_counts: None,
            back_counts: None,
            splits: None,
            consensi: None,
            mismatches: 0,
            insertions: 0,
            deletions: 0,
            overlap_correction: OC_NONE,
            string_consensus: Vec::new(),
        }
    }

    pub fn traceback_str(&self) -> &str {
        ascii_str(&self.traceback)
    }

    pub fn sequence_str(&self) -> &str {
        ascii_str(&self.sequence)
    }

    pub fn string_consensus_str(&self) -> &str {
        ascii_str(&self.string_consensus)
    }

    /// Consensi entry `i` as `&str`; panics if `consensi` is `None` or too
    /// short (C++ `consensi->at(i)` would crash or throw).
    pub fn consensi_str(&self, i: usize) -> &str {
        ascii_str(&self.consensi.as_ref().expect("consensi is null")[i])
    }

    /// C++ `logo[i][j]`.
    #[inline]
    pub fn logo_at(&self, i: usize, j: usize) -> i32 {
        self.logo.as_ref().expect("logo is null")[i * NUM_SYMBOLS + j]
    }

    fn alloc_logo(&mut self) -> Vec<i32> {
        let period = self.repeat_period as i64 as usize;
        // Rows are zeroed; the unused tail is zero too (garbage in C++, never read).
        vec![0i32; period * (NUM_SYMBOLS + 1)]
    }

    /// C++ `CreateLogo`: counts symbols per unit position by walking the
    /// matrix traceback forward from `window_start`, and accumulates
    /// `mismatches`, `insertions` and `deletions`.
    pub fn create_logo(&mut self, window: &SequenceWindow, matrix: &UMatrix) {
        let mut logo = self.alloc_logo();
        let period = self.repeat_period;
        let period_u = period as i64 as u64;
        let seq = window.seq();

        // `backset`, `insertionReset` and `deletionReset` are dropped: the
        // resets start at -1 and only decrement, so the branches that read
        // them never run and `backset` is never used.
        let mut endpoint: u64 = self.repeat_length.wrapping_sub(period_u);

        let mut fplace: u64 = period_u.wrapping_add(self.window_start);
        let mut bplace: u64 = self.window_start;

        let win_len = window.length as u64 + window.overlap as u64;
        if endpoint > win_len {
            endpoint = win_len;
        }

        let mut p: i32 = 0;
        let mut modp: i32 = 0;

        let mut i: u64 = 0;
        while i < endpoint {
            let frow = matrix.traceback[fplace as usize];
            let brow = matrix.traceback[bplace as usize];

            let fdesc = matrix.cell_descriptions[frow as usize];
            let mut bdesc = matrix.cell_descriptions[brow as usize];

            // unsigned long < int: the int converts to unsigned long.
            if i < period_u {
                bdesc = matrix.cell_descriptions[0];
            }

            if fdesc.cell_type == CT_INSERTION {
                self.insertions += fdesc.indel_number;
                i = i.wrapping_add(fdesc.indel_number as i64 as u64);
                fplace = fplace.wrapping_add(fdesc.indel_number as i64 as u64);
            } else if fdesc.cell_type == CT_DELETION {
                self.deletions += fdesc.indel_number;
                bplace = bplace.wrapping_add(fdesc.indel_number as i64 as u64);
                p += fdesc.indel_number;
                modp += fdesc.indel_number + 1;
            } else {
                modp += 1;
            }

            modp %= period;
            let _ = modp;

            if bdesc.cell_type == CT_INSERTION {
                bplace = bplace.wrapping_add(bdesc.indel_number as i64 as u64);
            } else if bdesc.cell_type == CT_DELETION {
                // The bases skipped here sit at the unit positions the
                // previous copy deleted. They are counted in those columns
                // (weight 1, as there is no base to compare them with) and
                // the phase moves past them. Upstream 1.2.2 does neither,
                // which shifts every later logo column by the size of the
                // deletion. Known limit: if the model calls a deletion in a
                // copy that is in fact full length, moving the phase leaves
                // it off instead (see TODO.md).
                for j in 0..bdesc.indel_number.max(0) {
                    let fp = fplace.wrapping_add(j as u64);
                    if fp < self.window_start.wrapping_add(self.repeat_length) {
                        let col = ((p + j) % period) as usize;
                        logo[col * NUM_SYMBOLS + seq[fp as usize] as usize] += 1;
                    }
                }
                fplace = fplace.wrapping_add(bdesc.indel_number as i64 as u64);
                i = i.wrapping_add(bdesc.indel_number as i64 as u64);
                // Known limit: if the model calls a deletion in a copy that
                // is in fact full length, this leaves the phase off instead
                // (see TODO.md).
                p += bdesc.indel_number;
            }

            if fplace >= self.window_start.wrapping_add(self.repeat_length) {
                break;
            }

            let s = seq[fplace as usize];
            let sb = seq[bplace as usize];

            p %= period;

            let idx = p as usize * NUM_SYMBOLS + s as usize;
            if s != sb {
                self.mismatches += 1;
                logo[idx] += 1;
            } else {
                logo[idx] += 4;
            }

            p += 1;
            fplace = fplace.wrapping_add(1);
            bplace = bplace.wrapping_add(1);
            i = i.wrapping_add(1);
        }

        self.logo = Some(logo);
    }

    /// C++ `CreateLogoWithoutMatrix`: builds the logo from `traceback` and
    /// `sequence`. Only reached from the dead `SubRepeat`.
    pub fn create_logo_without_matrix(&mut self) {
        if self.sequence.is_empty() || self.traceback.is_empty() {
            return;
        }

        let mut logo = self.alloc_logo();
        let period = self.repeat_period;

        let mut p: i32 = 0;
        let mut i: i32 = 0;
        // int < unsigned long: the int converts to unsigned long.
        while (i as i64 as u64) < self.repeat_length {
            let t = byte_at(&self.traceback, i as i64);
            if t == b'.' {
                let c = byte_at(&self.sequence, i as i64);
                logo[(p % period) as usize * NUM_SYMBOLS + symbol_for_char(c) as usize] += 1;
                p += 1;
            } else if t == b'I' {
                // Do nothing
            } else if t == b'D' {
                let mut dels = 1;
                let mut j = 1;
                while byte_at(&self.traceback, i as i64 + j) == b'D' {
                    dels += 1;
                    j += 1;
                }

                for _ in 0..dels {
                    let c = byte_at(&self.sequence, i as i64);
                    logo[((p + dels) % period) as usize * NUM_SYMBOLS
                        + symbol_for_char(c) as usize] += 1;
                    p += 1;
                }
            }
            i += 1;
        }

        self.logo = Some(logo);
    }

    /// C++ `CreateConsensusFromLogo`: the arg-max symbol of each logo row
    /// (lowest symbol on ties), then sets `string_consensus`.
    pub fn create_consensus_from_logo(&mut self) {
        let period = self.repeat_period.max(0) as usize;
        let logo = self.logo.as_ref().expect("logo is null");
        let mut consensus = vec![0 as Symbol; period];
        for (i, c) in consensus.iter_mut().enumerate() {
            let row = &logo[i * NUM_SYMBOLS..i * NUM_SYMBOLS + NUM_SYMBOLS];
            let mut s: usize = 0;
            for j in 1..NUM_SYMBOLS {
                if row[j] > row[s] {
                    s = j;
                }
            }
            *c = s as Symbol;
        }
        self.consensus = Some(consensus);

        self.string_consensus = self.get_consensus();
    }

    /// C++ `GetConsensus`: `consensus` as characters. Panics if `consensus`
    /// is `None` (the C++ dereferences NULL, e.g. on a joined region).
    pub fn get_consensus(&self) -> Vec<u8> {
        let consensus = self.consensus.as_ref().expect("consensus is null");
        (0..self.repeat_period.max(0) as usize)
            .map(|i| char_for_symbol(consensus[i]))
            .collect()
    }

    /// C++ `StoreSequence`: copies `repeat_length` symbols of the window,
    /// starting at `window_start`, as characters.
    pub fn store_sequence(&mut self, window: &SequenceWindow) {
        let seq = window.seq();
        let mut sequence = Vec::with_capacity(self.repeat_length as usize + 1);
        let mut i: i32 = 0;
        // int < unsigned long, and `i + windowStart` in unsigned long.
        while (i as i64 as u64) < self.repeat_length {
            let idx = (i as i64 as u64).wrapping_add(self.window_start);
            sequence.push(char_for_symbol(seq[idx as usize]));
            i += 1;
        }
        // The C++ then writes '\0' at index repeatLength, which is the
        // terminator already; nothing changes.
        self.sequence = sequence;
    }

    /// C++ `StoreTraceback`: appends (the C++ never clears it first)
    /// `repeat_period` dots, then one character per traceback column of the
    /// repeat, with each indel cell expanded to `indel_number` characters.
    pub fn store_traceback(&mut self, matrix: &UMatrix) {
        for _ in 0..self.repeat_period.max(0) {
            self.traceback.push(b'.');
        }

        let end: u64 = self.window_start.wrapping_add(self.repeat_length);
        let mut i: u64 = self
            .window_start
            .wrapping_add(self.repeat_period as i64 as u64);
        while i < end {
            let c = matrix.cell_descriptions[matrix.traceback[i as usize] as usize];

            if c.cell_type == CT_INSERTION {
                i = i.wrapping_add((c.indel_number - 1) as i64 as u64);
                for _ in 0..c.indel_number {
                    self.traceback.push(b'I');
                }
            } else if c.cell_type == CT_DELETION {
                i = i.wrapping_add((c.indel_number - 1) as i64 as u64);
                self.traceback.push(b'D');
                for _ in 1..c.indel_number {
                    self.traceback.push(b'D');
                }
            } else {
                self.traceback.push(b'.');
            }
            i = i.wrapping_add(1);
        }
    }

    /// C++ `StoreScores`: the score change into each column of the repeat.
    /// Reads matrix columns `window_start ..= window_start + repeat_length`.
    pub fn store_scores(&mut self, matrix: &UMatrix) {
        let mut scores = vec![0f32; self.repeat_length as usize];

        let mut i: u64 = 1;
        while i <= self.repeat_length {
            let p = i.wrapping_add(self.window_start) as usize;
            let mut cell = matrix.traceback[p];
            if matrix.cell_descriptions[cell as usize].cell_type != CT_MATCH {
                cell = matrix.cell_descriptions[cell as usize].parent_index;
            }
            let s = matrix.score_at(p, cell as usize);

            let mut prev_cell = matrix.traceback[p - 1];
            if matrix.cell_descriptions[prev_cell as usize].cell_type != CT_MATCH {
                // Upstream 1.2.2 takes the parent of `cell` here.
                prev_cell = matrix.cell_descriptions[prev_cell as usize].parent_index;
            }
            let ps = matrix.score_at(p - 1, prev_cell as usize);

            scores[(i - 1) as usize] = s - ps;
            i += 1;
        }

        self.scores = Some(scores);
    }

    /// C++ `LookBackDistance` (dead code). Requires `look_back` allocated.
    pub fn look_back_distance(&mut self) {
        let mut look_back = self.look_back.take().expect("lookBack is null");
        let period = self.repeat_period;
        let tb = &self.traceback;

        for i in 0..period.max(0) {
            look_back[i as usize] = 0;
        }

        let mut indel_num: i32 = 0;
        let mut in_period: i32 = 0;
        let mut del_period: i32 = 0;
        let mut indel_stop: i32 = 0;
        // `delSwitch` is declared and never used upstream.

        let mut i: i32 = period;
        while (i as i64 as u64) < self.repeat_length {
            let iu = i as usize;
            look_back[iu] = period;

            if indel_num == 0 {
                if byte_at(tb, i as i64) == b'D' {
                    indel_stop = i + period;

                    del_period = period - 1;
                    in_period = period + period - 1;

                    indel_num = 1;

                    let mut j = 1;
                    while byte_at(tb, i as i64 + j) == b'D' {
                        del_period -= 1;
                        in_period -= 1;
                        indel_num += 1;
                        j += 1;
                    }

                    look_back[iu] = del_period;
                } else if byte_at(tb, i as i64) == b'I' {
                    indel_stop = i + period;

                    indel_num = 1;
                    in_period = period + 1;

                    look_back[iu] = 0;
                }
            } else if byte_at(tb, i as i64) == b'I' {
                indel_stop += 1;

                in_period += 1;
                indel_num += 1;

                look_back[iu] = 0;
            } else if del_period > 0 {
                look_back[iu] = del_period;

                if i >= indel_stop {
                    look_back[iu] = period;
                    indel_num = 0;
                    in_period = 0;
                    del_period = 0;
                } else if i >= indel_stop - indel_num {
                    look_back[iu] = in_period;
                }
            } else if byte_at(tb, i as i64) == b'.' {
                if i > indel_stop {
                    indel_num = 0;
                    in_period = 0;
                    del_period = 0;
                } else if in_period > 0 {
                    look_back[iu] = in_period;
                }
            }
            i += 1;
        }

        self.look_back = Some(look_back);
    }

    /// C++ `LookForwardDistance` (dead code). Requires `look_forward`
    /// allocated.
    pub fn look_forward_distance(&mut self) {
        let mut look_forward = self.look_forward.take().expect("lookForward is null");
        let period = self.repeat_period;
        // int = unsigned long - int, narrowed.
        let start_tail = self.repeat_length.wrapping_sub(period as i64 as u64) as i32;
        let tb = &self.traceback;

        let mut i = start_tail;
        while (i as i64 as u64) < self.repeat_length {
            look_forward[i as usize] = 0;
            i += 1;
        }

        let mut indel_num: i32 = 0;
        let mut in_period: i32 = 0;
        let mut del_period: i32 = 0;
        let mut indel_stop: i32 = 0;

        let mut i: i32 = start_tail - 1;
        while i >= 0 {
            let iu = i as usize;
            look_forward[iu] = period;

            if indel_num == 0 {
                if byte_at(tb, i as i64) == b'D' {
                    del_period = period;
                    in_period = period + period;

                    indel_num = 0;
                    while byte_at(tb, i as i64) == b'D' {
                        del_period -= 1;
                        in_period -= 1;
                        indel_num += 1;
                        look_forward[i as usize] = period;
                        i -= 1;
                    }

                    indel_stop = i - period + 1;
                    // Upstream writes lookForward[-1] here when the run
                    // reaches index 0 (undefined); skipped in that case.
                    if i >= 0 {
                        look_forward[i as usize] = del_period;
                    }
                } else if byte_at(tb, i as i64) == b'I' {
                    indel_stop = i - period;
                    indel_num = 1;
                    in_period = period + 1;
                    look_forward[iu] = 0;
                }
            } else if del_period > 0 {
                look_forward[iu] = del_period;

                if i < indel_stop {
                    look_forward[iu] = period;
                    indel_num = 0;
                    in_period = 0;
                    del_period = 0;
                } else if i < indel_stop + indel_num {
                    look_forward[iu] = in_period;
                }
            } else if byte_at(tb, i as i64) == b'.' {
                if i < indel_stop {
                    indel_num = 0;
                    in_period = 0;
                    del_period = 0;
                } else if in_period > 0 {
                    look_forward[iu] = in_period;
                }
            } else if byte_at(tb, i as i64) == b'I' {
                indel_stop -= 1;
                in_period += 1;
                indel_num += 1;
                look_forward[iu] = 0;
            }
            i -= 1;
        }

        self.look_forward = Some(look_forward);
    }

    /// C++ `CountFromDistance` (dead code). `counts` are signed 8-bit and
    /// wrap like the C++ `char`.
    pub fn count_from_distance(&self, distance: &[i32], counts: &mut [i8], dir: i32, depth: i32) {
        let mut i: i32 = 0;
        while (i as i64 as u64) < self.repeat_length {
            let iu = i as usize;
            let mut count: i8 = 0;

            if distance[iu] == 0 {
                counts[iu] = -1;
                i += 1;
                continue;
            }

            let mut p = i;
            for _ in 0..depth {
                p += distance[p as usize] * dir;

                if p < 0 || (p as i64 as u64) >= self.repeat_length {
                    break;
                }

                if byte_at(&self.sequence, i as i64) == byte_at(&self.sequence, p as i64) {
                    count = count.wrapping_add(1);
                }

                if distance[p as usize] == 0 {
                    break;
                }
            }

            counts[iu] = count;
            i += 1;
        }
    }

    /// C++ `CountForwardBackwardConsensus` (dead code): allocates
    /// `look_back` / `look_forward` and fills them. `depth` is unused, as
    /// upstream.
    pub fn count_forward_backward_consensus(&mut self, _depth: i32) {
        let n = self.repeat_length as usize;
        // malloc'd and uninitialised in C++; every entry is written below.
        self.look_back = Some(vec![0; n]);
        self.look_forward = Some(vec![0; n]);

        self.look_back_distance();
        self.look_forward_distance();
    }

    /// C++ `SplitRepeats` (dead code). Upstream never allocates
    /// `forwardCounts` / `backCounts`, so the C++ dereferences NULL here; this
    /// panics in the same place.
    pub fn split_repeats(&mut self, depth: i32, cutoff: i32) -> Option<Vec<Box<RepeatRegion>>> {
        self.count_forward_backward_consensus(depth);

        let inv_cutoff = depth - cutoff;
        let mut rep_start: i32 = 0;
        let mut sub_reps: Option<Vec<Box<RepeatRegion>>> = None;

        let back_counts = self.back_counts.clone().expect("backCounts is null");
        let forward_counts = self.forward_counts.clone().expect("forwardCounts is null");
        let look_back = self.look_back.clone().expect("lookBack is null");
        let look_forward = self.look_forward.clone().expect("lookForward is null");

        let mut i: i32 = 0;
        while (i as i64 as u64) < self.repeat_length {
            let iu = i as usize;
            if back_counts[iu] < 0 || forward_counts[iu] < 0 || look_forward[iu] == 0 {
                i += 1;
                continue;
            }

            if back_counts[iu] as i32 > inv_cutoff && (forward_counts[iu] as i32) < cutoff {
                let f = i + look_forward[iu];
                let fu = f as usize;

                if back_counts[fu] < 0 || forward_counts[fu] < 0 || look_back[iu] == 0 {
                    i += 1;
                    continue;
                }

                if (back_counts[fu] as i32) < cutoff && forward_counts[fu] as i32 > inv_cutoff {
                    let sub = self.sub_repeat(rep_start, f - rep_start);
                    sub_reps.get_or_insert_with(Vec::new).push(sub);
                    rep_start = f;
                    i = f;
                }
            }
            i += 1;
        }

        if rep_start > 0 {
            let len = self.repeat_length.wrapping_sub(rep_start as i64 as u64) as i32;
            let sub = self.sub_repeat(rep_start, len);
            sub_reps.as_mut().expect("subReps is null").push(sub);
        }

        sub_reps
    }

    /// C++ `SubRepeat` (only reached from the dead `SplitRepeats`).
    pub fn sub_repeat(&self, start: i32, length: i32) -> Box<RepeatRegion> {
        let mut region = Box::new(RepeatRegion::new());

        region.sequence_name = self.sequence_name.clone();

        if !self.traceback.is_empty() {
            region.traceback = substr(
                &self.traceback,
                start as usize,
                length as i64 as u64 as usize,
            );
        }
        if !self.sequence.is_empty() {
            region.sequence = substr(
                &self.sequence,
                start as usize,
                length as i64 as u64 as usize,
            );
        }

        region.sequence_id = self.sequence_id;
        region.read_id = self.read_id;

        region.win_overlap_size = self.win_overlap_size;
        region.win_total_length = self.win_total_length;

        region.combined_repeat = self.combined_repeat;

        region.false_start = self.false_start;

        region.window_start = self.window_start.wrapping_add(start as i64 as u64);
        region.sequence_start = self.sequence_start.wrapping_add(start as i64 as u64);

        region.repeat_length = length as i64 as u64;

        region.repeat_period = self.repeat_period;

        let score_per_symbol = self.region_score / (self.repeat_length as f32);

        region.start_score = self.start_score + (score_per_symbol * (start as f32));
        region.end_score = region.start_score + (score_per_symbol * (length as f32));
        region.region_score = region.end_score - region.start_score;

        if !region.sequence.is_empty() && !region.traceback.is_empty() {
            region.create_logo_without_matrix();
            region.create_consensus_from_logo();
        }

        region
    }

    /// C++ `GetLogoNumbers`: the unit position of each repeat position (`-1`
    /// for insertions). Requires `traceback`.
    pub fn get_logo_numbers(&mut self) {
        let length: u64 = self.repeat_length;
        let period = self.repeat_period;
        // malloc(length + 1) in C++; entry `length` is never written.
        let mut logo_numbers = vec![0i32; length as usize + 1];
        let tb = &self.traceback;

        let mut num: i32 = 0;
        let mut i: u64 = 0;
        while i < length {
            let t = byte_at(tb, i as i64);
            if t == b'.' {
                logo_numbers[i as usize] = num % period;
                num += 1;
            }

            if t == b'I' {
                logo_numbers[i as usize] = -1;
            }

            let mut del_count = 0;
            let mut j: i64 = 0;
            while byte_at(tb, i as i64 + j) == b'D' {
                num += 1;
                del_count += 1;
                j += 1;
            }

            for _ in 0..del_count {
                logo_numbers[i as usize] = num % period;

                num += 1;
                i += 1;
            }
            if del_count > 0 {
                i -= 1;
            }
            i += 1;
        }

        self.logo_numbers = Some(logo_numbers);
    }

    /// C++ `PermutationForString`: `str` rotated left by `offset`.
    pub fn permutation_for_string(s: &[u8], offset: i32) -> Vec<u8> {
        let mut new_string = Vec::with_capacity(s.len());
        for i in 0..s.len() as i32 {
            // (int) % size_t: the int converts to size_t.
            let p = ((offset + i) as i64 as u64 % s.len() as u64) as usize;
            new_string.push(s[p]);
        }
        new_string
    }

    /// C++ `CompareStrPerm`: 0 if the two rotations are equal, 1 if rotation
    /// `perm2` sorts first, -1 if `perm1` does.
    pub fn compare_str_perm(s: &[u8], perm1: i32, perm2: i32) -> i32 {
        let size = s.len() as i32;
        for i in 0..size {
            let p1 = ((i + perm1) % size) as usize;
            let p2 = ((i + perm2) % size) as usize;

            // C++ `char` is signed; the strings are ASCII so this is the same.
            let (c1, c2) = (s[p1] as i8, s[p2] as i8);
            if c1 < c2 {
                return -1;
            }
            if c2 < c1 {
                return 1;
            }
        }
        0
    }

    /// The smallest rotation index of `s` as found by the C++ loop (which
    /// stops at the first rotation equal to the current best).
    fn best_rotation(s: &[u8]) -> i32 {
        let mut best_perm = 0;
        for i in 1..s.len() as i32 {
            let c = Self::compare_str_perm(s, best_perm, i);
            if c == 0 {
                break;
            }
            if c == 1 {
                best_perm = i;
            }
        }
        best_perm
    }

    /// C++ `SortConsensus`: rotates `consensi[index]` to its least rotation.
    pub fn sort_consensus(&mut self, index: usize) {
        let consensi = self.consensi.as_mut().expect("consensi is null");
        let best_perm = Self::best_rotation(&consensi[index]);
        if best_perm == 0 {
            return;
        }
        consensi[index] = Self::permutation_for_string(&consensi[index], best_perm);
    }

    /// C++ `SortConsensi`: rotates `string_consensus` and every `consensi`
    /// entry to its least rotation.
    pub fn sort_consensi(&mut self) {
        if self.string_consensus.len() > 1 {
            let best_perm = Self::best_rotation(&self.string_consensus);
            self.string_consensus = Self::permutation_for_string(&self.string_consensus, best_perm);
        }

        let n = match &self.consensi {
            None => return,
            Some(c) if c.is_empty() => return,
            Some(c) => c.len(),
        };

        for i in 0..n {
            self.sort_consensus(i);
        }
    }
}

/// C++ `GetNextRepeat`. Scans `matrix.traceback` from column `*pos` (0 is
/// treated as 1) for the next non-background column and returns the repeat
/// starting there, with its logo and consensus built, or `None`. Sets `*pos`
/// to the column after the repeat (or to `overlap + length` when none is
/// found). Requires `matrix.calculate_traceback` to have run on `window`.
pub fn get_next_repeat(
    window: &SequenceWindow,
    matrix: &UMatrix,
    pos: &mut i32,
) -> Option<Box<RepeatRegion>> {
    if *pos == 0 {
        *pos = 1;
    }

    let mut found_repeat = false;
    let seq_length: u64 = window.length as u64 + window.overlap as u64;

    let mut region: Option<Box<RepeatRegion>> = None;

    // int < unsigned long comparisons below convert the int to unsigned long.
    let mut i: i32 = *pos;
    while (i as i64 as u64) < seq_length {
        if matrix.traceback[i as usize] > 0 {
            found_repeat = true;
            break;
        }
        i += 1;
    }

    if found_repeat {
        let mut row = matrix.traceback[i as usize];
        let desc = matrix.cell_descriptions[row as usize];

        let mut r = Box::new(RepeatRegion::new());
        r.win_overlap_size = window.overlap as u64;
        r.win_total_length = window.length as u64 + window.overlap as u64;
        r.sequence_name = window.sequence_name.clone();
        r.sequence_id = window.seq_id;
        r.read_id = window.read_id;
        r.repeat_period = desc.order;

        // int arithmetic, then int -> unsigned long.
        r.window_start = (i - r.repeat_period - 1) as i64 as u64;
        // window->start + (i - window->overlap - repeatPeriod) - 1, all in
        // unsigned long with wraparound.
        r.sequence_start = window
            .start
            .wrapping_add(
                (i as i64 as u64)
                    .wrapping_sub(window.overlap as u64)
                    .wrapping_sub(r.repeat_period as i64 as u64),
            )
            .wrapping_sub(1);

        r.start_score = matrix.score_at(i as usize, row as usize);

        let mut length: i32 = r.repeat_period;

        // `<=`: one extra column is calculated.
        while (i as i64 as u64) <= seq_length {
            if matrix.traceback[i as usize] != 0 {
                length += 1;
            } else {
                i += 1; // always make i = end + 1
                break;
            }
            i += 1;
        }

        r.repeat_length = length as i64 as u64;

        row = matrix.traceback[(i - 1) as usize];

        r.end_score = matrix.score_at((i - 1) as usize, row as usize);
        r.region_score = r.end_score - r.start_score;
        region = Some(r);
    }

    *pos = i;

    if let Some(r) = region.as_mut() {
        r.read_id = window.read_id;
        r.create_logo(window, matrix);
        r.create_consensus_from_logo();
    }

    region
}

/// C++ `repeats_overlap`: same sequence (and, if `require_same_period`, same
/// period) and `r1` ends at or after `r2` starts. Ultra calls it with
/// `require_same_period = true` (the C++ default argument).
pub fn repeats_overlap(r1: &RepeatRegion, r2: &RepeatRegion, require_same_period: bool) -> bool {
    if r1.sequence_id == r2.sequence_id {
        if require_same_period {
            if r1.repeat_period == r2.repeat_period {
                return r1.sequence_start.wrapping_add(r1.repeat_length) >= r2.sequence_start;
            } else {
                return false;
            }
        }
        return r1.sequence_start.wrapping_add(r1.repeat_length) >= r2.sequence_start;
    }
    false
}

/// C++ `joint_repeat_region`: a new region covering `r1` and `r2`, or `None`
/// if they do not overlap (checked without requiring equal periods). Assumes
/// `r2` starts at or after `r1`. Does not modify or consume its inputs.
///
/// As upstream, the result has no logo, consensus, traceback-derived counts
/// (mismatches/insertions/deletions are 0), start/end scores or window
/// start; `r2`'s split positions are appended without shifting them.
pub fn joint_repeat_region(r1: &RepeatRegion, r2: &RepeatRegion) -> Option<Box<RepeatRegion>> {
    if !repeats_overlap(r1, r2, false) {
        return None;
    }

    let mut joint = Box::new(RepeatRegion::new());

    joint.string_consensus = r1.string_consensus.clone();
    joint.sequence_start = r1.sequence_start;

    // The joint region is the first `s1_seq_len` positions of r1 followed by
    // r2 without its first `r2_skip` positions. If r2 ends inside r1 the
    // joint region is all of r1 and none of r2. Upstream 1.2.2 always takes
    // r1 up to the start of r2 and then all of r2, which cuts r1 short in
    // that case.
    let r2_inside = r2.sequence_start.wrapping_add(r2.repeat_length)
        <= r1.sequence_start.wrapping_add(r1.repeat_length);
    // unsigned long narrowed to int.
    let mut s1_seq_len: i32 = r2.sequence_start.wrapping_sub(r1.sequence_start) as i32;
    let mut r2_skip: i32 = 0;
    if r2_inside {
        s1_seq_len = r1.repeat_length as i32;
        r2_skip = r2.repeat_length as i32;
    }
    joint.repeat_length =
        (s1_seq_len.wrapping_add((r2.repeat_length as i32).wrapping_sub(r2_skip))) as i64 as u64;
    let overlap: i32 = r1
        .sequence_start
        .wrapping_add(r1.repeat_length)
        .wrapping_sub(r2.sequence_start) as i32;

    joint.sequence_name = r1.sequence_name.clone();
    joint.sequence_id = r1.sequence_id;
    joint.read_id = r1.read_id;
    joint.win_overlap_size = r1.win_overlap_size;
    joint.win_total_length = r1.win_total_length;

    joint.combined_repeat = true;

    joint.repeat_period = r1.repeat_period;

    joint.overlap_correction = r1.overlap_correction + 1;

    // substr(0, n) with the int converted to size_t (clamped to the length).
    let s1_count = s1_seq_len as i64 as u64 as usize;

    if !r1.sequence.is_empty() && !r2.sequence.is_empty() {
        let mut s = substr(&r1.sequence, 0, s1_count);
        s.extend_from_slice(&r2.sequence[(r2_skip as usize).min(r2.sequence.len())..]);
        joint.sequence = s;
    }

    if !r1.traceback.is_empty() && !r2.traceback.is_empty() {
        let mut t = substr(&r1.traceback, 0, s1_count);
        t.extend_from_slice(&r2.traceback[(r2_skip as usize).min(r2.traceback.len())..]);
        joint.traceback = t;
    }

    if let (Some(sc1), Some(sc2)) = (&r1.scores, &r2.scores) {
        let mut scores = vec![0f32; joint.repeat_length as usize];
        let mut i: usize = 0;
        while (i as i32) < s1_seq_len {
            scores[i] = sc1[i];
            i += 1;
        }
        for j in r2_skip as usize..r2.repeat_length as usize {
            scores[i] = sc2[j];
            i += 1;
        }
        joint.scores = Some(scores);
    }

    if let (Some(ln1), Some(ln2)) = (&r1.logo_numbers, &r2.logo_numbers) {
        // malloc(repeatLength): no extra slot, unlike GetLogoNumbers.
        let mut ln = vec![0i32; joint.repeat_length as usize];
        let mut i: usize = 0;
        while (i as i32) < s1_seq_len {
            ln[i] = ln1[i];
            i += 1;
        }
        for j in r2_skip as usize..r2.repeat_length as usize {
            ln[i] = ln2[j];
            i += 1;
        }
        joint.logo_numbers = Some(ln);
    }

    // Weighted average of the scores. `1.0 - pct_seq1` is computed in double.
    let pct_seq1: f32 = s1_seq_len as f32 / joint.repeat_length as f32;
    let pct_seq2: f32 = (1.0f64 - pct_seq1 as f64) as f32;

    joint.region_score = (r1.region_score * pct_seq1) + (r2.region_score * pct_seq2);

    // r2 adds nothing, so the splits of r1 stand.
    if r2_inside {
        if let Some(s) = &r1.splits {
            joint.splits = Some(s.clone());
            joint.consensi = r1.consensi.clone();
        }
        return Some(joint);
    }

    if r1.repeat_period != r2.repeat_period {
        return Some(joint);
    }
    let (Some(r1_splits), Some(r2_splits)) = (&r1.splits, &r2.splits) else {
        return Some(joint);
    };

    let mut splits = r1_splits.clone();
    let mut consensi = r1.consensi.clone().expect("consensi is null");
    let r2_consensi = r2.consensi.as_ref().expect("consensi is null");

    let mut split_start: i32 = -1;
    for (i, &s) in r2_splits.iter().enumerate() {
        if s < overlap {
            continue;
        }
        split_start = i as i32;
        break;
    }

    if split_start != -1 {
        for i in split_start as usize..r2_splits.len() {
            // Splits are positions within their own repeat, so r2's move to
            // where r2 sits in the joint region (-1 marks a removed split).
            // Upstream 1.2.2 appends them unshifted.
            let s = r2_splits[i];
            splits.push(if s < 0 { s } else { s + s1_seq_len });
            consensi.push(r2_consensi[i + 1].clone());
        }
    }

    joint.splits = Some(splits);
    joint.consensi = Some(consensi);
    Some(joint)
}

/// C++ `dead_func`: joins a region with itself. Unused upstream.
pub fn dead_func(r1: &RepeatRegion) -> Option<Box<RepeatRegion>> {
    joint_repeat_region(r1, r1)
}

#[cfg(test)]
mod tests {
    //! Expected values come from a C++ harness linked against the upstream
    //! sources (`repeat.cpp`, `RepeatSplitter.cpp`, `umatrix.cpp`,
    //! `umodel.cpp`, ...) at commit 7947d0d (1.2.2 plus the fixes in
    //! TravisWheelerLab/ULTRA pull request 95), since the port makes the same
    //! fixes. The harness splits a sequence into windows of
    //! `ws` new symbols with `ov` overlap, runs the `AnalyzeSequenceWindow`
    //! repeat loop (store + split with the default split settings), then the
    //! join / `ValidateSplits` / `SortConsensi` steps of
    //! `OutputRepeats(flush)`, and dumps every field. The dump format is
    //! reproduced by `dump` below; floats are printed as bit patterns.
    use super::*;
    use crate::repeat_splitter::{validate_splits, SplitWindow};
    use crate::std_sort::std_sort;
    use crate::umodel::UModel;
    use std::fmt::Write as _;

    /// `upstream/ULTRA/examples/example_1.fa`.
    const T1: &[&str] = &[
        "GCTTAGTTCGAGATGAAATGTCCGGGGCCAAAGACAACCAGCATCTCGCGTCTTGCCTAACCCCCCTACATGTTGTTATAGACAATCAGTGGAAACCCAG",
        "TGCCAGACGATGGAATGACCTTAAGTCAGGCAGGAAATTAAAGGGAACGTATATTCAACGCAATGAAGCTGGAGGATTGGCGTGGGAATCGTGCTTCTGT",
        "CTAAGCAAGAATGGGTATGAGGTGGCAACCGTCCCCCTAGCGTACAGGGTGCACTTTGTAACGATTTGGGAGTCCAGAGACTCGCTGTTTTCGAAATTTC",
        "CCTTTAAGCGCGGGTATTGAACCAGGCTTATGCCCAAGATCGTAGCAAGCAGACTCAAACAAAATATATTTTGCCCGCGTTACAGACGAAACTAGTTGGA",
        "GGTTATGGAGCATACTAACACGTGGACGGCCACTGGTGAGTTGCTACATATATATATATATATATATATATATATATATATATATATATATATATATATA",
        "TATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCCCCTGCGGCAA",
        "CGTTGAAGCTCCTAAATTACTCTGGCTGGATCCAAGCCGTAACACCCGTCTCACTTCATAACCGTTTGCAATTCACGGCTTGATCTAGACTGGATTGCCA",
        "TTCTCTCAAAGTATTATGCAGGCCGGCGTACGCGTCCCATATAAACCTGTCATAACTTACCTGAGACTACTTGGAAATGTGGCTAGATCTTTGCCCACGC",
        "ACCTAATCGGTCCACGTTTGGTTTTTAGGACCCGATGATCTTCAAACACTGCAAGATTCCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCTCGCAGGC",
        "AGGCAGGCACGGCAGGCAGGTCAACCTGCTTTACAAAGCGCTGGATCCTACTCCAGCGGGATCTTTTATCTAAAGACGATGAGAGGAGTATTCGTCAGGC",
        "CACATAGCTTTCATGTTCTGATCGGAACGATCGTTGGCGCCCGACCCTCGGATTCCGTAGTGAGTTCTTTGTCCGAGCCATTGTATGCGAGATCGGTAGA",
        "CTGATAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGGGGATGCAGTATATCCCTGGATGCAATAGACGG",
        "ACAGGTTGGAATCCTAAGTGAAGTCGCGCGTCCGAACCCAGCTCTATTTTAGAGGTCACG",
    ];

    /// Two pieces of a synthetic sequence with indels and repeats that change
    /// unit part way (the test inputs' `mixed.fa`, seq2 [5200, 6000) and
    /// [9600, 10200)).
    const T2: &[&str] = &[
        "AGCTCTTAATCCCAATGGCCCAGCGATATGGGCATAAAAGTTAATACAGTTACACCGTTGTCAGCGTTTCCTATCCTAATTACGAGGTAATAGAGTTTGT",
        "TGGTCTGCGCTTCCAGGACACTAAGTCATATCGGACCAGACGAGCAAATTAAGAGTGTGAACAAGACTATGAGAAACTATTTTCGAGGAACAGGGCCCAG",
        "GGCCCAGGGCCCAGTGCCCAGGGCCCAGGGCCCAGGGGCCAGGGCCCAAGGGCCCAGGGCCCATGCCCCGGGCCCAGGGCCCAGGGCCCAGGACCCAGGG",
        "CCCAGGGCCCAGGGCCCTGTTTCATGTTTCATGTTTCATGTTTCATGTTTCATGTTTCATGTTTCATTGAGGTGCCCTATGGATACCGAGGCAAAGTAAC",
        "CTCAAGAGGCTTTGACTATGCAAATCATGGCCCTAACGCAACCATGATCCTATACGGTTTGACAATCTGGCCTTTAAGTCTTCTCGGCTATACGTAGATA",
        "CGGGTGATGCTCAAGGCCAATCTCCCGGATCCGTCAAGATCGCTACACGCTCCGGGCTCGGGCTCGCACTTTGTTTGCCACGTCAGGTCCACATGTCCGC",
        "GGAACCACCATAGAGGCCTAAGAAACAGCAGCCACAAAAGCCATTTCATCGGGTCGAGTGTCCCACCTGATTGTGTCTAGCGAGACACGCCGCTACGCCA",
        "ACTCAGGGATAGCGATTTAACGAGCGCCGCTTTCTGTCGCAGCATTCTAAACGTCTCATAAGATTTGGCGGGCGGCGAACGTGAGTCCGTAGGGCATACA",
        "AGTTTACCCGGGGTTGGCAATCGGGCCCGGCTCGGTATCTCATGGGGAGGATCTAAAATGAAAGGGTGTTTTAAGGTTAGCACTGACGTAGGCAATCTAA",
        "CTACACGTACCGCGAGACCGTTCACGCTGCGATAGCTCACACAATCCGACCTAACGGGGGCATTAGCTTAGCTTAGCTTGCTTAGCTTAGCCTAGCTTAG",
        "CTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCTTAGCCTTAGCTTAGCTTAGCTTAGCTTAGCTTGGCCGGACTGGGCGGA",
        "CTGGGCGGATGGGCGGCTGGGCGGACCTGGGCGGACTGGGCGGACTGGGCGGACTGGGCAGACTGGGCGGACTGGGCGGACTGGGCGGACTGGGCGGACT",
        "GGGCTGACTGGGCGGACTGGGCGGACTGGGCGGACTGGGTCAAAGGTTTACGGTTCAAGTCTCGGTCATAGACCAAAGCCATAGAAGCTGTCCTACGCTT",
        "ATGATAGAATCGGCCTTGCACAGTTCTTGCCTATACGTCACACCGTCACTCCACTATTTGATTCCACCCCGCTCGATATACAGGCGCAGACCATGTCTGA",
    ];

    fn fnv1a(data: &[u8]) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for &b in data {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
        h
    }

    fn dump(out: &mut String, tag: &str, r: &RepeatRegion) {
        writeln!(
            out,
            "{} name={} sid={} rid={} wov={} wtot={} comb={} fs={} ws={} ss={} len={} per={}",
            tag,
            r.sequence_name,
            r.sequence_id,
            r.read_id,
            r.win_overlap_size,
            r.win_total_length,
            r.combined_repeat as i32,
            r.false_start,
            r.window_start,
            r.sequence_start,
            r.repeat_length,
            r.repeat_period
        )
        .unwrap();
        writeln!(
            out,
            " sc {:08x} {:08x} {:08x}",
            r.start_score.to_bits(),
            r.end_score.to_bits(),
            r.region_score.to_bits()
        )
        .unwrap();
        writeln!(out, " tb {} {}", r.traceback.len(), r.traceback_str()).unwrap();
        writeln!(out, " seq {} {}", r.sequence.len(), r.sequence_str()).unwrap();
        writeln!(out, " scons {}", r.string_consensus_str()).unwrap();
        if r.logo.is_some() {
            out.push_str(" logo");
            for i in 0..r.repeat_period as usize {
                out.push_str(" |");
                for j in 0..NUM_SYMBOLS {
                    write!(out, " {}", r.logo_at(i, j)).unwrap();
                }
            }
            out.push('\n');
        } else {
            out.push_str(" logo null\n");
        }
        if let Some(c) = &r.consensus {
            out.push_str(" cons");
            for v in c.iter().take(r.repeat_period as usize) {
                write!(out, " {}", v).unwrap();
            }
            writeln!(out, "\n getcons {}", ascii_str(&r.get_consensus())).unwrap();
        } else {
            out.push_str(" cons null\n");
        }
        match &r.logo_numbers {
            Some(ln) => {
                out.push_str(" ln");
                for v in ln.iter().take(r.repeat_length as usize) {
                    write!(out, " {}", v).unwrap();
                }
                out.push('\n');
            }
            None => out.push_str(" ln null\n"),
        }
        match &r.scores {
            Some(s) => {
                out.push_str(" scores");
                for v in s.iter().take(r.repeat_length as usize) {
                    write!(out, " {:08x}", v.to_bits()).unwrap();
                }
                out.push('\n');
            }
            None => out.push_str(" scores null\n"),
        }
        writeln!(
            out,
            " mid {} {} {} oc {}",
            r.mismatches, r.insertions, r.deletions, r.overlap_correction
        )
        .unwrap();
        writeln!(
            out,
            " lb {} lf {} fc {} bc {}",
            r.look_back.is_some() as i32,
            r.look_forward.is_some() as i32,
            r.forward_counts.is_some() as i32,
            r.back_counts.is_some() as i32
        )
        .unwrap();
        match &r.splits {
            Some(s) => {
                write!(out, " splits {}:", s.len()).unwrap();
                for v in s {
                    write!(out, " {}", v).unwrap();
                }
                out.push('\n');
            }
            None => out.push_str(" splits null\n"),
        }
        match &r.consensi {
            Some(c) => {
                write!(out, " consensi {}:", c.len()).unwrap();
                for v in c {
                    write!(out, " {}", ascii_str(v)).unwrap();
                }
                out.push('\n');
            }
            None => out.push_str(" consensi null\n"),
        }
    }

    /// Runs the pipeline described in the module comment; returns the dump
    /// and the regions as they would be handed to the writers.
    fn run(
        seq: &[&str],
        ws: usize,
        ov: usize,
        per: i32,
        ins: i32,
        del: i32,
        split: bool,
    ) -> (String, Vec<Box<RepeatRegion>>) {
        let s: Vec<u8> = seq.concat().into_bytes();
        let max_split: i32 = if split { 10 } else { 0 };
        let split_threshold: f32 = 3.5;
        let split_depth: u64 = 5;
        let min_split_window: u64 = 15;

        let leng = (ws + ov + 2) as i32;
        let mut model = UModel::new(per + 1, ins, del, leng + 1);
        model.period_decay = 0.85;
        model.set_match_probabilities(0.7);
        model.set_atcg_probabilities(0.25, 0.25, 0.25, 0.25);
        model.tp_zero_to_match = 0.01;
        model.tp_match_to_zero = 0.05;
        model.tp_match_to_insertion = 0.02;
        model.tp_match_to_deletion = 0.02;
        model.tp_consecutive_insertion = 0.02;
        model.tp_consecutive_deletion = 0.02;
        model.calculate_scores();

        let mut splitter = SplitWindow::new();
        if split {
            splitter.allocate_split_window(4, per + 1);
        }

        let mut out = String::new();
        let mut all: Vec<Box<RepeatRegion>> = Vec::new();
        let mut window = SequenceWindow::new(ws, ov);
        window.sequence_name = "t".to_string();
        window.seq_id = 1;

        let mut k = 0usize;
        while k * ws < s.len() {
            let start = k * ws;
            let o = if k > 0 { ov } else { 0 };
            let len = ws.min(s.len() - start);
            window.start = start as u64;
            window.overlap = o;
            window.length = len;
            window.read_id = k as u64;
            for i in 0..o + len {
                window.seq_mem[i] = symbol_for_char(s[start - o + i]);
            }

            let sleng = window.length as i32 + window.overlap as i32;
            model.matrix.restart_matrix();
            for i in 0..sleng {
                model.calculate_current_column(&window, i, None);
            }
            model.calculate_current_column_without_emission(None);
            let pc = model.matrix.previous_column_index;
            model.matrix.calculate_traceback(pc);

            let mut i: i32 = 0;
            let mut r = get_next_repeat(&window, &model.matrix, &mut i);
            while let Some(mut rr) = r {
                rr.store_sequence(&window);
                rr.store_traceback(&model.matrix);
                rr.store_scores(&model.matrix);
                rr.get_logo_numbers();
                let split_window = std::cmp::max(
                    (rr.repeat_period as i64 as u64).wrapping_mul(split_depth),
                    min_split_window,
                ) as i32;
                if rr.repeat_period <= max_split {
                    if rr.repeat_length > (2 * split_window) as i64 as u64 {
                        let mut splits =
                            splitter.splits_for_region(&rr, split_window, split_threshold);
                        let mut consensi = splitter.consensi_for_split(&rr, &splits, 0.65);
                        validate_splits(&mut consensi, &mut splits, 0.85);
                        rr.splits = Some(splits);
                        rr.consensi = Some(consensi);
                    } else {
                        rr.splits = Some(Vec::new());
                        rr.consensi = Some(vec![rr.string_consensus.clone()]);
                    }
                }
                writeln!(out, "pos={}", i).unwrap();
                dump(&mut out, "R", &rr);
                all.push(rr);
                r = get_next_repeat(&window, &model.matrix, &mut i);
            }
            writeln!(out, "endpos={}", i).unwrap();
            k += 1;
        }

        std_sort(&mut all, |l, r| {
            if l.sequence_id != r.sequence_id {
                return l.sequence_id > r.sequence_id;
            }
            l.sequence_start > r.sequence_start
        });
        let mut written = Vec::new();
        while let Some(mut r) = all.pop() {
            while let Some(next) = all.last() {
                if !repeats_overlap(&r, next, true) {
                    break;
                }
                let mut j = joint_repeat_region(&r, next).unwrap();
                if let Some(splits) = j.splits.as_mut() {
                    validate_splits(j.consensi.as_mut().unwrap(), splits, 0.85);
                }
                dump(&mut out, "J", &j);
                all.pop();
                r = j;
            }
            r.sort_consensi();
            dump(&mut out, "O", &r);
            written.push(r);
        }
        (out, written)
    }

    #[test]
    fn example_1_single_window_matches_cpp() {
        let (out, written) = run(T1, 2000, 0, 100, 10, 10, true);
        assert_eq!(
            (out.len(), fnv1a(out.as_bytes())),
            (14978, 0x5038ce45a802c8e1)
        );

        let r = &written[0];
        assert_eq!(
            (
                r.window_start,
                r.sequence_start,
                r.repeat_length,
                r.repeat_period
            ),
            (447, 447, 145, 2)
        );
        assert_eq!(
            (
                r.start_score.to_bits(),
                r.end_score.to_bits(),
                r.region_score.to_bits()
            ),
            (0xc13b2a09, 0x4323f652, 0x432fa8f3)
        );
        assert_eq!((r.mismatches, r.insertions, r.deletions), (3, 0, 0));
        assert_eq!(r.splits, Some(vec![81]));
        assert_eq!(r.consensi, Some(vec![b"AT".to_vec(), b"CG".to_vec()]));
        assert_eq!(r.string_consensus, b"AT");
    }

    #[test]
    fn windows_with_joins_match_cpp() {
        let (out, written) = run(T2, 500, 150, 10, 3, 3, true);
        assert_eq!(
            (out.len(), fnv1a(out.as_bytes())),
            (18668, 0xde7570a479d0ee6a)
        );
        assert!(written.iter().any(|r| r.combined_repeat));

        let (out, written) = run(T2, 300, 120, 6, 6, 6, true);
        assert_eq!(
            (out.len(), fnv1a(out.as_bytes())),
            (17251, 0x7328a59618405f90)
        );
        // The two copies start together and the first is the longer (170
        // against 107), so the joined region is the first one whole. Upstream
        // 1.2.2 cuts it to 107 and keeps the split at 123.
        let j = written.iter().find(|r| r.combined_repeat).unwrap();
        assert_eq!(
            (j.sequence_start, j.repeat_length, j.repeat_period),
            (193, 170, 6)
        );
        assert_eq!(j.region_score.to_bits(), 0x42195a43);
        assert_eq!((j.start_score, j.end_score, j.mismatches), (0.0, 0.0, 0));
        assert_eq!(j.overlap_correction, 1);
        assert_eq!(j.splits, Some(vec![44, 123]));
        assert!(j.logo.is_none() && j.consensus.is_none());
        assert_eq!(j.scores.as_ref().unwrap().len(), 170);
        assert_eq!(j.logo_numbers.as_ref().unwrap().len(), 170);
        assert_eq!(j.sequence.len(), 170);
    }

    /// Twenty-four copies of `ACGGTCATTG` between random flanks, with one
    /// base deleted from copies 6, 11, 16 and 21.
    const FOUR_DELETIONS: &[&str] = &[
        "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG",
        "ACGGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGACGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGA",
        "CGGTCTTGACGGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGACGGCATTGACGGTCATTGACGGTCATTGACGGTCATTGACGGTCATTGACG",
        "GTCATGACGGTCATTGACGGTCATTGACGGTCATTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGA",
        "CAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTG",
    ];

    #[test]
    fn logo_phase_follows_deletions_in_the_previous_copy() {
        // Upstream 1.2.2 does not advance the logo phase when the previous
        // copy has a deletion, so each deletion shifts the later copies by
        // one column and this array comes out with the unit `GGGGGTTTTT`.
        for split in [true, false] {
            let (_, written) = run(FOUR_DELETIONS, 1000, 100, 20, 10, 10, split);
            let r = written
                .iter()
                .find(|r| r.repeat_period == 10)
                .expect("the period-10 array is found");
            assert_eq!((r.sequence_start, r.repeat_length), (97, 239));
            assert_eq!(r.deletions, 4);
            assert_eq!(r.string_consensus_str(), "ACGGTCATTG");
        }
    }

    /// A bare region for the join tests: `len` positions of `base` starting
    /// at `start`, period 3.
    fn plain_region(start: u64, len: usize, base: u8, score: f32, splits: Vec<i32>) -> RepeatRegion {
        let mut r = RepeatRegion::new();
        r.sequence_id = 1;
        r.sequence_start = start;
        r.repeat_length = len as u64;
        r.repeat_period = 3;
        r.region_score = score;
        r.sequence = vec![base; len];
        r.traceback = vec![b'.'; len];
        r.consensi = Some((0..=splits.len()).map(|i| vec![b'A' + i as u8; 3]).collect());
        r.splits = Some(splits);
        r
    }

    #[test]
    fn join_moves_the_second_repeats_splits_into_place() {
        // r1 covers 100..150 and r2 covers 130..190, so they overlap by 20.
        let r1 = plain_region(100, 50, b'a', 10.0, vec![20]);
        let r2 = plain_region(130, 60, b'b', 20.0, vec![10, 40, -1, 55]);
        let j = joint_repeat_region(&r1, &r2).unwrap();
        assert_eq!((j.sequence_start, j.repeat_length), (100, 90));
        assert_eq!(j.sequence, [vec![b'a'; 30], vec![b'b'; 60]].concat());
        // r2's split at 10 lies in the overlap and is dropped. Its splits at
        // 40 and 55 are positions 70 and 85 of the joint region, and the
        // removed split stays -1. Upstream 1.2.2 keeps 40 and 55.
        assert_eq!(j.splits, Some(vec![20, 70, -1, 85]));
        assert_eq!(j.consensi.as_ref().unwrap().len(), 5);
    }

    #[test]
    fn join_keeps_the_first_repeat_whole_when_the_second_ends_inside_it() {
        // Upstream 1.2.2 returns a region of length 30 here (same start) and
        // 50 below (later start), cutting r1 short.
        let r1 = plain_region(100, 80, b'a', 10.0, vec![25, 60]);
        for (start, len) in [(100, 30), (120, 30), (150, 30)] {
            let r2 = plain_region(start, len, b'b', 20.0, vec![5]);
            let j = joint_repeat_region(&r1, &r2).unwrap();
            assert_eq!((j.sequence_start, j.repeat_length), (100, 80));
            assert_eq!(j.sequence, r1.sequence);
            assert_eq!(j.traceback, r1.traceback);
            assert_eq!(j.region_score, 10.0);
            assert_eq!(j.splits, r1.splits);
            assert_eq!(j.consensi, r1.consensi);
            assert!(j.combined_repeat);
        }
    }

    #[test]
    fn no_split_matches_cpp() {
        let (out, written) = run(T2, 500, 150, 10, 3, 3, false);
        assert_eq!(
            (out.len(), fnv1a(out.as_bytes())),
            (18568, 0x3652a6ae83f8fc78)
        );
        assert!(written
            .iter()
            .all(|r| r.splits.is_none() && r.consensi.is_none()));
    }

    #[test]
    fn rotations() {
        assert_eq!(RepeatRegion::compare_str_perm(b"CAGG", 0, 1), 1);
        assert_eq!(RepeatRegion::compare_str_perm(b"ATAT", 0, 2), 0);
        assert_eq!(RepeatRegion::permutation_for_string(b"CAGG", 1), b"AGGC");

        let mut r = RepeatRegion::new();
        r.string_consensus = b"TTAG".to_vec();
        r.consensi = Some(vec![b"GCA".to_vec(), Vec::new(), b"ATAT".to_vec()]);
        r.sort_consensi();
        assert_eq!(r.string_consensus, b"AGTT");
        assert_eq!(
            r.consensi,
            Some(vec![b"AGC".to_vec(), Vec::new(), b"ATAT".to_vec()])
        );

        // A one-character consensus is left alone, and null consensi are fine.
        let mut r = RepeatRegion::new();
        r.string_consensus = b"G".to_vec();
        r.sort_consensi();
        assert_eq!(r.string_consensus, b"G");
    }

    #[test]
    fn overlap_rules() {
        let mut a = RepeatRegion::new();
        a.sequence_id = 1;
        a.sequence_start = 100;
        a.repeat_length = 50;
        a.repeat_period = 3;
        let mut b = a.clone();
        b.sequence_start = 150;
        assert!(repeats_overlap(&a, &b, true));
        b.repeat_period = 4;
        assert!(!repeats_overlap(&a, &b, true));
        assert!(repeats_overlap(&a, &b, false));
        b.sequence_start = 151;
        assert!(joint_repeat_region(&a, &b).is_none());
        b.sequence_id = 2;
        b.sequence_start = 120;
        assert!(!repeats_overlap(&a, &b, false));

        // Different periods: no splits joined even when both have them.
        let mut b = a.clone();
        b.sequence_start = 120;
        b.repeat_period = 4;
        a.splits = Some(vec![]);
        a.consensi = Some(vec![b"ACG".to_vec()]);
        b.splits = Some(vec![]);
        b.consensi = Some(vec![b"ACGT".to_vec()]);
        let j = joint_repeat_region(&a, &b).unwrap();
        assert_eq!((j.repeat_length, j.repeat_period), (70, 3));
        assert!(j.splits.is_none() && j.consensi.is_none());
    }
}
