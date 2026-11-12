//! Splitting repeats where the unit composition changes. Port of
//! `RepeatSplitter.hpp` / `RepeatSplitter.cpp`.
//!
//! [`SplitWindow`] slides two adjacent windows along a repeat, keeps per-unit
//! position symbol counts for each (`p` = left, `q` = right) and scores the
//! symmetric KL divergence between them. Positions where the score peaks
//! above a threshold become split points.
//!
//! The C++ row tables `cp[a]`, `cq[a]`, `logo[a]` point into the flat
//! `p_mat`, `q_mat`, `logo_mat` with stride `max_period`; here the flat
//! vectors are indexed as `a * max_period + i`.
//!
//! Every `log2` call here resolves to `log2f` in the C++ Release object
//! (`<math.h>` brings the float overloads into the global namespace), so all
//! of them are `f32::log2`.

use crate::repeat::RepeatRegion;
use crate::symbol::{char_for_symbol, symbol_for_char};

/// Per-thread splitter state (C++ `SplitWindow` / `s_splitWindow`).
///
/// `Default` matches `new SplitWindow()` (value-initialised: every scalar is
/// zero and every buffer is unallocated).
#[derive(Debug, Clone, Default)]
pub struct SplitWindow {
    pub np: f32,
    pub nq: f32,

    pub log_np: f32,
    pub log_nq: f32,

    /// `max_period` floats (C++ `logo_totals`).
    pub logo_totals: Vec<f32>,

    /// Flat `alphabet_size * max_period` counts; C++ `cp[a][i]`.
    pub p_mat: Vec<f32>,
    /// Flat `alphabet_size * max_period` counts; C++ `cq[a][i]`.
    pub q_mat: Vec<f32>,
    /// Flat `alphabet_size * max_period`; C++ `logo[a][i]`.
    pub logo_mat: Vec<f32>,

    pub plogp: f32,
    pub plogq: f32,
    pub qlogq: f32,
    pub qlogp: f32,

    pub alphabet_size: i32,
    pub max_period: i32,
    pub period: i32,

    pub place: i32,
    pub length: i32,
}

/// `SymbolForChar(c) - 1`: 0..3 for `TCAG`, -1 for anything else.
#[inline]
fn sym_index(c: u8) -> i32 {
    symbol_for_char(c) as i32 - 1
}

impl SplitWindow {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    fn idx(&self, a: i32, i: i32) -> usize {
        (a * self.max_period + i) as usize
    }

    /// C++ `AllocateSplitWindow`. The C++ buffers are malloc'd and left
    /// uninitialised; they are zeroed here (every entry that is read is
    /// written first).
    pub fn allocate_split_window(&mut self, alphabet: i32, max_period: i32) {
        self.alphabet_size = alphabet;
        self.max_period = max_period;

        let n = (max_period as i64 * alphabet as i64).max(0) as usize;
        self.logo_totals = vec![0.0; max_period.max(0) as usize];
        self.p_mat = vec![0.0; n];
        self.q_mat = vec![0.0; n];
        self.logo_mat = vec![0.0; n];
    }

    /// C++ `DeallocSplitWindow`: releases the buffers.
    pub fn dealloc_split_window(&mut self) {
        self.logo_totals = Vec::new();
        self.p_mat = Vec::new();
        self.q_mat = Vec::new();
        self.logo_mat = Vec::new();
    }

    /// C++ `StartSplitWindow`: sets every count of both windows to
    /// `starting_score / (alphabet_size * p)`.
    pub fn start_split_window(&mut self, p: i32, starting_score: f32) {
        self.period = p;
        self.place = 0;

        let d: f32 = starting_score / (self.alphabet_size as f32 * self.period as f32);

        self.plogp = d * d.log2() * self.period as f32 * self.alphabet_size as f32;
        self.plogq = self.plogp;
        self.qlogq = self.plogp;
        self.qlogp = self.plogp;

        self.np = starting_score;
        self.nq = starting_score;

        self.log_np = self.np.log2();
        self.log_nq = self.nq.log2();

        for a in 0..self.alphabet_size {
            for i in 0..self.period {
                let k = self.idx(a, i);
                self.p_mat[k] = d;
                self.q_mat[k] = d;
            }
        }
    }

    /// C++ `slow_KLD_pq` (dead code).
    pub fn slow_kld_pq(&self) -> f32 {
        let mut sum: f32 = 0.0;
        for a in 0..self.alphabet_size {
            for i in 0..self.period {
                let k = self.idx(a, i);
                let p = self.p_mat[k] / self.np;
                let q = self.q_mat[k] / self.nq;
                sum += p * (p / q).log2();
            }
        }
        sum
    }

    /// C++ `slow_KLD_qp` (dead code).
    pub fn slow_kld_qp(&self) -> f32 {
        let mut sum: f32 = 0.0;
        for a in 0..self.alphabet_size {
            for i in 0..self.period {
                let k = self.idx(a, i);
                let p = self.p_mat[k] / self.np;
                let q = self.q_mat[k] / self.nq;
                sum += q * (q / p).log2();
            }
        }
        sum
    }

    #[inline]
    pub fn kld_pq(&self) -> f32 {
        (self.plogp / self.np) - (self.plogq / self.np) + self.log_nq - self.log_np
    }

    #[inline]
    pub fn kld_qp(&self) -> f32 {
        (self.qlogq / self.nq) - (self.qlogp / self.np) + self.log_np - self.log_nq
    }

    /// C++ `score_func()` macro: `KLD_pq() + KLD_qp()`.
    #[inline]
    pub fn score_func(&self) -> f32 {
        self.kld_pq() + self.kld_qp()
    }

    /// C++ `AdjustMatrix`: adds `delta` to count `[symbol][i]` of the `p`
    /// (left) or `q` (right) window, spreading it over the alphabet when
    /// `symbol < 0`, and updates the running `xlogy` sums. Returns `delta`, or
    /// 0 when `i < 0` (an insertion position).
    pub fn adjust_matrix(&mut self, is_p: bool, i: i32, symbol: i32, delta: f32) -> f32 {
        if i < 0 {
            return 0.0;
        }

        let alphabet_size = self.alphabet_size;
        let mp = self.max_period;
        let SplitWindow {
            p_mat,
            q_mat,
            plogp,
            plogq,
            qlogq,
            qlogp,
            ..
        } = self;

        let (aloga, alogb, bloga, ca, cb): (
            &mut f32,
            &mut f32,
            &mut f32,
            &mut Vec<f32>,
            &Vec<f32>,
        ) = if is_p {
            (plogp, plogq, qlogp, p_mat, q_mat)
        } else {
            (qlogq, qlogp, plogq, q_mat, p_mat)
        };

        let mut update = |k: usize, d: f32| {
            *aloga -= ca[k] * ca[k].log2();
            *alogb -= ca[k] * cb[k].log2();
            *bloga -= cb[k] * ca[k].log2();

            ca[k] += d;

            *aloga += ca[k] * ca[k].log2();
            *alogb += ca[k] * cb[k].log2();
            *bloga += cb[k] * ca[k].log2();
        };

        if symbol < 0 {
            let d: f32 = delta / alphabet_size as f32;
            for a in 0..alphabet_size {
                update((a * mp + i) as usize, d);
            }
        } else {
            update((symbol * mp + i) as usize, delta);
        }

        delta
    }

    /// C++ `FillWindow`: counts the left window `[0, window_size)` into `p`
    /// and the right window `[window_size, 2 * window_size)` into `q`, on
    /// top of the pseudo-counts from
    /// [`start_split_window`](Self::start_split_window), then recomputes the
    /// sums from scratch. Requires `r.sequence` and `r.logo_numbers`.
    pub fn fill_window(&mut self, r: &RepeatRegion, window_size: i32) {
        self.start_split_window(r.repeat_period, 1.0);
        self.length = window_size;

        let logo_numbers = r.logo_numbers.as_ref().expect("logoNumbers is null");

        // `1.0 / alphabetSize` is a double, narrowed to float.
        let d: f32 = (1.0f64 / self.alphabet_size as f64) as f32;

        for i in 0..window_size {
            let symbol = sym_index(r.sequence[i as usize]);
            let indx = logo_numbers[i as usize];

            if indx < 0 {
                continue;
            }

            if symbol < 0 {
                for a in 0..self.alphabet_size {
                    let k = self.idx(a, indx);
                    self.p_mat[k] += d;
                }
            } else {
                let k = self.idx(symbol, indx);
                // float += 1.0 (double), narrowed.
                self.p_mat[k] = (self.p_mat[k] as f64 + 1.0) as f32;
            }
        }

        for i in self.length..2 * self.length {
            let symbol = sym_index(r.sequence[i as usize]);
            let indx = logo_numbers[i as usize];

            if indx < 0 {
                continue;
            }
            if symbol < 0 {
                for a in 0..self.alphabet_size {
                    let k = self.idx(a, indx);
                    self.q_mat[k] += d;
                }
            } else {
                let k = self.idx(symbol, indx);
                self.q_mat[k] = (self.q_mat[k] as f64 + 1.0) as f32;
            }
        }

        self.plogp = 0.0;
        self.plogq = 0.0;
        self.qlogq = 0.0;
        self.qlogp = 0.0;
        self.np = 0.0;
        self.nq = 0.0;

        for i in 0..self.period {
            for a in 0..self.alphabet_size {
                let k = self.idx(a, i);
                let cp = self.p_mat[k];
                let cq = self.q_mat[k];
                self.plogp += cp * cp.log2();
                self.plogq += cp * cq.log2();
                self.qlogp += cq * cp.log2();
                self.qlogq += cq * cq.log2();
                self.np += cp;
                self.nq += cq;
            }
        }

        self.log_np = self.np.log2();
        self.log_nq = self.nq.log2();
    }

    /// C++ `MoveWindowForward`: slides both windows one position right.
    pub fn move_window_forward(&mut self, r: &RepeatRegion) {
        let place = self.place as usize;
        let length = self.length as usize;
        let logo_numbers = r.logo_numbers.as_ref().expect("logoNumbers is null");

        let old_left_symbol = sym_index(r.sequence[place]);
        let old_right_symbol = sym_index(r.sequence[place + length]);
        let new_left_symbol = old_right_symbol;
        let new_right_symbol = sym_index(r.sequence[place + length + length]);

        let old_left_indx = logo_numbers[place];
        let old_right_indx = logo_numbers[place + length];
        let new_left_indx = old_right_indx;
        let new_right_indx = logo_numbers[place + length + length];

        let v = self.adjust_matrix(true, old_left_indx, old_left_symbol, -1.0);
        self.np += v;
        let v = self.adjust_matrix(false, old_right_indx, old_right_symbol, -1.0);
        self.nq += v;

        let v = self.adjust_matrix(true, new_left_indx, new_left_symbol, 1.0);
        self.np += v;
        let v = self.adjust_matrix(false, new_right_indx, new_right_symbol, 1.0);
        self.nq += v;

        self.log_np = self.np.log2();
        self.log_nq = self.nq.log2();

        self.place += 1;
    }

    /// C++ `RegionScores` (dead code): the score at every position, with
    /// `window_size` zeros at each end.
    pub fn region_scores(&mut self, r: &RepeatRegion, window_size: i32) -> Vec<f32> {
        let mut scores = Vec::new();

        self.fill_window(r, window_size);

        for _ in 0..self.length {
            scores.push(0.0);
        }
        scores.push(self.score_func());

        // int < size_t - int: computed in size_t.
        let limit = (r.sequence.len() as u64).wrapping_sub(self.length as i64 as u64);
        let mut i = self.length;
        while (i as i64 as u64) < limit {
            self.move_window_forward(r);
            scores.push(self.score_func());
            i += 1;
        }

        for _ in 0..self.length {
            scores.push(0.0);
        }
        scores
    }

    /// C++ `SplitsForRegion`: positions (relative to the repeat start) where
    /// the divergence between the two windows peaks above `threshold`.
    /// Requires `r.sequence` and `r.logo_numbers`, and
    /// `r.sequence.len() >= 2 * window_size`.
    pub fn splits_for_region(
        &mut self,
        r: &RepeatRegion,
        window_size: i32,
        threshold: f32,
    ) -> Vec<i32> {
        let mut splits = Vec::new();

        let mut split_pos: i32 = -1;
        let mut split_val: f32 = threshold;

        self.fill_window(r, window_size);

        let mut score = self.score_func();

        if score > split_val {
            split_val = score;
            split_pos = self.place + self.length;
        }

        // int < size_t - int - int: computed in size_t.
        let limit = (r.sequence.len() as u64)
            .wrapping_sub(self.length as i64 as u64)
            .wrapping_sub(self.length as i64 as u64);
        let mut i: i32 = 0;
        while (i as i64 as u64) < limit {
            self.move_window_forward(r);
            score = self.score_func();
            if score > split_val {
                split_val = score;
                split_pos = self.place + self.length;
            } else if split_pos >= 0 && self.place - split_pos > 0 {
                splits.push(split_pos);
                split_pos = -1;
                split_val = threshold;
            }
            i += 1;
        }

        if split_pos >= 0 {
            splits.push(split_pos);
        }

        splits
    }

    /// C++ `FillLogoForRegion`: counts of each symbol per unit position over
    /// `[start, end)`, with a `1 / alphabet_size` pseudo-count.
    pub fn fill_logo_for_region(&mut self, r: &RepeatRegion, start: i32, end: i32) {
        // `1.0 / (float)alphabetSize` is a double, narrowed to float.
        let d: f32 = (1.0f64 / (self.alphabet_size as f32) as f64) as f32;

        for i in 0..r.repeat_period {
            for a in 0..self.alphabet_size {
                let k = self.idx(a, i);
                self.logo_mat[k] = d;
            }
            self.logo_totals[i as usize] = 0.0;
        }

        let logo_numbers = r.logo_numbers.as_ref().expect("logoNumbers is null");
        for i in start..end {
            let symbol = sym_index(r.sequence[i as usize]);
            let indx = logo_numbers[i as usize];

            if indx < 0 || symbol < 0 {
                continue;
            }

            let k = self.idx(symbol, indx);
            // float += 1.0 (double), narrowed.
            self.logo_mat[k] = (self.logo_mat[k] as f64 + 1.0) as f32;
            self.logo_totals[indx as usize] = (self.logo_totals[indx as usize] as f64 + 1.0) as f32;
        }
    }

    /// C++ `ConsensusForRegion`: per unit position, the symbol whose share
    /// exceeds `consensus_threshold` (the largest one), or `'*'`.
    pub fn consensus_for_region(
        &mut self,
        r: &RepeatRegion,
        start: i32,
        end: i32,
        consensus_threshold: f32,
    ) -> Vec<u8> {
        let mut consensus = vec![b'*'; r.repeat_period.max(0) as usize];

        self.fill_logo_for_region(r, start, end);

        for i in 0..r.repeat_period {
            let mut best_val = consensus_threshold;
            for a in 0..self.alphabet_size {
                let mass = self.logo_mat[self.idx(a, i)] / self.logo_totals[i as usize];
                if mass > best_val {
                    best_val = mass;
                    consensus[i as usize] = char_for_symbol((a + 1) as u8);
                }
            }
        }

        consensus
    }

    /// C++ `ConsensiForSplit`: one consensus per segment, `splits.len() + 1`
    /// entries.
    pub fn consensi_for_split(
        &mut self,
        r: &RepeatRegion,
        splits: &[i32],
        consensus_threshold: f32,
    ) -> Vec<Vec<u8>> {
        let mut consensi = Vec::new();

        let mut start = 0;
        for &s in splits {
            consensi.push(self.consensus_for_region(r, start, s, consensus_threshold));
            start = s;
        }

        // size_t narrowed to the int parameter.
        consensi.push(self.consensus_for_region(
            r,
            start,
            r.sequence.len() as i32,
            consensus_threshold,
        ));

        consensi
    }
}

/// C++ `ConsensusSimilarity`: the fraction of positions where `c1` matches
/// `c2` rotated by `offset` (0.25 for a `'*'`).
///
/// Upstream quirk: the `'*'` test reads `c1[offi]` instead of `c2[offi]`.
pub fn consensus_similarity(c1: &[u8], c2: &[u8], offset: i32) -> f32 {
    let mut similarity: f32 = 0.0;
    for i in 0..c1.len() {
        // (int) % size_t: computed in size_t.
        let offi = ((i as i32 + offset) as i64 as u64 % c2.len() as u64) as usize;
        if c1[i] == c2[offi] {
            similarity = (similarity as f64 + 1.0) as f32;
        } else if c1[i] == b'*' || c1[offi] == b'*' {
            similarity = (similarity as f64 + 0.25) as f32;
        }
    }

    similarity / c1.len() as f32
}

/// C++ `ShouldJoinConsensus`: true if some rotation of `c2` is more than
/// `join_threshold` similar to `c1`.
pub fn should_join_consensus(c1: &[u8], c2: &[u8], join_threshold: f32) -> bool {
    for i in 0..c1.len() {
        if consensus_similarity(c1, c2, i as i32) > join_threshold {
            return true;
        }
    }
    false
}

/// C++ `ValidateSplits`: where consecutive consensi are similar, copies the
/// earlier consensus forward and marks the split between them as `-1`.
pub fn validate_splits(consensi: &mut [Vec<u8>], splits: &mut [i32], join_threshold: f32) {
    for i in 1..consensi.len() {
        if should_join_consensus(&consensi[i - 1], &consensi[i], join_threshold) {
            consensi[i] = consensi[i - 1].clone();
            splits[i - 1] = -1;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Expected bit patterns come from the upstream Release objects
    //! (`RepeatSplitter.cpp.o`). End-to-end splitting is pinned against the
    //! C++ in `repeat::tests`.
    use super::*;

    #[test]
    fn consensus_similarity_matches_cpp() {
        let cases: &[(&str, &str, &[u32])] = &[
            (
                "AC*G",
                "ACGG",
                &[0x3f500000, 0x3e000000, 0x3e000000, 0x3ea00000],
            ),
            ("ACGT", "CGTA", &[0, 0, 0, 0x3f800000]),
            // `c1[offi]` (not `c2[offi]`) decides the '*' credit.
            (
                "A*GT",
                "TTTT",
                &[0x3ea00000, 0x3ec00000, 0x3ea00000, 0x3ec00000],
            ),
            ("ACG", "AC", &[0x3f2aaaab, 0, 0x3f2aaaab]),
            ("**A", "A**", &[0x3ed55555, 0x3f800000, 0x3f000000]),
        ];
        let joins = [
            (false, true),
            (true, true),
            (false, false),
            (false, true),
            (true, true),
        ];
        for ((c1, c2, bits), join) in cases.iter().zip(joins) {
            for (o, &b) in bits.iter().enumerate() {
                assert_eq!(
                    consensus_similarity(c1.as_bytes(), c2.as_bytes(), o as i32).to_bits(),
                    b
                );
            }
            assert_eq!(
                (
                    should_join_consensus(c1.as_bytes(), c2.as_bytes(), 0.85),
                    should_join_consensus(c1.as_bytes(), c2.as_bytes(), 0.5)
                ),
                join
            );
        }
    }

    #[test]
    fn validate_splits_matches_cpp() {
        let mut c: Vec<Vec<u8>> = ["ACGT", "CGTA", "GGGG", "GGGC", "TTTT"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let mut s = vec![10, 20, 30, 40];
        validate_splits(&mut c, &mut s, 0.85);
        let c: Vec<&str> = c.iter().map(|v| std::str::from_utf8(v).unwrap()).collect();
        assert_eq!(c, ["ACGT", "ACGT", "GGGG", "GGGC", "TTTT"]);
        assert_eq!(s, [-1, 20, 30, 40]);
    }

    #[test]
    fn running_sums_match_cpp() {
        let mut w = SplitWindow::new();
        assert_eq!((w.np, w.place, w.alphabet_size), (0.0, 0, 0));
        w.allocate_split_window(4, 11);
        w.start_split_window(7, 1.0);
        let k = w.idx(3, 6);
        assert_eq!(
            [
                w.plogp.to_bits(),
                w.log_np.to_bits(),
                w.p_mat[k].to_bits(),
                w.score_func().to_bits()
            ],
            [0xc099d5da, 0, 0x3d124925, 0]
        );

        w.adjust_matrix(true, 2, 1, 1.0);
        w.adjust_matrix(false, 2, -1, 1.0);
        // Mirrors the C++ test driver: `Np += 1` (int promoted to float).
        w.np += 1.0;
        w.nq += 1.0;
        w.log_np = w.np.log2();
        w.log_nq = w.nq.log2();
        let (kp, kq) = (w.idx(1, 2), w.idx(0, 2));
        assert_eq!(
            [
                w.plogp.to_bits(),
                w.plogq.to_bits(),
                w.qlogq.to_bits(),
                w.qlogp.to_bits(),
                w.p_mat[kp].to_bits(),
                w.q_mat[kq].to_bits(),
                w.score_func().to_bits()
            ],
            [0xc092a9d1, 0xc0c5f4da, 0xc0c5f4d6, 0xc103a0a0, 0x3f849249, 0x3e924925, 0x3fe92ee6]
        );
        assert_eq!(w.adjust_matrix(true, -1, 0, 1.0), 0.0);
    }
}
