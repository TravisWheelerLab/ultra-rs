//! The ULTRA repeat model: transition and emission scores, and the column
//! recurrence that fills the [`UMatrix`]. Port of `umodel.hpp` /
//! `umodel.cpp`.
//!
//! Floating point mirrors the C++ Release build exactly. Every `log2` call in
//! `umodel.cpp` resolves to the C `double log2(double)` (the object file calls
//! `log2`, never `log2f`), so float arguments are widened, the log is taken in
//! double, and the result is narrowed when stored into a float.

use crate::sequence_window::SequenceWindow;
use crate::symbol::Symbol;
use crate::umatrix::{
    previous_column_of, Cell, UMatrix, CT_BACKGROUND, CT_DELETION, CT_INSERTION, CT_MATCH,
    NEG_INF_F32,
};

/// C `log2(double)` applied to a float argument, as `umodel.cpp` does.
#[inline(always)]
fn log2d(x: f64) -> f64 {
    x.log2()
}

/// C++ `EmissionScore` body.
#[inline(always)]
fn emission_score(
    seq: &[Symbol],
    mscore: &[[f32; 5]; 5],
    bscore: &[f32; 5],
    index: i32,
    order: i32,
) -> f32 {
    let s1 = seq[index as usize] as usize;
    let s2 = seq[(index - order) as usize] as usize;
    mscore[s2][s1] - bscore[s1]
}

/// C++ `PreviousEmissionScore` body. Like the C++, does not check that
/// `index - d - order >= 0` (here a negative index panics).
#[inline(always)]
fn previous_emission_score(
    seq: &[Symbol],
    mscore: &[[f32; 5]; 5],
    bscore: &[f32; 5],
    index: i32,
    order: i32,
    d: i32,
) -> f32 {
    let p = index - d;
    let b = bscore[seq[p as usize] as usize];
    let v = mscore[seq[(p - order) as usize] as usize][seq[p as usize] as usize];
    v - b
}

/// The score matrix split around the current column, so the current column
/// can be written while earlier and later columns are read.
struct ColumnSplit<'a> {
    before: &'a [f32],
    c: &'a mut [f32],
    after: &'a [f32],
    cur: usize,
    cpc: usize,
    length: usize,
}

impl<'a> ColumnSplit<'a> {
    fn new(score_matrix: &'a mut [f32], cur: usize, cpc: usize, length: usize) -> Self {
        let (before, rest) = score_matrix.split_at_mut(cur * cpc);
        let (c, after) = rest.split_at_mut(cpc);
        ColumnSplit {
            before,
            c,
            after,
            cur,
            cpc,
            length,
        }
    }

    /// Column `col` (must not be the current column).
    #[inline(always)]
    fn column(&self, col: usize) -> &'a [f32] {
        let (before, after): (&'a [f32], &'a [f32]) = (self.before, self.after);
        if col < self.cur {
            &before[col * self.cpc..(col + 1) * self.cpc]
        } else {
            let k = col - self.cur - 1;
            &after[k * self.cpc..(k + 1) * self.cpc]
        }
    }

    /// C++ `UMatrix::PreviousScore`, reading the current column live.
    #[inline(always)]
    fn previous_score(&self, row: i32, d: i32) -> f32 {
        let mut s = NEG_INF_F32;
        if (d as i64 as u64) < self.length as u64 {
            let col = previous_column_of(self.cur, self.length, d) as usize;
            let row = row as usize;
            s = if col < self.cur {
                self.before[col * self.cpc + row]
            } else if col == self.cur {
                self.c[row]
            } else {
                self.after[(col - self.cur - 1) * self.cpc + row]
            };
        }
        s
    }
}

/// The repeat model (C++ `UModel`). Owns its matrix.
#[derive(Debug, Clone)]
pub struct UModel {
    pub matrix: UMatrix,

    // (tp = transition probabilities)
    pub tp_zero_to_match: f32,
    pub tp_match_to_zero: f32,
    pub tp_match_to_insertion: f32,
    pub tp_match_to_deletion: f32,
    pub tp_consecutive_insertion: f32,
    pub tp_consecutive_deletion: f32,

    /// Before [`calculate_scores`](Self::calculate_scores) this is the decay
    /// probability; `calculate_scores` overwrites it with its log2. Callers
    /// that call `calculate_scores` again must set it again first.
    pub period_decay: f32,
    pub period_decay_offset: f32,

    /// If true then `tscore[BG][MATCH] = log2(tp_zero_to_match / max_period)`.
    pub adjust_for_maximum_period: bool,

    /// `[N][*]`, `[T][*]`, `[C][*]`, `[A][*]`, `[G][*]` (symbol order).
    pub match_probabilities: [[f32; 5]; 5],
    /// N, T, C, A, G.
    pub background_probabilties: [f32; 5],

    pub windowed_weight: f32,

    pub logmp: [[f32; 5]; 5],
    pub lbp: [[f32; 5]; 5],
    pub mscore: [[f32; 5]; 5],
    pub bscore: [f32; 5],
    /// `tscore[from][to]`, indexed by `CT_*`.
    pub tscore: [[f32; 4]; 4],

    pub immediate_transition_to_repeat: bool,
    pub adjust_match_mismatch_matrix: bool,
}

impl UModel {
    /// C++ `UModel(maxPeriod, maxInsertions, maxDeletions, matrixLength)`.
    /// Builds the matrix and calls [`calculate_scores`](Self::calculate_scores)
    /// with the default parameters.
    pub fn new(
        max_period: i32,
        max_insertions: i32,
        max_deletions: i32,
        matrix_length: i32,
    ) -> Self {
        let mut m = UModel {
            matrix: UMatrix::new(max_period, max_insertions, max_deletions, matrix_length),
            tp_zero_to_match: 0.01,
            tp_match_to_zero: 0.02,
            tp_match_to_insertion: 0.07,
            tp_match_to_deletion: 0.05,
            tp_consecutive_insertion: 0.1,
            tp_consecutive_deletion: 0.1,
            period_decay: 0.9,
            period_decay_offset: 0.0,
            adjust_for_maximum_period: false,
            match_probabilities: [
                [1.0, 0.3, 0.2, 0.3, 0.2],
                [1.00, 0.61, 0.13, 0.13, 0.13],
                [1.00, 0.13, 0.61, 0.13, 0.13],
                [1.00, 0.13, 0.13, 0.61, 0.13],
                [1.00, 0.13, 0.13, 0.13, 0.61],
            ],
            background_probabilties: [1.00, 0.3, 0.2, 0.3, 0.2],
            windowed_weight: 0.75,
            // Uninitialised in C++ until CalculateScores, which the
            // constructor calls.
            logmp: [[0.0; 5]; 5],
            lbp: [[0.0; 5]; 5],
            mscore: [[0.0; 5]; 5],
            bscore: [0.0; 5],
            tscore: [[0.0; 4]; 4],
            immediate_transition_to_repeat: false,
            adjust_match_mismatch_matrix: true,
        };
        m.calculate_scores();
        m
    }

    /// Sets the T/C/A/G match probabilities: `p` on the diagonal and
    /// `(1 - p) / 3` elsewhere (computed in double, stored as float).
    pub fn set_match_probabilities(&mut self, p: f32) {
        let mut ip = (1.0f64 - p as f64) as f32;
        ip = (ip as f64 / 3.0) as f32;

        for i in 1..5 {
            for j in 1..5 {
                if i == j {
                    self.match_probabilities[i][j] = p;
                } else {
                    self.match_probabilities[i][j] = ip;
                }
            }
        }
    }

    pub fn set_atcg_probabilities(&mut self, a: f32, t: f32, c: f32, g: f32) {
        use crate::symbol::{N_A, N_C, N_G, N_T};
        self.background_probabilties[N_A as usize] = a;
        self.background_probabilties[N_T as usize] = t;
        self.background_probabilties[N_C as usize] = c;
        self.background_probabilties[N_G as usize] = g;
    }

    /// Derives `tscore`, `logmp`, `mscore`, `lbp`, `bscore`,
    /// `period_decay_offset` and the log `period_decay`. Must be called before
    /// general use.
    pub fn calculate_scores(&mut self) {
        let bg = CT_BACKGROUND as usize;
        let ma = CT_MATCH as usize;
        let ins = CT_INSERTION as usize;
        let del = CT_DELETION as usize;

        self.tscore[bg][bg] = log2d(1.0 - self.tp_zero_to_match as f64) as f32;

        // float sum += pow(float, int): pow is evaluated in double.
        let mut sum: f32 = 0.0;
        let mut i: i32 = 0;
        while (i as i64 as u64) < self.matrix.max_period as u64 {
            sum = (sum as f64 + (self.period_decay as f64).powf(i as f64)) as f32;
            i += 1;
        }

        sum = (1.0f64 / sum as f64) as f32;
        self.period_decay_offset = log2d(sum as f64) as f32;
        self.period_decay = log2d(self.period_decay as f64) as f32;

        if self.adjust_for_maximum_period {
            self.tscore[bg][ma] =
                log2d((self.tp_zero_to_match / self.matrix.max_period as u64 as f32) as f64) as f32;
        } else {
            self.tscore[bg][ma] = log2d(self.tp_zero_to_match as f64) as f32;
        }

        self.tscore[bg][ins] = NEG_INF_F32;
        self.tscore[bg][del] = NEG_INF_F32;

        self.tscore[ma][bg] = log2d(self.tp_match_to_zero as f64) as f32;

        // The sum is in float, the subtraction from 1.0 in double.
        self.tscore[ma][ma] = log2d(
            1.0 - (self.tp_match_to_zero + self.tp_match_to_insertion + self.tp_match_to_deletion)
                as f64,
        ) as f32;
        self.tscore[ma][ins] = log2d(self.tp_match_to_insertion as f64) as f32;
        self.tscore[ma][del] = log2d(self.tp_match_to_deletion as f64) as f32;

        self.tscore[ins][bg] = NEG_INF_F32;
        self.tscore[ins][ma] = 0.0;
        self.tscore[ins][ins] = log2d(self.tp_consecutive_insertion as f64) as f32;
        self.tscore[ins][del] = NEG_INF_F32;

        self.tscore[del][bg] = NEG_INF_F32;
        self.tscore[del][ma] = 0.0;
        self.tscore[del][ins] = NEG_INF_F32;
        self.tscore[del][del] = log2d(self.tp_consecutive_deletion as f64) as f32;

        for i in 1..5 {
            self.match_probabilities[0][i] = self.background_probabilties[i];
        }

        for i in 0..5 {
            for j in 0..5 {
                self.logmp[i][j] = log2d(self.match_probabilities[i][j] as f64) as f32;
                self.mscore[i][j] = self.logmp[i][j];
                // Product in float, log in double.
                self.lbp[i][j] = log2d(
                    (self.background_probabilties[i] * self.background_probabilties[j]) as f64,
                ) as f32;
            }

            self.bscore[i] = log2d(self.background_probabilties[i] as f64) as f32;
        }
    }

    /// `mscore[seq[p - order]][seq[p]] - bscore[seq[p]]` with `p = index - d`.
    /// Does not check that `index - d - order >= 0`.
    #[inline]
    pub fn previous_emission_score(
        &self,
        seq: &SequenceWindow,
        index: i32,
        order: i32,
        d: i32,
    ) -> f32 {
        previous_emission_score(seq.seq(), &self.mscore, &self.bscore, index, order, d)
    }

    /// `mscore[seq[index - order]][seq[index]] - bscore[seq[index]]`.
    #[inline]
    pub fn emission_score(&self, seq: &SequenceWindow, index: i32, order: i32) -> f32 {
        emission_score(seq.seq(), &self.mscore, &self.bscore, index, order)
    }

    /// Fills the current column for symbol `nuc_index` of `sequence` from
    /// the previous column, then moves the matrix forward one column.
    /// `can_be_repetitive` is accepted and ignored, as upstream.
    pub fn calculate_current_column(
        &mut self,
        sequence: &SequenceWindow,
        nuc_index: i32,
        _can_be_repetitive: Option<&[bool]>,
    ) {
        let seq = sequence.seq();
        let mscore = &self.mscore;
        let bscore = &self.bscore;
        let tscore = &self.tscore;
        let period_decay = self.period_decay;
        let period_decay_offset = self.period_decay_offset;
        let immediate = self.immediate_transition_to_repeat;

        let bg = CT_BACKGROUND as usize;
        let ma = CT_MATCH as usize;
        let ins = CT_INSERTION as usize;
        let del = CT_DELETION as usize;

        let m = &mut self.matrix;
        let cpc = m.cells_per_column;
        let cur = m.current_column_index;
        let prev = m.previous_column_index;
        let matind = prev as i32;
        let tbr = m.max_period + 1;

        let desc: &[Cell] = &m.cell_descriptions[..cpc];
        let ct = &mut m.traceback_matrix[cur * tbr..(cur + 1) * tbr];
        let cols = ColumnSplit::new(&mut m.score_matrix, cur, cpc, m.length);
        // The previous column is never the current one (it is always cur - 1).
        let p: &[f32] = cols.column(prev);

        for (row, &d) in desc.iter().enumerate() {
            match d.cell_type {
                CT_BACKGROUND => {
                    cols.c[row] = p[row] + tscore[bg][bg];
                    ct[0] = row as i32;
                }

                CT_MATCH => {
                    let parent_index = d.parent_index as usize;
                    let order = d.order;
                    let mut score = p[row] + tscore[ma][bg];
                    score += order as f32 * period_decay;
                    score += period_decay_offset;

                    if score > cols.c[parent_index] {
                        cols.c[parent_index] = score;
                        ct[0] = row as i32;
                    }

                    if matind < order {
                        cols.c[row] = NEG_INF_F32;
                        ct[order as usize] = 0;
                        continue;
                    }

                    let mut zscore: f32;
                    if immediate {
                        // p[0] is the previous column; the loop has not
                        // written it.
                        zscore = p[0];
                    } else {
                        zscore = cols.previous_score(0, order - 1);
                        zscore += tscore[bg][bg] * order as f32;
                    }

                    let es = emission_score(seq, mscore, bscore, nuc_index, order);

                    score = zscore + tscore[bg][ma] + es;

                    cols.c[row] = score;
                    ct[order as usize] = 0;

                    score = p[row] + tscore[ma][ma] + es;
                    if score > cols.c[row] {
                        cols.c[row] = score;
                        ct[order as usize] = row as i32;
                    }
                }

                CT_INSERTION => {
                    let order = d.order;
                    let parent_index = d.parent_index;
                    let peridocity = order + d.indel_number;
                    let min_index = peridocity + order + 1;

                    if matind > min_index {
                        let mut score = p[row]
                            + tscore[ins][ma]
                            + emission_score(seq, mscore, bscore, nuc_index, order);
                        if score > cols.c[parent_index as usize] {
                            cols.c[parent_index as usize] = score;
                            ct[order as usize] = row as i32;
                        }

                        score = p[row];

                        score -= cols.previous_score(parent_index, peridocity + 1);
                        score += cols.previous_score(parent_index, peridocity);

                        score -= previous_emission_score(
                            seq, mscore, bscore, nuc_index, peridocity, order,
                        );
                        score += emission_score(seq, mscore, bscore, nuc_index, peridocity);

                        cols.c[row] = score;
                    } else if matind == min_index {
                        let mut score = cols.previous_score(parent_index, peridocity);
                        // float += 0.0 (double)
                        score = (score as f64 + 0.0) as f32;

                        score += tscore[ma][ins];
                        score += tscore[ins][ins] * (d.indel_number - 1) as f32;
                        // Each of the `order` bases after the indel is a
                        // match, so it pays the match-to-match transition, as
                        // in a match state. Upstream 1.2.2 charges
                        // log2(1 - tp_match_to_zero) instead, which is cheaper
                        // per base, so a period P - 1 (or k * P - 1) path with
                        // an insertion each copy outscored the true period P
                        // on long units.
                        score =
                            (score as f64 + tscore[ma][ma] as f64 * (order as f32) as f64) as f32;

                        for j in 0..order {
                            score += emission_score(seq, mscore, bscore, nuc_index - j, peridocity);
                        }

                        cols.c[row] = score;
                    } else {
                        cols.c[row] = NEG_INF_F32;
                    }
                }

                CT_DELETION => {
                    let parent_index = d.parent_index;
                    let order = d.order;
                    let indel_num = d.indel_number;
                    let min_index = 3 * order;

                    let del_order = order - indel_num;
                    let e_order = 2 * order - indel_num;

                    if matind > min_index {
                        let mut score = p[row]
                            + tscore[del][ma]
                            + emission_score(seq, mscore, bscore, nuc_index, order);
                        if score > cols.c[parent_index as usize] {
                            cols.c[parent_index as usize] = score;
                            ct[order as usize] = row as i32;
                        }

                        score = p[row];
                        score -= cols.previous_score(parent_index, order + 1);
                        score += cols.previous_score(parent_index, order);

                        score -= previous_emission_score(
                            seq, mscore, bscore, nuc_index, del_order, order,
                        );
                        score += previous_emission_score(
                            seq, mscore, bscore, nuc_index, del_order, indel_num,
                        );

                        score -= previous_emission_score(
                            seq, mscore, bscore, nuc_index, e_order, indel_num,
                        );
                        score += emission_score(seq, mscore, bscore, nuc_index, e_order);

                        cols.c[row] = score;
                    } else if matind == min_index {
                        let mut score = cols.previous_score(parent_index, order);

                        score += tscore[ma][del];
                        score += tscore[del][del] * (indel_num - 1) as f32;
                        // Match-to-match cost per base; see the insertion
                        // case of `calculate_current_column`.
                        score =
                            (score as f64 + tscore[ma][ma] as f64 * (order as f32) as f64) as f32;

                        for j in 0..(order - indel_num) {
                            let place = nuc_index - order + 1 + j;
                            score += emission_score(seq, mscore, bscore, place, del_order);
                        }

                        for j in 0..indel_num {
                            let place = nuc_index - j;
                            score += emission_score(seq, mscore, bscore, place, e_order);
                        }

                        cols.c[row] = score;
                    } else {
                        cols.c[row] = NEG_INF_F32;
                    }
                }

                _ => {}
            }
        }

        m.move_matrix_forward();
    }

    /// Like [`calculate_current_column`](Self::calculate_current_column) with
    /// every emission score taken as 0. ULTRA calls it once after the last
    /// symbol of a window.
    pub fn calculate_current_column_without_emission(
        &mut self,
        _can_be_repetitive: Option<&[bool]>,
    ) {
        let tscore = &self.tscore;
        let period_decay = self.period_decay;
        let period_decay_offset = self.period_decay_offset;
        let immediate = self.immediate_transition_to_repeat;

        let bg = CT_BACKGROUND as usize;
        let ma = CT_MATCH as usize;
        let ins = CT_INSERTION as usize;
        let del = CT_DELETION as usize;

        let m = &mut self.matrix;
        let cpc = m.cells_per_column;
        let cur = m.current_column_index;
        let prev = m.previous_column_index;
        let matind = prev as i32;
        let tbr = m.max_period + 1;

        let desc: &[Cell] = &m.cell_descriptions[..cpc];
        let ct = &mut m.traceback_matrix[cur * tbr..(cur + 1) * tbr];
        let cols = ColumnSplit::new(&mut m.score_matrix, cur, cpc, m.length);
        let p: &[f32] = cols.column(prev);

        for (row, &d) in desc.iter().enumerate() {
            match d.cell_type {
                CT_BACKGROUND => {
                    cols.c[row] = p[row] + tscore[bg][bg];
                    ct[0] = row as i32;
                }

                CT_MATCH => {
                    let parent_index = d.parent_index as usize;
                    let order = d.order;
                    let mut score = p[row] + tscore[ma][bg];
                    score += order as f32 * period_decay;
                    score += period_decay_offset;

                    if score > cols.c[parent_index] {
                        cols.c[parent_index] = score;
                        ct[0] = row as i32;
                    }

                    if matind < order {
                        cols.c[row] = NEG_INF_F32;
                        ct[order as usize] = 0;
                        continue;
                    }

                    let mut zscore: f32;
                    if immediate {
                        zscore = p[0];
                    } else {
                        zscore = cols.previous_score(0, order - 1);
                        zscore += tscore[bg][bg] * order as f32;
                    }

                    let es: f32 = 0.0;

                    score = zscore + tscore[bg][ma] + es;

                    cols.c[row] = score;
                    ct[order as usize] = 0;

                    score = p[row] + tscore[ma][ma] + es;
                    if score > cols.c[row] {
                        cols.c[row] = score;
                        ct[order as usize] = row as i32;
                    }
                }

                CT_INSERTION => {
                    let order = d.order;
                    let parent_index = d.parent_index;
                    let peridocity = order + d.indel_number;
                    let min_index = peridocity + order + 1;

                    if matind > min_index {
                        let mut score = p[row] + tscore[ins][ma];
                        if score > cols.c[parent_index as usize] {
                            cols.c[parent_index as usize] = score;
                            ct[order as usize] = row as i32;
                        }

                        score = p[row];

                        score -= cols.previous_score(parent_index, peridocity + 1);
                        score += cols.previous_score(parent_index, peridocity);

                        cols.c[row] = score;
                    } else if matind == min_index {
                        let mut score = cols.previous_score(parent_index, peridocity);
                        score = (score as f64 + 0.0) as f32;

                        score += tscore[ma][ins];
                        score += tscore[ins][ins] * (d.indel_number - 1) as f32;
                        // Match-to-match cost per base; see the insertion
                        // case of `calculate_current_column`.
                        score =
                            (score as f64 + tscore[ma][ma] as f64 * (order as f32) as f64) as f32;

                        cols.c[row] = score;
                    } else {
                        cols.c[row] = NEG_INF_F32;
                    }
                }

                CT_DELETION => {
                    let parent_index = d.parent_index;
                    let order = d.order;
                    let indel_num = d.indel_number;
                    let min_index = 3 * order;

                    if matind > min_index {
                        let mut score = p[row] + tscore[del][ma];
                        if score > cols.c[parent_index as usize] {
                            cols.c[parent_index as usize] = score;
                            ct[order as usize] = row as i32;
                        }

                        score = p[row];
                        score -= cols.previous_score(parent_index, order + 1);
                        score += cols.previous_score(parent_index, order);

                        cols.c[row] = score;
                    } else if matind == min_index {
                        let mut score = cols.previous_score(parent_index, order);

                        score += tscore[ma][del];
                        score += tscore[del][del] * (indel_num - 1) as f32;
                        // Match-to-match cost per base; see the insertion
                        // case of `calculate_current_column`.
                        score =
                            (score as f64 + tscore[ma][ma] as f64 * (order as f32) as f64) as f32;

                        cols.c[row] = score;
                    } else {
                        cols.c[row] = NEG_INF_F32;
                    }
                }

                _ => {}
            }
        }

        m.move_matrix_forward();
    }
}

#[cfg(test)]
mod tests {
    //! Expected values were dumped from the upstream C++ (`umodel.cpp`,
    //! `umatrix.cpp`, g++ -O3) driven the same way.

    use super::*;
    use crate::symbol::symbol_for_char;

    const SEQ_A: &str = "ACGTTGCANNNNGCAGCAGCAGCAGCATCAGCAGCAGCGCAGCAGCAGGTTACGATATATATANATATATATCGATGCTAGTACCAGTCAGTTAGCAGTCAGTTAGCAGTCATCAGTTAGCAGGTCAGTTAGCACCGTA";
    const SEQ_B: &str = "TTAGGGTTAGGGTTAGGTTAGGGTTAGGGTTAGGGATCGNNACACACACACAC";

    fn bits4(v: &[f32; 4]) -> [u32; 4] {
        v.map(f32::to_bits)
    }

    fn bits5(v: &[f32; 5]) -> [u32; 5] {
        v.map(f32::to_bits)
    }

    fn ints(s: &str) -> Vec<i32> {
        s.split_whitespace().map(|x| x.parse().unwrap()).collect()
    }

    /// Model configured as `Ultra`'s constructor does with default settings.
    fn ultra_model(
        mp: i32,
        ins: i32,
        del: i32,
        mlen: i32,
        adjmax: bool,
        immediate: bool,
    ) -> UModel {
        let mut m = UModel::new(mp, ins, del, mlen);
        m.period_decay = 0.85;
        m.set_match_probabilities(0.7);
        m.set_atcg_probabilities(0.25, 0.25, 0.25, 0.25);
        m.tp_zero_to_match = 0.01;
        m.tp_match_to_zero = 0.05;
        m.tp_match_to_insertion = 0.02;
        m.tp_match_to_deletion = 0.02;
        m.tp_consecutive_insertion = 0.02;
        m.tp_consecutive_deletion = 0.02;
        m.adjust_for_maximum_period = adjmax;
        m.immediate_transition_to_repeat = immediate;
        m.calculate_scores();
        m
    }

    fn fold(h: u64, v: u64) -> u64 {
        (h ^ v).wrapping_mul(0x100000001b3)
    }

    /// Runs one window as `Ultra::AnalyzeSequenceWindow` does and returns an
    /// FNV-style hash of every computed score column and traceback column.
    fn run_window(m: &mut UModel, w: &mut SequenceWindow, s: &str) -> u64 {
        w.length = s.len();
        w.overlap = 0;
        for (i, c) in s.bytes().enumerate() {
            w.seq_mem[i] = symbol_for_char(c);
        }
        let mut h: u64 = 0xcbf29ce484222325;
        let absorb = |m: &UModel, h: &mut u64| {
            let x = &m.matrix;
            let col = x.previous_column_index;
            for v in x.score_column(col) {
                *h = fold(*h, v.to_bits() as u64);
            }
            for v in &x.traceback_column(col)[..x.max_period] {
                *h = fold(*h, *v as u32 as u64);
            }
        };
        m.matrix.restart_matrix();
        for i in 0..s.len() as i32 {
            m.calculate_current_column(w, i, None);
            absorb(m, &mut h);
        }
        m.calculate_current_column_without_emission(None);
        absorb(m, &mut h);
        let pc = m.matrix.previous_column_index;
        m.matrix.calculate_traceback(pc);
        h
    }

    #[test]
    fn constructor_scores_match_cpp() {
        let m = UModel::new(11, 10, 10, 20);
        assert_eq!(m.period_decay.to_bits(), 0xbe1ba6b5);
        assert_eq!(m.period_decay_offset.to_bits(), 0xc031d4b1);
        assert_eq!(
            bits4(&m.tscore[0]),
            [0xbc6d8f9a, 0xc0d49a78, 0xd51184e7, 0xd51184e7]
        );
        assert_eq!(
            bits4(&m.tscore[1]),
            [0xc0b49a78, 0xbe5ed04a, 0xc075893d, 0xc08a4d3c]
        );
        assert_eq!(
            bits4(&m.tscore[2]),
            [0xd51184e7, 0x00000000, 0xc0549a78, 0xd51184e7]
        );
        assert_eq!(
            bits5(&m.logmp[1]),
            [0, 0xbf368ef4, 0xc03c60f0, 0xc03c60f0, 0xc03c60f0]
        );
        assert_eq!(
            bits5(&m.lbp[2]),
            [0xc0149a78, 0xc081e275, 0xc0949a78, 0xc081e275, 0xc0949a78]
        );
        assert_eq!(
            bits5(&m.bscore),
            [0, 0xbfde54e3, 0xc0149a78, 0xbfde54e3, 0xc0149a78]
        );
    }

    #[test]
    fn ultra_scores_match_cpp() {
        let m = ultra_model(11, 10, 10, 20, false, false);
        assert_eq!(m.period_decay.to_bits(), 0xbe7017a6);
        assert_eq!(m.period_decay_offset.to_bits(), 0xc01e41b1);
        assert_eq!(
            bits4(&m.tscore[0]),
            [0xbc6d8f9a, 0xc0d49a78, 0xd51184e7, 0xd51184e7]
        );
        assert_eq!(
            bits4(&m.tscore[1]),
            [0xc08a4d3c, 0xbe0b53b8, 0xc0b49a78, 0xc0b49a78]
        );
        assert_eq!(
            bits4(&m.tscore[3]),
            [0xd51184e7, 0x00000000, 0xd51184e7, 0xc0b49a78]
        );
        assert_eq!(
            bits5(&m.match_probabilities[1]),
            [0x3f800000, 0x3f333333, 0x3dcccccd, 0x3dcccccd, 0x3dcccccd]
        );
        assert_eq!(
            bits5(&m.mscore[2]),
            [0, 0xc0549a78, 0xbf03bb12, 0xc0549a78, 0xc0549a78]
        );
        assert_eq!(
            bits5(&m.lbp[1]),
            [0xc0000000, 0xc0800000, 0xc0800000, 0xc0800000, 0xc0800000]
        );
    }

    #[test]
    fn columns_and_tracebacks_match_cpp() {
        let mlen = SEQ_A.len() as i32 + 3;
        let mut m = ultra_model(11, 10, 10, mlen, false, false);
        assert_eq!(m.matrix.cells_per_column, 101);
        let mut w = SequenceWindow::new(mlen as usize, 0);

        assert_eq!(run_window(&mut m, &mut w, SEQ_A), 0x93e4fb65ae905247);
        assert_eq!(
            (
                m.matrix.current_column_index,
                m.matrix.previous_column_index
            ),
            (141, 140)
        );
        assert_eq!(m.matrix.score_at(140, 0).to_bits(), 0x41b732bf);
        assert_eq!(m.matrix.score_at(140, 7).to_bits(), 0x40c9e268);
        let trace = ints(
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 8 5 5 \
             5 5 5 5 5 5 5 5 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 \
             0 0 0 0 0 0 0 82 83 82 82 82 82 82 82 82 82 82 82 82 82 83 82 82 82 82 82 82 82 82 \
             82 82 82 82 82 82 82 82 82 82 84 82 82 82 82 82 82 82 82 82 82 82 82 82 82 82 82 82 \
             82 0 0 0",
        );
        assert_eq!(m.matrix.traceback, trace);

        let mut fwd = vec![-7; m.matrix.length + 2];
        assert!(m.matrix.forward_traceback(&mut fwd, SEQ_A.len() as i32, 0));
        let want_fwd = ints(
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 8 5 5 \
             5 5 5 5 5 5 5 5 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 \
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -7 -7 -7 -7",
        );
        assert_eq!(fwd, want_fwd);

        // A second, shorter window reuses the matrix.
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x29e911bc00278209);
        assert_eq!(m.matrix.traceback[15], 32);

        m.matrix.adjust_score_to_zero(5);
        assert_eq!(m.matrix.last_score_adjustment.to_bits(), 0xc1b8de78);
        assert_eq!(m.matrix.total_score_adjustment.to_bits(), 0xc1b8de78);
    }

    #[test]
    fn no_indel_model_matches_cpp() {
        let mlen = SEQ_A.len() as i32 + 3;
        let mut m = ultra_model(7, 0, 0, mlen, false, false);
        let mut w = SequenceWindow::new(mlen as usize, 0);
        assert_eq!(run_window(&mut m, &mut w, SEQ_A), 0x40b63a9eb7c6d539);
    }

    #[test]
    fn adjusted_and_immediate_models_match_cpp() {
        let mlen = SEQ_B.len() as i32 + 3;
        let mut m = ultra_model(5, 2, 1, mlen, true, false);
        let mut w = SequenceWindow::new(mlen as usize, 0);
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0xa5b59ec6f3742034);
        assert_eq!(
            run_window(&mut m, &mut w, "ACGNACGACGACG"),
            0x48116a8c4a543ec4
        );
        // Entries past the second window are stale from the first.
        let trace = ints(
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
             0 2 2 2 2 2 2 2 2 2 2 2 0 0 0",
        );
        assert_eq!(m.matrix.traceback, trace);

        let mut m = ultra_model(9, 3, 2, mlen, false, true);
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x3a614baf77750b0);
    }
}
