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

/// Emission scores by sequence position for the steady-state column. Every
/// emission the column scores is `diff[seq[a]][seq[b]]` (`diff[x][y] =
/// mscore[x][y] - bscore[y]`) with one of the two symbols fixed for a whole
/// run of the other, so with ten arrays, `against(y)[a] = diff[seq[a]][y]`
/// and `from(x)[a] = diff[x][seq[a]]`, each run is a contiguous slice.
///
/// A position's values depend only on its symbol and on `diff`, so the
/// cache keeps the symbol each position was filled from. Before every
/// column, [`fill`](Self::fill) checks the positions the column will read
/// against the sequence as it is now and refills from the first one that
/// differs, then fills up to the current position (never ahead of it). The
/// cached values are therefore always those the general path would compute
/// from the live sequence, whatever the caller does between columns
/// (writing symbols as it goes, or a new window without
/// `restart_matrix`). Everything is refilled when any bit of `mscore` or
/// `bscore` changes.
#[derive(Debug, Clone, Default)]
pub struct PositionEmissions {
    /// Bits of the `mscore` and `bscore` the values were computed from.
    mscore_bits: [[u32; 5]; 5],
    bscore_bits: [u32; 5],
    diff: [[f32; 5]; 5],
    /// Positions `0..filled` hold values; `symbols[a]` is the symbol
    /// position `a` was filled from.
    filled: usize,
    symbols: Vec<Symbol>,
    stride: usize,
    values: Vec<f32>,
}

impl PositionEmissions {
    /// `diff[seq[a]][y]` for every filled position `a`.
    #[inline(always)]
    fn against(&self, y: Symbol) -> &[f32] {
        let k = y as usize;
        &self.values[k * self.stride..k * self.stride + self.filled]
    }

    /// `diff[x][seq[a]]` for every filled position `a`.
    #[inline(always)]
    fn from(&self, x: Symbol) -> &[f32] {
        let k = 5 + x as usize;
        &self.values[k * self.stride..k * self.stride + self.filled]
    }

    /// Makes positions `upto - back ..= upto` of `seq` hold the values for
    /// their current symbols and these `mscore` and `bscore`; `back` is how
    /// far back the column reads.
    fn fill(
        &mut self,
        seq: &[Symbol],
        upto: usize,
        back: usize,
        mscore: &[[f32; 5]; 5],
        bscore: &[f32; 5],
    ) {
        // Any bit of mscore or bscore changed (branch-free: this runs for
        // every column).
        let mut changed = 0u32;
        for x in 0..5 {
            for y in 0..5 {
                changed |= mscore[x][y].to_bits() ^ self.mscore_bits[x][y];
            }
            changed |= bscore[x].to_bits() ^ self.bscore_bits[x];
        }
        if changed != 0 || self.values.is_empty() {
            self.mscore_bits = mscore.map(|row| row.map(f32::to_bits));
            self.bscore_bits = bscore.map(f32::to_bits);
            self.filled = 0;
            // diff[x][y] = mscore[x][y] - bscore[y]: the emission score of
            // symbol y against x, the same float as the general path.
            for x in 0..5 {
                for y in 0..5 {
                    self.diff[x][y] = mscore[x][y] - bscore[y];
                }
            }
        }
        if seq.len() > self.stride {
            self.stride = seq.len();
            self.values = vec![0.0; 10 * self.stride];
            self.symbols = vec![0; self.stride];
            self.filled = 0;
        }
        // Positions the column reads that are already filled: keep those up
        // to the first whose symbol has changed.
        let lo = upto.saturating_sub(back);
        let checked_to = self.filled.min(upto + 1);
        if lo < checked_to && self.symbols[lo..checked_to] != seq[lo..checked_to] {
            let first = (lo..checked_to).find(|&a| self.symbols[a] != seq[a]).unwrap_or(checked_to);
            self.filled = first;
        }
        if self.filled > upto {
            return;
        }
        let diff = self.diff;
        let stride = self.stride;
        let from = self.filled;
        let to = upto + 1;
        let symbols = &seq[from..to];
        self.symbols[from..to].copy_from_slice(symbols);
        for y in 0..5 {
            let against = &mut self.values[y * stride + from..y * stride + to];
            for (v, &s) in against.iter_mut().zip(symbols) {
                *v = diff[s as usize][y];
            }
        }
        for x in 0..5 {
            let k = 5 + x;
            let by = &mut self.values[k * stride + from..k * stride + to];
            for (v, &s) in by.iter_mut().zip(symbols) {
                *v = diff[x][s as usize];
            }
        }
        self.filled = to;
    }
}

/// The current column (written) and the previous column (read) of the
/// matrix, every row, plus the background and match rows of every column.
struct ColumnSplit<'a> {
    p: &'a [f32],
    c: &'a mut [f32],
    history: &'a [f32],
    desc: &'a [Cell],
    history_stride: usize,
    cur: usize,
    length: usize,
}

impl<'a> ColumnSplit<'a> {
    /// `full` is `UMatrix::full_columns`; the previous column is `cur - 1`.
    fn new(
        full: &'a mut [f32],
        history: &'a [f32],
        desc: &'a [Cell],
        cur: usize,
        history_stride: usize,
        length: usize,
    ) -> Self {
        let cpc = desc.len();
        let (slot0, slot1) = full.split_at_mut(cpc);
        let (c, p) = if cur & 1 == 0 { (slot0, &*slot1) } else { (slot1, &*slot0) };
        ColumnSplit { p, c, history, desc, history_stride, cur, length }
    }

    /// C++ `UMatrix::PreviousScore`, reading the current column live.
    /// `row` is the background row or a match row.
    #[inline(always)]
    fn previous_score(&self, row: i32, d: i32) -> f32 {
        let mut s = NEG_INF_F32;
        if (d as i64 as u64) < self.length as u64 {
            let col = previous_column_of(self.cur, self.length, d) as usize;
            let row = row as usize;
            s = if col == self.cur {
                self.c[row]
            } else {
                let cell = self.desc[row];
                debug_assert!(cell.cell_type == CT_MATCH || cell.cell_type == CT_BACKGROUND);
                self.history[cell.order as usize * self.history_stride + col]
            };
        }
        s
    }
}

/// `LANES` consecutive elements of `s` from `start`, as an array.
#[inline(always)]
fn block<T: Copy + Default>(s: &[T], start: usize) -> [T; LANES] {
    let mut a = [T::default(); LANES];
    a.copy_from_slice(&s[start..start + LANES]);
    a
}

/// Cells per block in the steady column's indel loops.
const LANES: usize = 4;

/// Whether a model with these indel counts fills its steady columns with
/// the blocked kernel, which also keeps [`POSITION_EMISSION_ARRAYS`] floats
/// per window position.
pub fn uses_blocked_kernel(max_insertions: usize, max_deletions: usize) -> bool {
    max_insertions / LANES + max_deletions / LANES >= 2
}

/// Floats per window position that the blocked kernel keeps (see
/// `PositionEmissions`).
pub const POSITION_EMISSION_ARRAYS: u64 = 10;

/// Updates `best` / `best_from` with the exits of consecutive rows
/// `first_row..` as the general path does, one at a time, keeping the first
/// strictly greater. The block's maximum only decides whether any exit needs
/// that check, so the result is the same.
#[inline(always)]
fn take_exits(exit: &[f32], best: &mut f32, best_from: &mut i32, first_row: usize) {
    let mut m = exit[0];
    for &e in &exit[1..] {
        m = m.max(e);
    }
    if m > *best {
        for (l, &e) in exit.iter().enumerate() {
            if e > *best {
                *best = e;
                *best_from = (first_row + l) as i32;
            }
        }
    }
}

/// Per-order scratch for the steady-state column.
#[derive(Debug, Clone, Default)]
struct OrderScratch {
    best: Vec<f32>,
    best_from: Vec<i32>,
    es: Vec<f32>,
}

impl OrderScratch {
    /// The three arrays, at least `n` long.
    fn for_orders(&mut self, n: usize) -> (&mut [f32], &mut [i32], &mut [f32]) {
        if self.best.len() < n {
            self.best.resize(n, 0.0);
            self.best_from.resize(n, 0);
            self.es.resize(n, 0.0);
        }
        (&mut self.best, &mut self.best_from, &mut self.es)
    }
}

/// One order's insertion and deletion cells of a steady column, filled by
/// [`SteadyOrder::fill`]. Kept out of line so that the loop over orders, which
/// most orders pass through with few or no indel cells, stays small.
struct SteadyOrder<'a> {
    seq: &'a [Symbol],
    i: usize,
    o: usize,
    /// The order's match row.
    r: usize,
    cur: usize,
    n_ins: usize,
    n_del: usize,
    /// The order's match emission, emission_score(i, o).
    es: f32,
    t_in_ma: f32,
    t_de_ma: f32,
    /// This order's history row.
    h_o: &'a [f32],
    pe: &'a PositionEmissions,
}

impl SteadyOrder<'_> {
    /// Writes the indel rows of the current column `c` from the previous
    /// column `p`, and takes their exits into `best_in` / `best_from_in`.
    #[inline(never)]
    fn fill(&self, p: &[f32], c: &mut [f32], best_in: &mut f32, best_from_in: &mut i32) {
        let SteadyOrder { seq, i, o, r, cur, n_ins, n_del, es, t_in_ma, t_de_ma, h_o, pe } = *self;
        let against = |y: Symbol| pe.against(y);
        let from = |x: Symbol| pe.from(x);
        let emit_i = against(seq[i]);
        let mut best = *best_in;
        let mut best_from = *best_from_in;
        // Every emission of this order's indel cells is a contiguous run
        // of one of the position arrays, and every look-back a
        // contiguous run of this order's history, so the indel loops run
        // on blocks of LANES cells. Within a cell the operations and
        // their order are those of the general path.
        // previous_emission_score(i, q, d) = diff[seq[i - d - q]][seq[i - d]].
        let at_o = against(seq[i - o]);
        let parent_o = h_o[cur - o];
        let parent_o1 = h_o[cur - o - 1];

        // Fewer cells than a block, and the cells after the last block,
        // are filled one at a time straight from the arrays.
        let ins_cell = |k: usize, pr: f32| -> f32 {
            let mut v = pr;
            v -= h_o[cur - o - k - 1];
            v += h_o[cur - o - k];
            v -= at_o[i - 2 * o - k];
            v += emit_i[i - o - k];
            v
        };
        if n_ins < LANES {
            for k in 1..=n_ins {
                let row = r + k;
                let pr = p[row];
                let exit = pr + t_in_ma + es;
                if exit > best {
                    best = exit;
                    best_from = row as i32;
                }
                c[row] = ins_cell(k, pr);
            }
        } else {
            let n = n_ins;
            // Index n - k (k = 1..=n) of each slice: insertion k's match
            // score o + k + 1 and o + k columns back, and its emissions,
            // i - 2o - k against i - o and i - o - k against i.
            let further = &h_o[cur - o - 1 - n..cur - o - 1];
            let nearer = &h_o[cur - o - n..cur - o];
            let e_back = &at_o[i - 2 * o - n..i - 2 * o];
            let e_here = &emit_i[i - o - n..i - o];
            let rows = r + 1..r + 1 + n;
            let pin = &p[rows.clone()];
            let cin = &mut c[rows];
            let one = |pr: f32, fu: f32, ne: f32, eb: f32, eh: f32| -> f32 {
                let mut v = pr;
                v -= fu;
                v += ne;
                v -= eb;
                v += eh;
                v
            };
            let mut j = 0;
            while j + LANES <= n {
                // Cells j..j + LANES; the reversed slices hold them at
                // n - LANES - j.. in reverse.
                let q = n - LANES - j;
                let (pv, fv, nv) = (block(pin, j), block(further, q), block(nearer, q));
                let (ebv, ehv) = (block(e_back, q), block(e_here, q));
                let mut v = [0f32; LANES];
                let mut exit = [0f32; LANES];
                for l in 0..LANES {
                    let x = LANES - 1 - l;
                    v[l] = one(pv[l], fv[x], nv[x], ebv[x], ehv[x]);
                    exit[l] = pv[l] + t_in_ma + es;
                }
                cin[j..j + LANES].copy_from_slice(&v);
                take_exits(&exit, &mut best, &mut best_from, r + 1 + j);
                j += LANES;
            }
            while j < n {
                let pr = pin[j];
                cin[j] = ins_cell(j + 1, pr);
                take_exits(&[pr + t_in_ma + es], &mut best, &mut best_from, r + 1 + j);
                j += 1;
            }
        }

        if n_del > 0 && n_del < LANES {
            let from_o = from(seq[i - o]);
            let from_2o = from(seq[i - 2 * o]);
            for k in 1..=n_del {
                let row = r + n_ins + k;
                let pr = p[row];
                let exit = pr + t_de_ma + es;
                if exit > best {
                    best = exit;
                    best_from = row as i32;
                }
                let mut v = pr;
                v -= parent_o1;
                v += parent_o;
                v -= at_o[i - 2 * o + k];
                v += from_o[i - k];
                v -= from_2o[i - k];
                v += emit_i[i - 2 * o + k];
                c[row] = v;
            }
        } else if n_del > 0 {
            let n = n_del;
            // Deletion k's emissions: i - 2o + k against i - o and i
            // (index k - 1), and i - k against i - o and i - 2o (index
            // n - k).
            let e_at = &at_o[i - 2 * o + 1..i - 2 * o + 1 + n];
            let e_up = &emit_i[i - 2 * o + 1..i - 2 * o + 1 + n];
            let e_by_o = &from(seq[i - o])[i - n..i];
            let e_by_2o = &from(seq[i - 2 * o])[i - n..i];
            let first = r + n_ins + 1;
            let rows = first..first + n;
            let pdel = &p[rows.clone()];
            let cdel = &mut c[rows];
            let one = |pr: f32, a: f32, b: f32, d: f32, u: f32| -> f32 {
                let mut v = pr;
                v -= parent_o1;
                v += parent_o;
                v -= a;
                v += b;
                v -= d;
                v += u;
                v
            };
            let mut j = 0;
            while j + LANES <= n {
                let q = n - LANES - j;
                let (pv, av, uv) = (block(pdel, j), block(e_at, j), block(e_up, j));
                let (bv, dv) = (block(e_by_o, q), block(e_by_2o, q));
                let mut v = [0f32; LANES];
                let mut exit = [0f32; LANES];
                for l in 0..LANES {
                    let x = LANES - 1 - l;
                    v[l] = one(pv[l], av[l], bv[x], dv[x], uv[l]);
                    exit[l] = pv[l] + t_de_ma + es;
                }
                cdel[j..j + LANES].copy_from_slice(&v);
                take_exits(&exit, &mut best, &mut best_from, first + j);
                j += LANES;
            }
            while j < n {
                let q = n - 1 - j;
                let pr = pdel[j];
                cdel[j] = one(pr, e_at[j], e_by_o[q], e_by_2o[q], e_up[j]);
                take_exits(&[pr + t_de_ma + es], &mut best, &mut best_from, first + j);
                j += 1;
            }
        }

        *best_in = best;
        *best_from_in = best_from;
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

    /// Charge each of the `order` bases that follow an indel
    /// `log2(1 - tp_match_to_zero)`, as the C++ does, instead of the
    /// match-to-match cost. That is cheaper per base than a match state, so
    /// on a perfect array of a long unit of length P, a path at period P - 1
    /// (or k * P - 1) with an insertion each copy outscores the true period.
    /// False by default; `Settings::upstream_compat` turns it on.
    pub upstream_indel_windows: bool,

    /// Fill columns past the warm-up with
    /// [`calculate_steady_column`](Self::calculate_steady_column), which
    /// gives the same scores and tracebacks bit for bit. True by default;
    /// tests turn it off to compare the two paths.
    pub fast_columns: bool,

    /// Emissions by position for the steady-state column.
    position_emissions: PositionEmissions,
    /// Per-order scratch for the steady-state column.
    order_scratch: OrderScratch,
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
            upstream_indel_windows: false,
            fast_columns: true,
            position_emissions: PositionEmissions::default(),
            order_scratch: OrderScratch::default(),
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

    /// Transition cost, in bits, of each of the `order` bases that follow an
    /// indel. Each is a match, so it pays the match-to-match transition as a
    /// base in a match state does. With `upstream_indel_windows` it pays
    /// `log2(1 - tp_match_to_zero)` instead, as in the C++.
    fn indel_window_cost(&self) -> f64 {
        if self.upstream_indel_windows {
            log2d(1.0 - self.tp_match_to_zero as f64)
        } else {
            self.tscore[CT_MATCH as usize][CT_MATCH as usize] as f64
        }
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
        let max_order = self.matrix.max_period.saturating_sub(1);
        if self.fast_columns && max_order > 0 && (nuc_index as usize) < seq.len() {
            // Only the blocked kernel reads emissions by position.
            if self.use_blocked_kernel() {
                let back = self.matrix.steady.max_sequence_back.max(0) as usize;
                self.position_emissions.fill(seq, nuc_index as usize, back, &self.mscore, &self.bscore);
            }
            if self.column_is_steady(seq.len(), nuc_index) {
                self.steady_column(seq, nuc_index as usize, max_order);
                self.matrix.move_matrix_forward();
                return;
            }
            // A column still in its warm-up: the low orders may already be
            // steady. Fill them the fast way and the rest cell by cell, in
            // row order, as the general path would.
            let n = self.steady_orders(seq.len(), nuc_index);
            if n > 0 {
                self.steady_column(seq, nuc_index as usize, n);
                let start_row = self.matrix.match_rows[n + 1];
                self.calculate_general_rows(seq, nuc_index, start_row);
                let cur = self.matrix.current_column_index;
                self.matrix.record_history(cur);
                self.matrix.move_matrix_forward();
                return;
            }
        }
        self.calculate_general_rows(seq, nuc_index, 0);
        let cur = self.matrix.current_column_index;
        self.matrix.record_history(cur);
        self.matrix.move_matrix_forward();
    }

    /// The general path of
    /// [`calculate_current_column`](Self::calculate_current_column) for rows
    /// `start_row..` of the current column. Rows before `start_row` must
    /// already hold their scores for this column. Does not record the
    /// column's history or move the matrix forward.
    fn calculate_general_rows(&mut self, seq: &[Symbol], nuc_index: i32, start_row: usize) {
        let mscore = &self.mscore;
        let bscore = &self.bscore;
        let tscore = &self.tscore;
        let period_decay = self.period_decay;
        let period_decay_offset = self.period_decay_offset;
        let immediate = self.immediate_transition_to_repeat;
        let window_cost = self.indel_window_cost();

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
        debug_assert_eq!(prev + 1, cur);
        let stride = m.length + 1;
        let cols = ColumnSplit::new(&mut m.full_columns, &m.history, desc, cur, stride, m.length);
        let p: &[f32] = cols.p;

        for (row, &d) in desc.iter().enumerate().skip(start_row) {
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
                        // The extensions of an insertion of more than one
                        // base. Skipped for one base: adding t * 0 changes
                        // nothing for a finite t, and with --ii 0 (t = -inf)
                        // it would make the score NaN, so that even
                        // single-base insertions could never be taken.
                        if d.indel_number > 1 {
                            score += tscore[ins][ins] * (d.indel_number - 1) as f32;
                        }
                        // The `order` bases after the indel; see
                        // `indel_window_cost`.
                        score = (score as f64 + window_cost * (order as f32) as f64) as f32;

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
                        // As for insertions (--dd 0).
                        if indel_num > 1 {
                            score += tscore[del][del] * (indel_num - 1) as f32;
                        }
                        // The `order` bases after the indel; see
                        // `indel_window_cost`.
                        score = (score as f64 + window_cost * (order as f32) as f64) as f32;

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

    }

    /// [`calculate_steady_column`](Self::calculate_steady_column) for models
    /// with few indel cells per order (see `use_blocked_kernel`): one pass
    /// over the orders, scalar, which is faster when the orders have fewer
    /// than two full blocks. The same operations in the same order, so the
    /// same result.
    fn calculate_steady_column_fused(&mut self, seq: &[Symbol], i: usize, n_orders: usize) {
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
        let t_bg_bg = tscore[bg][bg];
        let t_ma_bg = tscore[ma][bg];
        let t_bg_ma = tscore[bg][ma];
        let t_ma_ma = tscore[ma][ma];
        let t_in_ma = tscore[ins][ma];
        let t_de_ma = tscore[del][ma];

        let m = &mut self.matrix;
        let cpc = m.cells_per_column;
        let cur = m.current_column_index;
        let tbr = m.max_period + 1;
        let max_period = m.max_period;
        let max_ins = m.max_insertions;
        let max_del = m.max_deletions;

        let ct = &mut m.traceback_matrix[cur * tbr..(cur + 1) * tbr];
        let (slot0, slot1) = m.full_columns.split_at_mut(cpc);
        let (c, p) = if cur & 1 == 0 { (slot0, &*slot1) } else { (slot1, &*slot0) };
        let stride = m.length + 1;
        // Order `o` of column `col` is history[o * stride + col].
        let history = &mut m.history[..];
        // The background (`order` 0) or order-`order` match score of the
        // column `d >= 1` columns back.
        macro_rules! back {
            ($order:expr, $d:expr) => {
                history[($order) * stride + cur - ($d)]
            };
        }

        // diff[x][y] = mscore[x][y] - bscore[y]: the emission score of
        // symbol y against x, the same float as the general path computes.
        let mut diff = [[0f32; 5]; 5];
        for x in 0..5 {
            for y in 0..5 {
                diff[x][y] = mscore[x][y] - bscore[y];
            }
        }
        let s_i = seq[i] as usize;
        let emit = |q: usize| -> f32 { diff[seq[i - q] as usize][s_i] };
        // previous_emission_score(i, q, d): symbol i - d against i - d - q.
        let emit_back = |q: usize, d: usize| -> f32 {
            let at = i - d;
            diff[seq[at - q] as usize][seq[at] as usize]
        };

        c[0] = p[0] + t_bg_bg;
        ct[0] = 0;

        debug_assert!(n_orders < max_period);
        let mut r = 1usize;
        let mut o = 1usize;
        while o <= n_orders {
            let order = o as f32;

            // The match row: first its exit to the background row.
            let mut score = p[r] + t_ma_bg;
            score += order * period_decay;
            score += period_decay_offset;
            if score > c[0] {
                c[0] = score;
                ct[0] = r as i32;
            }

            let zscore = if immediate {
                p[0]
            } else {
                // Order 1 reads the current column's background row live.
                let z = if o == 1 { c[0] } else { back!(0, o - 1) };
                z + t_bg_bg * order
            };
            let es = emit(o);
            let mut best = zscore + t_bg_ma + es;
            let mut best_from = 0i32;
            let stay = p[r] + t_ma_ma + es;
            if stay > best {
                best = stay;
                best_from = r as i32;
            }

            let n_ins = max_ins.min(o - 1);
            let n_del = max_del.min(o - 1);

            // This order's match score o and o + 1 columns back. Insertion k
            // reads it o + k and o + k + 1 columns back, so each iteration
            // loads one new value and reuses the previous one.
            let parent_o = back!(o, o);
            let parent_o1 = back!(o, o + 1);
            let mut nearer = parent_o1;
            for k in 1..=n_ins {
                let row = r + k;
                let per = o + k;
                let pr = p[row];
                let exit = pr + t_in_ma + es;
                if exit > best {
                    best = exit;
                    best_from = row as i32;
                }
                let further = back!(o, per + 1);
                let mut v = pr;
                v -= further;
                v += nearer;
                v -= emit_back(per, o);
                v += emit(per);
                c[row] = v;
                nearer = further;
            }

            for k in 1..=n_del {
                let row = r + n_ins + k;
                let del_order = o - k;
                let e_order = 2 * o - k;
                let pr = p[row];
                let exit = pr + t_de_ma + es;
                if exit > best {
                    best = exit;
                    best_from = row as i32;
                }
                let mut v = pr;
                v -= parent_o1;
                v += parent_o;
                v -= emit_back(del_order, o);
                v += emit_back(del_order, k);
                v -= emit_back(e_order, k);
                v += emit(e_order);
                c[row] = v;
            }

            c[r] = best;
            history[o * stride + cur] = best;
            ct[o] = best_from;
            r += 1 + n_ins + n_del;
            o += 1;
        }
        debug_assert!(n_orders + 1 < max_period || r == cpc);
        history[cur] = c[0];
    }

    /// The steady-state column for orders `1..=n_orders`: the blocked kernel,
    /// or the fused scalar one for models with few indel cells per order.
    #[inline(always)]
    fn steady_column(&mut self, seq: &[Symbol], i: usize, n_orders: usize) {
        if self.use_blocked_kernel() {
            self.calculate_steady_column(seq, i, n_orders);
        } else {
            self.calculate_steady_column_fused(seq, i, n_orders);
        }
    }

    /// Whether the largest orders have at least two full blocks of indel
    /// cells; with fewer, the fused scalar kernel is faster.
    #[inline(always)]
    fn use_blocked_kernel(&self) -> bool {
        uses_blocked_kernel(self.matrix.max_insertions, self.matrix.max_deletions)
    }

    /// How many orders, from order 1 up, have every cell past its warm-up in
    /// the current column (see [`SteadyState::by_order`]). The bounds grow
    /// with the order, so these orders form a prefix.
    ///
    /// [`SteadyState::by_order`]: crate::umatrix::SteadyState::by_order
    fn steady_orders(&self, seq_len: usize, nuc_index: i32) -> usize {
        let m = &self.matrix;
        let prev = m.previous_column_index as i64;
        let cur = m.current_column_index as i64;
        let length = m.length as i64;
        let i = nuc_index as i64;
        if i >= seq_len as i64 {
            return 0;
        }
        let orders = &m.steady.by_order[1..];
        orders.partition_point(|b| {
            prev > b.min_previous_column
                && cur >= b.max_lookback
                && b.max_lookback < length
                && i >= b.max_sequence_back
        })
    }

    /// Whether every cell of the current column is past its warm-up (see
    /// [`SteadyState`](crate::umatrix::SteadyState)): no `matind` branch
    /// applies, every look-back column is at or before the current one
    /// without wrapping, and every sequence index is in range.
    #[inline(always)]
    fn column_is_steady(&self, seq_len: usize, nuc_index: i32) -> bool {
        let m = &self.matrix;
        let st = &m.steady;
        (m.previous_column_index as i64) > st.min_previous_column
            && (m.current_column_index as i64) >= st.max_lookback
            && (st.max_lookback as u64) < m.length as u64
            && (nuc_index as i64) >= st.max_sequence_back
            && (nuc_index as i64) < seq_len as i64
    }

    /// [`calculate_current_column`](Self::calculate_current_column) for a
    /// steady column. The rows are walked in layout order (match row of
    /// order `o`, then its insertions, then its deletions) without looking
    /// up each cell's description, and the look-back columns are addressed
    /// directly. Every cell does the same floating-point operations in the
    /// same order as the general path, and the max into a parent row keeps
    /// the first strictly greater score, so the result is identical.
    ///
    /// Fills the background row and orders `1..=n_orders`, which must all be
    /// steady (`n_orders` is at most `max_period - 1`). With every order this
    /// is the whole column; with fewer, the caller fills the remaining rows
    /// with [`calculate_general_rows`](Self::calculate_general_rows) and
    /// records the history. Does not move the matrix forward.
    fn calculate_steady_column(&mut self, seq: &[Symbol], i: usize, n_orders: usize) {
        let tscore = &self.tscore;
        let pe = &self.position_emissions;
        debug_assert!(pe.filled > i);
        let period_decay = self.period_decay;
        let period_decay_offset = self.period_decay_offset;
        let immediate = self.immediate_transition_to_repeat;

        let bg = CT_BACKGROUND as usize;
        let ma = CT_MATCH as usize;
        let ins = CT_INSERTION as usize;
        let del = CT_DELETION as usize;
        let t_bg_bg = tscore[bg][bg];
        let t_ma_bg = tscore[ma][bg];
        let t_bg_ma = tscore[bg][ma];
        let t_ma_ma = tscore[ma][ma];
        let t_in_ma = tscore[ins][ma];
        let t_de_ma = tscore[del][ma];

        let m = &mut self.matrix;
        let cpc = m.cells_per_column;
        let cur = m.current_column_index;
        let tbr = m.max_period + 1;
        let max_period = m.max_period;
        let max_ins = m.max_insertions;
        let max_del = m.max_deletions;

        let ct = &mut m.traceback_matrix[cur * tbr..(cur + 1) * tbr];
        let (slot0, slot1) = m.full_columns.split_at_mut(cpc);
        let (c, p) = if cur & 1 == 0 { (slot0, &*slot1) } else { (slot1, &*slot0) };
        let stride = m.length + 1;
        // Order `o` of column `col` is history[o * stride + col]. This
        // column's entries are written as its rows finish; only earlier
        // columns are read.
        let history = &mut m.history[..];
        // The background (`order` 0) or order-`order` match score of the
        // column `d >= 1` columns back.
        macro_rules! back {
            ($order:expr, $d:expr) => {
                history[($order) * stride + cur - ($d)]
            };
        }

        // Emissions by position (see PositionEmissions): against(y)[a] is
        // diff[seq[a]][y] and from(x)[a] is diff[x][seq[a]], the same floats
        // as the general path computes.
        let against = |y: Symbol| pe.against(y);
        // emission_score(i, q) = diff[seq[i - q]][seq[i]].
        let emit_i = against(seq[i]);

        // The background row's running max stays in locals and is written
        // once at the end, so the loop over orders does not go through
        // memory for it.
        let mut bg = p[0] + t_bg_bg;
        let mut bg_from = 0i32;

        debug_assert!(n_orders < max_period);
        let match_rows = &m.match_rows[..];
        let (best, best_from, order_es) = self.order_scratch.for_orders(n_orders + 1);

        // Pass 1: every order's match row, from the background row and from
        // itself, and its exit to the background row, in order.
        for o in 1..=n_orders {
            let r = match_rows[o];
            let order = o as f32;
            let mut score = p[r] + t_ma_bg;
            score += order * period_decay;
            score += period_decay_offset;
            if score > bg {
                bg = score;
                bg_from = r as i32;
            }

            let zscore = if immediate {
                p[0]
            } else {
                // Order 1 reads the current column's background row live.
                let z = if o == 1 { bg } else { back!(0, o - 1) };
                z + t_bg_bg * order
            };
            let es = emit_i[i - o];
            let mut b = zscore + t_bg_ma + es;
            let mut from = 0i32;
            let stay = p[r] + t_ma_ma + es;
            if stay > b {
                b = stay;
                from = r as i32;
            }
            best[o] = b;
            best_from[o] = from;
            order_es[o] = es;
        }

        // Pass 2: the insertion and deletion rows, whose exits can still
        // raise their match row (after its own scores, as in the general
        // path). Order 1 has none.
        if max_ins + max_del > 0 {
            for o in 2..=n_orders {
                let n_ins = max_ins.min(o - 1);
                let n_del = max_del.min(o - 1);
                let r = match_rows[o];
                let es = order_es[o];
                let h_o = &history[o * stride..(o + 1) * stride];
                if n_ins >= LANES || n_del >= LANES {
                    let order_cells = SteadyOrder {
                        seq,
                        i,
                        o,
                        r,
                        cur,
                        n_ins,
                        n_del,
                        es,
                        t_in_ma,
                        t_de_ma,
                        h_o,
                        pe,
                    };
                    order_cells.fill(p, c, &mut best[o], &mut best_from[o]);
                    continue;
                }
                // Fewer cells than a block: one at a time, inline.
                let at_o = against(seq[i - o]);
                let mut b = best[o];
                let mut from_row = best_from[o];
                for k in 1..=n_ins {
                    let row = r + k;
                    let pr = p[row];
                    let exit = pr + t_in_ma + es;
                    if exit > b {
                        b = exit;
                        from_row = row as i32;
                    }
                    let mut v = pr;
                    v -= h_o[cur - o - k - 1];
                    v += h_o[cur - o - k];
                    v -= at_o[i - 2 * o - k];
                    v += emit_i[i - o - k];
                    c[row] = v;
                }
                if n_del > 0 {
                    let parent_o = h_o[cur - o];
                    let parent_o1 = h_o[cur - o - 1];
                    let by_o = pe.from(seq[i - o]);
                    let by_2o = pe.from(seq[i - 2 * o]);
                    for k in 1..=n_del {
                        let row = r + n_ins + k;
                        let pr = p[row];
                        let exit = pr + t_de_ma + es;
                        if exit > b {
                            b = exit;
                            from_row = row as i32;
                        }
                        let mut v = pr;
                        v -= parent_o1;
                        v += parent_o;
                        v -= at_o[i - 2 * o + k];
                        v += by_o[i - k];
                        v -= by_2o[i - k];
                        v += emit_i[i - 2 * o + k];
                        c[row] = v;
                    }
                }
                best[o] = b;
                best_from[o] = from_row;
            }
        }

        // Pass 3: the match rows, their history and their tracebacks.
        for o in 1..=n_orders {
            c[match_rows[o]] = best[o];
            history[o * stride + cur] = best[o];
            ct[o] = best_from[o];
        }
        c[0] = bg;
        ct[0] = bg_from;
        history[cur] = bg;
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
        let window_cost = self.indel_window_cost();

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
        debug_assert_eq!(prev + 1, cur);
        let stride = m.length + 1;
        let cols = ColumnSplit::new(&mut m.full_columns, &m.history, desc, cur, stride, m.length);
        let p: &[f32] = cols.p;

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
                        // The extensions of an insertion of more than one
                        // base. Skipped for one base: adding t * 0 changes
                        // nothing for a finite t, and with --ii 0 (t = -inf)
                        // it would make the score NaN, so that even
                        // single-base insertions could never be taken.
                        if d.indel_number > 1 {
                            score += tscore[ins][ins] * (d.indel_number - 1) as f32;
                        }
                        // The `order` bases after the indel; see
                        // `indel_window_cost`.
                        score = (score as f64 + window_cost * (order as f32) as f64) as f32;

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
                        // As for insertions (--dd 0).
                        if indel_num > 1 {
                            score += tscore[del][del] * (indel_num - 1) as f32;
                        }
                        // The `order` bases after the indel; see
                        // `indel_window_cost`.
                        score = (score as f64 + window_cost * (order as f32) as f64) as f32;

                        cols.c[row] = score;
                    } else {
                        cols.c[row] = NEG_INF_F32;
                    }
                }

                _ => {}
            }
        }

        m.record_history(cur);
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
    /// The model of the C++ dumps the values below were taken from.
    fn ultra_model(
        mp: i32,
        ins: i32,
        del: i32,
        mlen: i32,
        adjmax: bool,
        immediate: bool,
    ) -> UModel {
        ultra_model_with(mp, ins, del, mlen, adjmax, immediate, true)
    }

    fn ultra_model_with(
        mp: i32,
        ins: i32,
        del: i32,
        mlen: i32,
        adjmax: bool,
        immediate: bool,
        upstream_indel_windows: bool,
    ) -> UModel {
        let mut m = UModel::new(mp, ins, del, mlen);
        m.upstream_indel_windows = upstream_indel_windows;
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

    /// Random sequence with tandem arrays (substitutions, insertions and
    /// deletions in the copies), lowercase stretches and runs of N.
    fn mixed_sequence(len: usize, max_unit: usize, seed: u64) -> String {
        let mut x = seed;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let bases = b"ACGT";
        let mut out = Vec::with_capacity(len);
        while out.len() < len {
            match next() % 4 {
                0 | 1 => {
                    for _ in 0..(next() % 200) {
                        out.push(bases[(next() % 4) as usize]);
                    }
                }
                2 => {
                    let unit_len = 1 + (next() as usize % max_unit);
                    let unit: Vec<u8> = (0..unit_len).map(|_| bases[(next() % 4) as usize]).collect();
                    let lower = next() % 3 == 0;
                    for _ in 0..(2 + next() % 8) {
                        for &b in &unit {
                            let b = if lower { b.to_ascii_lowercase() } else { b };
                            match next() % 20 {
                                0 => {}                                      // deletion
                                1 => { out.push(b); out.push(bases[(next() % 4) as usize]); } // insertion
                                2 => out.push(bases[(next() % 4) as usize]), // substitution
                                _ => out.push(b),
                            }
                        }
                    }
                }
                _ => {
                    for _ in 0..(1 + next() % 30) {
                        out.push(b'N');
                    }
                }
            }
        }
        out.truncate(len);
        String::from_utf8(out).unwrap()
    }

    /// --ii 0 / --dd 0 forbid extending an insertion or deletion, not making
    /// one: no score is NaN, and the traceback is that of a tiny positive
    /// probability (where an extension is never worth taking either).
    #[test]
    fn zero_extension_probabilities_allow_single_indels() {
        let s = mixed_sequence(800, 13, 21);
        let mlen = s.len() as i32 + 3;
        let run = |p: f32, fast: bool| {
            let mut m = ultra_model_with(14, 4, 4, mlen, false, false, false);
            m.tp_consecutive_insertion = p;
            m.tp_consecutive_deletion = p;
            m.calculate_scores();
            m.fast_columns = fast;
            let mut w = SequenceWindow::new(mlen as usize, 0);
            w.length = s.len();
            for (i, c) in s.bytes().enumerate() {
                w.seq_mem[i] = symbol_for_char(c);
            }
            m.matrix.restart_matrix();
            for i in 0..s.len() as i32 {
                m.calculate_current_column(&w, i, None);
                let col = m.matrix.previous_column_index;
                assert!(m.matrix.score_column(col).iter().all(|v| !v.is_nan()), "NaN in column {col}");
            }
            m.calculate_current_column_without_emission(None);
            let pc = m.matrix.previous_column_index;
            m.matrix.calculate_traceback(pc);
            m.matrix.traceback[..s.len() + 2].to_vec()
        };
        for fast in [false, true] {
            assert_eq!(run(0.0, fast), run(1e-30, fast));
        }
    }

    /// A negative indel count builds no indel cells, so every path must give
    /// the results of a count of 0 (the steady kernels once indexed rows that
    /// the layout does not have).
    #[test]
    fn negative_indel_counts_mean_none() {
        let s = mixed_sequence(500, 11, 7);
        let mlen = s.len() as i32 + 3;
        let run = |ins: i32, del: i32, fast: bool| {
            let mut m = ultra_model_with(12, ins, del, mlen, false, false, false);
            m.fast_columns = fast;
            let mut w = SequenceWindow::new(mlen as usize, 0);
            (run_window(&mut m, &mut w, &s), m.matrix.cells_per_column)
        };
        for (ins, del) in [(-1, 3), (4, i32::MIN), (-7, -7), (i32::MIN, 10)] {
            let want = run(ins.max(0), del.max(0), false);
            assert_eq!(run(ins, del, false), want, "general path, -i {ins} -d {del}");
            assert_eq!(run(ins, del, true), want, "fast paths, -i {ins} -d {del}");
        }
    }

    /// The blocked kernel's emission cache gives what the general path gives
    /// from the live sequence: for a caller that writes each symbol just
    /// before its column, and for a second window written into the same
    /// buffer (as the window pool reuses buffers).
    #[test]
    fn emission_cache_follows_the_live_sequence() {
        let a = mixed_sequence(600, 30, 11);
        let b = mixed_sequence(600, 30, 12);
        let mlen = 603;
        let hash_columns = |m: &mut UModel, w: &SequenceWindow, n: usize, h: &mut u64| {
            for i in 0..n as i32 {
                m.calculate_current_column(w, i, None);
                let x = &m.matrix;
                for v in x.score_column(x.previous_column_index) {
                    *h = fold(*h, v.to_bits() as u64);
                }
            }
        };
        for fast in [false, true] {
            let mut results = Vec::new();
            // Streaming: the window holds only the symbols written so far.
            let mut m = ultra_model_with(31, 10, 10, mlen, false, false, false);
            m.fast_columns = fast;
            let mut w = SequenceWindow::new(mlen as usize, 0);
            w.length = a.len();
            m.matrix.restart_matrix();
            let mut h = 0u64;
            for (i, c) in a.bytes().enumerate() {
                w.seq_mem[i] = symbol_for_char(c);
                m.calculate_current_column(&w, i as i32, None);
                let x = &m.matrix;
                for v in x.score_column(x.previous_column_index) {
                    h = fold(h, v.to_bits() as u64);
                }
            }
            // Then a second window in the same buffer.
            for (i, c) in b.bytes().enumerate() {
                w.seq_mem[i] = symbol_for_char(c);
            }
            m.matrix.restart_matrix();
            hash_columns(&mut m, &w, b.len(), &mut h);
            results.push(h);
            if !fast {
                continue;
            }
            // The same with every symbol written up front, general path.
            let mut g = ultra_model_with(31, 10, 10, mlen, false, false, false);
            g.fast_columns = false;
            let mut w2 = SequenceWindow::new(mlen as usize, 0);
            w2.length = a.len();
            for (i, c) in a.bytes().enumerate() {
                w2.seq_mem[i] = symbol_for_char(c);
            }
            g.matrix.restart_matrix();
            let mut h2 = 0u64;
            hash_columns(&mut g, &w2, a.len(), &mut h2);
            for (i, c) in b.bytes().enumerate() {
                w2.seq_mem[i] = symbol_for_char(c);
            }
            g.matrix.restart_matrix();
            hash_columns(&mut g, &w2, b.len(), &mut h2);
            assert_eq!(results[0], h2, "blocked kernel with a streamed / reused window");
        }
    }

    /// take_exits keeps the first strictly greater exit, as one comparison
    /// per cell does: ties keep the earlier row, and -0.0 does not beat 0.0.
    #[test]
    fn take_exits_keeps_the_first_strictly_greater() {
        let one_at_a_time = |exits: &[f32], mut best: f32, mut from: i32| {
            for (l, &e) in exits.iter().enumerate() {
                if e > best {
                    best = e;
                    from = 10 + l as i32;
                }
            }
            (best.to_bits(), from)
        };
        let cases: [(&[f32], f32); 6] = [
            (&[1.0, 3.0, 3.0, 2.0], 0.0),
            (&[3.0, 1.0, 3.0, 3.0], 0.0),
            (&[2.0, 2.0, 2.0, 2.0], 2.0),
            (&[0.0, -0.0, 0.0, -0.0], -0.0),
            (&[-0.0, 0.0, 1.0, 1.0], -1.0),
            (&[5.0], 4.0),
        ];
        for (exits, start) in cases {
            let (mut best, mut from) = (start, -1);
            take_exits(exits, &mut best, &mut from, 10);
            assert_eq!((best.to_bits(), from), one_at_a_time(exits, start, -1), "{exits:?} from {start}");
        }
    }

    /// The steady-state column gives the same scores and tracebacks, bit for
    /// bit, as the general path, for small and large periods, with indel
    /// counts below and above the order, with many states (13,796 for the
    /// largest shape), and with each model switch.
    #[test]
    fn steady_columns_match_general_path() {
        // (max period, insertions, deletions, sequence length)
        // Shapes on both sides of use_blocked_kernel, with exact blocks and
        // with remainders.
        let shapes = [
            (2, 0, 0, 300),
            (8, 3, 3, 400),
            (60, 3, 2, 500),
            (40, 5, 2, 400),
            (12, 20, 20, 400),
            (40, 4, 4, 400),
            (50, 9, 5, 500),
            (31, 10, 10, 700),
            (101, 10, 10, 1200),
            (300, 40, 8, 1400),
        ];
        for (si, &(mp, ins, del, n)) in shapes.iter().enumerate() {
            let s = mixed_sequence(n, (mp as usize - 1).max(1), 0x9e3779b97f4a7c15 ^ si as u64);
            let mlen = n as i32 + 3;
            for variant in 0..4 {
                let (adjmax, immediate, upstream) = match variant {
                    0 => (false, false, false),
                    1 => (false, false, true),
                    2 => (true, false, false),
                    _ => (false, true, false),
                };
                let mut hashes = Vec::new();
                for fast in [false, true] {
                    let mut m = ultra_model_with(mp, ins, del, mlen, adjmax, immediate, upstream);
                    m.fast_columns = fast;
                    assert!(
                        (m.matrix.steady.min_previous_column as usize) + 10 < n,
                        "shape {si} never reaches a steady column"
                    );
                    let mut w = SequenceWindow::new(mlen as usize, 0);
                    hashes.push((run_window(&mut m, &mut w, &s), m.matrix.traceback.clone()));
                }
                assert!(
                    hashes[0] == hashes[1],
                    "shape {:?} variant {variant}: steady columns differ from the general path",
                    shapes[si]
                );
            }
        }
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

    /// The input of `columns_and_tracebacks_match_cpp` under the default
    /// model, which charges the bases after an indel as matches. The C++
    /// model puts a repeat in columns 88 to 135 (rows 82 to 84 of `want_fwd`
    /// in that test); the default model leaves them background.
    #[test]
    fn default_model_charges_indel_window_bases_as_matches() {
        let mlen = SEQ_A.len() as i32 + 3;
        let mut m = ultra_model_with(11, 10, 10, mlen, false, false, false);
        let mut w = SequenceWindow::new(mlen as usize, 0);

        assert_eq!(run_window(&mut m, &mut w, SEQ_A), 0x93e4fb65ae905247);
        assert_eq!(m.matrix.score_at(140, 0).to_bits(), 0x41b732bf);
        assert_eq!(m.matrix.score_at(140, 7).to_bits(), 0x40c9e268);

        let mut fwd = vec![-7; m.matrix.length + 2];
        assert!(m.matrix.forward_traceback(&mut fwd, SEQ_A.len() as i32, 0));
        assert!(
            fwd[72..136].iter().all(|&row| row == 0),
            "{:?}",
            &fwd[72..136]
        );

        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x29e911bc00278209);
    }

    #[test]
    fn columns_and_tracebacks_match_cpp() {
        let mlen = SEQ_A.len() as i32 + 3;
        let mut m = ultra_model(11, 10, 10, mlen, false, false);
        assert_eq!(m.matrix.cells_per_column, 101);
        let mut w = SequenceWindow::new(mlen as usize, 0);

        assert_eq!(run_window(&mut m, &mut w, SEQ_A), 0x91a9ffebafb895bb);
        assert_eq!(
            (
                m.matrix.current_column_index,
                m.matrix.previous_column_index
            ),
            (141, 140)
        );
        assert_eq!(m.matrix.score_at(140, 0).to_bits(), 0x41bae9fe);
        assert_eq!(m.matrix.score_at(140, 7).to_bits(), 0x40de5c84);
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
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 8 5 5 5 \
             5 5 5 5 5 5 5 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 \
             0 0 0 0 0 0 82 83 82 82 82 82 82 82 82 82 82 82 82 82 83 82 82 82 82 82 82 82 82 82 \
             82 82 82 82 82 82 82 82 82 84 82 82 82 82 82 82 82 82 82 82 82 82 82 82 0 0 0 0 -7 \
             -7 -7 -7",
        );
        assert_eq!(fwd, want_fwd);

        // A second, shorter window reuses the matrix.
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x6b22681375298aa3);
        assert_eq!(m.matrix.traceback[15], 32);

        m.matrix.adjust_score_to_zero(5);
        assert_eq!(m.matrix.last_score_adjustment.to_bits(), 0xc1bbd913);
        assert_eq!(m.matrix.total_score_adjustment.to_bits(), 0xc1bbd913);
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
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x3f3cc269dba4402e);
        assert_eq!(
            run_window(&mut m, &mut w, "ACGNACGACGACG"),
            0x414aa4cda468176a
        );
        // Entries past the second window are stale from the first.
        let trace = ints(
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
             0 2 2 2 2 2 2 2 2 2 2 2 0 0 0",
        );
        assert_eq!(m.matrix.traceback, trace);

        let mut m = ultra_model(9, 3, 2, mlen, false, true);
        assert_eq!(run_window(&mut m, &mut w, SEQ_B), 0x6ebe7c72c1de7bb5);
    }
}
