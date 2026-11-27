//! Wraparound alignment of a tandem repeat to its consensus unit.
//!
//! The unit is treated as a circle, column 0 following the last column, so
//! one alignment runs through every copy of the repeat and may start in any
//! column (the phase). Each base is aligned to a column, or is an insertion
//! (a base with no column); a column can be skipped (a deletion: a column
//! with no base). The alignment is global in the repeat and free at both
//! ends in the unit: it starts and ends in whatever column scores best.
//!
//! [`RepeatRegion::realign_logo`](crate::repeat::RepeatRegion) builds the
//! consensus with this alignment, and the same alignment against the final
//! consensus describes the repeat in alignment output, so the two cannot
//! disagree.
//!
//! The traceback takes one byte per base and column. Above
//! [`SINGLE_PASS_CELLS`] cells the alignment keeps the score row at the start
//! of every block of rows and recomputes each block's traceback from it, so
//! memory stays bounded for long arrays; the result is the same alignment.

use crate::symbol::{symbol_for_char, Symbol, N_UNK};

/// Scores for [`align`]. An unknown base (`N_UNK`) scores 0 in every column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scores {
    pub matched: i32,
    pub mismatched: i32,
    /// For each inserted base and each skipped column.
    pub gap: i32,
}

/// The scores ULTRA builds its consensus with. In a simulation of 5,200
/// noisy tandem arrays, these gave a consensus that was a rotation of the
/// planted unit in 2,798; eight other settings gave between 2,791 and 2,798.
pub const ULTRA_SCORES: Scores = Scores {
    matched: 2,
    mismatched: -2,
    gap: -5,
};

/// Above this many cells (bases × unit length) the traceback is recomputed in
/// blocks instead of being kept whole (16 MiB of traceback at most).
pub const SINGLE_PASS_CELLS: usize = 1 << 24;

/// One kind of step along the alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// A base aligned to a column (a match or a mismatch).
    Aligned,
    /// A base with no column.
    Inserted,
    /// A column with no base.
    Skipped,
}

/// A wraparound alignment: the column it starts in and its steps as runs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    /// Length of the unit.
    pub period: usize,
    /// The first column the alignment uses (0-based), or 0 when it uses none.
    pub phase: usize,
    /// Runs of steps in sequence order; adjacent runs differ in [`Op`].
    pub runs: Vec<(Op, u64)>,
    pub score: i64,
}

impl Alignment {
    /// Bases the alignment covers (aligned and inserted).
    pub fn bases(&self) -> u64 {
        self.count(|op| op != Op::Skipped)
    }

    /// Columns the alignment uses (aligned and skipped), counting each column
    /// again on every pass around the unit.
    pub fn columns(&self) -> u64 {
        self.count(|op| op != Op::Inserted)
    }

    fn count(&self, keep: impl Fn(Op) -> bool) -> u64 {
        self.runs.iter().filter(|r| keep(r.0)).map(|r| r.1).sum()
    }

    /// Calls `f(base, column)` for every aligned base, in order. `base` is the
    /// offset in the aligned sequence.
    pub fn for_each_aligned(&self, mut f: impl FnMut(usize, usize)) {
        if self.period == 0 {
            return;
        }
        let (mut base, mut column) = (0usize, self.phase);
        for &(op, len) in &self.runs {
            let len = len as usize;
            match op {
                Op::Aligned => {
                    for _ in 0..len {
                        f(base, column);
                        base += 1;
                        column = (column + 1) % self.period;
                    }
                }
                Op::Inserted => base += len,
                Op::Skipped => column = (column + len) % self.period,
            }
        }
    }

    /// The steps as a CIGAR string with the unit as the reference, as ULTRA's
    /// traceback writes them: `M` aligned, `I` inserted base, `D` skipped
    /// column.
    pub fn cigar(&self) -> String {
        let mut s = String::new();
        for &(op, len) in &self.runs {
            let c = match op {
                Op::Aligned => 'M',
                Op::Inserted => 'I',
                Op::Skipped => 'D',
            };
            s.push_str(&len.to_string());
            s.push(c);
        }
        s
    }
}

const DIAGONAL: u8 = 0;
const INSERTED_BASE: u8 = 1;
const SKIPPED_COLUMN: u8 = 2;

/// Aligns `seq` to the circular `unit`, both as [`Symbol`]s. Columns are
/// compared by symbol, so a unit column that is `N_UNK` matches nothing.
///
/// Ties: a base is aligned rather than inserted, a column is skipped only for
/// a strictly better score, and the alignment ends in the lowest best column.
pub fn align(seq: &[Symbol], unit: &[Symbol], scores: Scores) -> Alignment {
    align_with_limit(seq, unit, scores, SINGLE_PASS_CELLS)
}

/// [`align`] for bytes; anything other than `ACGTacgt` is unknown.
pub fn align_bytes(seq: &[u8], unit: &[u8], scores: Scores) -> Alignment {
    let seq: Vec<Symbol> = seq.iter().map(|&c| symbol_for_char(c)).collect();
    let unit: Vec<Symbol> = unit.iter().map(|&c| symbol_for_char(c)).collect();
    align(&seq, &unit, scores)
}

/// A stretch of a repeat with one consensus, aligned to it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Part {
    /// The stretch `[start, end)` of the repeat, 0-based offsets into it.
    pub start: usize,
    pub end: usize,
    pub consensus: Vec<u8>,
    pub alignment: Alignment,
}

/// Cuts the repeat `seq` into parts and aligns each to its consensus with
/// [`ULTRA_SCORES`]. `parts` gives each part's start (an offset into `seq`)
/// and consensus, in order; a part ends where the next one starts, the first
/// starts at 0 whatever its offset, and an empty part is dropped. A start
/// before an earlier one is taken as that earlier start, so parts never
/// overlap. With no parts, nothing is aligned.
pub fn align_parts(seq: &[u8], parts: &[(usize, &[u8])]) -> Vec<Part> {
    let mut bounds = Vec::with_capacity(parts.len() + 1);
    let mut at = 0usize;
    for (k, &(start, _)) in parts.iter().enumerate() {
        if k > 0 {
            at = at.max(start.min(seq.len()));
        }
        bounds.push(at);
    }
    bounds.push(seq.len());

    let mut out = Vec::with_capacity(parts.len());
    for (k, &(_, consensus)) in parts.iter().enumerate() {
        let (start, end) = (bounds[k], bounds[k + 1]);
        if end <= start {
            continue;
        }
        out.push(Part {
            start,
            end,
            consensus: consensus.to_vec(),
            alignment: align_bytes(&seq[start..end], consensus, ULTRA_SCORES),
        });
    }
    out
}

/// [`align`] with the cell count above which the traceback is recomputed in
/// blocks.
pub fn align_with_limit(
    seq: &[Symbol],
    unit: &[Symbol],
    scores: Scores,
    single_pass_cells: usize,
) -> Alignment {
    let period = unit.len();
    let n = seq.len();
    if period == 0 || n == 0 {
        return Alignment {
            period,
            phase: 0,
            runs: Vec::new(),
            score: 0,
        };
    }

    // Rows per block of traceback, and the score row before each block.
    let block = (single_pass_cells / period).max(1).min(n);
    let n_blocks = n.div_ceil(block);
    let mut checkpoints: Vec<Vec<i32>> = Vec::with_capacity(n_blocks);

    let mut prev = vec![0i32; period];
    let mut cur = vec![0i32; period];
    let mut moves = vec![DIAGONAL; block * period];

    // Forward pass. With one block its traceback is kept from here.
    let keep = n_blocks == 1;
    for rows in seq.chunks(block) {
        checkpoints.push(prev.clone());
        for (r, &sym) in rows.iter().enumerate() {
            let row = if keep { r } else { 0 };
            fill_row(sym, unit, scores, &prev, &mut cur, &mut moves[row * period..(row + 1) * period]);
            std::mem::swap(&mut prev, &mut cur);
        }
    }

    // Best column for the last base (lowest column on ties).
    let mut j = 0;
    for t in 1..period {
        if prev[t] > prev[j] {
            j = t;
        }
    }
    let score = prev[j] as i64;

    // Traceback, last block first, recomputing each block's moves when they
    // were not kept.
    let mut rev: Vec<(Op, u64)> = Vec::new();
    let mut push = |op: Op| match rev.last_mut() {
        Some(last) if last.0 == op => last.1 += 1,
        _ => rev.push((op, 1)),
    };
    let mut first_column = None;
    let mut i = n;
    'blocks: for b in (0..n_blocks).rev() {
        let lo = b * block;
        if n_blocks > 1 {
            let mut p = checkpoints[b].clone();
            for (r, &sym) in seq[lo..(lo + block).min(n)].iter().enumerate() {
                fill_row(sym, unit, scores, &p, &mut cur, &mut moves[r * period..(r + 1) * period]);
                std::mem::swap(&mut p, &mut cur);
            }
        }
        let mut skipped_in_a_row = 0;
        while i > lo {
            match moves[(i - 1 - lo) * period + j] {
                DIAGONAL => {
                    push(Op::Aligned);
                    first_column = Some(j);
                    j = (j + period - 1) % period;
                    i -= 1;
                    skipped_in_a_row = 0;
                }
                INSERTED_BASE => {
                    push(Op::Inserted);
                    i -= 1;
                    skipped_in_a_row = 0;
                }
                _ => {
                    push(Op::Skipped);
                    first_column = Some(j);
                    j = (j + period - 1) % period;
                    skipped_in_a_row += 1;
                    // Skipping every column gains nothing, so the alignment
                    // never does it; this only guards the loop.
                    if skipped_in_a_row > period {
                        break 'blocks;
                    }
                }
            }
        }
    }
    rev.reverse();

    Alignment {
        period,
        phase: first_column.unwrap_or(0),
        runs: rev,
        score,
    }
}

/// One row of the recurrence: scores into `cur` and moves into `moves`, from
/// the previous row `prev`.
#[inline]
fn fill_row(sym: Symbol, unit: &[Symbol], scores: Scores, prev: &[i32], cur: &mut [i32], moves: &mut [u8]) {
    let period = unit.len();
    for j in 0..period {
        let before = if j == 0 { period - 1 } else { j - 1 };
        let diagonal = prev[before]
            + if sym == N_UNK {
                0
            } else if sym == unit[j] {
                scores.matched
            } else {
                scores.mismatched
            };
        let inserted = prev[j] + scores.gap;
        if diagonal >= inserted {
            cur[j] = diagonal;
            moves[j] = DIAGONAL;
        } else {
            cur[j] = inserted;
            moves[j] = INSERTED_BASE;
        }
    }
    // A run of skipped columns may pass from the last column to column 0, so
    // the row is swept twice.
    for _ in 0..2 {
        for j in 0..period {
            let before = if j == 0 { period - 1 } else { j - 1 };
            let skipped = cur[before] + scores.gap;
            if skipped > cur[j] {
                cur[j] = skipped;
                moves[j] = SKIPPED_COLUMN;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(seq: &str, unit: &str) -> Alignment {
        align_bytes(seq.as_bytes(), unit.as_bytes(), ULTRA_SCORES)
    }

    /// Score of `a` recomputed from its steps.
    fn rescore(a: &Alignment, seq: &[u8], unit: &[u8], s: Scores) -> i64 {
        let mut total = 0i64;
        let (mut base, mut col) = (0usize, a.phase);
        for &(op, len) in &a.runs {
            for _ in 0..len {
                match op {
                    Op::Aligned => {
                        let x = symbol_for_char(seq[base]);
                        total += if x == N_UNK {
                            0
                        } else if x == symbol_for_char(unit[col]) {
                            s.matched
                        } else {
                            s.mismatched
                        } as i64;
                        base += 1;
                        col = (col + 1) % unit.len();
                    }
                    Op::Inserted => {
                        total += s.gap as i64;
                        base += 1;
                    }
                    Op::Skipped => {
                        total += s.gap as i64;
                        col = (col + 1) % unit.len();
                    }
                }
            }
        }
        assert_eq!(base, seq.len(), "path covers every base");
        total
    }

    #[test]
    fn perfect_repeat_starts_in_its_phase() {
        let a = aln("CAGACAGACAGA", "ACAG");
        assert_eq!((a.phase, a.cigar(), a.columns()), (1, "12M".to_string(), 12));
        assert_eq!(a.score, 24);
    }

    #[test]
    fn insertion_and_deletion() {
        // Copy 2 has an extra T, copy 4 lacks its G.
        let a = aln("ACGTTGCA ACGTTTGCA ACGTTGCA ACGTTCA ACGTTGCA".replace(' ', "").as_str(), "ACGTTGCA");
        assert_eq!(a.phase, 0);
        assert_eq!(a.bases(), 40);
        assert_eq!(a.columns(), 40);
        let ins: u64 = a.runs.iter().filter(|r| r.0 == Op::Inserted).map(|r| r.1).sum();
        let del: u64 = a.runs.iter().filter(|r| r.0 == Op::Skipped).map(|r| r.1).sum();
        assert_eq!((ins, del), (1, 1), "{}", a.cigar());
    }

    #[test]
    fn deletion_runs_across_the_end_of_the_unit() {
        // Copy 2 lacks the last base of one copy and the first of the next.
        let a = aln("ACGTACGTACGACGTACGTACGT", "ACGTACGT");
        assert_eq!(a.bases(), 23);
        let s = rescore(&a, b"ACGTACGTACGACGTACGTACGT", b"ACGTACGT", ULTRA_SCORES);
        assert_eq!(s, a.score);
    }

    #[test]
    fn unknown_bases_align_for_free() {
        let a = aln("ACGNNNACGACG", "ACG");
        assert_eq!((a.cigar(), a.score), ("12M".to_string(), 18));
    }

    #[test]
    fn empty_inputs() {
        assert_eq!(aln("", "ACG").runs, vec![]);
        assert_eq!(aln("ACG", "").runs, vec![]);
    }

    /// Deterministic noisy arrays: path score equals the DP score, the path
    /// covers every base, and recomputing the traceback in blocks of any size
    /// gives the same alignment.
    #[test]
    fn blocks_give_the_same_alignment() {
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut rnd = |m: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % m
        };
        for case in 0..60 {
            let p = 1 + rnd(23) as usize;
            let unit: Vec<u8> = (0..p).map(|_| b"ACGT"[rnd(4) as usize]).collect();
            let mut seq = Vec::new();
            let copies = 2 + rnd(30);
            for _ in 0..copies {
                for &c in &unit {
                    match rnd(20) {
                        0 => {}                                          // deletion
                        1 => seq.extend([c, b"ACGTN"[rnd(5) as usize]]), // insertion
                        2 => seq.push(b"ACGTN"[rnd(5) as usize]),        // substitution
                        _ => seq.push(c),
                    }
                }
            }
            let s: Vec<Symbol> = seq.iter().map(|&c| symbol_for_char(c)).collect();
            let u: Vec<Symbol> = unit.iter().map(|&c| symbol_for_char(c)).collect();
            let whole = align(&s, &u, ULTRA_SCORES);
            assert_eq!(rescore(&whole, &seq, &unit, ULTRA_SCORES), whole.score, "case {case}");
            assert!(whole.phase < p);
            for limit in [1, p, 3 * p, 7 * p + 1, s.len() * p - 1] {
                let blocked = align_with_limit(&s, &u, ULTRA_SCORES, limit.max(1));
                assert_eq!(blocked, whole, "case {case} limit {limit}");
            }
        }
    }

    #[test]
    fn parts_cut_at_their_starts() {
        let seq = b"ACACACACACAGTAGTAGTAGT";
        let p = align_parts(seq, &[(0, b"AC"), (10, b"AGT"), (10, b"GG")]);
        assert_eq!(p.len(), 2, "the empty part is dropped");
        assert_eq!((p[0].start, p[0].end, p[0].alignment.cigar()), (0, 10, "10M".to_string()));
        assert_eq!((p[1].start, p[1].end, p[1].consensus.as_slice()), (10, 22, &b"GG"[..]));
        assert!(align_parts(seq, &[]).is_empty());

        // A start that goes backwards is held at the earlier start.
        let p = align_parts(seq, &[(0, b"AC"), (12, b"AGT"), (4, b"AC")]);
        let spans: Vec<_> = p.iter().map(|x| (x.start, x.end)).collect();
        assert_eq!(spans, [(0, 12), (12, 22)]);
    }

    #[test]
    fn for_each_aligned_visits_columns_in_order() {
        let a = aln("CAGTACAGACAG", "ACAG");
        let mut seen = Vec::new();
        a.for_each_aligned(|b, c| seen.push((b, c)));
        assert_eq!(seen.first(), Some(&(0, 1)));
        assert_eq!(seen.len() as u64, a.bases() - a.runs.iter().filter(|r| r.0 == Op::Inserted).map(|r| r.1).sum::<u64>());
    }
}
