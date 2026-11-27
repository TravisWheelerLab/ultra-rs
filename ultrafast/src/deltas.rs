//! Per-base score deltas for a tandem array, for Aurora's tandem-repeat rows
//! (the role of ULTRA's `--json --show_delta` `PositionScoreDeltas`).
//!
//! Both models work from the same path: the array aligned base by base to its
//! cyclic consensus (unit-cost wraparound alignment, as in `consensus.rs`).
//!
//! * `Ultra`: ULTRA 1.2.1's default model, in bits. Each base is compared with
//!   the base aligned to the same consensus column in the previous copy (ULTRA
//!   compares with the base one period back; following the column keeps that
//!   comparison in phase across indels). Constants as measured from ULTRA's own
//!   deltas: first copy -0.0145/base (background), first compared base -5.187
//!   (entering the repeat), then match +1.349 / mismatch -1.457 (log2(0.7/0.25),
//!   log2(0.1/0.25), plus log2(0.91) to stay in the repeat). An unknown base
//!   scores log2(0.91). Inserted bases and skipped consensus columns cost
//!   log2(0.02) = -5.644 each (ULTRA's indel open and extend probabilities; its
//!   own deltas around indels are path artefacts and are not reproduced).
//!   As in ULTRA, a substitution costs twice: against the copy before it and
//!   the copy after.
//! * `Te`: the same per-column score Aurora gives a TE alignment, in nats: for
//!   genome base t aligned to consensus base q, ln(ratio[t][q] / f[t]), with
//!   ratio[t][q] = exp(lambda*S[q][t]) * fj[t] and f the genome's base
//!   frequencies in a 61-bp window centred on the base (Aurora's
//!   `score_with_background`, `windowed_scores.rs`); a non-ACGT genome base
//!   scores lambda*S[q][N]. A run of k skipped consensus columns (a gap in the
//!   genome) costs lambda*(open + (k-1)*ext) on the genome base before it; a run
//!   of k inserted genome bases costs lambda*(open + (k-1)*ext)/k on each.
//!
//! In both models the deltas sum to the reported score.

use crate::consensus::code;
use anyhow::{anyhow, bail, Context, Result};

pub use crate::wrap::Op;

fn align_window_path(s: &[u8], u: &[u8], cm: i32, ci: i32, out: &mut Vec<Op>) {
    crate::wrap::align(s, u, cm, ci, &mut crate::wrap::Scratch::new(), out);
}

/// Alignment path of `array` to the cyclic `unit` (both ASCII), minimising
/// cm per substitution + ci per inserted base or skipped column; pick cm:ci to
/// match the scoring model (`ULTRA_COSTS`, `TE_COSTS`) so the path is the one
/// that model would prefer. Every array base appears exactly once, as `M` or `I`,
/// and the columns of successive `M`/`D` steps are consecutive around the unit.
///
/// Long arrays are aligned in overlapping windows, each with a free starting
/// phase, and stitched at a base both windows put in the same column (`stitch`).
pub fn path(array: &[u8], unit: &[u8], max_cells: usize, (cm, ci): (i32, i32)) -> Vec<Op> {
    let s: Vec<u8> = array.iter().map(|&c| code(c)).collect();
    let u: Vec<u8> = unit.iter().map(|&c| code(c)).filter(|&c| c < 4).collect();
    let mut out = Vec::with_capacity(s.len() + s.len() / 8);
    if u.is_empty() {
        out.extend(std::iter::repeat(Op::I).take(s.len()));
        return out;
    }
    // windows of at least three units, up to 200 kb; but a unit so long that three of
    // them pass max_cells (above ~9.5 kb by default) gets windows of max_cells cells,
    // shorter than three units, so memory stays bounded (2 bits per cell)
    let l = u.len();
    let w = if 3 * l * l <= max_cells { (max_cells / l).clamp(3 * l, 200_000.max(3 * l)) } else { (max_cells / l).max(1) };
    let w = w.min(s.len().max(1));
    if s.len() <= w {
        align_window_path(&s, &u, cm, ci, &mut out);
        return out;
    }
    // each window runs `ov` bases into the next one; a last piece no longer than
    // that is covered by the window before it
    let ov = l.max(64).min(w / 4);
    let mut spans: Vec<(usize, usize)> = Vec::new();
    let mut b = 0;
    while spans.last().is_none_or(|&(_, e)| e < s.len()) {
        spans.push((b, (b + w + ov).min(s.len())));
        b += w;
    }
    // windows are independent: align them in parallel, then stitch in order
    use rayon::prelude::*;
    let parts: Vec<Vec<Op>> = spans
        .par_iter()
        .with_max_len(1)
        .map(|&(b, e)| {
            let mut o = Vec::with_capacity(e - b + (e - b) / 8);
            align_window_path(&s[b..e], &u, cm, ci, &mut o);
            o
        })
        .collect();
    let mut end = 0;
    for (p, &(b, e)) in parts.iter().zip(&spans) {
        if b == 0 {
            out.extend_from_slice(p);
        } else {
            stitch(&mut out, p, end - b, l);
        }
        end = e;
    }
    out
}

/// Append window `next` to `out`, whose last `ov` bases are also the first `ov`
/// bases of `next`. The cut goes after the overlap base nearest the middle that both
/// put in the same column, so the path stays one alignment. If they agree nowhere, cut
/// in the middle and pay for the phase change: skip columns, or turn the last aligned
/// bases before the cut into insertions, whichever moves fewer columns.
fn stitch(out: &mut Vec<Op>, next: &[Op], ov: usize, l: usize) {
    let is_base = |o: &Op| !matches!(o, Op::D(_));
    // op index of each overlap base, in `out` and in `next`
    let mut ia = Vec::with_capacity(ov);
    let mut i = out.len();
    while ia.len() < ov {
        i -= 1;
        if is_base(&out[i]) {
            ia.push(i);
        }
    }
    ia.reverse();
    let ib: Vec<usize> = next.iter().enumerate().filter(|(_, o)| is_base(o)).map(|(i, _)| i).take(ov).collect();
    let mid = ov / 2;
    let agree = (0..ov)
        .filter(|&k| matches!(out[ia[k]], Op::M(_)) && out[ia[k]] == next[ib[k]])
        .min_by_key(|&k| k.abs_diff(mid));
    if let Some(k) = agree {
        out.truncate(ia[k] + 1);
        out.extend_from_slice(&next[ib[k] + 1..]);
        return;
    }
    let tail = if mid < ov { &next[ib[mid]..] } else { &next[next.len()..] };
    if mid < ov {
        out.truncate(ia[mid]);
    }
    let col = |o: &Op| match o {
        Op::M(c) | Op::D(c) => Some(*c as usize),
        Op::I => None,
    };
    if let (Some(prev), Some(to)) = (out.iter().rev().find_map(col), tail.iter().find_map(col)) {
        let from = (prev + 1) % l;
        let fwd = (to + l - from) % l;
        if fwd <= l / 2 {
            out.extend((0..fwd).map(|d| Op::D(((from + d) % l) as u32)));
        } else {
            // step back l - fwd columns: aligned bases become insertions, skipped
            // columns are no longer skipped
            let (mut back, mut ins) = (l - fwd, 0);
            while back > 0 {
                match out.pop() {
                    Some(Op::M(_)) => {
                        back -= 1;
                        ins += 1;
                    }
                    Some(Op::D(_)) => back -= 1,
                    Some(Op::I) => ins += 1,
                    None => break,
                }
            }
            out.extend(std::iter::repeat_n(Op::I, ins));
        }
    }
    out.extend_from_slice(tail);
}

// ---------------------------------------------------------------- ULTRA model

/// Path costs for the ULTRA model: a substitution loses MATCH - MISMATCH = 2.8
/// bits, an indel about 5.6 + 1.35 = 7.0 (gap plus a lost match): ~2:5.
pub const ULTRA_COSTS: (i32, i32) = (2, 5);
/// Path costs for the TE model (RepeatMasker matrices: a substitution ~25 below
/// a match, gap open 35 plus a lost match ~45): ~5:9.
pub const TE_COSTS: (i32, i32) = (5, 9);

pub const U_BG: f64 = -0.0145;
pub const U_ENTRY: f64 = -5.187428;
pub const U_MATCH: f64 = 1.349365; // log2(0.7/0.25) + log2(0.91)
pub const U_MISMATCH: f64 = -1.457990; // log2(0.1/0.25) + log2(0.91)
pub const U_STAY: f64 = -0.136062; // log2(0.91): unknown base inside the repeat
pub const U_INDEL: f64 = -5.643856; // log2(0.02)

/// ULTRA-style deltas (bits), one per array base.
pub fn ultra_deltas(array: &[u8], unit_len: usize, ops: &[Op]) -> Vec<f64> {
    let mut out = Vec::with_capacity(array.len());
    let mut last: Vec<u8> = vec![255; unit_len.max(1)];
    let mut entered = false;
    let mut pending = 0.0; // deletions before the next base
    let mut i = 0;
    for op in ops {
        match *op {
            Op::D(_) => {
                if entered {
                    pending += U_INDEL;
                }
            }
            Op::I => {
                out.push(if entered { U_INDEL } else { U_BG } + pending);
                pending = 0.0;
                i += 1;
            }
            Op::M(col) => {
                let col = col as usize;
                let a = code(array[i]);
                let prev = last[col];
                let v = if prev == 255 && !entered {
                    U_BG
                } else if !entered {
                    entered = true;
                    U_ENTRY
                } else if a == 4 || prev == 4 || prev == 255 {
                    U_STAY
                } else if a == prev {
                    U_MATCH
                } else {
                    U_MISMATCH
                };
                out.push(v + pending);
                pending = 0.0;
                last[col] = a;
                i += 1;
            }
        }
    }
    if pending != 0.0 {
        if let Some(x) = out.last_mut() {
            *x += pending;
        }
    }
    out
}

// ------------------------------------------------------------------- TE model

/// One Aurora-format substitution matrix (`#matrix`, `#gap-open`, `#gap-ext`,
/// `#lambda`, `#fi`, `#fj`, then a header row of letters and the score rows;
/// rows = consensus base, columns = genome base, as in RepeatMasker's
/// cross_match matrices).
#[derive(Clone, Debug)]
pub struct Matrix {
    pub name: String,
    pub lambda: f64,
    pub gap_open: f64,
    pub gap_ext: f64,
    /// ln(ratio[t][q]) for genome t, consensus q in A,C,G,T.
    ln_ratio: [[f64; 4]; 4],
    /// lambda * S[q][N] for an unknown genome base.
    n_score: [f64; 4],
}

pub fn read_matrices(path: &str) -> Result<Vec<Matrix>> {
    let txt = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut out = Vec::new();
    let mut lines = txt.lines().map(str::trim).filter(|l| !l.is_empty()).peekable();
    while let Some(l) = lines.next() {
        let Some(name) = l.strip_prefix("#matrix") else { continue };
        let name = name.trim().trim_end_matches(".matrix").to_string();
        let (mut lambda, mut go, mut ge, mut fj) = (0.0, 0.0, 0.0, [0.25; 4]);
        while let Some(h) = lines.peek().filter(|h| h.starts_with('#')).copied() {
            lines.next();
            let t: Vec<&str> = h.split_whitespace().collect();
            let f = |i: usize| -> Result<f64> { t.get(i).ok_or_else(|| anyhow!("bad line {h}"))?.parse().map_err(Into::into) };
            match t[0] {
                "#lambda" => lambda = f(1)?,
                "#gap-open" => go = f(1)?,
                "#gap-ext" => ge = f(1)?,
                "#fj" => fj = [f(1)?, f(2)?, f(3)?, f(4)?],
                _ => {}
            }
        }
        let header: Vec<u8> = lines.next().ok_or_else(|| anyhow!("{name}: no letters"))?
            .split_whitespace().map(|c| c.as_bytes()[0]).collect();
        let mut s = vec![vec![0.0f64; header.len()]; header.len()];
        for row in s.iter_mut() {
            let r: Vec<f64> = lines.next().ok_or_else(|| anyhow!("{name}: short matrix"))?
                .split_whitespace().map(|v| v.parse()).collect::<Result<_, _>>()?;
            if r.len() != header.len() {
                bail!("{name}: row length");
            }
            row.copy_from_slice(&r);
        }
        let idx = |c: u8| header.iter().position(|&h| h == c).ok_or_else(|| anyhow!("{name}: no {}", c as char));
        let acgt = [idx(b'A')?, idx(b'C')?, idx(b'G')?, idx(b'T')?];
        let ni = idx(b'N')?;
        let mut ln_ratio = [[0.0; 4]; 4];
        let mut n_score = [0.0; 4];
        for q in 0..4 {
            for t in 0..4 {
                ln_ratio[t][q] = lambda * s[acgt[q]][acgt[t]] + fj[t].ln();
            }
            n_score[q] = lambda * s[acgt[q]][ni];
        }
        out.push(Matrix { name, lambda, gap_open: go, gap_ext: ge, ln_ratio, n_score });
    }
    if out.is_empty() {
        bail!("{path}: no matrices");
    }
    Ok(out)
}

/// Pick the matrix whose name `<div>p<gc>g` is nearest to the array's
/// divergence (percent) and the sequence's GC (percent); `div` is weighted
/// first. Names that don't follow the pattern are never picked automatically.
pub fn pick<'a>(ms: &'a [Matrix], div_pct: f64, gc_pct: f64) -> Option<&'a Matrix> {
    let parse = |n: &str| -> Option<(f64, f64)> {
        let (d, g) = n.split_once('p')?;
        Some((d.parse().ok()?, g.strip_suffix('g')?.parse().ok()?))
    };
    ms.iter()
        .filter_map(|m| parse(&m.name).map(|(d, g)| (m, d, g)))
        .min_by(|a, b| {
            let ka = ((a.1 - div_pct).abs(), (a.2 - gc_pct).abs());
            let kb = ((b.1 - div_pct).abs(), (b.2 - gc_pct).abs());
            ka.partial_cmp(&kb).unwrap()
        })
        .map(|x| x.0)
}

/// Base frequencies in a `win`-bp window centred on each position of
/// seq[beg..end] (window clipped to the sequence; unknown bases count 1/4 each,
/// as in Aurora's NUCLEOTIDE_WEIGHTS for N).
pub fn local_freqs(seq: &[u8], beg: usize, end: usize, win: usize) -> Vec<[f64; 4]> {
    let h = (win - 1) / 2;
    let add = |c: &mut [f64; 4], b: u8, w: f64| match code(b) {
        4 => c.iter_mut().for_each(|x| *x += 0.25 * w),
        k => c[k as usize] += w,
    };
    let mut cnt = [0.0; 4];
    let (mut lo, mut hi) = (beg.saturating_sub(h), beg.saturating_sub(h)); // window [lo, hi)
    let mut out = Vec::with_capacity(end - beg);
    for p in beg..end {
        let (nl, nh) = (p.saturating_sub(h), (p + h + 1).min(seq.len()));
        while hi < nh {
            add(&mut cnt, seq[hi], 1.0);
            hi += 1;
        }
        while lo < nl {
            add(&mut cnt, seq[lo], -1.0);
            lo += 1;
        }
        let n = (hi - lo) as f64;
        out.push([cnt[0] / n, cnt[1] / n, cnt[2] / n, cnt[3] / n]);
    }
    out
}

/// TE-scale deltas (nats), one per array base. `seq` is the whole sequence and
/// the array is seq[beg..beg+len].
pub fn te_deltas(seq: &[u8], beg: usize, unit: &[u8], ops: &[Op], m: &Matrix, win: usize) -> Vec<f64> {
    let n = ops.iter().filter(|o| !matches!(o, Op::D(_))).count();
    let f = local_freqs(seq, beg, beg + n, win);
    let u: Vec<u8> = unit.iter().map(|&c| code(c)).filter(|&c| c < 4).collect();
    let gap = |k: usize| m.lambda * (m.gap_open + (k as f64 - 1.0) * m.gap_ext);
    let mut out = vec![0.0; n];
    let (mut i, mut j) = (0usize, 0usize);
    while j < ops.len() {
        match ops[j] {
            Op::M(col) => {
                let t = code(seq[beg + i]);
                let q = u[col as usize] as usize;
                out[i] += if t == 4 { m.n_score[q] } else { m.ln_ratio[t as usize][q] - f[i][t as usize].max(1e-9).ln() };
                i += 1;
                j += 1;
            }
            Op::I => {
                let s = j;
                while j < ops.len() && ops[j] == Op::I {
                    j += 1;
                }
                let k = j - s;
                for x in 0..k {
                    out[i + x] += gap(k) / k as f64;
                }
                i += k;
            }
            Op::D(_) => {
                let s = j;
                while j < ops.len() && matches!(ops[j], Op::D(_)) {
                    j += 1;
                }
                let k = j - s;
                let at = if i > 0 { i - 1 } else { 0 };
                if n > 0 {
                    out[at] += gap(k);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_windows_for_units_too_long_for_three_per_window() {
        // three units would pass max_cells: windows shorter than a unit, no panic, and
        // every array base appears once
        let unit: Vec<u8> = (0..200).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect();
        let arr = [unit.clone(), unit.clone(), unit[..50].to_vec()].concat();
        for max_cells in [1_000, 30_000, 1 << 28] {
            let ops = path(&arr, &unit, max_cells, ULTRA_COSTS);
            assert_eq!(ops.iter().filter(|o| !matches!(o, Op::D(_))).count(), arr.len(), "max_cells {max_cells}");
        }
    }

    /// Columns of successive M/D steps are consecutive around the unit.
    fn continuous(ops: &[Op], l: usize) -> bool {
        let mut exp: Option<u32> = None;
        ops.iter().all(|o| match o {
            Op::M(c) | Op::D(c) => {
                let ok = exp.map_or(true, |e| e == *c);
                exp = Some((*c + 1) % l as u32);
                ok
            }
            Op::I => true,
        })
    }

    #[test]
    fn windowed_path_is_one_alignment() {
        // a period-20 array over several windows, with a base deleted or inserted at
        // each offset around the first window seam (max_cells 4000 -> windows of 200)
        let unit = b"ACGTTGCAAGCTTAGGCATC";
        let arr: Vec<u8> = unit.iter().cycle().take(1000).copied().collect();
        let whole = path(&arr, unit, usize::MAX, ULTRA_COSTS);
        assert_eq!(whole.iter().filter(|o| matches!(o, Op::M(_))).count(), 1000);
        for at in 190..=260 {
            for del in [true, false] {
                let mut a = arr.clone();
                if del {
                    a.remove(at);
                } else {
                    a.insert(at, if arr[at] == b'A' { b'C' } else { b'A' });
                }
                let ops = path(&a, unit, 4000, ULTRA_COSTS);
                assert_eq!(ops.iter().filter(|o| !matches!(o, Op::D(_))).count(), a.len(), "at {at}");
                assert!(continuous(&ops, unit.len()), "at {at} del {del}");
                // one indel, nothing else
                let (i, d) = (ops.iter().filter(|o| **o == Op::I).count(), ops.iter().filter(|o| matches!(o, Op::D(_))).count());
                assert_eq!((i, d), if del { (0, 1) } else { (1, 0) }, "at {at} del {del}");
                let want = ultra_deltas(&a, unit.len(), &path(&a, unit, usize::MAX, ULTRA_COSTS)).iter().sum::<f64>();
                let got = ultra_deltas(&a, unit.len(), &ops).iter().sum::<f64>();
                assert!((want - got).abs() < 1e-6, "at {at} del {del}: {got} vs {want}");
            }
        }
    }

    #[test]
    fn stitch_pays_for_a_phase_change() {
        // no base agrees in the overlap: forward by 2 columns, then back by 2
        let l = 10;
        let a: Vec<Op> = (0..20).map(|i| Op::M(i % 10)).collect();
        for shift in [2u32, 8] {
            let next: Vec<Op> = (0..20).map(|i| Op::M((i + 6 + shift) % 10)).collect();
            let mut out = a.clone();
            stitch(&mut out, &next, 4, l);
            assert!(continuous(&out, l), "shift {shift}: {out:?}");
            assert_eq!(out.iter().filter(|o| !matches!(o, Op::D(_))).count(), 36);
            let (i, d) = (out.iter().filter(|o| **o == Op::I).count(), out.iter().filter(|o| matches!(o, Op::D(_))).count());
            assert_eq!((i, d), if shift == 2 { (0, 2) } else { (2, 0) });
        }
    }

    #[test]
    fn path_covers_every_base() {
        let unit = b"ACGTTGCA";
        let mut arr = unit.repeat(6);
        arr.remove(20); // deletion
        arr.insert(33, b'G'); // insertion
        let ops = path(&arr, unit, 1 << 20, ULTRA_COSTS);
        assert_eq!(ops.iter().filter(|o| !matches!(o, Op::D(_))).count(), arr.len());
        assert!(ops.iter().any(|o| matches!(o, Op::D(_))));
        assert!(ops.contains(&Op::I));
    }

    #[test]
    fn ultra_matches_ultra_on_perfect_array() {
        let unit = b"AACGTTGC";
        let arr = unit.repeat(8);
        let ops = path(&arr, unit, 1 << 20, ULTRA_COSTS);
        let d = ultra_deltas(&arr, unit.len(), &ops);
        assert_eq!(d.len(), 64);
        assert!(d[..8].iter().all(|&x| x == U_BG));
        assert_eq!(d[8], U_ENTRY);
        assert!(d[9..].iter().all(|&x| x == U_MATCH));
        // ULTRA 1.2.1 on this array: sum of its deltas 68.91
        let s: f64 = d.iter().sum();
        assert!((s - 68.91).abs() < 0.05, "{s}");
    }

    #[test]
    fn ultra_substitution_counts_twice() {
        let unit = b"ACGTTGCA";
        let mut arr = unit.repeat(8);
        arr[30] = if arr[30] == b'A' { b'C' } else { b'A' };
        let ops = path(&arr, unit, 1 << 20, ULTRA_COSTS);
        let d = ultra_deltas(&arr, unit.len(), &ops);
        assert_eq!(d.iter().filter(|&&x| x == U_MISMATCH).count(), 2);
    }
}
