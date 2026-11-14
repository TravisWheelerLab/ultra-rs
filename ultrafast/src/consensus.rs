//! TRF-style iterative consensus for one tandem array.
//!
//! Starting from a seed unit, repeat up to `rounds` times:
//!  1. wraparound-align the whole array to the cyclic current consensus
//!     (unit-cost edit distance, free starting and ending phase), with
//!     traceback; long arrays are cut into windows, each aligned on its own
//!     to the same consensus so the columns stay consistent;
//!  2. build a profile: per consensus column, counts of A/C/G/T/deletion;
//!     per slot between columns, how often it is crossed and how often an
//!     insertion occurs there (plus the first inserted letter);
//!  3. vote with TRF's unweighted rule: a column takes the plurality of
//!     A, C, G, T, '-' (the column is dropped when '-' wins, ties included);
//!     a slot gets the most common inserted letter when insertions occur in
//!     at least half of its crossings;
//!  4. stop when the consensus no longer changes or its total edit distance
//!     stops improving (the best consensus seen is kept).

use rayon::prelude::*;

const BM: u8 = 1;
const BI: u8 = 2;
const BD: u8 = 3;

/// Base code: A C G T -> 0..3, anything else -> 4.
#[inline]
pub fn code(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

const LET: &[u8; 5] = b"ACGTN";

#[derive(Clone, Debug)]
pub struct Params {
    /// Maximum refinement rounds.
    pub rounds: usize,
    /// Window length in bp; 0 = max(20 kb, 40 * unit), limited by `max_cells`.
    pub window: usize,
    /// Upper bound on traceback cells (window * unit) per alignment.
    pub max_cells: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params { rounds: 6, window: 0, max_cells: 1 << 28 }
    }
}

#[derive(Clone, Debug)]
pub struct Consensus {
    /// Consensus unit (A/C/G/T).
    pub unit: Vec<u8>,
    /// Total edit distance of the array to the cyclic consensus.
    pub ed: usize,
    /// Rounds run.
    pub rounds: usize,
    /// Times the alignment wrapped past column 0 (about the copy number).
    pub wraps: usize,
}

struct Profile {
    cnt: Vec<u64>,   // L*5: A C G T del per column
    cross: Vec<u64>, // L: crossings of slot k (before column k)
    ins: Vec<u64>,   // L: crossings with >= 1 insertion
    ilet: Vec<u64>,  // L*4: first inserted letter per slot
    wraps: u64,
}

impl Profile {
    fn new(l: usize) -> Profile {
        Profile { cnt: vec![0; l * 5], cross: vec![0; l], ins: vec![0; l], ilet: vec![0; l * 4], wraps: 0 }
    }
    fn add(&mut self, o: &Profile) {
        let add = |a: &mut Vec<u64>, b: &Vec<u64>| a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
        add(&mut self.cnt, &o.cnt);
        add(&mut self.cross, &o.cross);
        add(&mut self.ins, &o.ins);
        add(&mut self.ilet, &o.ilet);
        self.wraps += o.wraps;
    }
}

impl Scratch {
    fn new() -> Scratch {
        Scratch { bt: Vec::new(), d: Vec::new(), e: Vec::new() }
    }
}

struct Scratch {
    bt: Vec<u8>,
    d: Vec<i32>,
    e: Vec<i32>,
}

/// Align s to cyclic u; add to the profile; return the edit distance.
fn align_window(s: &[u8], u: &[u8], sc: &mut Scratch, p: &mut Profile, mut cols: Option<&mut [i32]>) -> u64 {
    let l = u.len();
    let n = s.len();
    let need = (n + 1) * l;
    if sc.bt.len() < need {
        sc.bt.resize(need, 0);
    }
    sc.d.clear();
    sc.d.resize(l, 0);
    sc.e.clear();
    sc.e.resize(l, 0);
    let (mut d, mut e) = (std::mem::take(&mut sc.d), std::mem::take(&mut sc.e));
    let bt = &mut sc.bt;
    for i in 1..=n {
        let a = s[i - 1];
        let b = &mut bt[i * l..(i + 1) * l];
        for k in 0..l {
            let pk = if k > 0 { k - 1 } else { l - 1 };
            let m = d[pk] + if a == 4 || u[pk] != a { 1 } else { 0 };
            let ins = d[k] + 1;
            if m <= ins {
                e[k] = m;
                b[k] = BM;
            } else {
                e[k] = ins;
                b[k] = BI;
            }
        }
        // deletions wrap around the cycle: two passes settle them
        let mut j = 0;
        for _ in 0..2 * l {
            let nj = if j + 1 == l { 0 } else { j + 1 };
            if e[j] + 1 < e[nj] {
                e[nj] = e[j] + 1;
                b[nj] = BD;
            }
            j = nj;
        }
        std::mem::swap(&mut d, &mut e);
    }
    let mut best = 0;
    for j in 1..l {
        if d[j] < d[best] {
            best = j;
        }
    }
    let ed = d[best] as u64;

    // traceback: events come in reverse; an insertion run at slot k ends
    // (going backwards) when slot k is left by a match or deletion
    let (mut i, mut k) = (n, best);
    let (mut pending, mut first) = (false, 4u8);
    while i > 0 {
        let b = bt[i * l + k];
        let pk = if k > 0 { k - 1 } else { l - 1 };
        if b == BI {
            pending = true;
            first = s[i - 1];
            if let Some(c) = cols.as_deref_mut() {
                c[i - 1] = -1;
            }
            i -= 1;
        } else {
            p.cross[k] += 1;
            if pending {
                p.ins[k] += 1;
                if first < 4 {
                    p.ilet[k * 4 + first as usize] += 1;
                }
                pending = false;
            }
            if b == BM {
                let c = s[i - 1];
                if c < 4 {
                    p.cnt[pk * 5 + c as usize] += 1;
                }
                if let Some(cc) = cols.as_deref_mut() {
                    cc[i - 1] = pk as i32;
                }
                i -= 1;
            } else {
                p.cnt[pk * 5 + 4] += 1;
            }
            if pk == 0 {
                p.wraps += 1;
            }
            k = pk;
        }
    }
    if pending {
        p.cross[k] += 1;
        p.ins[k] += 1;
        if first < 4 {
            p.ilet[k * 4 + first as usize] += 1;
        }
    }
    sc.d = d;
    sc.e = e;
    ed
}

fn window_len(n: usize, l: usize, par: &Params) -> usize {
    let mut w = if par.window > 0 { par.window } else { (40 * l).max(20_000) };
    let cap = (par.max_cells / l.max(1)).max(3 * l);
    w = w.min(cap);
    w.min(n).max(1)
}

/// Align s to cyclic u window by window, adding to the profile. Windows are
/// independent, so a long array's windows run in parallel (rayon); their profiles
/// are integer counts and are summed, so the result doesn't depend on threading.
fn align_all(s: &[u8], u: &[u8], par: &Params, sc: &mut Scratch, p: &mut Profile) -> u64 {
    let w = window_len(s.len(), u.len(), par);
    if s.len() <= w {
        return align_window(s, u, sc, p, None);
    }
    let parts: Vec<(u64, Profile)> = s
        .par_chunks(w)
        .with_max_len(1)
        .map_init(Scratch::new, |sc, c| {
            let mut q = Profile::new(u.len());
            let e = align_window(c, u, sc, &mut q, None);
            (e, q)
        })
        .collect();
    parts.iter().map(|(e, q)| {
        p.add(q);
        e
    }).sum()
}

/// As `align_all`, also returning each position's consensus column (-1 = insertion).
fn align_cols(s: &[u8], u: &[u8], par: &Params, sc: &mut Scratch) -> (u64, Vec<i32>) {
    let w = window_len(s.len(), u.len(), par);
    let mut cols = vec![-1i32; s.len()];
    if s.len() <= w {
        let mut p = Profile::new(u.len());
        let ed = align_window(s, u, sc, &mut p, Some(&mut cols));
        return (ed, cols);
    }
    let ed = s
        .par_chunks(w)
        .zip(cols.par_chunks_mut(w))
        .with_max_len(1)
        .map_init(Scratch::new, |sc, (c, cc)| {
            let mut p = Profile::new(u.len());
            align_window(c, u, sc, &mut p, Some(cc))
        })
        .sum();
    (ed, cols)
}

fn vote(p: &Profile, l: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(2 * l);
    for k in 0..l {
        if p.cross[k] > 0 && 2 * p.ins[k] >= p.cross[k] {
            let il = &p.ilet[k * 4..k * 4 + 4];
            let mut bl = 0;
            for c in 1..4 {
                if il[c] > il[bl] {
                    bl = c;
                }
            }
            out.push(bl as u8);
        }
        let c = &p.cnt[k * 5..k * 5 + 5];
        let (mut bc, mut bv) = (4, c[4]); // '-' wins ties
        for x in 0..4 {
            if c[x] > bv {
                bv = c[x];
                bc = x;
            }
        }
        if bc != 4 {
            out.push(bc as u8);
        }
    }
    out
}

/// Consensus of `array` (ASCII) starting from `seed` (ASCII).
pub fn refine(array: &[u8], seed: &[u8], par: &Params) -> Consensus {
    let s: Vec<u8> = array.iter().map(|&c| code(c)).collect();
    let u: Vec<u8> = seed.iter().map(|&c| code(c)).filter(|&c| c < 4).collect();
    if u.is_empty() || s.is_empty() {
        return Consensus { unit: seed.to_vec(), ed: s.len(), rounds: 0, wraps: 0 };
    }
    let mut sc = Scratch::new();
    let (bu, ed, rounds, wraps) = refine_codes(&s, u, par, &mut sc);
    Consensus { unit: bu.iter().map(|&c| LET[c as usize]).collect(), ed: ed as usize, rounds, wraps: wraps as usize }
}

/// `refine` on base codes: (unit, edit distance, rounds, wraps).
fn refine_codes(s: &[u8], mut u: Vec<u8>, par: &Params, sc: &mut Scratch) -> (Vec<u8>, u64, usize, u64) {
    let mut best: Option<(Vec<u8>, u64, u64)> = None;
    let mut rounds = 0;
    for r in 0..par.rounds {
        let mut p = Profile::new(u.len());
        let ed = align_all(s, &u, par, sc, &mut p);
        rounds = r + 1;
        match &best {
            Some((_, be, _)) if ed >= *be => break, // no improvement: keep the best
            _ => best = Some((u.clone(), ed, p.wraps)),
        }
        let nu = vote(&p, u.len());
        if nu.is_empty() || nu == u {
            break;
        }
        u = nu;
    }
    let (bu, ed, wraps) = best.unwrap();
    (bu, ed, rounds, wraps)
}

// ------------------------------------------------------------- splitting

/// Splitting an array where its unit changes (A...AB...B with the same period).
#[derive(Clone, Debug)]
pub struct SplitParams {
    /// Each part has at least this many copies of the unit ...
    pub min_copies: usize,
    /// ... and at least this many bp.
    pub min_len: usize,
    /// A split must save at least this many edits (parts aligned to their own
    /// re-voted units vs both aligned to the parent unit) ...
    pub min_gain: f64,
    /// ... and at least this fraction of the parent unit's edits in the two parts.
    pub min_gain_frac: f64,
    /// Per-column G statistic above which a column counts towards the change-point score.
    pub col_g: f64,
    /// At most this many parts per array.
    pub max_parts: usize,
}

impl Default for SplitParams {
    fn default() -> Self {
        SplitParams { min_copies: 3, min_len: 24, min_gain: 8.0, min_gain_frac: 0.2, col_g: 8.0, max_parts: 64 }
    }
}

#[inline]
fn xlx(x: u32) -> f64 {
    if x > 1 { let x = x as f64; x * x.ln() } else { 0.0 }
}

/// G statistic of one column's base counts, left part vs right part.
#[inline]
fn col_g(l: &[u32], t: &[u32]) -> f64 {
    let (mut a, mut nl, mut nt) = (0.0, 0u32, 0u32);
    for b in 0..4 {
        a += xlx(l[b]) + xlx(t[b] - l[b]) - xlx(t[b]);
        nl += l[b];
        nt += t[b];
    }
    2.0 * (a - xlx(nl) - xlx(nt - nl) + xlx(nt))
}

/// Most likely change point in s (codes) given each position's consensus column:
/// the position maximising the sum over columns of max(0, G_c - col_g), with at
/// least `minseg` bp on each side. Returns (position, score).
fn change_point(s: &[u8], cols: &[i32], l: usize, minseg: usize, cg: f64) -> Option<(usize, f64)> {
    let n = s.len();
    if n < 2 * minseg {
        return None;
    }
    let mut tot = vec![0u32; l * 4];
    for (i, &c) in cols.iter().enumerate() {
        if c >= 0 && s[i] < 4 {
            tot[c as usize * 4 + s[i] as usize] += 1;
        }
    }
    let mut left = vec![0u32; l * 4];
    let mut gc = vec![0.0f64; l];
    let mut score = 0.0;
    let mut best: Option<(usize, f64)> = None;
    for i in 0..n - minseg {
        // position i moves to the left part; split point is i + 1
        let c = cols[i];
        if c >= 0 && s[i] < 4 {
            let c = c as usize;
            left[c * 4 + s[i] as usize] += 1;
            let g = col_g(&left[c * 4..c * 4 + 4], &tot[c * 4..c * 4 + 4]);
            score += (g - cg).max(0.0) - (gc[c] - cg).max(0.0);
            gc[c] = g;
        }
        if i + 1 >= minseg && best.is_none_or(|(_, b)| score > b + 1e-9) {
            best = Some((i + 1, score));
        }
    }
    best.filter(|&(_, b)| b > 0.0)
}

fn cyclic_eq(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && (a.is_empty() || (0..a.len()).any(|r| a[r..].iter().chain(&a[..r]).eq(b.iter())))
}

/// Split `array` (ASCII) where its unit changes, starting from its consensus
/// `unit` (ASCII). Returns the parts as (offset in array, re-voted unit); one part
/// (offset 0, `unit`) when there's no split. Recursive binary segmentation: the best
/// change point by per-column G statistic is accepted when the two sides' re-voted
/// units differ and they save enough edits over the parent unit.
pub fn split(array: &[u8], unit: &[u8], par: &Params, sp: &SplitParams) -> Vec<(usize, Vec<u8>)> {
    let s: Vec<u8> = array.iter().map(|&c| code(c)).collect();
    let u: Vec<u8> = unit.iter().map(|&c| code(c)).filter(|&c| c < 4).collect();
    if u.is_empty() {
        return vec![(0, unit.to_vec())];
    }
    let mut sc = Scratch::new();
    let mut out = Vec::new();
    split_rec(&s, 0, u, par, sp, &mut sc, &mut out);
    if out.len() == 1 {
        return vec![(0, unit.to_vec())];
    }
    out.into_iter().map(|(o, u)| (o, u.iter().map(|&c| LET[c as usize]).collect())).collect()
}

fn split_rec(s: &[u8], off: usize, u: Vec<u8>, par: &Params, sp: &SplitParams, sc: &mut Scratch,
             out: &mut Vec<(usize, Vec<u8>)>) {
    let l = u.len();
    let minseg = (sp.min_copies * l).max(sp.min_len);
    // remaining budget: this part plus the ones still queued can't exceed max_parts
    if s.len() < 2 * minseg || out.len() + 2 > sp.max_parts {
        out.push((off, u));
        return;
    }
    let (ed, cols) = align_cols(s, &u, par, sc);
    let Some((k, _)) = change_point(s, &cols, l, minseg, sp.col_g) else {
        out.push((off, u));
        return;
    };
    // cheap screen before re-voting both sides: substitutions saved if each side
    // took its own plurality base per column (ignores indels; half the thresholds)
    let mut cnt = vec![[0u32; 4]; 2 * l];
    for (i, &c) in cols.iter().enumerate() {
        if c >= 0 && s[i] < 4 {
            cnt[(i >= k) as usize * l + c as usize][s[i] as usize] += 1;
        }
    }
    let saved: u32 = cnt.iter().enumerate().map(|(j, x)| x.iter().max().unwrap() - x[u[j % l] as usize]).sum();
    if (saved as f64) < 0.5 * sp.min_gain || (saved as f64) < 0.5 * sp.min_gain_frac * ed as f64 {
        out.push((off, u));
        return;
    }
    let (sl, sr) = s.split_at(k);
    let (ul, el, _, _) = refine_codes(sl, u.clone(), par, sc);
    let (ur, er, _, _) = refine_codes(sr, u.clone(), par, sc);
    if cyclic_eq(&ul, &ur) {
        out.push((off, u));
        return;
    }
    let mut p = Profile::new(l);
    let base = align_all(sl, &u, par, sc, &mut p) + align_all(sr, &u, par, sc, &mut p);
    let gain = base as f64 - (el + er) as f64;
    let ok = gain >= sp.min_gain && gain >= sp.min_gain_frac * base as f64;
    if std::env::var_os("FASTANRS_SPLIT_DEBUG").is_some() {
        eprintln!("splitcand\t{}\t{}\t{}\t{}\t{}\t{}", off + k, s.len(), l, base, gain, ok as u8);
    }
    if !ok {
        out.push((off, u));
        return;
    }
    split_rec(sl, off, ul, par, sp, sc, out);
    split_rec(sr, off + k, ur, par, sp, sc, out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_unit_with_indels() {
        let unit = b"ACGTTGCAAGGCTTAACCGA";
        let mut arr = Vec::new();
        for i in 0..20 {
            let mut c = unit.to_vec();
            if i % 4 == 1 {
                c.remove(5);
            }
            if i % 5 == 2 {
                c.insert(12, b'T');
            }
            if i % 3 == 0 {
                c[17] = b'A';
            }
            arr.extend(c);
        }
        let seed = &arr[40..60];
        let r = refine(&arr, seed, &Params::default());
        let doubled: Vec<u8> = [unit.as_slice(), unit.as_slice()].concat();
        assert_eq!(r.unit.len(), unit.len());
        assert!(doubled.windows(unit.len()).any(|w| w == r.unit.as_slice()), "{}", String::from_utf8_lossy(&r.unit));
    }

    #[test]
    fn splits_adjacent_same_period_arrays() {
        let sp = SplitParams::default();
        let par = Params::default();
        // (ATG)20 (ATC)20: one split at the junction, each part with its own unit
        let arr = [b"ATG".repeat(20), b"ATC".repeat(20)].concat();
        let r = refine(&arr, b"ATG", &par);
        let p = split(&arr, &r.unit, &par, &sp);
        assert_eq!(p.len(), 2, "{:?}", p);
        assert!((p[1].0 as i64 - 60).abs() <= 3, "{:?}", p);
        assert!(cyclic_eq(&p[0].1, b"ATG") && cyclic_eq(&p[1].1, b"ATC"), "{:?}", p);
        // a 20-bp unit with 2 differences, 10 copies each side
        let a = b"ACGTTGCAAGGCTTAACCGA".to_vec();
        let mut b = a.clone();
        b[3] = b'A';
        b[11] = b'T';
        let arr = [a.repeat(10), b.repeat(10)].concat();
        let r = refine(&arr, &a, &par);
        assert_eq!(split(&arr, &r.unit, &par, &sp).len(), 2);
        // a single array is left alone
        let arr = b"ATG".repeat(40);
        assert_eq!(split(&arr, b"ATG", &par, &sp), vec![(0, b"ATG".to_vec())]);
    }
}
