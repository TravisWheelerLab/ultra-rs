//! Tandem-array finder: a Rust port of FasTAN (Gene Myers, v0.8, `FasTAN.c`
//! `spectrum_block` and `compress_thread`) and of the wave aligner it calls
//! (`align.c` `Local_Alignment`, `forward_wave`, `reverse_wave`, from
//! daligner).
//!
//! Each contig (a maximal run of A/C/G/T) is scanned in 32 kb blocks that
//! start 24 kb apart.  In a block, consecutive occurrences of the same 8-mer
//! give seeds (position, offset d < 30,000); offsets holding more than
//! max(1, d/64) seeds are tried, most seeds first, by a local alignment of
//! the sequence against itself along diagonal d (band +-0.2 d, unit-cost
//! edit distance, extended while the last 60 columns stay near 70% identity).
//! An alignment spanning at least 1.95 d is a tandem array of period about d;
//! its bases are set to N so it is not found again.
//!
//! The original code carries this notice:
//!
//!   Copyright (c) 2020, Dr. Eugene W. Myers (EWM). All rights reserved.
//!
//!   Redistribution and use in source and binary forms, with or without
//!   modification, are permitted provided that the following conditions are met:
//!    - Redistributions of source code or derivative codes must retain the
//!      above copyright notice, this list of conditions and the following
//!      disclaimer.
//!    - Redistributions in binary form must reproduce the above copyright
//!      notice, this list of conditions and the following disclaimer in the
//!      documentation and/or other materials provided with the distribution.
//!    - The name of EWM may not be used to endorse or promote products derived
//!      from this software without specific prior written permission.
//!
//!   THIS SOFTWARE IS PROVIDED BY EWM "AS IS" AND ANY EXPRESS OR IMPLIED
//!   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
//!   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
//!   EVENT SHALL EWM BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//!   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//!   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//!   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//!   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use rayon::prelude::*;
use std::collections::VecDeque;

// FasTAN.c
const TSPACE: i32 = 100;
const DIAG_MAX: usize = 30000;
const DIAG_MIN: usize = 2;
const BLOCK: usize = 0x8000;
const BLOCK_STEP: usize = 0x6000;
const AVE_CORR: f64 = 0.7;

// align.c
const TRIM_LEN: u32 = 15;
const DUB_TRIM: i32 = 45;
const PATH_LEN: i32 = 60;
const PATH_TOP: u64 = 1 << 60;
const PATH_INT: u64 = PATH_TOP - 1;
const TRIM_MASK: u64 = 0x7fff;
const TRIM_MLAG: i32 = 250;
const WAVE_LAG: i32 = 70;
const FRACTION: f64 = 1000.0;
const BIAS_FACTOR: [f64; 10] = [0.690, 0.690, 0.690, 0.690, 0.780, 0.850, 0.900, 0.933, 0.966, 1.000];

/// One tandem array, in sequence coordinates: the alignment of
/// `[abeg, aend)` (the later copies) to `[bbeg, bend)` (the earlier ones);
/// the array spans `[bbeg, aend)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TandemArray {
    pub seq: usize,
    pub bbeg: usize,
    pub aend: usize,
    pub abeg: usize,
    pub bend: usize,
    pub period: usize,
    pub diffs: usize,
}

/// Path tip tables and the identity floor (`New_Align_Spec`).
pub struct Spec {
    ave_path: i32,
    score: Vec<i16>,
    table: Vec<i16>,
}

impl Spec {
    pub fn new(ave_corr: f64, freq: [f32; 4]) -> Spec {
        let mut m = (freq[0] + freq[3]) as f64;
        if m.is_nan() {
            m = 0.5;
        }
        if m > 0.5 {
            m = 1.0 - m;
        }
        let bias = if m < 0.2 { 3 } else { ((m + 0.025) * 20.0 - 1.0) as usize };
        let bf = BIAS_FACTOR[bias.min(9)];
        let ave_path = (PATH_LEN as f64 * (1.0 - bf * (1.0 - ave_corr))) as i32;
        let mscore = (FRACTION * bf * (1.0 - ave_corr)) as i32;
        let dscore = FRACTION as i32 - mscore;
        let n = (TRIM_MASK + 1) as usize;
        let (mut score, mut table) = (vec![0i16; n], vec![0i16; n]);
        set_table(0, 0, 0, 0, mscore, dscore, &mut table, &mut score);
        Spec { ave_path, score, table }
    }
}

/// table[x] = minimum over prefixes of the trim score (relative to the best
/// prefix), score[x] = total score, for the last TRIM_LEN columns x
/// (1 = match, 0 = difference).
fn set_table(bit: u32, prefix: usize, score: i32, max: i32, ms: i32, ds: i32, table: &mut [i16], sc: &mut [i16]) {
    if bit >= TRIM_LEN {
        table[prefix] = (score - max) as i16;
        sc[prefix] = score as i16;
    } else {
        let max = if score > max { score } else { max };
        set_table(bit + 1, prefix << 1, score - ds, max, ms, ds, table, sc);
        set_table(bit + 1, (prefix << 1) | 1, score + ms, max, ms, ds, table, sc);
    }
}

#[derive(Clone, Copy)]
struct Pebble {
    ptr: i32,
    diag: i32,
    diff: i32,
    mark: i32,
}

#[derive(Default, Debug)]
struct Path {
    abpos: i32,
    bbpos: i32,
    aepos: i32,
    bepos: i32,
    diffs: i32,
    /// (differences, b advance) per trace interval of TSPACE a-bases
    trace: VecDeque<(u16, u16)>,
}

/// Wave vectors over diagonals [base, base + len).
struct Work {
    base: i32,
    v: Vec<i32>,
    m: Vec<i32>,
    ha: Vec<i32>,
    na: Vec<i32>,
    t: Vec<u64>,
    cells: Vec<Pebble>,
}

impl Work {
    fn new() -> Work {
        Work { base: 0, v: vec![], m: vec![], ha: vec![], na: vec![], t: vec![], cells: vec![] }
    }
    /// Room for diagonals minp - 2 ..= maxp + 2.  Entries are not cleared:
    /// the waves write every diagonal before reading it, as in align.c.
    fn reset(&mut self, minp: i32, maxp: i32) {
        self.base = minp - 2;
        let n = (maxp - minp + 5) as usize;
        if self.v.len() < n {
            for x in [&mut self.v, &mut self.m, &mut self.ha, &mut self.na] {
                x.resize(n, 0);
            }
            self.t.resize(n, 0);
        }
    }
}

/// `buf` holds the contig with a 4 (N) on either side: base i is buf[i + 1].
fn forward_wave(w: &mut Work, spec: &Spec, buf: &[u8], path: &mut Path, mind: &mut i32, maxd: i32, mida: i32, minp: i32, maxp: i32) {
    let a_at = |x: i32| buf[(x + 1) as usize];
    let b_at = |x: i32, k: i32| buf[(x - k + 1) as usize];
    let ix = |k: i32, base: i32| (k - base) as usize;
    let base = w.base;
    let (v, m_, ha_, na_, t_) = (&mut w.v, &mut w.m, &mut w.ha, &mut w.na, &mut w.t);
    let cells = &mut w.cells;
    cells.clear();

    let mut hgh = maxd;
    let mut low = *mind;
    let mut dif = 0;
    let mut more = true;
    let mut aclip = i32::MAX;
    let mut bclip = -i32::MAX;
    let (mut besta, mut trima, mut lasta) = (mida, mida, mida);
    let (mut bestx, mut trimx) = ((mida + hgh) >> 1, (mida + hgh) >> 1);
    let mut trimd = 0;
    let mut trimha: i32 = 0;

    for k in (low..=hgh).rev() {
        let mut x = (mida + k) >> 1;
        let mut na = ((x + TSPACE) / TSPACE - 1) * TSPACE;
        cells.push(Pebble { ptr: -1, diag: k, diff: 0, mark: na });
        let mut ha = cells.len() as i32 - 1;
        na += TSPACE;
        loop {
            let c = b_at(x, k);
            if c == 4 {
                more = false;
                if bclip < k {
                    bclip = k;
                }
                break;
            }
            let d = a_at(x);
            if c != d {
                if d == 4 {
                    more = false;
                    aclip = k;
                }
                break;
            }
            x += 1;
        }
        let c = (x << 1) - k;
        while x >= na {
            cells.push(Pebble { ptr: ha, diag: k, diff: 0, mark: na });
            ha = cells.len() as i32 - 1;
            na += TSPACE;
        }
        if c > besta {
            besta = c;
            trima = c;
            lasta = c;
            bestx = x;
            trimx = x;
            trimha = ha;
        }
        let i = ix(k, base);
        v[i] = c;
        t_[i] = PATH_INT;
        m_[i] = PATH_LEN;
        ha_[i] = ha;
        na_[i] = na;
    }

    if !more {
        if b_at(bestx, bestx - (besta - bestx)) != 4 && a_at(bestx) != 4 {
            more = true;
        }
        if hgh >= aclip {
            hgh = aclip - 1;
        }
        if low <= bclip {
            low = bclip + 1;
        }
        aclip = i32::MAX;
        bclip = -i32::MAX;
    }

    while more && lasta >= besta - TRIM_MLAG {
        low -= 1;
        hgh += 1;
        if low >= minp {
            na_[ix(low, base)] = na_[ix(low + 1, base)];
            v[ix(low, base)] = -1;
        } else {
            low += 1;
        }
        let mut am;
        if hgh <= maxp {
            na_[ix(hgh, base)] = na_[ix(hgh - 1, base)];
            v[ix(hgh, base)] = -1;
            am = -1;
        } else {
            hgh -= 1;
            am = v[ix(hgh, base)];
        }
        dif += 1;

        v[ix(hgh + 1, base)] = -1;
        v[ix(low - 1, base)] = -1;
        let mut ac = -1;
        let mut t = PATH_INT;
        let mut n = PATH_LEN;
        let mut ua = -1;
        for k in (low..=hgh).rev() {
            let ap = ac;
            ac = am;
            let dd = k - 1;
            am = v[ix(dd, base)];
            let ik = ix(k, base);
            let id = ix(dd, base);
            let (mut c, mut m, mut b, mut ha);
            if ac < am {
                if am < ap {
                    c = ap + 1;
                    m = n;
                    b = t;
                    ha = ua;
                } else {
                    c = am + 1;
                    m = m_[id];
                    b = t_[id];
                    ha = ha_[id];
                }
            } else if ac < ap {
                c = ap + 1;
                m = n;
                b = t;
                ha = ua;
            } else {
                c = ac + 2;
                m = m_[ik];
                b = t_[ik];
                ha = ha_[ik];
            }
            if b & PATH_TOP != 0 {
                m -= 1;
            }
            b <<= 1;

            let mut x = (c + k) >> 1;
            loop {
                let cc = b_at(x, k);
                if cc == 4 {
                    more = false;
                    if bclip < k {
                        bclip = k;
                    }
                    break;
                }
                let d = a_at(x);
                if cc != d {
                    if d == 4 {
                        more = false;
                        aclip = k;
                    }
                    break;
                }
                x += 1;
                if b & PATH_TOP == 0 {
                    m += 1;
                }
                b = (b << 1) | 1;
            }
            c = (x << 1) - k;

            while x >= na_[ik] {
                if cells[ha as usize].mark < na_[ik] {
                    cells.push(Pebble { ptr: ha, diag: k, diff: dif, mark: na_[ik] });
                    ha = cells.len() as i32 - 1;
                }
                na_[ik] += TSPACE;
            }

            if c > besta {
                besta = c;
                bestx = x;
                if m >= spec.ave_path {
                    lasta = c;
                    if spec.table[(b & TRIM_MASK) as usize] >= 0
                        && spec.table[((b >> TRIM_LEN) & TRIM_MASK) as usize] as i32 + spec.score[(b & TRIM_MASK) as usize] as i32 >= 0
                    {
                        trima = c;
                        trimx = x;
                        trimd = dif;
                        trimha = ha;
                    }
                }
            }

            t = t_[ik];
            n = m_[ik];
            ua = ha_[ik];
            v[ik] = c;
            t_[ik] = b;
            m_[ik] = m;
            ha_[ik] = ha;
        }

        if !more {
            if b_at(bestx, bestx - (besta - bestx)) != 4 && a_at(bestx) != 4 {
                more = true;
            }
            if hgh >= aclip {
                hgh = aclip - 1;
            }
            if low <= bclip {
                low = bclip + 1;
            }
            aclip = i32::MAX;
            bclip = -i32::MAX;
        }

        let n = besta - WAVE_LAG;
        while hgh >= low {
            if v[ix(hgh, base)] < n {
                hgh -= 1;
            } else {
                while v[ix(low, base)] < n {
                    low += 1;
                }
                break;
            }
        }
    }

    // trace points along the trimmed path
    let trimy = trima - trimx;
    let mut trace = VecDeque::new();
    let mut a = -1;
    let mut h = trimha;
    while h >= 0 {
        let nb = cells[h as usize].ptr;
        cells[h as usize].ptr = a;
        a = h;
        h = nb;
    }
    h = a;
    let mut k = cells[h as usize].diag;
    let mut b = (mida - k) >> 1;
    let mut e = 0;
    low = k;
    h = cells[h as usize].ptr;
    while h >= 0 {
        let c = cells[h as usize];
        k = c.diag;
        let a = c.mark - k;
        let d = c.diff;
        trace.push_back(((d - e) as u16, (a - b) as u16));
        b = a;
        e = d;
        h = c.ptr;
    }
    if b + k != trimx {
        trace.push_back(((trimd - e) as u16, (trimy - b) as u16));
    } else if b != trimy {
        match trace.back_mut() {
            Some(p) => {
                p.0 = p.0.wrapping_add((trimd - e) as u16);
                p.1 = p.1.wrapping_add((trimy - b) as u16);
            }
            None => trace.push_back(((trimd - e) as u16, (trimy - b) as u16)),
        }
    }
    path.aepos = trimx;
    path.bepos = trimy;
    path.diffs = trimd;
    path.trace = trace;
    *mind = low;
}

fn reverse_wave(w: &mut Work, spec: &Spec, buf: &[u8], path: &mut Path, mind: i32, maxd: i32, mida: i32, minp: i32, maxp: i32) {
    // aseq - 1 and bseq - 1 in the C: base x - 1 is buf[x]
    let a_at = |x: i32| buf[x as usize];
    let b_at = |x: i32, k: i32| buf[(x - k) as usize];
    let ix = |k: i32, base: i32| (k - base) as usize;
    let base = w.base;
    let (v, m_, ha_, na_, t_) = (&mut w.v, &mut w.m, &mut w.ha, &mut w.na, &mut w.t);
    let cells = &mut w.cells;
    cells.clear();

    let mut hgh = maxd;
    let mut low = mind;
    let mut dif = 0;
    let mut more = true;
    let mut aclip = -i32::MAX;
    let mut bclip = i32::MAX;
    let (mut besta, mut trima, mut lasta) = (mida, mida, mida);
    let (mut bestx, mut trimx) = ((mida + hgh) >> 1, (mida + hgh) >> 1);
    let mut trimd = 0;
    let mut trimha: i32 = 0;

    for k in low..=hgh {
        let mut x = (mida + k) >> 1;
        let mut na = ((x + TSPACE - 1) / TSPACE - 1) * TSPACE;
        cells.push(Pebble { ptr: -1, diag: k, diff: 0, mark: x });
        let mut ha = cells.len() as i32 - 1;
        loop {
            let c = b_at(x, k);
            if c == 4 {
                more = false;
                if bclip > k {
                    bclip = k;
                }
                break;
            }
            let d = a_at(x);
            if c != d {
                if d == 4 {
                    more = false;
                    aclip = k;
                }
                break;
            }
            x -= 1;
        }
        let c = (x << 1) - k;
        while x <= na {
            cells.push(Pebble { ptr: ha, diag: k, diff: 0, mark: na });
            ha = cells.len() as i32 - 1;
            na -= TSPACE;
        }
        if c < besta {
            besta = c;
            trima = c;
            lasta = c;
            bestx = x;
            trimx = x;
            trimha = ha;
        }
        let i = ix(k, base);
        v[i] = c;
        t_[i] = PATH_INT;
        m_[i] = PATH_LEN;
        ha_[i] = ha;
        na_[i] = na;
    }

    if !more {
        if b_at(bestx, bestx - (besta - bestx)) != 4 && a_at(bestx) != 4 {
            more = true;
        }
        if low <= aclip {
            low = aclip + 1;
        }
        if hgh >= bclip {
            hgh = bclip - 1;
        }
        aclip = -i32::MAX;
        bclip = i32::MAX;
    }

    while more && lasta <= besta + TRIM_MLAG {
        low -= 1;
        hgh += 1;
        let mut ap;
        if low >= minp {
            na_[ix(low, base)] = na_[ix(low + 1, base)];
            v[ix(low, base)] = i32::MAX;
            ap = i32::MAX;
        } else {
            low += 1;
            ap = v[ix(low, base)];
        }
        if hgh <= maxp {
            na_[ix(hgh, base)] = na_[ix(hgh - 1, base)];
            v[ix(hgh, base)] = i32::MAX;
        } else {
            hgh -= 1;
        }
        dif += 1;

        v[ix(hgh + 1, base)] = i32::MAX;
        v[ix(low - 1, base)] = i32::MAX;
        let mut ac = i32::MAX;
        let mut t = PATH_INT;
        let mut n = PATH_LEN;
        let mut ua = -1;
        for k in low..=hgh {
            let am = ac;
            ac = ap;
            let dd = k + 1;
            ap = v[ix(dd, base)];
            let ik = ix(k, base);
            let id = ix(dd, base);
            let (mut c, mut m, mut b, mut ha);
            if ac > ap {
                if ap > am {
                    c = am - 1;
                    m = n;
                    b = t;
                    ha = ua;
                } else {
                    c = ap - 1;
                    m = m_[id];
                    b = t_[id];
                    ha = ha_[id];
                }
            } else if ac > am {
                c = am - 1;
                m = n;
                b = t;
                ha = ua;
            } else {
                c = ac - 2;
                m = m_[ik];
                b = t_[ik];
                ha = ha_[ik];
            }
            if b & PATH_TOP != 0 {
                m -= 1;
            }
            b <<= 1;

            let mut x = (c + k) >> 1;
            loop {
                let cc = b_at(x, k);
                if cc == 4 {
                    more = false;
                    if bclip > k {
                        bclip = k;
                    }
                    break;
                }
                let d = a_at(x);
                if cc != d {
                    if d == 4 {
                        more = false;
                        aclip = k;
                    }
                    break;
                }
                x -= 1;
                if b & PATH_TOP == 0 {
                    m += 1;
                }
                b = (b << 1) | 1;
            }
            c = (x << 1) - k;

            while x <= na_[ik] {
                if cells[ha as usize].mark > na_[ik] {
                    cells.push(Pebble { ptr: ha, diag: k, diff: dif, mark: na_[ik] });
                    ha = cells.len() as i32 - 1;
                }
                na_[ik] -= TSPACE;
            }

            if c < besta {
                besta = c;
                bestx = x;
                if m >= spec.ave_path {
                    lasta = c;
                    if spec.table[(b & TRIM_MASK) as usize] >= 0
                        && spec.table[((b >> TRIM_LEN) & TRIM_MASK) as usize] as i32 + spec.score[(b & TRIM_MASK) as usize] as i32 >= 0
                    {
                        trima = c;
                        trimx = x;
                        trimd = dif;
                        trimha = ha;
                    }
                }
            }

            t = t_[ik];
            n = m_[ik];
            ua = ha_[ik];
            v[ik] = c;
            t_[ik] = b;
            m_[ik] = m;
            ha_[ik] = ha;
        }

        if !more {
            if b_at(bestx, bestx - (besta - bestx)) != 4 && a_at(bestx) != 4 {
                more = true;
            }
            if low <= aclip {
                low = aclip + 1;
            }
            if hgh >= bclip {
                hgh = bclip - 1;
            }
            aclip = -i32::MAX;
            bclip = i32::MAX;
        }

        let n = besta + WAVE_LAG;
        while hgh >= low {
            if v[ix(hgh, base)] > n {
                hgh -= 1;
            } else {
                while v[ix(low, base)] > n {
                    low += 1;
                }
                break;
            }
        }
    }

    // trace points, prepended to the forward ones
    let trimy = trima - trimx;
    let fwd_empty = path.trace.is_empty();
    let trace = &mut path.trace;
    let add_front = |trace: &mut VecDeque<(u16, u16)>, d: i32, db: i32| match trace.front_mut() {
        Some(p) => {
            p.0 = p.0.wrapping_add(d as u16);
            p.1 = p.1.wrapping_add(db as u16);
        }
        None => trace.push_front((d as u16, db as u16)),
    };
    let mut a = -1;
    let mut h = trimha;
    while h >= 0 {
        let nb = cells[h as usize].ptr;
        cells[h as usize].ptr = a;
        a = h;
        h = nb;
    }
    h = a;
    let mut k = cells[h as usize].diag;
    let mut b = cells[h as usize].mark - k;
    let mut e = 0;
    if (b + k) % TSPACE != 0 {
        h = cells[h as usize].ptr;
        let (a, d) = if h < 0 {
            (trimy, trimd)
        } else {
            k = cells[h as usize].diag;
            (cells[h as usize].mark - k, cells[h as usize].diff)
        };
        if fwd_empty {
            trace.push_front(((d - e) as u16, (b - a) as u16));
        } else {
            add_front(trace, d - e, b - a);
        }
        b = a;
        e = d;
    }
    if h >= 0 {
        h = cells[h as usize].ptr;
        while h >= 0 {
            let c = cells[h as usize];
            k = c.diag;
            let a = c.mark - k;
            let d = c.diff;
            trace.push_front(((d - e) as u16, (b - a) as u16));
            b = a;
            e = d;
            h = c.ptr;
        }
        if b + k != trimx {
            trace.push_front(((trimd - e) as u16, (b - trimy) as u16));
        } else if b != trimy {
            add_front(trace, trimd - e, b - trimy);
        }
    }
    path.abpos = trimx;
    path.bbpos = trimy;
    path.diffs += trimd;
}

/// Longest local alignment through anti-diagonal `anti` on diagonal `d`
/// (a - b), kept within d +- `wide` (`Local_Alignment` with low = hgh = d).
fn local_alignment(w: &mut Work, spec: &Spec, buf: &[u8], d: i32, anti: i32, wide: i32) -> Path {
    let (minp, maxp) = (d - wide, d + wide);
    w.reset(minp, maxp);
    let mut path = Path::default();
    let mut low = d;
    forward_wave(w, spec, buf, &mut path, &mut low, d, anti, minp, maxp);
    let fshort = (path.aepos + path.bepos) - anti < DUB_TRIM;
    reverse_wave(w, spec, buf, &mut path, low, low, anti, minp, maxp);
    let rshort = anti - (path.abpos + path.bbpos) < DUB_TRIM;
    if fshort {
        if rshort {
            path.abpos = (path.abpos + path.aepos) >> 1;
            path.aepos = path.abpos;
            path.bbpos = (path.bbpos + path.bepos) >> 1;
            path.bepos = path.bbpos;
            path.trace.clear();
        } else {
            let mut low = path.abpos - path.bbpos;
            let anti = path.abpos + path.bbpos;
            path.trace.clear();
            let hgh = low;
            forward_wave(w, spec, buf, &mut path, &mut low, hgh, anti, minp, maxp);
        }
    } else if rshort {
        let low = path.aepos - path.bepos;
        let anti = path.aepos + path.bepos;
        path.trace.clear();
        path.diffs = 0;
        reverse_wave(w, spec, buf, &mut path, low, low, anti, minp, maxp);
    }
    path
}

/// Mean diagonal over the trace points: the unit length (`ave_tp_diag`).
fn ave_tp_diag(path: &Path) -> usize {
    let tlen = 2 * path.trace.len() as i64 - 2;
    let mut ab = path.abpos as i64;
    let mut bb = path.bbpos as i64;
    let mut ave = ab - bb;
    ab = (ab / TSPACE as i64) * TSPACE as i64;
    for (_, db) in path.trace.iter().take(((tlen.max(0)) / 2) as usize) {
        ab += TSPACE as i64;
        bb += *db as i64;
        ave += ab - bb;
    }
    ave += (path.aepos - path.bepos) as i64;
    let den = (tlen / 2) as f64 + 2.0;
    ((ave as f64 / den) + 0.5) as usize
}

struct Block {
    count: Vec<u16>,
    index: Vec<u16>,
    diags: Vec<u16>,
    post: Vec<(u16, u16)>, // (diag, ibeg)
    hits: Vec<(u16, u16)>,
    work: Work,
}

impl Block {
    fn new() -> Block {
        Block {
            count: vec![0; 0x10000],
            index: vec![0; BLOCK + 1],
            diags: vec![0; BLOCK],
            post: vec![(0, 0); BLOCK],
            hits: vec![(0, 0); BLOCK],
            work: Work::new(),
        }
    }
}

/// One block of `len` bases at contig offset `off` (`spectrum_block`);
/// returns the furthest array end, or -1 for a block too short for an 8-mer.
fn spectrum_block(bk: &mut Block, spec: &Spec, buf: &mut [u8], off: usize, len: usize, out: &mut Vec<(usize, usize, usize, usize, usize, usize)>) -> i64 {
    if len < 8 {
        return -1;
    }
    let s = |buf: &[u8], i: usize| buf[off + i + 1];
    let (count, index, diags, post, hits) = (&mut bk.count, &mut bk.index, &mut bk.diags, &mut bk.post, &mut bk.hits);
    count.iter_mut().for_each(|x| *x = 0);
    let l7 = len - 7;

    // count each 8-mer
    let mut kmer: u16 = s(buf, 0) as u16;
    for i in 1..7 {
        kmer = (kmer << 2) | s(buf, i) as u16;
    }
    for i in 0..l7 {
        kmer = (kmer << 2) | s(buf, i + 7) as u16;
        count[kmer as usize] = count[kmer as usize].wrapping_add(1);
    }
    let mut p: u16 = 0;
    for c in count.iter_mut() {
        let x = *c;
        *c = p;
        p = p.wrapping_add(x);
    }
    // positions in 8-mer order
    let mut kmer: u16 = s(buf, 0) as u16;
    for i in 1..7 {
        kmer = (kmer << 2) | s(buf, i) as u16;
    }
    for i in 0..l7 {
        kmer = (kmer << 2) | s(buf, i + 7) as u16;
        let c = &mut count[kmer as usize];
        index[*c as usize] = i as u16;
        *c = c.wrapping_add(1);
    }
    index[l7] = 0;
    index[0] |= 0x8000;
    for i in 0..0xffff {
        index[count[i] as usize] |= 0x8000;
        count[i] = 0;
    }
    count[0xffff] = 0;

    // seeds: consecutive positions of one 8-mer within DIAG_MAX
    let mut e = (index[0] & 0x7fff) as usize;
    for i in 1..l7 {
        let f = index[i] as usize;
        if f < 0x8000 {
            if f - e < DIAG_MAX {
                count[e] += 1;
            }
            e = f;
        } else {
            e = f & 0x7fff;
        }
    }
    let mut p: u16 = 0;
    for c in count[..0x8000].iter_mut() {
        let x = *c;
        *c = p;
        p += x;
    }
    diags[..DIAG_MAX].iter_mut().for_each(|x| *x = 0);
    let mut e = (index[0] & 0x7fff) as usize;
    for i in 1..l7 {
        let f = index[i] as usize;
        if f < 0x8000 {
            let d = f - e;
            if d < DIAG_MAX {
                let c = count[e] as usize;
                count[e] += 1;
                post[c] = (d as u16, e as u16);
                diags[d] += 1;
            }
            e = f;
        } else {
            e = f & 0x7fff;
        }
    }
    let mut p: u16 = 0;
    for c in diags[..DIAG_MAX].iter_mut() {
        let x = *c;
        *c = p;
        p += x;
    }
    for i in 0..p as usize {
        let c = post[i].0 as usize;
        hits[diags[c] as usize] = post[i];
        diags[c] += 1;
    }

    // offsets with enough seeds, most first
    let mut hist: Vec<(u16, usize)> = Vec::new();
    let mut p = diags[DIAG_MIN - 1];
    for i in DIAG_MIN..DIAG_MAX {
        let f = diags[i];
        let n = f - p;
        if n > 1 && n as usize > (i >> 6) {
            hist.push((n, i));
        }
        p = f;
    }
    hist.sort_by(|a, b| b.0.cmp(&a.0));

    let mut outhit: i64 = 0;
    for &(_, d) in &hist {
        let mut last: i64 = -1;
        for x in diags[d - 1] as usize + 1..diags[d] as usize {
            let p = hits[x].1 as usize;
            if (p as i64) < last || p - hits[x - 1].1 as usize > d || s(buf, p) >= 4 || s(buf, p + 7) >= 4 {
                continue;
            }
            let wide = ((0.2 * d as f64) as i32).max(1);
            let anti = 2 * (off + p) as i32 + d as i32;
            let path = local_alignment(&mut bk.work, spec, buf, d as i32, anti, wide);
            let end = path.aepos as i64 - off as i64;
            let beg = path.bbpos as i64 - off as i64;
            if beg > p as i64 || end <= p as i64 {
                continue;
            }
            if end > last {
                last = end;
            }
            if ((end - beg) as f64) < 1.95 * d as f64 {
                continue;
            }
            let unit = ave_tp_diag(&path);
            out.push((path.bbpos as usize, path.aepos as usize, path.abpos as usize, path.bepos as usize, unit, path.diffs as usize));
            let mut mask = |from: i64, to: i64| {
                for f in from..to {
                    buf[(off as i64 + f + 1) as usize] = 4;
                }
            };
            if path.bepos < path.abpos {
                mask(path.abpos as i64 - off as i64, end);
                mask(beg, path.bepos as i64 - off as i64);
            } else {
                mask(beg, end);
            }
            if path.aepos as i64 > outhit {
                outhit = path.aepos as i64;
            }
        }
    }
    outhit
}

/// All arrays in one contig (`compress_thread`), in contig coordinates.
fn scan_contig(spec: &Spec, contig: &[u8], bk: &mut Block) -> Vec<(usize, usize, usize, usize, usize, usize)> {
    let clen = contig.len();
    let mut buf = Vec::with_capacity(clen + 2);
    buf.push(4u8);
    buf.extend(contig.iter().map(|&c| match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        _ => 3,
    }));
    buf.push(4);
    let mut out = Vec::new();
    if clen <= BLOCK {
        spectrum_block(bk, spec, &mut buf, 0, clen, &mut out);
    } else {
        let mut p = 0usize;
        while p + 7 < clen {
            if p + BLOCK >= clen {
                spectrum_block(bk, spec, &mut buf, p, clen - p, &mut out);
                break;
            }
            let last = spectrum_block(bk, spec, &mut buf, p, BLOCK, &mut out);
            if last >= (p + BLOCK_STEP) as i64 {
                p = last as usize - BLOCK_STEP;
            }
            p += BLOCK_STEP;
        }
    }
    out
}

#[inline]
fn is_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Tandem arrays in every sequence of `text` (sequence k is
/// text[starts[k]..starts[k+1]]), in sequence coordinates, sorted.
pub fn tandem_arrays(text: &[u8], starts: &[usize]) -> Vec<TandemArray> {
    // contigs: maximal A/C/G/T runs
    let mut contigs: Vec<(usize, usize, usize)> = Vec::new(); // (seq, start, end) in text
    let mut fq = [0u64; 4];
    for k in 0..starts.len() {
        let end = if k + 1 < starts.len() { starts[k + 1] } else { text.len() };
        let mut i = starts[k];
        while i < end {
            while i < end && !is_acgt(text[i]) {
                i += 1;
            }
            let s = i;
            while i < end && is_acgt(text[i]) {
                fq[match text[i] {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    _ => 3,
                }] += 1;
                i += 1;
            }
            if i > s {
                contigs.push((k, s, i));
            }
        }
    }
    let tot = fq.iter().sum::<u64>().max(1) as f32;
    let spec = Spec::new(AVE_CORR, [fq[0] as f32 / tot, fq[1] as f32 / tot, fq[2] as f32 / tot, fq[3] as f32 / tot]);
    // Long contigs (e.g. gap-free T2T chromosomes) are scanned as overlapping windows so they use
    // all threads; calls cut by a window edge are merged afterwards (merge_windows).
    // tasks: (seq, start, end, window id); window ids are unique across the whole input
    // FASTANRS_WINDOW overrides the window length (0 = scan each contig whole, as C FasTAN does)
    let window = std::env::var("FASTANRS_WINDOW").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(WINDOW);
    let mut tasks: Vec<(usize, usize, usize, usize)> = Vec::new();
    for &(k, s, e) in &contigs {
        if window == 0 || e - s <= window + WINDOW_OVERLAP {
            tasks.push((k, s, e, tasks.len()));
        } else {
            let mut ws = s;
            loop {
                let we = (ws + window + WINDOW_OVERLAP).min(e);
                tasks.push((k, ws, we, tasks.len()));
                if we == e {
                    break;
                }
                ws += window;
            }
        }
    }
    // longest first, so one chromosome does not finish alone
    let mut order: Vec<usize> = (0..tasks.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(tasks[i].2 - tasks[i].1));
    let mut all: Vec<(TandemArray, usize)> = order
        .par_iter()
        .with_max_len(1) // one window per task: sizes are very uneven
        .map_init(Block::new, |bk, &i| {
            let (k, s, e, w) = tasks[i];
            let off = s - starts[k];
            scan_contig(&spec, &text[s..e], bk)
                .into_iter()
                .map(|(bb, ae, ab, be, unit, diffs)| {
                    (TandemArray {
                        seq: k,
                        bbeg: bb + off,
                        aend: ae + off,
                        abeg: ab + off,
                        bend: be + off,
                        period: unit,
                        diffs,
                    }, w)
                })
                .collect::<Vec<_>>()
        })
        .flatten()
        .collect();
    all.sort_unstable_by_key(|(t, w)| (t.seq, t.bbeg, t.aend, *w));
    merge_windows(all)
}

/// Window length and overlap for scanning long contigs in parallel.
pub const WINDOW: usize = 2_000_000;
pub const WINDOW_OVERLAP: usize = 100_000;

/// Merge calls of the same array found in different (overlapping) windows: two calls from
/// different windows that overlap and have similar periods (within 25%) become one call
/// spanning both. Calls from the same window are never merged, so nested arrays FasTAN
/// reports in one window are kept as they are. Input sorted by (seq, bbeg).
fn merge_windows(all: Vec<(TandemArray, usize)>) -> Vec<TandemArray> {
    let mut out: Vec<(TandemArray, Vec<usize>)> = Vec::with_capacity(all.len());
    for (a, w) in all {
        let mut merged = false;
        for (b, ws) in out.iter_mut().rev().take(64) {
            if b.seq != a.seq {
                break;
            }
            if b.aend <= a.bbeg || ws.contains(&w) {
                continue;
            }
            let (p, q) = (a.period.max(1) as f64, b.period.max(1) as f64);
            if p / q > 1.25 || q / p > 1.25 {
                continue;
            }
            let (la, lb) = ((a.aend - a.bbeg) as f64, (b.aend - b.bbeg) as f64);
            let rate = (a.diffs as f64 + b.diffs as f64) / (la + lb).max(1.0);
            let longer = if la > lb { a } else { *b };
            b.bbeg = b.bbeg.min(a.bbeg);
            b.aend = b.aend.max(a.aend);
            b.abeg = b.abeg.min(a.abeg);
            b.bend = b.bend.max(a.bend);
            b.period = longer.period;
            b.diffs = (rate * (b.aend - b.bbeg) as f64).round() as usize;
            ws.push(w);
            merged = true;
            break;
        }
        if !merged {
            out.push((a, vec![w]));
        }
    }
    let mut v: Vec<TandemArray> = out.into_iter().map(|(t, _)| t).collect();
    v.sort_unstable_by_key(|t| (t.seq, t.bbeg, t.aend));
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_seq(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    fn mutate(s: &[u8], every: usize, mut x: u64) -> Vec<u8> {
        s.iter()
            .enumerate()
            .map(|(i, &c)| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                if i % every == (x as usize) % every { b"ACGT"[((x >> 8) % 4) as usize] } else { c }
            })
            .collect()
    }

    #[test]
    fn finds_satellite_and_unit() {
        let unit = rand_seq(171, 7);
        let mut t = rand_seq(20_000, 1);
        for i in 0..30 {
            t.extend(mutate(&unit, 12, 100 + i));
        }
        t.extend(rand_seq(20_000, 2));
        let arr = tandem_arrays(&t, &[0]);
        assert!(!arr.is_empty());
        let a = arr[0];
        assert!((165..=177).contains(&a.period), "{arr:?}");
        assert!(a.bbeg <= 20_000 + 171 && a.aend >= 20_000 + 171 * 29, "{arr:?}");
    }

    #[test]
    fn finds_long_unit_array() {
        let unit = rand_seq(4_000, 11);
        let mut t = rand_seq(30_000, 3);
        for i in 0..4 {
            t.extend(mutate(&unit, 10, 200 + i));
        }
        t.extend(rand_seq(30_000, 4));
        let arr = tandem_arrays(&t, &[0]);
        assert!(arr.iter().any(|a| (3_900..=4_100).contains(&a.period) && a.aend - a.bbeg >= 12_000), "{arr:?}");
    }

    #[test]
    fn ignores_ltr_pair() {
        let ltr = rand_seq(500, 9);
        let mut e = rand_seq(10_000, 3);
        e.extend(&ltr);
        e.extend(rand_seq(2_500, 4));
        e.extend(&ltr);
        e.extend(rand_seq(10_000, 5));
        assert!(tandem_arrays(&e, &[0]).is_empty());
    }
}
