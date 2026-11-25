//! Wraparound alignment of an array to its cyclic unit: the one aligner behind the
//! consensus (`consensus.rs`) and the scoring paths (`deltas.rs`).
//!
//! The alignment is global in the array and free at both ends in the unit: every base
//! may sit in any column (cells = bases x unit). It minimises `cm` per substitution plus
//! `ci` per inserted base or skipped column.

/// One step of an array's alignment to its cyclic consensus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Op {
    /// Array base aligned to consensus column `col`.
    M(u32),
    /// Array base not in the consensus (insertion).
    I,
    /// Consensus column `col` skipped (deletion from the array).
    D(u32),
}

const BM: u8 = 1;
const BI: u8 = 2;
const BD: u8 = 3;

/// Buffers reused between alignments.
pub(crate) struct Scratch {
    /// Traceback, 2 bits per cell (`pack_row`).
    bt: Vec<u8>,
    /// The current row's traceback, a byte per cell.
    row: Vec<u8>,
    d: Vec<i32>,
    e: Vec<i32>,
}

impl Scratch {
    pub(crate) fn new() -> Scratch {
        Scratch { bt: Vec::new(), row: Vec::new(), d: Vec::new(), e: Vec::new() }
    }
}

/// Bytes per traceback row of `l` cells at 2 bits per cell.
#[inline]
fn row_bytes(l: usize) -> usize {
    l.div_ceil(4)
}

/// Pack a row of traceback codes (0-3) into 2 bits per cell.
#[inline]
fn pack_row(row: &[u8], dst: &mut [u8]) {
    let mut chunks = row.chunks_exact(4);
    for (d, ch) in dst.iter_mut().zip(&mut chunks) {
        // four codes, one per byte, to 2 bits each
        let x = u32::from_le_bytes([ch[0], ch[1], ch[2], ch[3]]);
        *d = ((x & 3) | ((x >> 6) & 0xC) | ((x >> 12) & 0x30) | ((x >> 18) & 0xC0)) as u8;
    }
    let rest = chunks.remainder();
    if !rest.is_empty() {
        let mut x = 0u8;
        for (t, &b) in rest.iter().enumerate() {
            x |= b << (2 * t);
        }
        dst[row.len() / 4] = x;
    }
}

/// Traceback code of cell `k` in the packed row starting at byte `off`.
#[inline]
fn bt_at(bt: &[u8], off: usize, k: usize) -> u8 {
    (bt[off + k / 4] >> (2 * (k % 4))) & 3
}

/// Align `s` to the cyclic `u` (base codes: 0-3 = ACGT, 4 = other; `u` holds 0-3 only
/// and is not empty). Appends the steps to `out` in array order and returns the cost and
/// the slot the alignment ends in (the column after the last one used).
pub(crate) fn align(s: &[u8], u: &[u8], cm: i32, ci: i32, sc: &mut Scratch, out: &mut Vec<Op>) -> (u64, usize) {
    let l = u.len();
    let n = s.len();
    let rb = row_bytes(l);
    let need = (n + 1) * rb;
    if sc.bt.len() < need {
        sc.bt.resize(need, 0);
    }
    sc.row.clear();
    sc.row.resize(l, 0);
    sc.d.clear();
    sc.d.resize(l, 0);
    sc.e.clear();
    sc.e.resize(l, 0);
    let (mut d, mut e) = (std::mem::take(&mut sc.d), std::mem::take(&mut sc.e));
    let bt = &mut sc.bt;
    let b: &mut [u8] = &mut sc.row[..l];
    for i in 1..=n {
        let a = s[i - 1];
        for k in 0..l {
            let pk = if k > 0 { k - 1 } else { l - 1 };
            let m = d[pk] + if a == 4 || u[pk] != a { cm } else { 0 };
            let ins = d[k] + ci;
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
            if e[j] + ci < e[nj] {
                e[nj] = e[j] + ci;
                b[nj] = BD;
            }
            j = nj;
        }
        pack_row(b, &mut bt[i * rb..(i + 1) * rb]);
        std::mem::swap(&mut d, &mut e);
    }
    let mut best = 0;
    for j in 1..l {
        if d[j] < d[best] {
            best = j;
        }
    }
    let cost = d[best] as u64;
    let (mut i, mut k) = (n, best);
    let start = out.len();
    while i > 0 {
        let pk = if k > 0 { k - 1 } else { l - 1 };
        match bt_at(bt, i * rb, k) {
            BI => {
                out.push(Op::I);
                i -= 1;
            }
            BM => {
                out.push(Op::M(pk as u32));
                i -= 1;
                k = pk;
            }
            _ => {
                out.push(Op::D(pk as u32));
                k = pk;
            }
        }
    }
    out[start..].reverse();
    sc.d = d;
    sc.e = e;
    (cost, best)
}
