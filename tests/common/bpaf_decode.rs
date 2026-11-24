//! A reader for the BPAF files `ultra --bpaf` and `ultrafast --bpaf` write,
//! made from the format's specification (BPAF version 0, §4-§6) and not from
//! the writer, for the tests of both tools.

pub struct Rec {
    pub query: String,
    /// 1-based.
    pub q_start: u64,
    pub q_len: u64,
    pub unit_id: String,
    pub unit: Vec<u8>,
    /// 1-based column of the unit.
    pub t_start: u64,
    pub t_len: u64,
    pub score: f32,
    pub scoring: String,
    pub match0: u64,
    pub pairs: Vec<(i64, u64)>,
}

impl Rec {
    pub fn period(&self) -> usize {
        self.unit.len()
    }
    pub fn copies(&self) -> f64 {
        self.t_len as f64 / self.unit.len() as f64
    }
    /// Bases with no unit column, and unit columns with no base.
    pub fn gaps(&self) -> (u64, u64) {
        let ins = self
            .pairs
            .iter()
            .filter(|p| p.0 < 0)
            .map(|p| -p.0 as u64)
            .sum();
        let del = self
            .pairs
            .iter()
            .filter(|p| p.0 > 0)
            .map(|p| p.0 as u64)
            .sum();
        (ins, del)
    }
    /// Walks the alignment over `seq` (the whole sequence the repeat is on)
    /// and the unit laid end to end, checking that the CIGAR uses exactly the
    /// record's spans. Returns (aligned columns, those where the letters agree).
    pub fn walk(&self, seq: &[u8]) -> (u64, u64) {
        let p = self.unit.len() as u64;
        assert!(
            self.t_start >= 1 && self.t_start <= p,
            "t_start {} of {p}",
            self.t_start
        );
        assert!(self.q_len >= 1 && self.t_len >= 1);
        let q = &seq[self.q_start as usize - 1..(self.q_start + self.q_len) as usize - 1];
        let (mut qi, mut x) = (0u64, self.t_start);
        let (mut aligned, mut same) = (0u64, 0u64);
        let mut run = |n: u64, qi: &mut u64, x: &mut u64| {
            for _ in 0..n {
                let (a, b) = (q[*qi as usize], self.unit[((*x - 1) % p) as usize]);
                aligned += 1;
                same += a.eq_ignore_ascii_case(&b) as u64;
                *qi += 1;
                *x += 1;
            }
        };
        run(self.match0, &mut qi, &mut x);
        for &(indel, m) in &self.pairs {
            assert!(indel != 0);
            if indel > 0 {
                x += indel as u64;
            } else {
                qi += (-indel) as u64;
            }
            run(m, &mut qi, &mut x);
        }
        assert_eq!(
            (qi, x - self.t_start),
            (self.q_len, self.t_len),
            "CIGAR and spans"
        );
        // The writer's rules (§4.4.1).
        if let (Some(first), Some(last)) = (self.pairs.first(), self.pairs.last()) {
            assert!(
                !(self.match0 == 0 && first.0 > 0),
                "starts with skipped columns"
            );
            assert!(!(last.1 == 0 && last.0 > 0), "ends with skipped columns");
        }
        (aligned, same)
    }
}

struct Cur<'a>(&'a [u8], usize);

impl Cur<'_> {
    fn uleb(&mut self) -> u64 {
        let (mut v, mut s) = (0u64, 0);
        loop {
            let b = self.0[self.1];
            self.1 += 1;
            v |= ((b & 0x7F) as u64) << s;
            if b & 0x80 == 0 {
                return v;
            }
            s += 7;
        }
    }
    fn bytes(&mut self, n: usize) -> &[u8] {
        self.1 += n;
        &self.0[self.1 - n..self.1]
    }
    fn strings(&mut self) -> Vec<String> {
        (0..self.uleb())
            .map(|_| {
                let n = self.uleb() as usize;
                String::from_utf8(self.bytes(n).to_vec()).unwrap()
            })
            .collect()
    }
}

/// §6.2.3.
fn sequence(e: &[u8]) -> Vec<u8> {
    let u = |i: usize| u32::from_le_bytes(e[i..i + 4].try_into().unwrap()) as usize;
    let n = u(0);
    let ends: Vec<usize> = (0..n).map(|i| u(4 + 4 * i)).collect();
    let dends: Vec<usize> = (0..n).map(|i| u(4 + 4 * n + 4 * i)).collect();
    let mut at = 4 + 8 * n;
    let m = u(at);
    let masks: Vec<(usize, usize)> = (0..m)
        .map(|i| (u(at + 4 + 4 * i), u(at + 4 + 4 * m + 4 * i)))
        .collect();
    at += 4 + 8 * m;
    let total = |parity| {
        (0..n)
            .rev()
            .find(|i| i % 2 == parity)
            .map_or(0, |i| dends[i])
    };
    let two = &e[at..at + total(0).div_ceil(4)];
    let four = &e[at + two.len()..];
    assert_eq!(four.len(), total(1).div_ceil(2), "sequence entry length");
    let (mut out, mut pos) = (Vec::new(), 0);
    for i in 0..n {
        let start = if i >= 2 { dends[i - 2] } else { 0 };
        for k in 0..ends[i] - pos {
            let j = start + k % (dends[i] - start);
            out.push(if i % 2 == 0 {
                b"ACGT"[(two[j / 4] >> (6 - 2 * (j % 4)) & 3) as usize]
            } else {
                b"ACGTBDHKMNRSVWXY"[(four[j / 2] >> (4 - 4 * (j % 2)) & 15) as usize]
            });
        }
        pos = ends[i];
    }
    for (s, z) in masks {
        out[s..s + z].make_ascii_lowercase();
    }
    out
}

/// Every record of a tandem BPAF file, checking the layout on the way: the
/// header and footer, `EXT_SEQUENCES` with no query letters and every unit
/// stored, distinct unit ids, and `TARGET_CIRCULAR` with a score on each record.
pub fn read(b: &[u8]) -> Vec<Rec> {
    assert_eq!(&b[..7], b"BPAF\x01\x00\x01", "header");
    assert_eq!(&b[b.len() - 5..], b"BPAF\x01", "footer");
    let tables_start =
        u64::from_le_bytes(b[b.len() - 13..b.len() - 5].try_into().unwrap()) as usize;
    let mut t = Cur(b, tables_start);
    let (qids, tids, mats) = (t.strings(), t.strings(), t.strings());
    assert_eq!(t.uleb() as usize, qids.len(), "qseqs count");
    for _ in &qids {
        assert_eq!(t.uleb(), 0, "a query sequence is stored");
    }
    assert_eq!(t.uleb() as usize, tids.len(), "tseqs count");
    let units: Vec<Vec<u8>> = tids
        .iter()
        .map(|_| {
            let n = t.uleb() as usize;
            assert!(n > 0, "a unit is not stored");
            sequence(t.bytes(n))
        })
        .collect();
    assert_eq!(t.1, b.len() - 13, "tables end at the footer");
    let mut ids = tids.clone();
    ids.sort();
    ids.dedup();
    assert_eq!(ids.len(), tids.len(), "unit ids repeat");

    let mut c = Cur(b, 7);
    let mut recs = Vec::new();
    while c.1 < tables_start {
        let len = c.uleb() as usize;
        let end = c.1 + len;
        let (qi, ti) = (c.uleb() as usize, c.uleb() as usize);
        let (q_start, q_len, t_start, t_len) = (c.uleb(), c.uleb(), c.uleb(), c.uleb());
        let mi = c.uleb() as usize;
        assert_eq!(c.bytes(1), [0x22], "flags: TARGET_CIRCULAR | HAS_SCORE");
        let score = f32::from_le_bytes(c.bytes(4).try_into().unwrap());
        let (npairs, match0) = (c.uleb(), c.uleb());
        let pairs = (0..npairs)
            .map(|_| {
                let z = c.uleb();
                ((z >> 1) as i64 ^ -((z & 1) as i64), c.uleb())
            })
            .collect();
        assert_eq!(c.1, end, "record length");
        recs.push(Rec {
            query: qids[qi].clone(),
            q_start,
            q_len,
            unit_id: tids[ti].clone(),
            unit: units[ti].clone(),
            t_start,
            t_len,
            score,
            scoring: mats[mi].clone(),
            match0,
            pairs,
        });
    }
    assert_eq!(c.1, tables_start, "records end at the tables");
    recs
}
