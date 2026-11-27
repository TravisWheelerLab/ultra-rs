//! BPAF (Binary Pairwise Alignment Format, version 0) output of tandem
//! repeats aligned to their consensus units.
//!
//! Each repeat is one record with `TARGET_CIRCULAR` set (spec §4.4.1): the
//! repeat on the sequence is the query, the consensus unit is the target,
//! `t_start` is the unit column the alignment starts in (1-based) and `t_len`
//! is the number of unit columns it passes through, so the copy number is
//! `t_len / period`. The units are stored in the file (`EXT_SEQUENCES`,
//! `tseqs`); the sequences the repeats are on are not (`qseqs` entries of
//! length 0), so a reader needs the original FASTA only to rebuild the query
//! row.
//!
//! The file is written front to back with no seeking: header, one record per
//! call to [`TandemWriter::write_record`], then the tables and footer from
//! [`TandemWriter::finish`].

use std::collections::{HashMap, HashSet};
use std::io::{self, Write};

use crate::wrap_align::{Alignment, Op};

pub const MAGIC: &[u8; 5] = b"BPAF\x01";
pub const VERSION: u8 = 0;
/// Header `extensions` bit: the `qseqs` and `tseqs` tables follow `mats`.
pub const EXT_SEQUENCES: u8 = 0x01;
/// Record flag: a `score` field follows the flags byte.
pub const HAS_SCORE: u8 = 0x02;
/// Record flag: target coordinates wrap around the target.
pub const TARGET_CIRCULAR: u8 = 0x20;

/// Units up to this long are named `(UNIT)n`; longer ones `TR<period>_<hash>`.
pub const SHORT_UNIT_NAME_MAX: usize = 32;

/// A run of one repeated non-`ACGT` letter at least this long gets a 4-bit
/// run of its own, stored as a single letter. Shorter ones stay in the run
/// around them: a separate run costs 16 bytes of run bookkeeping.
const REPEATED_LETTER_RUN_MIN: usize = 64;

fn put_uleb(mut x: u64, out: &mut Vec<u8>) {
    while x >= 0x80 {
        out.push((x & 0x7F) as u8 | 0x80);
        x >>= 7;
    }
    out.push(x as u8);
}

fn zigzag(n: i64) -> u64 {
    ((n << 1) ^ (n >> 63)) as u64
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// The record coordinates and CIGAR of a wraparound alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Encoded {
    /// Bases of the repeat the alignment covers.
    pub q_len: u64,
    /// The unit column the alignment starts in, 1-based.
    pub t_start: u64,
    /// Unit columns the alignment passes through.
    pub t_len: u64,
    pub match0: u64,
    /// `(indel, match_after)`: `indel < 0` is that many bases with no column,
    /// `indel > 0` that many columns with no base.
    pub pairs: Vec<(i64, u64)>,
}

/// `a` as BPAF coordinates and CIGAR, or `None` when it has no bases or uses
/// no columns (BPAF has no zero-length spans). Skipped columns at either end
/// are dropped, moving the start column for leading ones.
pub fn encode_alignment(a: &Alignment) -> Option<Encoded> {
    if a.period == 0 {
        return None;
    }
    let mut runs: &[(Op, u64)] = &a.runs;
    let mut phase = a.phase as u64 % a.period as u64;
    while let Some((&(Op::Skipped, n), rest)) = runs.split_first() {
        phase = (phase + n % a.period as u64) % a.period as u64;
        runs = rest;
    }
    while let Some((&(Op::Skipped, _), rest)) = runs.split_last() {
        runs = rest;
    }

    let (mut q_len, mut t_len, mut match0) = (0u64, 0u64, 0u64);
    let mut pairs: Vec<(i64, u64)> = Vec::new();
    for &(op, n) in runs {
        if n == 0 {
            continue;
        }
        match op {
            Op::Aligned => {
                q_len += n;
                t_len += n;
                match pairs.last_mut() {
                    Some(p) => p.1 += n,
                    None => match0 += n,
                }
            }
            Op::Inserted => {
                q_len += n;
                pairs.push((-(n as i64), 0));
            }
            Op::Skipped => {
                t_len += n;
                pairs.push((n as i64, 0));
            }
        }
    }
    if q_len == 0 || t_len == 0 {
        return None;
    }
    Some(Encoded {
        q_len,
        t_start: phase + 1,
        t_len,
        match0,
        pairs,
    })
}

/// The 4-bit code of an upper-case letter (spec §6.2.2), or `None`.
fn four_bit(c: u8) -> Option<u8> {
    b"ACGTBDHKMNRSVWXY"
        .iter()
        .position(|&x| x == c)
        .map(|i| i as u8)
}

fn is_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T')
}

/// Length of the stretch of one repeated letter starting at `seq[i]`.
fn same_letter_len(seq: &[u8], i: usize) -> usize {
    seq[i..].iter().take_while(|&&c| c == seq[i]).count()
}

/// A sequence as an `EXT_SEQUENCES` entry (spec §6.2), without its length
/// prefix. Lower-case letters become soft-mask intervals. Fails on a letter
/// the format cannot store or a sequence of 2³² letters or more.
pub fn encode_sequence(seq: &[u8]) -> io::Result<Vec<u8>> {
    if seq.len() > u32::MAX as usize {
        return Err(invalid(format!(
            "BPAF: a stored sequence is limited to {} letters, this one has {}",
            u32::MAX,
            seq.len()
        )));
    }
    let upper: Vec<u8> = seq.iter().map(|c| c.to_ascii_uppercase()).collect();
    if let Some(&c) = upper.iter().find(|&&c| four_bit(c).is_none()) {
        return Err(invalid(format!(
            "BPAF: letter {:?} cannot be stored in a sequence",
            c as char
        )));
    }

    let mut run_ends: Vec<u32> = Vec::new();
    let mut run_data_ends: Vec<u32> = Vec::new();
    let mut two: Vec<u8> = Vec::new(); // one letter code per element
    let mut four: Vec<u8> = Vec::new();
    let n = upper.len();
    let mut i = 0usize;
    while i < n {
        // Even run: 2-bit letters, possibly none.
        let start = i;
        while i < n && is_acgt(upper[i]) {
            i += 1;
        }
        two.extend(upper[start..i].iter().map(|&c| four_bit(c).expect("ACGT")));
        run_ends.push(i as u32);
        run_data_ends.push(two.len() as u32);
        if i == n {
            break;
        }
        // Odd run: 4-bit letters, at least one. A long stretch of one letter
        // is a run of its own, stored once.
        let start = i;
        if same_letter_len(&upper, i) >= REPEATED_LETTER_RUN_MIN {
            i += same_letter_len(&upper, i);
            four.push(four_bit(upper[start]).expect("checked"));
        } else {
            while i < n
                && !is_acgt(upper[i])
                && same_letter_len(&upper, i) < REPEATED_LETTER_RUN_MIN
            {
                i += 1;
            }
            four.extend(
                upper[start..i]
                    .iter()
                    .map(|&c| four_bit(c).expect("checked")),
            );
        }
        run_ends.push(i as u32);
        run_data_ends.push(four.len() as u32);
    }

    let mut masks: Vec<(u32, u32)> = Vec::new();
    let mut i = 0usize;
    while i < n {
        if seq[i].is_ascii_lowercase() {
            let start = i;
            while i < n && seq[i].is_ascii_lowercase() {
                i += 1;
            }
            masks.push((start as u32, (i - start) as u32));
        } else {
            i += 1;
        }
    }

    let mut out = Vec::with_capacity(8 + 8 * run_ends.len() + 8 * masks.len() + n / 4 + 2);
    let put32 = |out: &mut Vec<u8>, x: u32| out.extend_from_slice(&x.to_le_bytes());
    put32(&mut out, run_ends.len() as u32);
    run_ends.iter().for_each(|&x| put32(&mut out, x));
    run_data_ends.iter().for_each(|&x| put32(&mut out, x));
    put32(&mut out, masks.len() as u32);
    masks.iter().for_each(|m| put32(&mut out, m.0));
    masks.iter().for_each(|m| put32(&mut out, m.1));
    // First letter in the most significant bits; unused bits 0.
    for chunk in two.chunks(4) {
        let mut b = 0u8;
        for (k, &code) in chunk.iter().enumerate() {
            b |= code << (6 - 2 * k);
        }
        out.push(b);
    }
    for chunk in four.chunks(2) {
        out.push(chunk[0] << 4 | chunk.get(1).copied().unwrap_or(0));
    }
    Ok(out)
}

fn fnv1a(data: &[u8]) -> u64 {
    data.iter().fold(0xcbf29ce484222325u64, |h, &b| {
        (h ^ b as u64).wrapping_mul(0x100000001b3)
    })
}

/// The target id of a unit: `(UNIT)n` up to [`SHORT_UNIT_NAME_MAX`] bases,
/// else `TR<period>_<hash>` with a 64-bit hash of the unit, so the same unit
/// has the same id in every file.
pub fn unit_name(unit: &[u8]) -> String {
    if unit.len() <= SHORT_UNIT_NAME_MAX {
        format!("({})n", String::from_utf8_lossy(unit))
    } else {
        format!("TR{}_{:016x}", unit.len(), fnv1a(unit))
    }
}

/// A unit as it is stored: upper-case, with every letter other than `ACGT`
/// as `N`. ULTRA writes `*` for a consensus column it has no base for, and
/// the aligner treats every such column as unknown.
pub fn stored_unit(unit: &[u8]) -> Vec<u8> {
    unit.iter()
        .map(|c| match c.to_ascii_uppercase() {
            c if is_acgt(c) => c,
            _ => b'N',
        })
        .collect()
}

/// One repeat (or part of one) to write.
pub struct TandemRecord<'a> {
    /// Name of the sequence the repeat is on.
    pub query: &'a str,
    /// First base of the aligned stretch on that sequence, 1-based.
    pub q_start: u64,
    /// The consensus unit, as the alignment's columns.
    pub unit: &'a [u8],
    pub alignment: &'a Alignment,
    pub score: Option<f32>,
    /// What `score` is measured in (the `mats` table entry).
    pub scoring: &'a str,
}

/// Writes a BPAF file of tandem repeat records.
#[derive(Default)]
pub struct TandemWriter {
    /// Bytes written so far.
    pos: u64,
    qids: Vec<String>,
    qid_index: HashMap<String, u64>,
    /// Target id and unit, in table order.
    units: Vec<(String, Vec<u8>)>,
    unit_index: HashMap<Vec<u8>, u64>,
    unit_names: HashSet<String>,
    mats: Vec<String>,
    buf: Vec<u8>,
}

impl TandemWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Records written so far would be followed by tables at this offset.
    pub fn position(&self) -> u64 {
        self.pos
    }

    fn emit(&mut self, out: &mut dyn Write, bytes: &[u8]) -> io::Result<()> {
        out.write_all(bytes)?;
        self.pos += bytes.len() as u64;
        Ok(())
    }

    /// The 7-byte header. Call once, first.
    pub fn write_header(&mut self, out: &mut dyn Write) -> io::Result<()> {
        let mut h = MAGIC.to_vec();
        h.push(VERSION);
        h.push(EXT_SEQUENCES);
        self.emit(out, &h)
    }

    fn qidx(&mut self, name: &str) -> u64 {
        if let Some(&i) = self.qid_index.get(name) {
            return i;
        }
        let i = self.qids.len() as u64;
        self.qids.push(name.to_string());
        self.qid_index.insert(name.to_string(), i);
        i
    }

    fn tidx(&mut self, unit: &[u8]) -> u64 {
        if let Some(&i) = self.unit_index.get(unit) {
            return i;
        }
        // Two different units must not share an id, whatever the hash says.
        let base = unit_name(unit);
        let mut name = base.clone();
        let mut k = 2;
        while !self.unit_names.insert(name.clone()) {
            name = format!("{base}.{k}");
            k += 1;
        }
        let i = self.units.len() as u64;
        self.units.push((name, unit.to_vec()));
        self.unit_index.insert(unit.to_vec(), i);
        i
    }

    fn midx(&mut self, scoring: &str) -> u64 {
        match self.mats.iter().position(|m| m == scoring) {
            Some(i) => i as u64,
            None => {
                self.mats.push(scoring.to_string());
                self.mats.len() as u64 - 1
            }
        }
    }

    /// Writes one record. Returns `false`, writing nothing, when the
    /// alignment has no bases or uses no unit columns. The unit is stored as
    /// [`stored_unit`] gives it.
    pub fn write_record(&mut self, out: &mut dyn Write, r: &TandemRecord) -> io::Result<bool> {
        let Some(e) = encode_alignment(r.alignment) else {
            return Ok(false);
        };
        if r.unit.len() != r.alignment.period {
            return Err(invalid(format!(
                "BPAF: unit of {} bases for an alignment of period {}",
                r.unit.len(),
                r.alignment.period
            )));
        }
        let unit = stored_unit(r.unit);
        let (qidx, tidx, midx) = (self.qidx(r.query), self.tidx(&unit), self.midx(r.scoring));

        let mut body = std::mem::take(&mut self.buf);
        body.clear();
        for x in [qidx, tidx, r.q_start, e.q_len, e.t_start, e.t_len, midx] {
            put_uleb(x, &mut body);
        }
        match r.score {
            Some(s) => {
                body.push(TARGET_CIRCULAR | HAS_SCORE);
                body.extend_from_slice(&s.to_le_bytes());
            }
            None => body.push(TARGET_CIRCULAR),
        }
        put_uleb(e.pairs.len() as u64, &mut body);
        put_uleb(e.match0, &mut body);
        for &(indel, m) in &e.pairs {
            put_uleb(zigzag(indel), &mut body);
            put_uleb(m, &mut body);
        }

        let mut prefix = Vec::with_capacity(4);
        put_uleb(body.len() as u64, &mut prefix);
        self.emit(out, &prefix)?;
        let res = self.emit(out, &body);
        self.buf = body;
        res.map(|()| true)
    }

    /// The tables and footer. Call once, last.
    pub fn finish(&mut self, out: &mut dyn Write) -> io::Result<()> {
        let tables_start = self.pos;
        let mut t: Vec<u8> = Vec::new();
        let put_str = |t: &mut Vec<u8>, s: &str| {
            put_uleb(s.len() as u64, t);
            t.extend_from_slice(s.as_bytes());
        };
        put_uleb(self.qids.len() as u64, &mut t);
        self.qids.iter().for_each(|s| put_str(&mut t, s));
        put_uleb(self.units.len() as u64, &mut t);
        self.units
            .iter()
            .for_each(|(name, _)| put_str(&mut t, name));
        put_uleb(self.mats.len() as u64, &mut t);
        self.mats.iter().for_each(|s| put_str(&mut t, s));
        // qseqs: no query letters are stored (entries of length 0).
        put_uleb(self.qids.len() as u64, &mut t);
        t.resize(t.len() + self.qids.len(), 0);
        // tseqs: every unit.
        put_uleb(self.units.len() as u64, &mut t);
        for (_, unit) in &self.units {
            let entry = encode_sequence(unit)?;
            put_uleb(entry.len() as u64, &mut t);
            t.extend_from_slice(&entry);
        }
        t.extend_from_slice(&tables_start.to_le_bytes());
        t.extend_from_slice(MAGIC);
        self.emit(out, &t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment(period: usize, phase: usize, runs: &[(Op, u64)]) -> Alignment {
        Alignment {
            period,
            phase,
            runs: runs.to_vec(),
            score: 0,
        }
    }

    /// Spec §6.2.3, for checking the encoder.
    fn decode_sequence(e: &[u8]) -> Vec<u8> {
        let u32_at = |i: usize| u32::from_le_bytes(e[i..i + 4].try_into().unwrap()) as usize;
        let n = u32_at(0);
        let ends: Vec<usize> = (0..n).map(|i| u32_at(4 + 4 * i)).collect();
        let data_ends: Vec<usize> = (0..n).map(|i| u32_at(4 + 4 * n + 4 * i)).collect();
        let mut at = 4 + 8 * n;
        let m = u32_at(at);
        let starts: Vec<usize> = (0..m).map(|i| u32_at(at + 4 + 4 * i)).collect();
        let sizes: Vec<usize> = (0..m).map(|i| u32_at(at + 4 + 4 * m + 4 * i)).collect();
        at += 4 + 8 * m;
        let last = |parity: usize| {
            (0..n)
                .rev()
                .find(|i| i % 2 == parity)
                .map_or(0, |i| data_ends[i])
        };
        let two = &e[at..at + last(0).div_ceil(4)];
        let four = &e[at + two.len()..];
        assert_eq!(four.len(), last(1).div_ceil(2), "entry length");
        let mut out = Vec::new();
        let mut pos = 0;
        for i in 0..n {
            let data_start = if i >= 2 { data_ends[i - 2] } else { 0 };
            let data_len = data_ends[i] - data_start;
            for k in 0..ends[i] - pos {
                let j = data_start + k % data_len;
                out.push(if i % 2 == 0 {
                    b"ACGT"[(two[j / 4] >> (6 - 2 * (j % 4)) & 3) as usize]
                } else {
                    b"ACGTBDHKMNRSVWXY"[(four[j / 2] >> (4 - 4 * (j % 2)) & 15) as usize]
                });
            }
            pos = ends[i];
        }
        for (s, z) in starts.into_iter().zip(sizes) {
            out[s..s + z].make_ascii_lowercase();
        }
        out
    }

    #[test]
    fn decoder_reads_the_spec_example() {
        // Spec §6.2.4.
        let entry: Vec<u8> = [
            &[5u8, 0, 0, 0][..],
            &[
                4, 0, 0, 0, 12, 0, 0, 0, 12, 0, 0, 0, 14, 0, 0, 0, 18, 0, 0, 0,
            ],
            &[4, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 0, 8, 0, 0, 0],
            &[1, 0, 0, 0, 14, 0, 0, 0, 4, 0, 0, 0],
            &[0x1b, 0x1b, 0x9a, 0xf0],
        ]
        .concat();
        assert_eq!(entry.len(), 60);
        assert_eq!(decode_sequence(&entry), b"ACGTNNNNNNNNRYacgt");
    }

    #[test]
    fn sequences_round_trip() {
        let long_n = [b"AC".to_vec(), vec![b'N'; 200], b"RYNNT".to_vec()].concat();
        let cases: Vec<Vec<u8>> = vec![
            b"".to_vec(),
            b"A".to_vec(),
            b"ACG".to_vec(),
            b"ACGTACGTA".to_vec(),
            b"ACGTNNNNNNNNRYacgt".to_vec(),
            b"NACGT".to_vec(),
            b"ACGTN".to_vec(),
            b"nnnn".to_vec(),
            b"acgtRYKMswbdhvxN".to_vec(),
            vec![b'N'; 64],
            vec![b'N'; 65],
            long_n,
        ];
        for seq in cases {
            let e = encode_sequence(&seq).unwrap();
            assert_eq!(
                decode_sequence(&e),
                seq,
                "{}",
                String::from_utf8_lossy(&seq)
            );
        }
        // Plain ACGT: one run, no masks, four letters to a byte.
        assert_eq!(
            encode_sequence(b"ACG").unwrap(),
            [
                &[1u8, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0][..],
                &[0x18]
            ]
            .concat()
        );
        // 200 Ns are stored as one letter.
        let e = encode_sequence(&[b'N'; 200]).unwrap();
        assert_eq!(e.len(), 4 + 16 + 4 + 1);
        assert!(encode_sequence(b"ACGU").is_err());
        assert!(encode_sequence(b"AC-G").is_err());
    }

    #[test]
    fn alignments_become_coordinates_and_cigar() {
        use Op::*;
        // Spec §4.4.2: unit ACG, query CGACGTACAC.
        let a = alignment(
            3,
            1,
            &[
                (Aligned, 5),
                (Inserted, 1),
                (Aligned, 2),
                (Skipped, 1),
                (Aligned, 2),
            ],
        );
        assert_eq!(
            encode_alignment(&a),
            Some(Encoded {
                q_len: 10,
                t_start: 2,
                t_len: 10,
                match0: 5,
                pairs: vec![(-1, 2), (1, 2)],
            })
        );
        // Gaps at the ends and next to each other; skipped columns at the
        // ends are dropped and move the start column around the unit.
        let a = alignment(
            4,
            3,
            &[
                (Skipped, 2),
                (Inserted, 3),
                (Skipped, 1),
                (Aligned, 6),
                (Inserted, 1),
                (Skipped, 5),
            ],
        );
        assert_eq!(
            encode_alignment(&a),
            Some(Encoded {
                q_len: 10,
                t_start: 2,
                t_len: 7,
                match0: 0,
                pairs: vec![(-3, 0), (1, 6), (-1, 0)],
            })
        );
        assert_eq!(encode_alignment(&alignment(3, 0, &[(Inserted, 4)])), None);
        assert_eq!(encode_alignment(&alignment(3, 0, &[(Skipped, 4)])), None);
        assert_eq!(encode_alignment(&alignment(3, 0, &[])), None);
        assert_eq!(encode_alignment(&alignment(0, 0, &[(Aligned, 4)])), None);
    }

    #[test]
    fn file_bytes() {
        use Op::*;
        let a = alignment(
            3,
            1,
            &[
                (Aligned, 5),
                (Inserted, 1),
                (Aligned, 2),
                (Skipped, 1),
                (Aligned, 2),
            ],
        );
        let rec = |score| TandemRecord {
            query: "chr1",
            q_start: 101,
            unit: b"ACG",
            alignment: &a,
            score,
            scoring: "ultra",
        };
        let mut w = TandemWriter::new();
        let mut out: Vec<u8> = Vec::new();
        w.write_header(&mut out).unwrap();
        assert!(w.write_record(&mut out, &rec(None)).unwrap());
        assert!(w.write_record(&mut out, &rec(Some(1.5))).unwrap());
        let empty = alignment(3, 0, &[(Inserted, 2)]);
        assert!(!w
            .write_record(
                &mut out,
                &TandemRecord {
                    alignment: &empty,
                    ..rec(None)
                }
            )
            .unwrap());
        w.finish(&mut out).unwrap();

        // The record of spec §4.4.2, byte for byte.
        let spec = [
            0x0e, 0, 0, 0x65, 0x0a, 0x02, 0x0a, 0, 0x20, 2, 5, 1, 2, 2, 2,
        ];
        let scored = [
            &[0x12, 0, 0, 0x65, 0x0a, 0x02, 0x0a, 0, 0x22][..],
            &1.5f32.to_le_bytes(),
            &[2, 5, 1, 2, 2, 2],
        ]
        .concat();
        let tables = [
            &[1u8, 4][..],
            b"chr1",
            &[1, 6],
            b"(ACG)n",
            &[1, 5],
            b"ultra",
            &[1, 0], // qseqs: one entry, not stored
            &[1, 17, 1, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0x18],
        ]
        .concat();
        let tables_start = 7 + spec.len() + scored.len();
        let want = [
            &b"BPAF\x01\x00\x01"[..],
            &spec,
            &scored,
            &tables,
            &(tables_start as u64).to_le_bytes(),
            b"BPAF\x01",
        ]
        .concat();
        assert_eq!(out, want);
    }

    #[test]
    fn unit_ids() {
        assert_eq!(unit_name(b"CA"), "(CA)n");
        let long = vec![b'A'; 33];
        assert_eq!(unit_name(&long), format!("TR33_{:016x}", fnv1a(&long)));
        // The same unit is one target; a different unit is another.
        let a = alignment(2, 0, &[(Op::Aligned, 6)]);
        let mut w = TandemWriter::new();
        let mut out: Vec<u8> = Vec::new();
        for unit in [&b"CA"[..], b"TG", b"CA"] {
            let r = TandemRecord {
                query: "s",
                q_start: 1,
                unit,
                alignment: &a,
                score: None,
                scoring: "ultra",
            };
            w.write_record(&mut out, &r).unwrap();
        }
        assert_eq!(
            w.units.iter().map(|u| u.0.as_str()).collect::<Vec<_>>(),
            ["(CA)n", "(TG)n"]
        );
        // A column with no consensus base is stored as N.
        let r = TandemRecord {
            query: "s",
            q_start: 1,
            unit: b"c*",
            alignment: &a,
            score: None,
            scoring: "ultra",
        };
        assert!(w.write_record(&mut out, &r).unwrap());
        assert_eq!(w.units[2], ("(CN)n".to_string(), b"CN".to_vec()));
    }
}
