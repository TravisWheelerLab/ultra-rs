//! `ultra --bpaf`: the file is read back with a reader made from the BPAF
//! specification and checked against the input and the BED output.

mod common {
    pub mod bpaf_decode;
}
use common::bpaf_decode;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn ultra(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_ultra"))
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}

/// xorshift64*
fn bases(state: &mut u64, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *state ^= *state >> 12;
            *state ^= *state << 25;
            *state ^= *state >> 27;
            b"ACGT"[(state.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 62) as usize]
        })
        .collect()
}

fn is_rotation(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && [a, a].concat().windows(b.len().max(1)).any(|w| w == b)
}

#[test]
fn one_array_with_an_insertion_a_deletion_and_a_substitution() {
    let d = dir("bpaf_one");
    let unit = b"ACGTTGCA";
    let mut copies: Vec<&[u8]> = vec![unit; 12];
    copies[4] = b"ACGTGCA"; // a T deleted
    copies[7] = b"ACGTTAGCA"; // an A inserted
    copies[9] = b"ACCTTGCA"; // G -> C
    let mut st = 7u64;
    let seq = [bases(&mut st, 300), copies.concat(), bases(&mut st, 300)].concat();
    fs::write(
        d.join("in.fa"),
        [b">chrT one array\n", &seq[..], b"\n"].concat(),
    )
    .unwrap();

    let out = ultra(&d, &["--bpaf", "--hs", "-o", "out.bpaf", "in.fa"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let recs = bpaf_decode::read(&fs::read(d.join("out.bpaf")).unwrap());
    let arr: Vec<_> = recs.iter().filter(|r| r.period() == 8).collect();
    assert_eq!(arr.len(), 1, "one period-8 record");
    let r = arr[0];
    assert_eq!(r.query, "chrT", "the id is the first word of the header");
    assert!(
        is_rotation(&r.unit, unit),
        "{}",
        String::from_utf8_lossy(&r.unit)
    );
    assert_eq!(
        r.unit_id,
        format!("({})n", String::from_utf8_lossy(&r.unit))
    );
    assert_eq!(r.scoring, "ultra");
    assert!(r.score > 50.0, "{}", r.score);
    // the array is bases 301..=396
    assert!(
        r.q_start.abs_diff(301) <= 8 && (r.q_start + r.q_len).abs_diff(397) <= 8,
        "{} +{}",
        r.q_start,
        r.q_len
    );
    assert!((r.copies() - 12.0).abs() <= 1.0, "{} copies", r.copies());
    assert_eq!(
        r.gaps(),
        (1, 1),
        "one inserted base, one skipped column: {:?}",
        r.pairs
    );
    let (aligned, same) = r.walk(&seq);
    assert!(aligned - same <= 3, "{same} of {aligned} columns agree");
    for r in &recs {
        r.walk(&seq);
    }
}

#[test]
fn records_cover_the_repeats_of_the_bed_output() {
    let d = dir("bpaf_bed");
    // random sequence with noisy arrays of several periods, and an array
    // whose unit changes half way (sub-repeats)
    let mut st = 0x1234_5678_9abc_def1u64;
    let mut seqs: Vec<(String, Vec<u8>)> = Vec::new();
    for s in 0..3 {
        let mut seq = Vec::new();
        for period in [2usize, 5, 9, 14, 3, 7] {
            seq.extend(bases(&mut st, 400));
            let unit = bases(&mut st, period);
            let mut array = unit.repeat(30);
            for i in (0..array.len()).step_by(23) {
                array[i] = bases(&mut st, 1)[0];
            }
            seq.extend(array);
        }
        seq.extend(bases(&mut st, 400));
        let a = bases(&mut st, 10);
        let mut b = a.clone();
        for i in [1, 3, 6, 8] {
            b[i] = if b[i] == b'A' { b'C' } else { b'A' };
        }
        seq.extend(a.repeat(40));
        seq.extend(b.repeat(40));
        seq.extend(bases(&mut st, 400));
        seqs.push((format!("seq{s}"), seq));
    }
    let mut fa = Vec::new();
    for (name, seq) in &seqs {
        fa.extend(format!(">{name} extra words\n").bytes());
        for line in seq.chunks(70) {
            fa.extend(line);
            fa.push(b'\n');
        }
    }
    fs::write(d.join("in.fa"), fa).unwrap();

    let out = ultra(&d, &["--bpaf", "--bed", "-t", "2", "-o", "out", "in.fa"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let recs = bpaf_decode::read(&fs::read(d.join("out.bpaf")).unwrap());
    let seq_of = |name: &str| &seqs.iter().find(|s| s.0 == name).expect("known sequence").1;
    let (mut aligned, mut same) = (0, 0);
    for r in &recs {
        let (a, s) = r.walk(seq_of(&r.query));
        aligned += a;
        same += s;
    }
    assert!(
        same as f64 > 0.9 * aligned as f64,
        "{same} of {aligned} columns agree"
    );

    // Each BED line is one repeat; its records, one per sub-repeat, follow
    // one another and cover it exactly.
    let bed = fs::read_to_string(d.join("out.bed")).unwrap();
    let mut it = recs.iter().peekable();
    let (mut repeats, mut split) = (0, 0);
    for line in bed.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (name, start, end): (&str, u64, u64) =
            (f[0], f[1].parse().unwrap(), f[2].parse().unwrap());
        let mut at = start;
        let mut parts = 0;
        while at < end {
            let r = it.next().unwrap_or_else(|| panic!("no record for {line}"));
            assert_eq!((r.query.as_str(), r.q_start - 1), (name, at), "{line}");
            at += r.q_len;
            parts += 1;
        }
        assert_eq!(at, end, "{line}");
        repeats += 1;
        split += (parts > 1) as usize;
    }
    assert!(it.next().is_none(), "records beyond the BED lines");
    assert!(repeats >= 18, "{repeats} repeats");
    assert!(split >= 1, "no repeat was written as several records");
}

#[test]
fn bpaf_needs_an_output_file() {
    let d = dir("bpaf_noout");
    fs::write(d.join("in.fa"), b">s\nACACACACACACACACACACACACACACACAC\n").unwrap();
    // as for ULTRA's other option errors: the message, and nothing is run
    let out = ultra(&d, &["--bpaf", "in.fa"]);
    let text = String::from_utf8_lossy(&out.stderr);
    assert!(text.contains("--bpaf requires an output file"), "{text}");
    assert!(out.stdout.is_empty(), "output without --out");
    // on its own with -o: the one output file, plus the settings
    let out = ultra(&d, &["--bpaf", "-o", "only.bpaf", "in.fa"]);
    assert!(out.status.success());
    assert!(fs::read(d.join("only.bpaf"))
        .unwrap()
        .starts_with(b"BPAF\x01"));
}
