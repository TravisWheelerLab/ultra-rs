//! Inputs from the pre-merge review of the branch, and the cases its tests did not
//! cover: several records, N runs, --filler, -p 0, SubRepeats.

use serde_json::Value;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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

fn dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn write_fa(path: &Path, recs: &[(&str, Vec<u8>)], eol: &str) {
    let mut s = Vec::new();
    for (name, seq) in recs {
        s.extend(format!(">{name}{eol}").bytes());
        for line in seq.chunks(80) {
            s.extend(line);
            s.extend(eol.bytes());
        }
    }
    fs::write(path, s).unwrap();
}

fn run(args: &[&str], fa: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_ultrafast")).args(args).arg(fa).output().unwrap()
}

/// Run with JSON on stdout; the records.
fn records(args: &[&str], fa: &Path) -> Vec<Value> {
    let out = run(args, fa);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let j: Value = serde_json::from_slice(&out.stdout).unwrap();
    j["Repeats"].as_array().unwrap().clone()
}

fn span(r: &Value) -> (u64, u64, u64) {
    let s = r["Start"].as_u64().unwrap();
    (s, s + r["Length"].as_u64().unwrap(), r["Period"].as_u64().unwrap())
}

fn of<'a>(recs: &'a [Value], name: &str) -> Vec<&'a Value> {
    recs.iter().filter(|r| r["SequenceName"] == name).collect()
}

#[test]
fn units_just_above_ultras_cap_keep_their_period() {
    // ULTRA at -p 20 calls these at period 20 with shifted ends; FasTAN's call must win
    let d = dir("cap_edge");
    let mut st = 77u64;
    let recs: Vec<(String, Vec<u8>)> = [19usize, 20, 21, 22, 23, 25]
        .iter()
        .map(|&p| {
            let unit = bases(&mut st, p);
            let mut s = bases(&mut st, 2000);
            s.extend(unit.repeat(30));
            s.extend(bases(&mut st, 2000));
            (format!("p{p}"), s)
        })
        .collect();
    let fa = d.join("edge.fa");
    write_fa(&fa, &recs.iter().map(|(n, s)| (n.as_str(), s.clone())).collect::<Vec<_>>(), "\n");
    let out = records(&["-t", "2"], &fa);
    for p in [19u64, 20, 21, 22, 23, 25] {
        let end = 2000 + 30 * p;
        let r = of(&out, &format!("p{p}"));
        assert!(
            r.iter().map(|x| span(x)).any(|(s, e, q)| q == p && s <= 2005 && e + 5 >= end),
            "p{p}: {:?}",
            r.iter().map(|x| span(x)).collect::<Vec<_>>()
        );
        assert!(
            !r.iter().map(|x| span(x)).any(|(s, e, q)| q != p && s < end && e > 2000 && (e.min(end) - s.max(2000)) * 2 > e - s),
            "p{p}: a call at another period inside the array: {:?}",
            r.iter().map(|x| span(x)).collect::<Vec<_>>()
        );
    }
}

#[test]
fn satellite_with_a_microsatellite_in_each_copy() {
    // 40 copies of a 60 bp unit holding (CAG)x8: one period-60 call, no CAG calls inside it
    let d = dir("nested");
    let mut st = 13u64;
    let mut unit = bases(&mut st, 18);
    unit.extend(b"CAG".repeat(8));
    unit.extend(bases(&mut st, 18));
    let mut s = bases(&mut st, 2000);
    s.extend(unit.repeat(40));
    s.extend(bases(&mut st, 2000));
    let fa = d.join("nested.fa");
    write_fa(&fa, &[("nested60", s)], "\n");
    let out = records(&["-t", "2"], &fa);
    let spans: Vec<_> = out.iter().map(span).collect();
    assert!(spans.iter().any(|&(s, e, p)| p == 60 && s <= 2005 && e + 5 >= 4400), "{spans:?}");
    assert!(!spans.iter().any(|&(s, e, p)| p <= 20 && s >= 2000 && e <= 4400), "{spans:?}");
    // FasTAN alone reports it too
    let out0 = records(&["-p", "0"], &fa);
    assert!(out0.iter().map(span).any(|(s, e, p)| p == 60 && s <= 2005 && e + 5 >= 4400));
}

/// `s` with each base replaced by another with probability `pct`/100.
fn substituted(state: &mut u64, s: &[u8], pct: u64) -> Vec<u8> {
    let draws = bases(state, 5 * s.len());
    s.iter()
        .zip(draws.chunks(5))
        .map(|(&c, d)| {
            // four draws make a number from 0 to 255; the fifth is the replacement
            let x = d[..4].iter().fold(0u64, |a, &b| 4 * a + b"ACGT".iter().position(|&y| y == b).unwrap() as u64);
            if x * 100 < pct * 256 && d[4] != c { d[4] } else { c }
        })
        .collect()
}

#[test]
fn microsatellite_inside_part_of_a_satellite_is_kept() {
    // a period-40 satellite with substitutions that holds 300 bp of (AGAT)n: FasTAN calls
    // the whole span at period 40 and scores badly on the AGAT bases, so ULTRA's AGAT call
    // must stay
    let d = dir("micro_in_sat");
    let mut st = 5u64;
    let unit = bases(&mut st, 40);
    let mut s = bases(&mut st, 2000);
    for _ in 0..50 {
        s.extend(substituted(&mut st, &unit, 8));
    }
    s.extend(b"AGAT".repeat(75));
    for _ in 0..24 {
        s.extend(substituted(&mut st, &unit, 8));
    }
    s.extend(bases(&mut st, 2000));
    let fa = d.join("in.fa");
    write_fa(&fa, &[("x", s)], "\n");
    let spans: Vec<_> = records(&["-t", "2"], &fa).iter().map(span).collect();
    assert!(spans.iter().any(|&(s, e, p)| p == 4 && s <= 4005 && e >= 4290), "{spans:?}");
    assert!(spans.iter().any(|&(s, e, p)| p == 40 && s <= 2100 && e >= 5200), "{spans:?}");
}

#[test]
fn divergent_array_is_not_reported_at_twice_its_period() {
    // 32 copies of a 12 bp unit with 15% substitutions (simulated). FasTAN calls it at
    // period 24, which scores a fraction of a bit below ULTRA's period-12 call; ULTRA's stays
    let d = dir("divergent");
    let s = b"TGGGGCTTATATTTGCAGGGAGCGGCTTTTCCATTCACTCGAATCTGATAATCGGGGACACGCGTAGACGATGCTAGTGGCTCCGTGCTAC\
CCCCCAGGTTAAACTGTACCATAAACGGTACAATATACTATACCATAAACTGTACCATAAGCTGTACCATAGTCTGTACCATAAACTGTATCACATACTG\
CACCATATGCTGTATCATAAACTGCACCATAAACTGTACCATGAACTGTACCATAAATGGTACCATGAACTGAATCATAAACTGTAGCATAAACTGTACC\
ATAAACTGCACCATACACTGTCGCCTAAACTGTAGCATAAGCTGTACCATACACTGTGCCATTGACTTTACCCAAAACTCTTCTCTTAACTGTACCATAA\
ACTGTACCATTAACTGGACCATAAACTGTACCATTAACTGGACAATAATCTGTCCCATAAAATGTTCCAGAACCTGTACCATAAACTATACCCCGGCACT\
CCGCCTGGGACTAAATTCACTGCAGACCCCTAGGGGGCACTGATTAGAAATGATTATGGCCGCAGACCTAACAAAAATGCTGGCGAATGTTAG"
        .to_vec();
    let fa = d.join("in.fa");
    write_fa(&fa, &[("d", s)], "\n");
    let spans: Vec<_> = records(&["-t", "2"], &fa).iter().map(span).collect();
    assert!(spans.iter().any(|&(s, e, p)| p == 12 && s <= 110 && e >= 470), "{spans:?}");
    assert!(!spans.iter().any(|&(_, _, p)| p == 24), "{spans:?}");
    // FasTAN alone does call it at period 24
    let alone: Vec<_> = records(&["-p", "0"], &fa).iter().map(span).collect();
    assert!(alone.iter().any(|&(_, _, p)| p == 24), "{alone:?}");
}

#[test]
fn higher_order_repeat_keeps_its_long_period() {
    // two 12 bp monomers that differ at three bases alternate: the unit is 24 bp. ULTRA at
    // -p 20 calls period 12 with mismatches in every copy; FasTAN's period-24 call must win
    let d = dir("hor");
    let mut st = 41u64;
    let a = bases(&mut st, 12);
    let mut b = a.clone();
    for i in [2, 6, 10] {
        b[i] = if a[i] == b'A' { b'C' } else { b'A' };
    }
    let mut s = bases(&mut st, 1000);
    for _ in 0..30 {
        s.extend(&a);
        s.extend(&b);
    }
    s.extend(bases(&mut st, 1000));
    let fa = d.join("in.fa");
    write_fa(&fa, &[("h", s)], "\n");
    let spans: Vec<_> = records(&["-t", "2"], &fa).iter().map(span).collect();
    assert!(spans.iter().any(|&(s, e, p)| p == 24 && s <= 1005 && e >= 1715), "{spans:?}");
    assert!(!spans.iter().any(|&(s, e, p)| p == 12 && s < 1700 && e > 1020), "{spans:?}");
}

#[test]
fn crlf_and_awkward_names() {
    let d = dir("names");
    let mut st = 5u64;
    let mut s = bases(&mut st, 1000);
    s.extend(b"CAG".repeat(40));
    s.extend(bases(&mut st, 1000));
    let crlf = d.join("crlf.fa");
    write_fa(&crlf, &[("seq1 a description", s.clone())], "\r\n");
    let r = records(&[], &crlf);
    assert!(r.iter().any(|x| x["SequenceName"] == "seq1" && span(x).2 == 3), "{r:?}");
    let odd = d.join("name.fa");
    write_fa(&odd, &[("we\"ird\\name\tx", s)], "\n");
    let r = records(&[], &odd);
    assert!(r.iter().any(|x| x["SequenceName"] == "we\"ird\\name" && span(x).2 == 3), "{r:?}");
}

#[test]
fn duplicate_names_are_an_error() {
    let d = dir("dup");
    let mut st = 5u64;
    let fa = d.join("dup.fa");
    write_fa(&fa, &[("seq1", bases(&mut st, 500)), ("seq1 other", bases(&mut st, 500))], "\n");
    let out = run(&[], &fa);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("more than once"));
}

#[test]
fn several_records_n_runs_and_subrepeats() {
    let d = dir("multi");
    let mut st = 99u64;
    // r1: random only; r2: N run, then (CAG)x40 after it; r3: A x15 then B x15 at period 40
    let r1 = bases(&mut st, 3000);
    let mut r2 = bases(&mut st, 1000);
    r2.extend(vec![b'N'; 500]);
    r2.extend(bases(&mut st, 500));
    r2.extend(b"CAG".repeat(40)); // 2000..2120
    r2.extend(bases(&mut st, 1000));
    let a = bases(&mut st, 40);
    let mut b = a.clone();
    for i in [3, 9, 15, 22, 28, 35] {
        b[i] = if b[i] == b'A' { b'C' } else { b'A' };
    }
    let mut r3 = bases(&mut st, 1500);
    r3.extend(a.repeat(15)); // 1500..2100
    r3.extend(b.repeat(15)); // 2100..2700
    r3.extend(bases(&mut st, 1500));
    let fa = d.join("multi.fa");
    write_fa(&fa, &[("r1", r1), ("r2", r2), ("r3", r3)], "\n");
    let out = records(&["-t", "3"], &fa);
    // records come out in input order, sorted by start, each within its sequence
    let names: Vec<&str> = out.iter().map(|r| r["SequenceName"].as_str().unwrap()).collect();
    let mut sorted = names.clone();
    sorted.sort_by_key(|n| ["r1", "r2", "r3"].iter().position(|x| x == n));
    assert_eq!(names, sorted);
    assert!(of(&out, "r2").iter().map(|x| span(x)).any(|(s, e, p)| p == 3 && s <= 2002 && e >= 2118), "{out:?}");
    let sat: Vec<_> = of(&out, "r3").into_iter().filter(|x| span(x).2 == 40).collect();
    assert_eq!(sat.len(), 1, "{out:?}");
    let subs = sat[0]["SubRepeats"].as_array().expect("SubRepeats");
    assert_eq!(subs.len(), 2, "{subs:?}");
    let j = subs[1]["Start"].as_u64().unwrap() + span(sat[0]).0;
    assert!(j.abs_diff(2100) <= 40, "split at {j}");
    for r in &out {
        assert_eq!(r["Length"].as_u64().unwrap() as usize, r["PositionScoreDeltas"].as_array().unwrap().len());
    }
    // -p 0: FasTAN only, which misses the short CAG array but keeps the satellite
    let out0 = records(&["-p", "0"], &fa);
    assert!(of(&out0, "r3").iter().any(|x| span(x).2 == 40));
}

#[test]
fn filler_file_and_unknown_sequences() {
    let d = dir("filler");
    let mut st = 21u64;
    let mut s = bases(&mut st, 1000);
    s.extend(b"CAG".repeat(40));
    s.extend(bases(&mut st, 1000));
    let fa = d.join("in.fa");
    write_fa(&fa, &[("s1 desc", s)], "\n");
    let deltas: Vec<f64> = (0..120).map(|i| if i < 3 { -0.0145 } else { 1.349365 }).collect();
    let filler = serde_json::json!({"Repeats": [
        {"SequenceName": "s1 desc", "Start": 1000, "Length": 120, "Period": 3, "Score": 155.0,
         "Consensus": "CAG", "PositionScoreDeltas": deltas},
        {"SequenceName": "elsewhere", "Start": 0, "Length": 30, "Period": 3, "Score": 20.0, "Consensus": "CAG"}]});
    let fj = d.join("filler.json");
    fs::write(&fj, filler.to_string()).unwrap();
    let out = run(&["--filler", fj.to_str().unwrap()], &fa);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stderr).contains("1 --filler calls name sequences not in"));
    let j: Value = serde_json::from_slice(&out.stdout).unwrap();
    let r = j["Repeats"].as_array().unwrap();
    let f: Vec<_> = r.iter().filter(|x| span(x) == (1000, 1120, 3)).collect();
    assert_eq!(f.len(), 1, "{r:?}");
    assert_eq!(f[0]["Score"], 155.0, "the filler's own Score is kept");
    assert_eq!(f[0]["SequenceName"], "s1");
    // -p with --filler is a usage error
    assert!(!run(&["--filler", fj.to_str().unwrap(), "-p", "10"], &fa).status.success());
}

#[test]
fn output_to_a_device_and_option_checks() {
    let d = dir("misc");
    let mut st = 3u64;
    let mut s = bases(&mut st, 500);
    s.extend(b"AT".repeat(30));
    s.extend(bases(&mut st, 500));
    let fa = d.join("in.fa");
    write_fa(&fa, &[("s", s)], "\n");
    assert!(run(&["-o", "/dev/null"], &fa).status.success());
    assert!(!run(&["--rounds", "0"], &fa).status.success());
    assert!(!run(&["--gc", "40"], &fa).status.success(), "--gc without --delta-model te");
    assert!(!run(&["--no-split", "--split-min-gain", "3"], &fa).status.success());
    assert!(!run(&["--rescore-only"], &fa).status.success(), "--rescore-only without --filler");
    assert!(!run(&["--rescore-only", "-p", "0"], &fa).status.success(), "--rescore-only without --filler");
    assert!(!run(&["--delta-scale", "nan"], &fa).status.success());
    assert!(!run(&["--filler-min-score", "nan"], &fa).status.success());
    // nothing left behind next to the output
    let o = d.join("o.json");
    assert!(run(&["-o", o.to_str().unwrap()], &fa).status.success());
    let left: Vec<_> = fs::read_dir(&d).unwrap().map(|e| e.unwrap().file_name()).collect();
    assert_eq!(left.len(), 2, "{left:?}");
}

/// Gaps, '*' and ';' comment lines are not sequence (as ULTRA counts it), so
/// FasTAN's and ULTRA's coordinates agree and an array is reported once.
#[test]
fn non_sequence_bytes_do_not_shift_coordinates() {
    let d = dir("non_sequence_bytes");
    let mut st = 0x5eed_0001u64;
    let mut seq = bases(&mut st, 200);
    seq.extend(b"CA".repeat(40));
    seq.extend(bases(&mut st, 200));
    let mut text = b">s1\n".to_vec();
    text.extend(b"-".repeat(50));
    // A comment long enough that counting it as sequence would shift FasTAN
    // off ULTRA's call.
    text.extend(b"\n;");
    text.extend(b"comment".repeat(12));
    text.extend(b"\n");
    text.extend(&seq);
    text.extend(b"\n");
    let fa = d.join("gaps.fa");
    fs::write(&fa, &text).unwrap();
    let clean = d.join("clean.fa");
    write_fa(&clean, &[("s1", seq)], "\n");
    let got: Vec<_> = records(&["-t", "1"], &fa).iter().map(span).collect();
    let want: Vec<_> = records(&["-t", "1"], &clean).iter().map(span).collect();
    assert_eq!(got, want);
    assert_eq!(got.len(), 1, "{got:?}");
}

/// Sequence before the first header is a record named "", as a header with
/// nothing after '>' is. Two such records holding sequence are refused with
/// a message naming them; a record without sequence may share a name, and
/// ULTRA's calls go to the record that has the sequence.
#[test]
fn empty_and_shared_names() {
    let d = dir("empty_names");
    let mut st = 0x5eed_0005u64;
    let mut seq = bases(&mut st, 1000);
    seq.extend(b"CAG".repeat(40));
    seq.extend(bases(&mut st, 1000));
    let other = bases(&mut st, 500);
    let fa = |name: &str, parts: &[&[u8]]| {
        let p = d.join(name);
        fs::write(&p, parts.concat()).unwrap();
        p
    };
    for f in [
        fa("before_and_bare.fa", &[&other, b"\n>\n", &seq, b"\n"]),
        fa("two_bare.fa", &[b">\n", &other, b"\n> \n", &seq, b"\n"]),
    ] {
        let out = run(&["-t", "1"], &f);
        assert!(!out.status.success());
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(err.contains("\"\" occurs more than once (records 1 and 2)") && err.contains("empty name"), "{err}");
    }
    // no sequence after the bare header, or in the first or the last "seq1"
    for f in [
        fa("before_then_empty.fa", &[&seq, b"\n>\n"]),
        fa("empty_first.fa", &[b">seq1\n>seq1 x\n", &seq, b"\n"]),
        fa("empty_last.fa", &[b">seq1\n", &seq, b"\n>seq1 x\n"]),
    ] {
        let r = records(&["-t", "1"], &f);
        let cag: Vec<_> = r.iter().filter(|x| span(x).2 == 3).collect();
        assert_eq!(cag.len(), 1, "{f:?}: {r:?}");
        let (s, e, _) = span(cag[0]);
        assert!(s <= 1010 && e >= 1110, "{f:?}: {r:?}");
    }
}

/// --revote-filler keeps the period when the middle copy of an ULTRA call is
/// partly N (it seeded the vote with a shorter copy before).
#[test]
fn revote_seed_skips_n() {
    let d = dir("revote_n");
    let mut st = 0x5eed_0002u64;
    let unit = b"AGGTCCATGTTAGC";
    let mut seq = bases(&mut st, 100);
    for _ in 0..6 {
        seq.extend_from_slice(unit);
    }
    seq.extend_from_slice(&unit[..3]);
    seq.extend(b"N".repeat(7));
    seq.extend_from_slice(&unit[10..]);
    for _ in 0..6 {
        seq.extend_from_slice(unit);
    }
    seq.extend(bases(&mut st, 100));
    let fa = d.join("nmid.fa");
    write_fa(&fa, &[("n1", seq)], "\n");
    let recs = records(&["-t", "1", "--revote-filler"], &fa);
    let r = recs.iter().find(|r| r["Length"].as_u64().unwrap() > 100).expect("the array");
    assert_eq!(r["Period"], 14, "{r}");
    let c = r["Consensus"].as_str().unwrap();
    assert!(format!("{0}{0}", unit.iter().map(|&b| b as char).collect::<String>()).contains(c), "{c}");
}

/// Input ULTRA cannot read a second time (a pipe) is refused rather than
/// giving ULTRA an empty file.
#[test]
fn pipe_input_is_refused() {
    use std::io::Write;
    use std::process::Stdio;
    let mut child = Command::new(env!("CARGO_BIN_EXE_ultrafast"))
        .args(["-t", "1", "/dev/stdin"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let _ = child.stdin.take().unwrap().write_all(b">s\nACGTACGT\n");
    let out = child.wait_with_output().unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("not a pipe"), "{}", String::from_utf8_lossy(&out.stderr));
}

/// A header that is not UTF-8 is accepted; a file name starting with '-'
/// reaches ULTRA.
#[test]
fn latin1_headers_and_dash_file_names() {
    let d = dir("latin1_dash");
    let mut st = 0x5eed_0003u64;
    let mut seq = bases(&mut st, 200);
    seq.extend(b"CA".repeat(40));
    seq.extend(bases(&mut st, 200));
    let mut text = b">s1 caf\xe9\n".to_vec();
    text.extend(&seq);
    text.extend(b"\n");
    let fa = d.join("latin1.fa");
    fs::write(&fa, &text).unwrap();
    let recs = records(&["-t", "1"], &fa);
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0]["SequenceName"], "s1");

    fs::write(d.join("-x.fa"), &text).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_ultrafast"))
        .current_dir(&d)
        .args(["-t", "1", "--", "-x.fa"])
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let j: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j["Repeats"].as_array().unwrap().len(), 1);
}
