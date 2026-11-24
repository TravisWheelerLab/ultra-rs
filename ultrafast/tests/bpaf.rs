//! `ultrafast --bpaf`: the file is read back with a reader made from the BPAF
//! specification and checked against the input and the JSON output.

#[path = "../../tests/common/bpaf_decode.rs"]
mod bpaf_decode;

use serde_json::Value;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

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

#[test]
fn bpaf_records_match_the_json_calls() {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("uf_bpaf");
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    let mut st = 99u64;
    // r1: random only; r2: an N run, then (CAG)x40; r3: a period-40 satellite whose
    // unit changes half way (SubRepeats), then a period-90 array with a few edits
    let r1 = bases(&mut st, 3000);
    let mut r2 = bases(&mut st, 1000);
    r2.extend(vec![b'N'; 500]);
    r2.extend(bases(&mut st, 500));
    r2.extend(b"CAG".repeat(40));
    r2.extend(bases(&mut st, 1000));
    let a = bases(&mut st, 40);
    let mut b = a.clone();
    for i in [3, 9, 15, 22, 28, 35] {
        b[i] = if b[i] == b'A' { b'C' } else { b'A' };
    }
    let mut r3 = bases(&mut st, 1500);
    r3.extend(a.repeat(15));
    r3.extend(b.repeat(15));
    r3.extend(bases(&mut st, 1500));
    let unit = bases(&mut st, 90);
    let mut long = unit.repeat(10);
    long.remove(200);
    long.insert(500, b'A');
    long[700] = if long[700] == b'A' { b'C' } else { b'A' };
    r3.extend(long);
    r3.extend(bases(&mut st, 1500));
    let seqs = [("r1", r1), ("r2", r2), ("r3", r3)];
    let fa = d.join("in.fa");
    let mut text = Vec::new();
    for (name, seq) in &seqs {
        text.extend(format!(">{name} more words\n").bytes());
        for line in seq.chunks(80) {
            text.extend(line);
            text.push(b'\n');
        }
    }
    fs::write(&fa, text).unwrap();

    let (json, bpaf) = (d.join("out.json"), d.join("out.bpaf"));
    let out = Command::new(env!("CARGO_BIN_EXE_ultrafast"))
        .args(["-t", "3", "-o"])
        .arg(&json)
        .arg("--bpaf")
        .arg(&bpaf)
        .arg(&fa)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let bytes = fs::read(&bpaf).unwrap();
    let recs = bpaf_decode::read(&bytes);

    // every alignment is consistent with the sequence it is on
    let seq_of = |name: &str| &seqs.iter().find(|s| s.0 == name).expect("known sequence").1;
    let (mut aligned, mut same) = (0, 0);
    for r in &recs {
        let (x, y) = r.walk(seq_of(&r.query));
        aligned += x;
        same += y;
        assert_eq!(r.scoring, "ultra");
    }
    assert!(same as f64 > 0.97 * aligned as f64, "{same} of {aligned} columns agree");

    // one record per call, or per SubRepeat, in the JSON's order
    let j: Value = serde_json::from_slice(&fs::read(&json).unwrap()).unwrap();
    let mut want: Vec<(String, u64, u64, String, f64)> = Vec::new();
    for c in j["Repeats"].as_array().unwrap() {
        let (start, len) = (c["Start"].as_u64().unwrap(), c["Length"].as_u64().unwrap());
        let (name, score) = (c["SequenceName"].as_str().unwrap().to_string(), c["Score"].as_f64().unwrap());
        match c.get("SubRepeats").and_then(|s| s.as_array()) {
            Some(subs) => {
                for (k, s) in subs.iter().enumerate() {
                    let from = s["Start"].as_u64().unwrap();
                    let to = subs.get(k + 1).map_or(len, |n| n["Start"].as_u64().unwrap());
                    want.push((name.clone(), start + from, to - from, s["Consensus"].as_str().unwrap().into(), score));
                }
            }
            None => want.push((name, start, len, c["Consensus"].as_str().unwrap().into(), score)),
        }
    }
    let got: Vec<(String, u64, u64, String, f64)> = recs
        .iter()
        .map(|r| (r.query.clone(), r.q_start - 1, r.q_len, String::from_utf8(r.unit.clone()).unwrap(), r.score as f64))
        .collect();
    assert_eq!(got.len(), want.len(), "{got:?}\n{want:?}");
    for (g, w) in got.iter().zip(&want) {
        assert_eq!((&g.0, g.1, g.2, &g.3), (&w.0, w.1, w.2, &w.3));
        assert!((g.4 - w.4).abs() <= 1e-3 * w.4.abs().max(1.0), "score {} for {}", g.4, w.4);
    }

    // the cases this input was built for
    let cag = recs.iter().find(|r| r.query == "r2" && r.period() == 3).expect("the CAG array");
    assert_eq!(cag.unit_id, format!("({})n", String::from_utf8_lossy(&cag.unit)));
    assert!((cag.copies() - 40.0).abs() <= 2.0, "{} copies", cag.copies());
    let sat: Vec<_> = recs.iter().filter(|r| r.query == "r3" && r.period() == 40).collect();
    assert_eq!(sat.len(), 2, "the satellite's two parts");
    assert_eq!(sat[0].q_start + sat[0].q_len, sat[1].q_start, "parts follow one another");
    assert_ne!(sat[0].unit, sat[1].unit);
    assert!(sat.iter().all(|r| r.unit_id.starts_with("TR40_")), "{} {}", sat[0].unit_id, sat[1].unit_id);
    let long = recs.iter().find(|r| r.period() == 90).expect("the period-90 array");
    assert_eq!(long.gaps(), (1, 1), "{:?}", long.pairs);
    assert!((long.copies() - 10.0).abs() <= 1.0, "{} copies", long.copies());

    // the library writes the same bytes
    let opt = ultrafast::HybridOptions { threads: 3, ..Default::default() };
    let h = ultrafast::hybrid(fa.to_str().unwrap(), &opt).unwrap();
    let mut lib = Vec::new();
    h.write_bpaf(&mut lib).unwrap();
    assert_eq!(lib, bytes, "library and command disagree");
}
