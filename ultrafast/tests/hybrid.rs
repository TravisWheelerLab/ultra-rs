//! ultrafast runs ULTRA in-process: a short array comes back as an ULTRA call with
//! ULTRA's own unit, and a long-period array is called by FasTAN.

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
fn hybrid_calls_short_and_long_arrays() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("hybrid");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut st = 0x1234_5678_9abc_def1u64;
    let mut seq = bases(&mut st, 2000);
    seq.extend(b"CAG".repeat(40)); // 2000..2120, period 3: the filler's
    seq.extend(bases(&mut st, 2000));
    let unit = bases(&mut st, 150);
    for _ in 0..8 {
        seq.extend(&unit); // 4120..5320, period 150: FasTAN's
    }
    seq.extend(bases(&mut st, 2000));
    let fa = dir.join("in.fa");
    fs::write(&fa, [b">s1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let js = dir.join("out.json");
    let out = Command::new(env!("CARGO_BIN_EXE_ultrafast"))
        .args(["-t", "2", "-o"])
        .arg(&js)
        .arg(&fa)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(!dir.join("out.json.ultra.partial").exists(), "ULTRA's temporary JSON was left behind");

    let j: serde_json::Value = serde_json::from_str(&fs::read_to_string(&js).unwrap()).unwrap();
    let reps = j["Repeats"].as_array().unwrap();
    let hit = |lo: u64, hi: u64, period: u64| {
        reps.iter().any(|r| {
            let s = r["Start"].as_u64().unwrap();
            let e = s + r["Length"].as_u64().unwrap();
            r["Period"].as_u64() == Some(period) && s <= lo + 10 && e + 10 >= hi
        })
    };
    assert!(hit(2000, 2120, 3), "CAG array not reported at period 3: {reps:?}");
    assert!(hit(4120, 5320, 150), "150 bp array not reported at period 150");

    // the same records on stdout
    let out = Command::new(env!("CARGO_BIN_EXE_ultrafast")).args(["-t", "2"]).arg(&fa).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let j2: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j, j2, "stdout JSON differs from -o JSON");
}

#[test]
fn library_matches_command() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("hybrid_lib");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut st = 0x0fed_cba9_8765_4321u64;
    let mut seq = bases(&mut st, 1500);
    seq.extend(b"AT".repeat(50));
    seq.extend(bases(&mut st, 1500));
    let unit = bases(&mut st, 90);
    for _ in 0..10 {
        seq.extend(&unit);
    }
    seq.extend(bases(&mut st, 1500));
    let fa = dir.join("in.fa");
    fs::write(&fa, [b">s1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let opt = ultrafast::HybridOptions { threads: 2, ..Default::default() };
    let h = ultrafast::hybrid(fa.to_str().unwrap(), &opt).unwrap();
    assert!(h.n_ultra >= 1 && h.repeats.len() > h.n_ultra, "expected ULTRA and FasTAN calls: {:?}", h.repeats);
    let mut lib = Vec::new();
    h.write_json(&mut lib).unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_ultrafast")).args(["-t", "2"]).arg(&fa).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let a: serde_json::Value = serde_json::from_slice(&lib).unwrap();
    let b: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(a, b, "library and command disagree");
}
