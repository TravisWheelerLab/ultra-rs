//! `-t N` must write exactly what `-t 1` writes. Needs no C++ binary.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// xorshift64*: a small deterministic generator for the test input.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
}

/// Random background with noisy tandem repeats, several records.
fn input_fasta() -> Vec<u8> {
    let mut rng = Rng(0x5EED_1234_ABCD_0001);
    let mut out = Vec::new();
    for s in 0..3 {
        let mut seq = Vec::new();
        while seq.len() < 30_000 {
            for _ in 0..(50 + rng.below(800)) {
                seq.push(rng.base());
            }
            let period = 1 + rng.below(12) as usize;
            let unit: Vec<u8> = (0..period).map(|_| rng.base()).collect();
            for k in 0..period * (4 + rng.below(30) as usize) {
                let b = if rng.below(100) < 4 {
                    rng.base()
                } else {
                    unit[k % period]
                };
                seq.push(b);
            }
        }
        out.extend_from_slice(format!(">seq{s}\n").as_bytes());
        for chunk in seq.chunks(70) {
            out.extend_from_slice(chunk);
            out.push(b'\n');
        }
    }
    out
}

fn work_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("determinism")
}

/// Runs ultra in a fresh directory and returns its output files (the
/// `.settings` file is left out: it echoes the thread count).
fn run(dir: &Path, input: &Path, threads: u32, extra: &str) -> Vec<(String, Vec<u8>)> {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .current_dir(dir)
        .arg(input)
        .args(["-t", &threads.to_string()])
        .args(extra.split_whitespace())
        .args([
            "--json",
            "--bed",
            "--tsv",
            "-o",
            "out",
            "--mask",
            "masked.fa",
        ])
        .output()
        .unwrap();
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let mut files: Vec<(String, Vec<u8>)> = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap())
        .map(|e| {
            (
                e.file_name().to_string_lossy().into_owned(),
                fs::read(e.path()).unwrap(),
            )
        })
        .filter(|(name, _)| !name.ends_with(".settings"))
        .collect();
    files.sort();
    files
}

#[test]
fn thread_count_does_not_change_output() {
    let root = work_dir();
    fs::create_dir_all(&root).unwrap();
    let input = root.join("in.fa");
    fs::write(&input, input_fasta()).unwrap();

    for (label, extra) in [
        // Many windows per batch.
        ("batch", "--win_size 1000 --overlap 150 -p 12 --windows 24"),
        // As many windows as threads: many batches.
        ("pool", "--win_size 800 --overlap 100 -p 8 --windows 8"),
    ] {
        let reference = run(&root.join(format!("{label}-t1")), &input, 1, extra);
        assert!(reference
            .iter()
            .any(|(n, d)| n == "out.tsv" && d.len() > 1000));
        for threads in [2, 3, 8] {
            for rep in 0..2 {
                let got = run(
                    &root.join(format!("{label}-t{threads}-{rep}")),
                    &input,
                    threads,
                    extra,
                );
                assert!(
                    got == reference,
                    "{label}: -t {threads} (run {rep}) output differs from -t 1"
                );
            }
        }
    }
}
