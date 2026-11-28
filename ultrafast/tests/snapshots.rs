//! Snapshot tests: ultrafast's exit status, stdout, stderr and the files it
//! writes must equal the files in `ultrafast/tests/snapshots/` byte for byte.
//! After an intended change in output, run
//! `UPDATE_SNAPSHOTS=1 cargo test -p ultrafast --test snapshots` to rewrite
//! the files and review the change with `git diff`.
//!
//! The inputs are those of `hybrid.rs` and of
//! `cases.rs::several_records_n_runs_and_subrepeats`. Each case runs with the
//! default thread count, and with `-t 1` and `-t 3`, which must give the same
//! bytes.

#[path = "../../tests/common/snapshot.rs"]
mod snapshot;
use snapshot::Outcome;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
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

fn fasta(recs: &[(&str, Vec<u8>)]) -> Vec<u8> {
    let mut s = Vec::new();
    for (name, seq) in recs {
        s.extend(format!(">{name}\n").bytes());
        for line in seq.chunks(80) {
            s.extend(line);
            s.push(b'\n');
        }
    }
    s
}

/// A period-3 array and a period-150 array in one record (`hybrid.rs`).
fn hybrid_input() -> Vec<u8> {
    let mut st = 0x1234_5678_9abc_def1u64;
    let mut seq = bases(&mut st, 2000);
    seq.extend(b"CAG".repeat(40));
    seq.extend(bases(&mut st, 2000));
    let unit = bases(&mut st, 150);
    for _ in 0..8 {
        seq.extend(&unit);
    }
    seq.extend(bases(&mut st, 2000));
    fasta(&[("s1", seq)])
}

/// Three records: random only; an N run then (CAG)x40; a period-40 array
/// whose unit changes half way (`cases.rs`).
fn multi_input() -> Vec<u8> {
    let mut st = 99u64;
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
    fasta(&[("r1", r1), ("r2", r2), ("r3", r3)])
}

const OPTIONS: [&str; 5] = ["", "-p 20", "--no-split", "--no-join", "-p 0"];

struct Case {
    name: String,
    input: fn() -> Vec<u8>,
    args: Vec<String>,
}

fn cases() -> Vec<Case> {
    let mut cases = Vec::new();
    for (input_name, input) in [("hybrid", hybrid_input as fn() -> Vec<u8>), ("multi", multi_input)] {
        for opts in OPTIONS {
            let slug: String = opts.chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
            cases.push(Case {
                name: format!("{input_name}__{}", if opts.is_empty() { "default" } else { slug.trim_matches('_') }),
                input,
                args: opts.split_whitespace().map(str::to_owned).collect(),
            });
        }
    }
    cases
}

fn snapshot_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/snapshots")
}

fn run(dir: &Path, input: &[u8], args: &[String]) -> Outcome {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    fs::write(dir.join("in.fa"), input).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_ultrafast"))
        .current_dir(dir)
        .args(args)
        .arg("in.fa")
        .output()
        .unwrap();
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().to_string_lossy().into_owned();
        if name != "in.fa" {
            files.insert(name, fs::read(entry.path()).unwrap());
        }
    }
    Outcome { code: output.status.code(), stdout: output.stdout, stderr: output.stderr, files }
}

/// The JSON is one line; a backslash and a line break after each record make
/// the stored file readable in a diff. JSON text cannot hold a backslash
/// followed by a raw line break, so the stored file still determines the
/// output exactly.
fn one_record_per_line(text: &str) -> String {
    text.replace("},{\"SequenceName\"", "},\\\n{\"SequenceName\"")
}

fn same(a: &Outcome, b: &Outcome) -> bool {
    a.code == b.code && a.stdout == b.stdout && a.stderr == b.stderr && a.files == b.files
}

#[test]
fn snapshot_ultrafast() {
    let cases = cases();
    let mut failures = Vec::new();
    for case in &cases {
        let base = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("snapshots").join(&case.name);
        let input = (case.input)();
        let command = format!("ultrafast {}in.fa", case.args.iter().map(|a| format!("{a} ")).collect::<String>());
        let outcome = run(&base.join("default"), &input, &case.args);
        let rendered = snapshot::render(&command, &outcome, &|t: &str| one_record_per_line(&snapshot::without_version(t)));
        if let Err(e) = snapshot::compare(&snapshot_dir(), &case.name, &command, &rendered) {
            failures.push(e);
        }
        for t in ["1", "3"] {
            let mut args = vec!["-t".to_owned(), t.to_owned()];
            args.extend(case.args.iter().cloned());
            if !same(&run(&base.join(format!("t{t}")), &input, &args), &outcome) {
                failures.push(format!("[{}] {command}\n    -t {t} gives different output from the default thread count", case.name));
            }
        }
    }
    snapshot::report(failures, cases.len());
}

#[test]
fn snapshot_files_all_have_a_case() {
    let names: BTreeSet<String> = cases().into_iter().map(|c| c.name).collect();
    snapshot::check_no_stale(&snapshot_dir(), &names);
}
