//! Differential tests: the Rust binary must reproduce the C++ reference
//! binary byte for byte (exit code, stdout, stderr and every output file).
//!
//! The reference is the C++ ULTRA at commit 7947d0d, which is 1.2.2 plus the
//! bug fixes this port also makes (TravisWheelerLab/ULTRA pull request 95).
//! The binary is taken from `$ULTRA_REF_BIN`, or from
//! `upstream/build-fixed/ultra` under the crate root. When neither exists
//! the tests print a notice and pass, so the suite still runs on a machine
//! without a C++ build. Pointing `$ULTRA_REF_BIN` at a build of stock 1.2.2
//! shows which outputs the fixes change.
//!
//! To build the reference from the crate root:
//!
//! ```text
//! git clone https://github.com/TravisWheelerLab/ULTRA upstream/ULTRA-fixed
//! git -C upstream/ULTRA-fixed checkout 7947d0d
//! cmake -S upstream/ULTRA-fixed -B upstream/build-fixed -DCMAKE_BUILD_TYPE=Release
//! cmake --build upstream/build-fixed
//! ```
//!
//! The Rust binary builds its consensus differently by default (indel test
//! and realignment, see `RepeatRegion::create_logo`), so every run here sets
//! `ULTRA_UPSTREAM_CONSENSUS=1`, which makes it build the consensus as the
//! C++ does. The C++ binary does not read that variable. The unit tests in
//! `src/repeat.rs` test the default consensus, and `tests/determinism.rs`
//! checks that it is the same at every thread count.
//!
//! Each case runs both binaries in their own scratch directory with the input
//! linked in as `in.fa` and `argv[0]` forced to `ultra`, because the settings
//! block that ULTRA writes into its output echoes the command line.
//!
//! The C++ binary is only deterministic with one thread, so thread-count cases
//! compare the Rust binary at `-t N` against the reference at `-t 1`, looking
//! at output files only.

use std::collections::BTreeMap;
use std::fs;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};

// ---------------------------------------------------------------------------
// Locations
// ---------------------------------------------------------------------------

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The binary under test. `$ULTRA_TEST_BIN` overrides the freshly built Rust
/// binary, which lets the harness be checked by running the reference against
/// itself.
fn rust_bin() -> PathBuf {
    match std::env::var_os("ULTRA_TEST_BIN") {
        Some(p) => PathBuf::from(p),
        None => PathBuf::from(env!("CARGO_BIN_EXE_ultra")),
    }
}

fn ref_bin() -> Option<PathBuf> {
    static REF: OnceLock<Option<PathBuf>> = OnceLock::new();
    REF.get_or_init(|| {
        let candidate = match std::env::var_os("ULTRA_REF_BIN") {
            Some(p) => PathBuf::from(p),
            None => crate_root().join("upstream/build-fixed/ultra"),
        };
        if candidate.is_file() {
            Some(candidate)
        } else {
            eprintln!(
                "differential tests skipped: no reference binary at {} (set ULTRA_REF_BIN)",
                candidate.display()
            );
            None
        }
    })
    .clone()
}

fn work_root() -> PathBuf {
    PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("differential")
}

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

/// Deterministic generator for synthetic inputs (xorshift64*).
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
    fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.below(hi - lo + 1)
    }
    fn chance(&mut self, percent: u64) -> bool {
        self.below(100) < percent
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
}

/// Random background with tandem repeats of assorted periods spliced in. The
/// repeats carry substitutions and indels, some are lower-case, some
/// backgrounds end in a run of N, and some repeats change unit part-way
/// through so that repeat splitting has something to find.
fn synthetic_fasta(seed: u64, n_seq: usize, seq_len: usize, line_widths: &[usize]) -> Vec<u8> {
    const PERIODS: [u64; 16] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 21, 33, 47, 60, 95];
    let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut out = Vec::new();

    for s in 0..n_seq {
        let mut seq: Vec<u8> = Vec::with_capacity(seq_len + 4096);
        while seq.len() < seq_len {
            for _ in 0..rng.range(40, 2500) {
                let b = rng.base();
                seq.push(b);
            }
            if rng.chance(6) {
                for _ in 0..rng.range(1, 150) {
                    seq.push(b'N');
                }
            }

            let lower = rng.chance(20);
            let segments = if rng.chance(25) { 2 } else { 1 };
            for _ in 0..segments {
                let period = PERIODS[rng.below(PERIODS.len() as u64) as usize] as usize;
                let unit: Vec<u8> = (0..period).map(|_| rng.base()).collect();
                let copies = rng.range(3, 40) as usize;
                for k in 0..period * copies {
                    let mut b = unit[k % period];
                    let roll = rng.below(100);
                    if roll < 3 {
                        b = rng.base();
                    } else if roll < 4 {
                        continue; // deletion
                    } else if roll < 5 {
                        let extra = rng.base();
                        seq.push(if lower { extra.to_ascii_lowercase() } else { extra });
                    }
                    seq.push(if lower { b.to_ascii_lowercase() } else { b });
                }
            }
        }

        out.extend_from_slice(format!(">seq{s} synthetic sequence {s}\n").as_bytes());
        let width = line_widths[s % line_widths.len()];
        for chunk in seq.chunks(width) {
            out.extend_from_slice(chunk);
            out.push(b'\n');
        }
    }
    out
}

/// Hand-written input exercising FASTA parsing corners: CRLF line endings,
/// blank lines, comment lines, an empty record, a record too short to report,
/// duplicate names, mixed case, IUPAC codes, digits and spaces inside
/// sequence lines, and no newline at end of file.
fn edge_case_fasta() -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b">crlf record with windows line endings\r\n");
    for _ in 0..6 {
        out.extend_from_slice(b"ACGTTGCAGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATTACGGA\r\n");
    }
    out.extend_from_slice(b">empty\n");
    out.extend_from_slice(b">tiny\nACG\n");
    out.extend_from_slice(b">blank lines and comments\n\n");
    out.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n\n");
    out.extend_from_slice(b";a comment line in the middle of a record\n");
    out.extend_from_slice(b"ttagggttagggttagggttagggttagggttagggttagggttaggg\n");
    out.extend_from_slice(b">dup\nATATATATATATATATATATATATATATATATATCGCGCGCGCGCGCGCGCGCGCGCG\n");
    out.extend_from_slice(b">dup\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCC\n");
    out.extend_from_slice(b">iupac and junk\n");
    out.extend_from_slice(b"ACGTRYKMSWBDHVNACGT ACGT 1234 ACGT-ACGT*ACGT\n");
    out.extend_from_slice(b"CAGCAGCAGNNNNNNNNNNCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG\n");
    out.extend_from_slice(b"TGTGTGTGTGTGTGTGTGTGTGTGRRRRRTGTGTGTGTGTGTGTGTGTGTGTGTG\n");
    out.extend_from_slice(b">one very long line\n");
    let mut rng = Rng(0xC0FFEE);
    for i in 0..9000 {
        let b = if (2000..2600).contains(&i) {
            b"AAGGC"[i % 5]
        } else if (5000..5900).contains(&i) {
            b"ctt"[i % 3]
        } else {
            rng.base()
        };
        out.push(b);
    }
    out.extend_from_slice(b"\n>no trailing newline\n");
    out.extend_from_slice(b"GGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCC");
    out
}

/// Input for the masked-FASTA corners that upstream 1.2.2 gets wrong: a
/// record whose repeat runs to its last base followed by records with no
/// repeats, a repeat after IUPAC codes, digits and a comment line, and two
/// repeats separated only by a run of N.
fn mask_corner_fasta() -> Vec<u8> {
    let mut rng = Rng(0xBADC0DE);
    let mut random = |n: usize| -> Vec<u8> { (0..n).map(|_| rng.base()).collect() };
    let mut out = Vec::new();

    out.extend_from_slice(b">ends in a repeat\n");
    out.extend_from_slice(&random(120));
    out.push(b'\n');
    for _ in 0..3 {
        out.extend_from_slice(b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG\n");
    }
    out.extend_from_slice(b">no repeats here\n");
    for _ in 0..8 {
        out.extend_from_slice(&random(50));
        out.push(b'\n');
    }
    out.extend_from_slice(b">empty record\n>iupac digits and a comment before the repeat\n");
    out.extend_from_slice(&random(40));
    out.extend_from_slice(b"RYKMSWBDHV 0123456789 -*.\n;acgt acgt acgt acgt acgt\n");
    out.extend_from_slice(&random(30));
    out.push(b'\n');
    for _ in 0..3 {
        out.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n");
    }
    out.extend_from_slice(&random(80));
    out.extend_from_slice(b"\n>two repeats around an N run\n");
    out.extend_from_slice(&random(60));
    out.push(b'\n');
    out.extend_from_slice(b"TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\n");
    out.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n");
    out.extend_from_slice(b"AAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGC\n");
    out.extend_from_slice(&random(60));
    out.push(b'\n');
    out
}

/// Returns the path of a named input, creating it on first use.
fn input(key: &str) -> PathBuf {
    static LOCK: Mutex<()> = Mutex::new(());
    let _guard = LOCK.lock().unwrap();

    let dir = work_root().join("inputs");
    fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("{key}.fa"));
    if path.exists() {
        return path;
    }

    let data = crate_root().join("examples");
    let bytes = match key {
        "ex1" => fs::read(data.join("example_1.fa")).unwrap(),
        "ex2" => fs::read(data.join("example_2.fa")).unwrap(),
        "ex3" => fs::read(data.join("example_3.fa")).unwrap(),
        "edge" => edge_case_fasta(),
        "maskcorners" => mask_corner_fasta(),
        // Several sequences, mixed line widths.
        "mixed" => synthetic_fasta(1, 4, 40_000, &[60, 70, 80, 200]),
        // One sequence long enough to span many small windows.
        "long" => synthetic_fasta(2, 1, 150_000, &[80]),
        // Enough sequence to span several default-size windows.
        "big" => synthetic_fasta(3, 3, 400_000, &[60, 100]),
        other => panic!("unknown input {other}"),
    };
    // Write then rename so a concurrent reader never sees a partial file.
    let tmp = dir.join(format!("{key}.fa.tmp"));
    fs::write(&tmp, bytes).unwrap();
    fs::rename(&tmp, &path).unwrap();
    path
}

// ---------------------------------------------------------------------------
// Running and comparing
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum Compare {
    /// Exit code, stdout, stderr and output files.
    Everything,
    /// Exit code and output files only (used when the two command lines
    /// differ, e.g. thread-count cases).
    FilesOnly,
    /// Exit code plus stdout/stderr/files with every number blanked out (used
    /// for the shuffled-sequence paths, which are random by design).
    Shape,
}

struct Case {
    name: String,
    input: &'static str,
    /// Arguments after the input file, for the reference binary.
    ref_args: Vec<String>,
    /// Arguments after the input file, for the Rust binary.
    rust_args: Vec<String>,
    compare: Compare,
    /// Pass the input file argument at all (false for bare invocations).
    with_input: bool,
}

fn args(s: &str) -> Vec<String> {
    s.split_whitespace().map(str::to_owned).collect()
}

fn slug(s: &str) -> String {
    let cleaned: String = s
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    cleaned.trim_matches('_').to_owned()
}

/// A case where both binaries get the same arguments.
fn case(input: &'static str, arg_str: &str) -> Case {
    Case {
        name: format!("{input}__{}", slug(arg_str)),
        input,
        ref_args: args(arg_str),
        rust_args: args(arg_str),
        compare: Compare::Everything,
        with_input: true,
    }
}

struct Outcome {
    code: Option<i32>,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    files: BTreeMap<String, Vec<u8>>,
}

fn run(bin: &Path, dir: &Path, input_path: &Path, case_args: &[String], with_input: bool) -> Outcome {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    std::os::unix::fs::symlink(input_path, dir.join("in.fa")).unwrap();
    // A parameter file for `--tune_file` cases.
    fs::copy(crate_root().join("examples/tune_file"), dir.join("tune_file")).unwrap();

    let mut cmd = Command::new(bin);
    cmd.arg0("ultra")
        .current_dir(dir)
        .env("ULTRA_UPSTREAM_CONSENSUS", "1");
    if with_input {
        cmd.arg("in.fa");
    }
    cmd.args(case_args);
    let output = cmd.output().unwrap_or_else(|e| panic!("failed to run {}: {e}", bin.display()));

    let mut files = BTreeMap::new();
    for entry in fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().to_string_lossy().into_owned();
        if name == "in.fa" || name == "tune_file" {
            continue;
        }
        files.insert(name, fs::read(entry.path()).unwrap());
    }

    Outcome {
        code: output.status.code(),
        stdout: output.stdout,
        stderr: output.stderr,
        files,
    }
}

/// Blanks every run of digits and dots, and everything after the label of
/// lines whose value is random.
fn shape(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    let mut out = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("Estimated false discovery rate:") {
            let _ = rest;
            out.push_str("Estimated false discovery rate: #\n");
            continue;
        }
        if line.starts_with("Best coverage within FDR limit")
            || line.starts_with("No parameters found within FDR limit")
        {
            // Which of the two is printed depends on the shuffle.
            out.push_str("<tune verdict>\n");
            continue;
        }
        let mut in_number = false;
        for c in line.chars() {
            if c.is_ascii_digit() || (in_number && c == '.') {
                if !in_number {
                    out.push('#');
                    in_number = true;
                }
            } else {
                in_number = false;
                out.push(c);
            }
        }
        out.push('\n');
    }
    out.into_bytes()
}

fn describe_difference(what: &str, expected: &[u8], actual: &[u8]) -> String {
    let exp = String::from_utf8_lossy(expected);
    let act = String::from_utf8_lossy(actual);
    let clip = |s: &str| -> String {
        if s.len() > 300 {
            let mut end = 300;
            while !s.is_char_boundary(end) {
                end -= 1;
            }
            format!("{}... ({} bytes)", &s[..end], s.len())
        } else {
            s.to_owned()
        }
    };
    let mut exp_lines = exp.split('\n');
    let mut act_lines = act.split('\n');
    let mut line_no = 1;
    loop {
        match (exp_lines.next(), act_lines.next()) {
            (Some(e), Some(a)) if e == a => line_no += 1,
            (e, a) => {
                return format!(
                    "{what} differs at line {line_no} ({} vs {} bytes total)\n      c++ : {}\n      rust: {}",
                    expected.len(),
                    actual.len(),
                    e.map(clip).unwrap_or_else(|| "<end of output>".into()),
                    a.map(clip).unwrap_or_else(|| "<end of output>".into()),
                );
            }
        }
    }
}

/// Runs one case against both binaries. Returns a description of every
/// difference found, empty when the outputs agree.
fn check(case: &Case, reference: &Path) -> Vec<String> {
    let input_path = input(case.input);
    let base = work_root().join("runs").join(&case.name);
    let ref_dir = base.join("cpp");
    let rust_dir = base.join("rust");

    let expected = run(reference, &ref_dir, &input_path, &case.ref_args, case.with_input);
    let actual = run(&rust_bin(), &rust_dir, &input_path, &case.rust_args, case.with_input);

    let mut problems = Vec::new();
    if expected.code != actual.code {
        problems.push(format!(
            "exit code differs: c++ {:?}, rust {:?}",
            expected.code, actual.code
        ));
    }

    let normalise = |b: &[u8]| -> Vec<u8> {
        if case.compare == Compare::Shape {
            shape(b)
        } else {
            b.to_vec()
        }
    };

    if case.compare != Compare::FilesOnly {
        let (e, a) = (normalise(&expected.stdout), normalise(&actual.stdout));
        if e != a {
            problems.push(describe_difference("stdout", &e, &a));
        }
        let (e, a) = (normalise(&expected.stderr), normalise(&actual.stderr));
        if e != a {
            problems.push(describe_difference("stderr", &e, &a));
        }
    }

    let names: Vec<&String> = expected.files.keys().chain(actual.files.keys()).collect();
    let mut seen = std::collections::BTreeSet::new();
    for name in names {
        if !seen.insert(name.clone()) {
            continue;
        }
        match (expected.files.get(name), actual.files.get(name)) {
            (Some(e), Some(a)) => {
                let (e, a) = (normalise(e), normalise(a));
                if e != a {
                    problems.push(describe_difference(&format!("file {name}"), &e, &a));
                }
            }
            (Some(_), None) => problems.push(format!("file {name} written by c++ only")),
            (None, Some(_)) => problems.push(format!("file {name} written by rust only")),
            (None, None) => unreachable!(),
        }
    }

    if problems.is_empty() {
        let _ = fs::remove_dir_all(&base);
    } else {
        problems.insert(
            0,
            format!(
                "c++ : ultra {}{}\n    rust: ultra {}{}\n    kept: {}",
                if case.with_input { "in.fa " } else { "" },
                case.ref_args.join(" "),
                if case.with_input { "in.fa " } else { "" },
                case.rust_args.join(" "),
                base.display()
            ),
        );
    }
    problems
}

/// Runs a batch of cases a few at a time and fails with a report of every
/// case that differed.
fn run_cases(cases: Vec<Case>) {
    let Some(reference) = ref_bin() else {
        return;
    };

    let next = AtomicUsize::new(0);
    let failures: Mutex<Vec<(String, Vec<String>)>> = Mutex::new(Vec::new());
    let workers = cases.len().min(8).max(1);

    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                let Some(case) = cases.get(i) else { break };
                let problems = check(case, &reference);
                if !problems.is_empty() {
                    failures.lock().unwrap().push((case.name.clone(), problems));
                }
            });
        }
    });

    let mut failures = failures.into_inner().unwrap();
    failures.sort();
    if !failures.is_empty() {
        let mut report = format!("{} of {} cases differ from the C++ reference\n", failures.len(), cases.len());
        for (name, problems) in &failures {
            report.push_str(&format!("\n[{name}]\n"));
            for p in problems {
                report.push_str(&format!("    {p}\n"));
            }
        }
        panic!("{report}");
    }
}

/// Builds the cross product of inputs and argument strings.
fn matrix(inputs: &[&'static str], arg_strs: &[&str]) -> Vec<Case> {
    let mut cases = Vec::new();
    for input in inputs {
        for a in arg_strs {
            cases.push(case(input, a));
        }
    }
    cases
}

// ---------------------------------------------------------------------------
// Cases
// ---------------------------------------------------------------------------


#[test]
fn default_run() {
    run_cases(matrix(&["ex1", "ex2", "ex3", "edge", "mixed"], &[""]));
}

#[test]
fn output_formats() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "--tsv",
            "--json",
            "--bed",
            "-o out.txt",
            "--tsv -o out.tsv",
            "--json -o out.json",
            "--bed -o out.bed",
            "--json --bed -o out",
            "--json --bed --tsv -o out",
            "--json --tsv",
            "--hs",
            "--hs --json",
            "--hs --bed -o out.bed",
            "--suppress",
            "--suppress --json -o out.json",
            "--disable_summary",
            "--disable_summary --hs --json",
            "--disable_streaming_out",
            "--disable_streaming_out --json -o out.json",
        ],
    ));
}

#[test]
fn output_detail_flags() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "-c",
            "--show_counts --bed",
            "--show_seq",
            "--show_seq --json",
            "--show_seq --bed",
            "--json --show_delta",
            "--json --show_trace",
            "--json --show_wid",
            "--json --show_logo",
            "--json -c --show_seq --show_delta --show_trace --show_wid --show_logo",
            "-c --show_seq --show_delta --show_trace --show_wid --show_logo",
            "--max_consensus 5",
            "--max_consensus 5 --json",
            "--max_consensus 0 --bed",
            "--max_consensus 0",
            "--max_consensus 0 -c --show_seq",
            "--max_consensus -1",
            "--pval",
            "--pval --bed",
            "--pval --json",
            "--pval --bed --pval_loc 3.5 --pval_scale 2.25",
        ],
    ));
}

#[test]
fn masking() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed", "maskcorners"],
        &[
            "--mask masked.fa",
            "--mask masked.fa --nmask",
            "--mask masked.fa --hs --suppress",
            "--mask masked.fa -p 10 -s 30",
            "--mask masked.fa --json -o out.json",
            "--mask masked.fa --no_split --min_length 40",
        ],
    ));
}

#[test]
fn filters() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "-s 20",
            "--score 75.5 --json",
            "-s -5",
            "--min_length 30",
            "--min_length 1",
            "--min_unit 5",
            "--min_unit 1 --min_length 2",
            "-s 15 --min_length 25 --min_unit 4 --json",
        ],
    ));
}

#[test]
fn model_size() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "-p 1",
            "-p 2",
            "-p 3 --json",
            "-p 10",
            "-p 30 -i 3 -d 5",
            "-p 200",
            "-p 300 -i 2 -d 2 --json",
            "-i 0 -d 0",
            "-i 0 -d 4",
            "-i 4 -d 0",
            "-i 1 -d 1 --json",
            "--period 50 --inserts 20 --deletes 20",
            "-p 6 -i 6 -d 6",
            "-p 4 -i 10 -d 10",
        ],
    ));
}

#[test]
fn probabilities() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "--at 0.3",
            "--at 0.72 --json",
            "--acgt 0.3 0.2 0.2 0.3",
            "--acgt 0.1 0.4 0.4 0.1 --json",
            "-m 0.8",
            "--match 0.55",
            "-m 0.95 --json",
            "--decay 0.7",
            "--decay 0.99",
            "--nr 0.001",
            "--rn 0.2",
            "--nr 0.05 --rn 0.01",
            "--ri 0.001 --rd 0.001",
            "--ri 0.1 --rd 0.08 --ii 0.3 --dd 0.25",
            "--ii 0.001 --dd 0.5",
            "-m 0.6 --at 0.6 --rn 0.001 --nr 0.005",
        ],
    ));
}

#[test]
fn splitting() {
    run_cases(matrix(
        &["ex1", "ex2", "ex3", "edge", "mixed"],
        &[
            "--no_split",
            "--no_split --json",
            "--max_split 3",
            "--max_split 50 --json",
            "--max_split 0",
            "--split_threshold 2",
            "--split_threshold 0.5 --json",
            "--split_threshold 9",
            "--split_depth 2",
            "--split_depth 12 --json",
            "--min_split_window 5",
            "--min_split_window 60",
            "--max_split 20 --split_threshold 1.5 --split_depth 3 --min_split_window 8 --json --show_seq",
        ],
    ));
}

#[test]
fn windowing() {
    let mut cases = matrix(
        &["ex2", "ex3", "edge", "mixed", "long"],
        &[
            "-r",
            "--read_all --json",
            "--windows 1",
            "--windows 2 --json",
            "--windows 3 --win_size 3000 --overlap 400 -p 20",
            "--win_size 2000 --overlap 300 -p 20",
            "--win_size 5000 --overlap 1000 -p 40 --json",
            "--win_size 1000 --overlap 200 -p 10 -i 2 -d 2",
            "--win_size 1500 --overlap 500 -p 12 --json --show_seq --show_wid",
            "--win_size 800 --overlap 150 -p 6 --mask masked.fa",
            "--win_size 20000 --overlap 2000",
            "--win_size 4000 --overlap 600 -p 30 --disable_streaming_out --json",
            "--win_size 4000 --overlap 600 -p 30 -r --bed",
            "--overlap 50 -p 10",
            "--win_size 500 -p 10",
        ],
    );
    // A few default-size windows end to end.
    cases.extend(matrix(&["big"], &["-p 30 -i 3 -d 3", "-p 30 -i 3 -d 3 --json --mask masked.fa"]));
    run_cases(cases);
}

#[test]
fn informational_and_errors() {
    let mut cases = matrix(
        &["ex1"],
        &[
            "--mem",
            "--mem -p 500 -i 20 -d 20 -t 8",
            "--mem -p 5 -i 0 -d 0 --windows 16",
            "--mem --win_size 5000 --overlap 500",
            "--cite",
        ],
    );

    for (name, a, with_input) in [
        ("no_arguments", "", false),
        ("missing_file", "does_not_exist.fa", false),
        ("missing_file_json", "does_not_exist.fa --json -o out.json", false),
    ] {
        cases.push(Case {
            name: name.to_owned(),
            input: "ex1",
            ref_args: args(a),
            rust_args: args(a),
            compare: Compare::Everything,
            with_input,
        });
    }

    // Help text and rejected or unusual arguments. The parser mirrors CLI11,
    // so usage and error messages are compared in full as well.
    for a in [
        "--help",
        "-h",
        "--no_such_flag",
        "-p notanumber",
        "-p",
        "--acgt 0.3 0.2",
        "-t 0",
        "--at 1.5",
        "-m 0",
        "-p 0",
        "--win_size 10 --overlap 100",
        "extra_positional.fa",
        "-p 5 -p 6",
        "--pval_loc",
        "--pval_loc=3.5 --pval",
        "--acgt 0.25 0.25 0.25 0.25 --json",
        "--acgt 0.1 0.1 0.1 0.1",
        "--acgt 0.4 0.4 0.4 0.4",
        "--acgt 0.25 0.25 0.25 0.2495",
        "--mask no_such_directory/masked.fa",
        "-p5",
        "-t 2 --windows 1",
        "--tune_file missing_tune_file --tune_only",
        "--mask",
        "-o",
        "-s abc",
        "--min_length -3",
    ] {
        cases.push(case("ex1", a));
    }
    run_cases(cases);
}

#[test]
fn threads_match_single_threaded_reference() {
    let mut cases = Vec::new();
    for input in ["ex2", "ex3", "edge", "mixed", "long"] {
        for rest in [
            "--hs -o out.txt",
            "--hs --json -o out.json",
            "--hs --bed --mask masked.fa -o out.bed",
            "--hs --win_size 2000 --overlap 300 -p 20 --json --show_seq -o out.json",
            "--hs --win_size 1000 --overlap 200 -p 10 --windows 16 -o out.txt",
            "--hs --disable_streaming_out --win_size 3000 --overlap 500 -p 25 -o out.txt",
        ] {
            for threads in [2, 4, 16] {
                cases.push(Case {
                    name: format!("{input}__t{threads}__{}", slug(rest)),
                    input,
                    ref_args: args(&format!("-t 1 {rest}")),
                    rust_args: args(&format!("-t {threads} {rest}")),
                    compare: Compare::FilesOnly,
                    with_input: true,
                });
            }
        }
    }
    run_cases(cases);
}

#[test]
fn shuffled_paths_have_the_same_shape() {
    let mut cases = Vec::new();
    for (input, a) in [
        ("ex1", "--fdr --hs"),
        ("ex2", "--fdr"),
        ("ex1", "--tune --tune_only --hs"),
        ("ex1", "--tune_file tune_file --tune_only --hs"),
        ("ex1", "--tune_file tune_file --tune_indel --tune_only --hs --tune_fdr 0.5"),
    ] {
        let mut c = case(input, a);
        c.compare = Compare::Shape;
        cases.push(c);
    }
    run_cases(cases);
}
