//! End-to-end checks of bugs fixed in this port (the C++ has them): joins
//! across window boundaries, masking before the first header, tuning runs
//! emptying the -o file, a closed stdout, and non-UTF-8 arguments. Needs no
//! C++ binary.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// xorshift64*, as in `determinism.rs`.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(self.next() % 4) as usize]).collect()
    }
    /// `copies` copies of `unit` with about 1 in 12 bases substituted,
    /// dropped or doubled.
    fn noisy(&mut self, unit: &[u8], copies: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for _ in 0..copies {
            for &b in unit {
                match self.next() % 36 {
                    0 => {}
                    1 => {
                        out.push(b);
                        out.push(b"ACGT"[(self.next() % 4) as usize]);
                    }
                    2 => out.push(b"ACGT"[(self.next() % 4) as usize]),
                    _ => out.push(b),
                }
            }
        }
        out
    }
}

fn dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("output_fixes").join(name);
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// The TSV rows (`-c --hs`), split into fields.
fn rows(input: &Path, args: &[&str], compat: bool) -> Vec<Vec<String>> {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_ultra"));
    cmd.args(["-t", "1", "--tsv", "-c", "--hs"]).args(args).arg(input);
    if compat {
        cmd.env("ULTRA_UPSTREAM_COMPAT", "1");
    } else {
        cmd.env_remove("ULTRA_UPSTREAM_COMPAT");
    }
    let out = cmd.output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("s\t"))
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

/// A noisy array that crosses a window boundary is joined into one record
/// with the counts and nearly the score it has inside one window.
#[test]
fn a_repeat_across_a_window_boundary_keeps_its_counts_and_score() {
    let d = dir("boundary");
    let mut rng = Rng(0x5eed_1234_9876_0001);
    let unit = rng.bases(10);
    let mut s = rng.bases(400);
    s.extend(rng.noisy(&unit, 150));
    s.extend(rng.bases(400));
    let input = d.join("cross.fa");
    fs::write(&input, [b">s\n".as_slice(), &s, b"\n"].concat()).unwrap();

    let whole = rows(&input, &[], false);
    let split = rows(&input, &["--win_size", "1000", "--overlap", "100"], false);
    assert_eq!(whole.len(), 1, "{whole:?}");
    assert_eq!(split.len(), 1, "{split:?}");
    // Start, end, period and the substitution / insertion / deletion counts.
    for col in [1, 2, 3, 7, 8, 9] {
        assert_eq!(split[0][col], whole[0][col], "column {col}: {split:?} vs {whole:?}");
    }
    let (a, b): (f64, f64) = (whole[0][4].parse().unwrap(), split[0][4].parse().unwrap());
    assert!((a - b).abs() < 0.05 * a, "scores {a} and {b}");

    // The C++ gives the joined record no counts and about half the score.
    let compat = rows(&input, &["--win_size", "1000", "--overlap", "100"], true);
    assert_eq!(compat.len(), 1);
    assert_eq!(&compat[0][7..10], ["0", "0", "0"]);
}

/// A repeat of another period inside a long-period repeat sorts between its
/// two pieces from neighbouring windows; the pieces are still joined.
#[test]
fn pieces_of_a_repeat_join_past_a_repeat_of_another_period() {
    let d = dir("other_period");
    let mut rng = Rng(0x5eed_1234_9876_0002);
    // A period-51 unit with a CAG run in its first copy region.
    let mut unit = rng.bases(15);
    unit.extend_from_slice(b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG");
    unit.extend(rng.bases(3));
    let mut s = Vec::new();
    for _ in 0..6 {
        s.extend_from_slice(&unit);
    }
    let input = d.join("sat.fa");
    fs::write(&input, [b">s\n".as_slice(), &s, b"\n"].concat()).unwrap();

    let windowed = rows(&input, &["--win_size", "200", "--overlap", "100"], false);
    let long: Vec<_> = windowed.iter().filter(|r| r[3] == "51").collect();
    assert_eq!(long.len(), 1, "one period-51 record: {windowed:?}");
    let whole = rows(&input, &[], false);
    let whole_long: Vec<_> = whole.iter().filter(|r| r[3] == "51").collect();
    assert_eq!(long[0][1..3], whole_long[0][1..3], "{windowed:?} vs {whole:?}");
}

/// Sequence before the first header is read as record 0 and its repeats
/// are reported, so --mask masks them.
#[test]
fn sequence_before_the_first_header_is_masked() {
    let d = dir("no_header");
    let mut rng = Rng(0x5eed_1234_9876_0003);
    let mut s = rng.bases(200);
    for _ in 0..30 {
        s.extend_from_slice(b"ACGTT");
    }
    s.extend(rng.bases(200));
    let input = d.join("nohdr.fa");
    fs::write(&input, [s.as_slice(), b"\n"].concat()).unwrap();
    let masked = d.join("masked.fa");
    let out = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .args(["-t", "1", "--hs", "--mask"])
        .arg(&masked)
        .arg(&input)
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    let row = text.lines().find(|l| l.starts_with('\t')).expect("a repeat in the headerless record");
    let f: Vec<&str> = row.split('\t').collect();
    let (start, end): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
    let m = fs::read(&masked).unwrap();
    let lower: Vec<usize> = (0..s.len()).filter(|&i| m[i].is_ascii_lowercase()).collect();
    assert_eq!(lower, (start..end).collect::<Vec<_>>());
}

/// Tuning runs must not open the -o file: with --tune_only it is not written
/// at all, so what was there stays.
#[test]
fn tuning_leaves_the_output_file_alone() {
    let d = dir("tune_only");
    let out = d.join("res.tsv");
    fs::write(&out, "precious\n").unwrap();
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/example_1.fa");
    let status = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .args(["-t", "1", "--tune_only", "-o"])
        .arg(&out)
        .arg(&input)
        .output()
        .unwrap();
    assert!(status.status.success());
    assert_eq!(fs::read_to_string(&out).unwrap(), "precious\n");
}

/// A reader that goes away (`ultra ... | head -1`) ends ultra quietly with
/// status 141, as SIGPIPE ends the C++, instead of a panic.
#[test]
fn a_closed_stdout_ends_quietly() {
    use std::io::Read;
    use std::process::Stdio;
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/example_2.fa");
    let mut child = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .args(["-t", "1", "--hs"])
        .arg(&input)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut first = [0u8; 1];
    child.stdout.take().unwrap().read_exact(&mut first).unwrap();
    // The pipe's read end is dropped here.
    let out = child.wait_with_output().unwrap();
    assert_eq!(out.status.code(), Some(141), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.stderr.is_empty(), "{}", String::from_utf8_lossy(&out.stderr));
}

/// A non-UTF-8 argument is an error, not silently replaced characters.
#[test]
fn non_utf8_arguments_are_rejected() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;
    let out = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .arg(OsStr::from_bytes(b"seq\xff.fa"))
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(255));
    assert!(String::from_utf8_lossy(&out.stderr).contains("is not valid UTF-8"));
}
