//! End-to-end checks of bugs fixed in this port (the C++ has them): joins
//! across window boundaries, output that depended on when repeats were
//! written, masking before the first header, tuning runs emptying the -o
//! file, a closed stdout, and non-UTF-8 arguments. Needs no C++ binary.

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
fn rows(input: &Path, args: &[&str]) -> Vec<Vec<String>> {
    let out = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .args(["-t", "1", "--tsv", "-c", "--hs"])
        .args(args)
        .arg(input)
        .output()
        .unwrap();
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

    let whole = rows(&input, &[]);
    let split = rows(&input, &["--win_size", "1000", "--overlap", "100"]);
    assert_eq!(whole.len(), 1, "{whole:?}");
    assert_eq!(split.len(), 1, "{split:?}");
    // Start, end, period and the substitution / insertion / deletion counts.
    for col in [1, 2, 3, 7, 8, 9] {
        assert_eq!(split[0][col], whole[0][col], "column {col}: {split:?} vs {whole:?}");
    }
    let (a, b): (f64, f64) = (whole[0][4].parse().unwrap(), split[0][4].parse().unwrap());
    assert!((a - b).abs() < 0.05 * a, "scores {a} and {b}");
}

/// What ultra writes to stdout for `args` and `input`.
fn ultra_stdout(args: &[&str], input: &Path) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .args(args)
        .arg(input)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// A repeat that starts in a window overlap is found by both windows, with
/// one start. The two pieces tie in the order the join takes repeats in, and
/// the join treated its first and second piece differently, so the record
/// depended on how the sort happened to arrange them. That differed between
/// a run that writes repeats as it goes and one that writes them all at the
/// end. The test compares four run orders.
/// `-r` was refused as needing too much memory: the estimate multiplied by
/// its window count, which has no limit. The repeats are those of a streamed
/// run; only the memory estimate printed above them differs.
#[test]
fn read_all_finds_the_same_repeats_as_streaming() {
    fn repeats(args: &[&str]) -> String {
        let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/example_1.fa");
        let out = ultra_stdout(args, &input);
        let (_, rows) = out.split_once("============================\n").expect("no memory estimate");
        rows.to_string()
    }
    let streamed = repeats(&["-t", "1", "--hs"]);
    assert!(streamed.lines().count() > 2, "no repeats in:\n{streamed}");
    assert_eq!(repeats(&["-t", "1", "--hs", "-r"]), streamed);
    assert_eq!(repeats(&["-t", "1", "--hs", "--read_all", "--win_size", "2000"]), repeats(&["-t", "1", "--hs", "--win_size", "2000"]));
}

#[test]
fn output_does_not_depend_on_when_repeats_are_written() {
    let d = dir("run_order");
    let mut rng = Rng(0x5eed_1234_9876_0004);
    // Windows of 250 bases with an overlap of 40. Each array starts 30 bases
    // before a window boundary and ends 60 after it.
    let mut s = Vec::new();
    for boundary in (250..=750_000).step_by(250) {
        let filler = rng.bases(boundary - 30 - s.len());
        s.extend(filler);
        let unit_len = 2 + (rng.next() % 4) as usize;
        let unit = rng.bases(unit_len);
        for _ in 0..90 / unit_len {
            s.extend_from_slice(&unit);
        }
    }
    s.extend(rng.bases(100));
    let input = d.join("starts_in_overlaps.fa");
    fs::write(&input, [b">s\n".as_slice(), &s, b"\n"].concat()).unwrap();

    // The record lines (the lines before them name the thread count).
    let run = |order: &[&str]| -> Vec<String> {
        let fixed = ["--tsv", "-c", "--hs", "--win_size", "250", "--overlap", "40", "-p", "6"];
        let text = ultra_stdout(&[order, &fixed].concat(), &input);
        text.lines().filter(|l| l.starts_with("s\t")).map(str::to_string).collect()
    };
    let streamed = run(&["-t", "1"]);
    // A streaming run starts writing once it holds more than 2,000
    // repeats, so this one writes some before it has read the last window.
    assert!(streamed.len() > 2500, "{} records", streamed.len());
    for order in [
        &["-t", "1", "--disable_streaming_out"][..],
        &["-t", "4", "--windows", "4"],
        &["-t", "7", "--windows", "9", "--disable_streaming_out"],
    ] {
        assert!(run(order) == streamed, "output with {order:?} differs from -t 1");
    }
}

/// One record of `--json -c --show_trace --show_delta`.
struct JsonRecord {
    start: usize,
    length: usize,
    period: usize,
    score: f64,
    insertions: usize,
    deletions: usize,
    traceback: String,
    deltas: Vec<f64>,
}

impl JsonRecord {
    /// Score minus the sum of the score changes.
    fn score_beyond_deltas(&self) -> f64 {
        self.score - self.deltas.iter().sum::<f64>()
    }
}

fn json_records(input: &Path, args: &[&str]) -> Vec<JsonRecord> {
    let fixed = ["-t", "1", "--json", "-c", "--hs", "--show_trace", "--show_delta"];
    let text = ultra_stdout(&[&fixed, args].concat(), input);
    let mut records = Vec::new();
    for record in text.split("{\"SequenceName\"").skip(1) {
        // The value on the first line that starts with the key. The
        // subrepeats' keys do not start a line.
        let field = |key: &str| -> &str {
            let tag = format!("\n\"{key}\": ");
            let from = record.find(&tag).unwrap_or_else(|| panic!("no {key} in {record}")) + tag.len();
            let line = record[from..].lines().next().unwrap();
            line.trim_end_matches(|c| ",}]\"".contains(c)).trim_start_matches(|c| "[\"".contains(c))
        };
        records.push(JsonRecord {
            start: field("Start").parse().unwrap(),
            length: field("Length").parse().unwrap(),
            period: field("Period").parse().unwrap(),
            score: field("Score").parse().unwrap(),
            insertions: field("Insertions").parse().unwrap(),
            deletions: field("Deletions").parse().unwrap(),
            traceback: field("Traceback").to_string(),
            deltas: field("PositionScoreDeltas").split(',').map(|x| x.parse().unwrap()).collect(),
        });
    }
    records
}

/// A record joined across a window boundary describes one alignment: its
/// traceback shows the insertions and deletions it counts, and its score
/// differs from the sum of its score changes by as much as it does when one
/// window holds the whole repeat. The join used to stitch the counts and the
/// score at the end of the first piece, and the traceback and the score
/// changes at the start of the second.
#[test]
fn a_joined_record_describes_one_alignment() {
    let d = dir("one_alignment");
    let mut rng = Rng(0x5eed_1234_9876_0005);
    let mut s = rng.bases(300);
    for _ in 0..24 {
        let unit_len = 7 + (rng.next() % 9) as usize;
        let unit = rng.bases(unit_len);
        s.extend(rng.noisy(&unit, 600 / unit_len));
        s.extend(rng.bases(300));
    }
    let input = d.join("arrays.fa");
    fs::write(&input, [b">s\n".as_slice(), &s, b"\n"].concat()).unwrap();

    // The default window is longer than the sequence.
    let whole = json_records(&input, &[]);
    let split = json_records(&input, &["--win_size", "1000", "--overlap", "100"]);
    let mut joined_and_compared = 0;
    for r in &split {
        let at = format!("the record at {}", r.start);
        assert_eq!(r.traceback.len(), r.length, "{at}");
        assert_eq!(r.deltas.len(), r.length, "{at}");
        assert_eq!(r.traceback.matches('I').count(), r.insertions, "{at}: {}", r.traceback);
        assert_eq!(r.traceback.matches('D').count(), r.deletions, "{at}: {}", r.traceback);

        // Windows start every 1,000 bases and reach 100 back, so a record
        // that starts more than 100 before such a position and ends after
        // it was joined.
        let end = r.start + r.length;
        let joined = (1..).map(|k| k * 1000).take_while(|&b| b < end).any(|b| r.start + 100 < b);
        let same = |w: &&JsonRecord| (w.start, w.length, w.period) == (r.start, r.length, r.period);
        if let (true, Some(w)) = (joined, whole.iter().find(same)) {
            let (a, b) = (r.score_beyond_deltas(), w.score_beyond_deltas());
            assert!((a - b).abs() < 0.05, "{at}: {a} joined, {b} in one window");
            joined_and_compared += 1;
        }
    }
    assert!(joined_and_compared >= 8, "{joined_and_compared} joined records compared");
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

    let windowed = rows(&input, &["--win_size", "200", "--overlap", "100"]);
    let long: Vec<_> = windowed.iter().filter(|r| r[3] == "51").collect();
    assert_eq!(long.len(), 1, "one period-51 record: {windowed:?}");
    let whole = rows(&input, &[]);
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
