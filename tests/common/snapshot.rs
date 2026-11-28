//! Expected-output files for the snapshot tests (`tests/snapshots.rs` and
//! `ultrafast/tests/snapshots.rs`).
//!
//! One case is one run of a program. Its exit status, stdout, stderr and the
//! files it wrote are rendered as one text file with a labelled section each,
//! and the rendering must equal the stored file. With `UPDATE_SNAPSHOTS=1`
//! the stored file is rewritten instead.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

/// What one run produced.
pub struct Outcome {
    pub code: Option<i32>,
    pub stdout: Vec<u8>,
    pub stderr: Vec<u8>,
    /// Files the run wrote, by name.
    pub files: BTreeMap<String, Vec<u8>>,
}

/// Placeholder for the crate version in stored files, so that a version bump
/// does not change them.
pub const VERSION_MARK: &str = "@VERSION@";

/// Replaces the crate version in the places the programs print it: the
/// `--help` banner (`Version X`) and the JSON header (`"Version": "X"`).
pub fn without_version(text: &str) -> String {
    let v = env!("CARGO_PKG_VERSION");
    text.replace(&format!("Version {v}\n"), &format!("Version {VERSION_MARK}\n"))
        .replace(&format!("\"Version\": \"{v}\""), &format!("\"Version\": \"{VERSION_MARK}\""))
}

fn utf8<'a>(what: &str, bytes: &'a [u8]) -> &'a str {
    std::str::from_utf8(bytes).unwrap_or_else(|e| panic!("{what} is not UTF-8: {e}"))
}

fn section(out: &mut String, label: &str, text: &str) {
    out.push_str(&format!("==== {label} ====\n"));
    out.push_str(text);
    if !text.is_empty() && !text.ends_with('\n') {
        out.push_str("\n\\ no newline at end\n");
    }
}

/// Renders an outcome as the text stored for it. `normalise` is applied to
/// stdout, stderr and every file.
pub fn render(command: &str, o: &Outcome, normalise: &dyn Fn(&str) -> String) -> String {
    let status = o.code.map_or_else(|| "killed by a signal".to_owned(), |c| c.to_string());
    let mut out = format!("command: {command}\nstatus: {status}\n");
    section(&mut out, "stdout", &normalise(utf8("stdout", &o.stdout)));
    section(&mut out, "stderr", &normalise(utf8("stderr", &o.stderr)));
    for (name, bytes) in &o.files {
        section(&mut out, &format!("file {name}"), &normalise(utf8(name, bytes)));
    }
    out
}

pub fn updating() -> bool {
    std::env::var_os("UPDATE_SNAPSHOTS").is_some_and(|v| !v.is_empty() && v != "0")
}

/// Compares `actual` with the stored file `dir/name.txt` (or rewrites the file
/// when updating). Returns a description of the first difference.
pub fn compare(dir: &Path, name: &str, command: &str, actual: &str) -> Result<(), String> {
    let path = dir.join(format!("{name}.txt"));
    if updating() {
        fs::create_dir_all(dir).unwrap();
        if fs::read_to_string(&path).ok().as_deref() != Some(actual) {
            fs::write(&path, actual).unwrap();
        }
        return Ok(());
    }
    let expected = match fs::read_to_string(&path) {
        Ok(s) => s,
        Err(e) => {
            return Err(format!(
                "[{name}] {command}\n    no expected output at {} ({e}); run with UPDATE_SNAPSHOTS=1 to write it",
                path.display()
            ))
        }
    };
    if expected == actual {
        return Ok(());
    }
    let mut exp = expected.split('\n');
    let mut act = actual.split('\n');
    let mut line = 1;
    loop {
        match (exp.next(), act.next()) {
            (Some(e), Some(a)) if e == a => line += 1,
            (e, a) => {
                return Err(format!(
                    "[{name}] {command}\n    {} differs at line {line}\n    expected: {}\n    actual:   {}",
                    path.display(),
                    e.unwrap_or("<end of file>"),
                    a.unwrap_or("<end of file>"),
                ))
            }
        }
    }
}

/// Fails on stored files that belong to no case (deletes them when updating).
pub fn check_no_stale(dir: &Path, names: &BTreeSet<String>) {
    let mut stale = Vec::new();
    if updating() {
        fs::create_dir_all(dir).unwrap();
    }
    for entry in fs::read_dir(dir).unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display())) {
        let path = entry.unwrap().path();
        let file = path.file_name().unwrap().to_string_lossy().into_owned();
        let known = file.strip_suffix(".txt").is_some_and(|n| names.contains(n));
        if !known {
            if updating() {
                fs::remove_file(&path).unwrap();
            } else {
                stale.push(file);
            }
        }
    }
    assert!(
        stale.is_empty(),
        "files in {} with no case (delete them, or run with UPDATE_SNAPSHOTS=1): {stale:?}",
        dir.display()
    );
}

/// Panics with every failure, sorted.
pub fn report(mut failures: Vec<String>, total: usize) {
    failures.sort();
    if !failures.is_empty() {
        panic!(
            "{} of {total} cases differ from their expected output\n\n{}",
            failures.len(),
            failures.join("\n\n")
        );
    }
}
