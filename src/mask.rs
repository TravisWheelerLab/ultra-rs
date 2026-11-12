//! Masked FASTA output. Port of `mask.hpp` / `mask.cpp`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

use crate::std_sort::std_sort;

/// A half-open interval `[start, end)` of sequence positions to mask (C++
/// `mregion`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MRegion {
    pub start: u64,
    pub end: u64,
}

/// Ordering used to sort mask regions: by start, then by end (C++
/// `mregion::operator()`).
#[inline]
pub fn mregion_less(m1: &MRegion, m2: &MRegion) -> bool {
    if m1.start != m2.start {
        m1.start < m2.start
    } else {
        m1.end < m2.end
    }
}

/// Sorts `mask` in place and returns a copy with overlapping and touching
/// regions merged. Returns `None` for an empty input, where the C++ returns
/// a null pointer.
pub fn cleaned_masks(mask: &mut [MRegion]) -> Option<Vec<MRegion>> {
    if mask.is_empty() {
        return None;
    }

    std_sort(mask, mregion_less);

    // Combine overlapped regions.
    let mut cmask: Vec<MRegion> = vec![mask[0]];
    for m in &mask[1..] {
        let back = cmask.last_mut().expect("cmask starts non-empty");
        if m.start <= back.end {
            if m.end > back.end {
                back.end = m.end;
            }
        } else {
            cmask.push(*m);
        }
    }

    Some(cmask)
}

/// Copies the FASTA at `in_path` to `out_file`, masking the regions in
/// `masks` (keyed by sequence id; the first record in the file has id 1).
/// Masked `ACGT` become lower case, or `n` when `n_mask` is set. Other
/// letters inside a region (`N`, IUPAC codes) are left as they are.
///
/// The vectors in `masks` are sorted in place, as in the C++, where the map
/// holds pointers to them.
///
/// Sequence positions are counted the way the reader counts them
/// (`SequenceWindow::read_line`): every ASCII alphanumeric byte is one
/// position, anything else is skipped, and `;` comment lines are not
/// sequence. Each record is masked with its own regions only.
///
/// This departs from upstream 1.2.2 in three places, all upstream bugs:
/// - Upstream advanced the position only on `ACGTacgt` and `Nn` and treated
///   `;` lines as sequence, so masks drifted after any other letter, digit
///   or comment line.
/// - Upstream kept the previous record's regions when a record had none of
///   its own, so a record with no repeats could be masked with its
///   predecessor's coordinates.
/// - Upstream did not test a base against the next region when the position
///   had passed the end of the current one, which left the first base of a
///   region unmasked when `N`s separated it from the previous region.
///
/// Output lines always end in `\n`, so a missing final newline is added;
/// `\r` is passed through.
pub fn output_masked_fasta<W: Write>(
    in_path: &str,
    out_file: &mut W,
    masks: &mut HashMap<u64, Vec<MRegion>>,
    n_mask: bool,
) -> io::Result<()> {
    let file = match File::open(in_path) {
        Ok(f) => f,
        Err(_) => {
            eprintln!("Unable to open {in_path}");
            return Ok(());
        }
    };
    let mut in_file = BufReader::new(file);

    let mut line: Vec<u8> = Vec::new();
    let mut seq_id: u64 = 0;
    let mut seq_pos: u64 = 0;
    let mut mask_i: usize = 0;
    // Merged regions of the current record; empty when it has none.
    let mut cmask: Vec<MRegion> = Vec::new();

    loop {
        // `std::getline`: split on '\n' only; a read error ends the loop
        // like end of file.
        line.clear();
        match in_file.read_until(b'\n', &mut line) {
            Ok(0) | Err(_) => break,
            Ok(_) => {}
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }

        if line.first() == Some(&b'>') {
            seq_pos = 0;
            seq_id = seq_id.wrapping_add(1);
            mask_i = 0;
            cmask = masks.get_mut(&seq_id).and_then(|regions| cleaned_masks(regions)).unwrap_or_default();
        } else if line.first() != Some(&b';') && mask_i < cmask.len() {
            for c in line.iter_mut() {
                if !c.is_ascii_alphanumeric() {
                    continue;
                }
                // Regions are sorted and disjoint, so those that end at or
                // before this position are finished for good.
                while mask_i < cmask.len() && cmask[mask_i].end <= seq_pos {
                    mask_i += 1;
                }
                if mask_i < cmask.len()
                    && cmask[mask_i].start <= seq_pos
                    && matches!(c.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                {
                    *c = if n_mask { b'n' } else { c.to_ascii_lowercase() };
                }
                seq_pos += 1;
            }
        }

        out_file.write_all(&line)?;
        out_file.write_all(b"\n")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn region(start: u64, end: u64) -> MRegion {
        MRegion { start, end }
    }

    struct TempFasta(PathBuf);

    impl TempFasta {
        fn new(tag: &str, contents: &[u8]) -> Self {
            let path = std::env::temp_dir().join(format!("ultra_mask_test_{}_{tag}.fa", std::process::id()));
            std::fs::write(&path, contents).unwrap();
            TempFasta(path)
        }
        fn path(&self) -> &str {
            self.0.to_str().unwrap()
        }
    }

    impl Drop for TempFasta {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn mask(contents: &[u8], masks: &[(u64, &[(u64, u64)])], n_mask: bool, tag: &str) -> String {
        let input = TempFasta::new(tag, contents);
        let mut map: HashMap<u64, Vec<MRegion>> = HashMap::new();
        for (id, regions) in masks {
            map.insert(*id, regions.iter().map(|&(s, e)| region(s, e)).collect());
        }
        let mut out = Vec::new();
        output_masked_fasta(input.path(), &mut out, &mut map, n_mask).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn cleaned_masks_sorts_and_merges() {
        let mut m = vec![region(50, 60), region(10, 20), region(15, 30), region(30, 35), region(10, 12)];
        let cleaned = cleaned_masks(&mut m).unwrap();
        assert_eq!(cleaned, vec![region(10, 35), region(50, 60)]);
        // The input is left sorted, as the C++ sorts through the pointer.
        assert_eq!(m[0], region(10, 12));
        assert_eq!(m[4], region(50, 60));
        assert_eq!(cleaned_masks(&mut []), None);
    }

    #[test]
    fn masks_by_case_and_by_n() {
        let fasta = b">s1 first\nACGTACGTAC\nGTACGTACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(2, 5), (8, 13)])];
        assert_eq!(mask(fasta, masks, false, "case"), ">s1 first\nACgtaCGTac\ngtaCGTACGT\n");
        assert_eq!(mask(fasta, masks, true, "n"), ">s1 first\nACnnnCGTnn\nnnnCGTACGT\n");
    }

    #[test]
    fn sequence_ids_start_at_one_and_unmasked_records_pass_through() {
        let fasta = b">a\nACGTACGT\n>b\nACGTACGT\n>c\nACGTACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(2, &[(0, 3)])];
        assert_eq!(mask(fasta, masks, false, "ids"), ">a\nACGTACGT\n>b\nacgTACGT\n>c\nACGTACGT\n");
    }

    #[test]
    fn line_endings_and_missing_final_newline() {
        let fasta = b">a\r\nACGT\r\n\r\nACGT";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(3, 5)])];
        assert_eq!(mask(fasta, masks, false, "eol"), ">a\r\nACGt\r\n\r\naCGT\n");
    }

    #[test]
    fn every_alphanumeric_byte_is_one_position() {
        // R, 7 and N each take a position, as in the reader; the space and
        // the dash do not. Only ACGT are changed inside the region.
        let fasta = b">a\nARC7G NT-ACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(1, 8)])];
        assert_eq!(mask(fasta, masks, false, "pos"), ">a\nARc7g Nt-aCGT\n");
        assert_eq!(mask(fasta, masks, true, "posn"), ">a\nARn7n Nn-nCGT\n");
    }

    #[test]
    fn comment_lines_are_not_sequence() {
        let fasta = b">a\nACGT\n;acgt comment\nACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(3, 6)])];
        assert_eq!(mask(fasta, masks, true, "comment"), ">a\nACGn\n;acgt comment\nnnGT\n");
    }

    #[test]
    fn masks_do_not_leak_into_the_next_record() {
        // Record 1's region runs to the end of its sequence and record 2 has
        // no masks. Upstream 1.2.2 masked record 2 with record 1's region.
        let fasta = b">a\nACGTACGT\n>b\nACGTACGTACGT\n>c\nACGTACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(4, 8)]), (3, &[(0, 2)])];
        assert_eq!(mask(fasta, masks, false, "leak"), ">a\nACGTacgt\n>b\nACGTACGTACGT\n>c\nacGTACGT\n");
    }

    #[test]
    fn region_starting_right_after_an_n_run_is_masked_from_its_first_base() {
        // Upstream 1.2.2 left position 8 upper case here.
        let fasta = b">a\nACGTANNNACGTACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(0, 5), (8, 12)])];
        assert_eq!(mask(fasta, masks, false, "nrun"), ">a\nacgtaNNNacgtACGT\n");
    }

    #[test]
    fn region_past_the_end_of_the_sequence_is_harmless() {
        let fasta = b">a\nACGT\n>b\nACGT\n";
        let masks: &[(u64, &[(u64, u64)])] = &[(1, &[(2, 50), (60, 70)])];
        assert_eq!(mask(fasta, masks, false, "past"), ">a\nACgt\n>b\nACGT\n");
    }

    #[test]
    fn missing_input_writes_nothing() {
        let mut map = HashMap::new();
        let mut out = Vec::new();
        output_masked_fasta("/nonexistent/ultra/input.fa", &mut out, &mut map, false).unwrap();
        assert!(out.is_empty());
    }
}
