//! Perfect tandem arrays of long units must be reported at their own period.
//! Needs no C++ binary.
//!
//! Upstream 1.2.2 charged the bases of an indel window less than the same
//! bases in a match state, so for long units a path at period P - 1 (or
//! k * P - 1) with an insertion each copy outscored the true period P. These
//! arrays came out at periods 119 (or 239 and 479), 170 (or 341), 249 and
//! 499.

use std::fs;
use std::path::PathBuf;
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
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }
}

/// (unit length, copies) of each perfect array.
const ARRAYS: [(usize, usize); 5] = [(120, 5), (120, 100), (171, 20), (250, 5), (500, 5)];
const FLANK: usize = 1000;

#[test]
fn perfect_long_unit_arrays_get_their_own_period() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("long_period");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();

    let mut rng = Rng(0x5EED_0000_0000_0120);
    let mut fasta = Vec::new();
    for (i, &(unit_len, copies)) in ARRAYS.iter().enumerate() {
        let unit = rng.bases(unit_len);
        let mut seq = rng.bases(FLANK);
        for _ in 0..copies {
            seq.extend_from_slice(&unit);
        }
        seq.extend(rng.bases(FLANK));
        fasta.extend_from_slice(format!(">a{i}\n").as_bytes());
        for chunk in seq.chunks(70) {
            fasta.extend_from_slice(chunk);
            fasta.push(b'\n');
        }
    }
    let input = dir.join("in.fa");
    fs::write(&input, fasta).unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_ultra"))
        .current_dir(&dir)
        .arg(&input)
        .args(["-p", "500", "-t", "4", "-o", "out.tsv"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let tsv = fs::read_to_string(dir.join("out.tsv")).unwrap();
    for (i, &(unit_len, copies)) in ARRAYS.iter().enumerate() {
        let name = format!("a{i}");
        let array_len = unit_len * copies;
        // The call covering the most of the array.
        let best = tsv
            .lines()
            .skip(1)
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .filter(|f| f[0] == name)
            .map(|f| {
                let start: usize = f[1].parse().unwrap();
                let end: usize = f[2].parse().unwrap();
                let period: usize = f[3].parse().unwrap();
                let overlap = end.min(FLANK + array_len).saturating_sub(start.max(FLANK));
                (overlap, period)
            })
            .max()
            .unwrap_or_else(|| panic!("{name}: no repeat reported"));
        assert_eq!(
            best.1, unit_len,
            "{name}: {copies} copies of a {unit_len} bp unit reported at period {}",
            best.1
        );
        assert!(
            best.0 * 100 >= array_len * 95,
            "{name}: call covers {} of {array_len} bp",
            best.0
        );
    }
}
