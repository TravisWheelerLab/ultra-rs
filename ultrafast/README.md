# ultrafast

A hybrid tandem-repeat annotator built on this repository's ULTRA library. It writes ULTRA's JSON format.

- **FasTAN** (a Rust port of Gene Myers' tandem-array finder, `src/fastan.rs`) finds arrays of any period
  (up to about 30 kb) quickly. Each array gets a TRF-style iterative consensus unit (`src/consensus.rs`)
  and is split where the unit changes.
- **ULTRA** runs in-process at a short maximum period (`-p`, default 20), alongside FasTAN.
  At short periods its sensitivity matters most and its cost is lowest.
- The two call sets are merged (see [Merging](#merging)).

ULTRA calls keep ULTRA's own units, `Score` and `PositionScoreDeltas`. FasTAN calls are scored per base on
the same scale (`--delta-model ultra`, bits), so the output works anywhere ULTRA's
`--json --show_delta` output does, e.g. as RepeatMasker Aurora's `-U` input.

## Building
From the repository root:
```
cargo build --release                  # builds ultra and ultrafast: target/release/ultrafast
cargo install --path ultrafast         # or install it
```
`ultrafast` is a separate crate in the workspace, so the `ultra` library keeps no dependencies.

## Usage
```
ultrafast -t 16 -o genome.json genome.fa                 # hybrid, ULTRA at periods <= 20
ultrafast -t 16 genome.fa > genome.json                  # the same, JSON on stdout
ultrafast -p 0 genome.fa > fastan_only.json              # FasTAN only (stdin allowed: '-')
ultrafast --tsv arrays.tsv -o genome.json genome.fa      # also FasTAN's arrays as TSV
ultrafast --filler ultra_p20.json -o genome.json genome.fa   # use a separate ULTRA run's calls
```
`-t N` is the total thread count (default: all cores), split between ULTRA and FasTAN, which run at the
same time; `-t 1` runs them one after the other. `ultrafast -h` lists every option. The main groups are:
- **Consensus:** `--rounds`, `--max-cons-period`, `--no-split`, `--split-min-gain`, `--split-min-frac`.
- **Score deltas:** `--delta-model te --matrices all.mat` scores FasTAN calls with Aurora's TE-alignment
  column scores (nats; experimental). `--delta-scale` multiplies every delta and `Score`.
- **Merging:**
  - `--filler-min-score` sets a score cutoff for ULTRA calls.
  - `--revote-filler` re-votes ULTRA's units with the TRF-style caller.
  - `--rescore-filler` recomputes ULTRA's deltas.
  - `--rescore-only` puts a `--filler` file on the `te` scale without running FasTAN.

Output is `{"Repeats": [...]}` with ULTRA's fields: `SequenceName`, `Start` (0-based), `Length`, `Period`,
`Score`, `Consensus`, `PositionScoreDeltas`, and `SubRepeats` for arrays split where the unit changes.
`--delta-model te` adds `Matrix`. Numbers have at most six decimals, as ULTRA writes them.
- `SequenceName` is the first word of the FASTA header (RepeatMasker's Aurora matches on it); `ultra
  --json` writes the whole header line. The names must be unique: a repeated name is an error.
- `Consensus` starts at the phase of the array's first copy, as ULTRA's units and the `SubRepeats` do.
- ULTRA's calls are taken from the library in memory; no intermediate files are written. ULTRA reads
  the FASTA itself, so stdin works only with `-p 0` or `--filler`. Input must be uncompressed, as for
  `ultra`.
- Any byte other than A/C/G/T (an N, or an IUPAC code such as R) ends a stretch of sequence for FasTAN,
  so a FasTAN array is reported in two parts on either side of one, as FasTAN does.

### Merging
- A FasTAN call with a period up to `-p` is dropped when ULTRA calls cover more than half of it.
- A FasTAN call with a longer period, which ULTRA cannot report, competes with the ULTRA calls it
  overlaps: its score under ULTRA's model is compared with the ULTRA calls' per-base deltas summed over
  the same bases. If it scores higher it is kept, and ULTRA calls it covers more than half of are dropped.
  Otherwise the first rule applies. This keeps a 21-23 bp unit that ULTRA would call at period 20 with an
  indel in every copy, and a satellite whose copies each hold a microsatellite (reported as the
  satellite). It also lets ULTRA win where FasTAN called a short, divergent array at a multiple of its
  period.
- A FasTAN call that is a fragment of one larger ULTRA call (the ULTRA call covers more than half of it,
  and it covers no more than half of the ULTRA call) never wins this way.
- Calls can still overlap in part. `--prefer-fastan` keeps every FasTAN call and drops ULTRA calls covered
  more than half by FasTAN calls. With `--filler`, the cap is the longest period in the file.

### Memory
The consensus and score-delta alignments keep a traceback of 2 bits per cell, in windows of at most 2^28
cells (64 MB) each; windows of one array run in parallel, so peak memory is about 64 MB per busy thread
plus the sequence. Units so long that three copies would pass 2^28 cells (about 9.5 kb) are aligned in
windows shorter than three units.

### Environment
- `ULTRA_UPSTREAM_COMPAT`: as for `ultra`, anything but empty or `0` runs ULTRA with the C++-compatible
  scoring.
- `ULTRAFAST_WINDOW` (FasTAN's scan window in bp, default 2,000,000; 0 = each stretch whole, as C FasTAN)
  and `ULTRAFAST_SPLIT_DEBUG` (prints split decisions to stderr) are for testing.

## As a library
The crate is also a library: the command is a thin wrapper around `ultrafast::hybrid`. It isn't on crates.io,
so depend on it through git or a path:
```toml
[dependencies]
ultrafast = { git = "https://github.com/TravisWheelerLab/ultra-rs", branch = "master" }
```
```rust
let opt = ultrafast::HybridOptions { ultra_period: 20, threads: 8, ..Default::default() };
let h = ultrafast::hybrid("genome.fa", &opt)?;
for r in &h.repeats {   // ULTRA's JSON records, sorted by sequence and start
    println!("{}\t{}\t{}\t{}\t{}", r.sequence_name, r.start, r.start + r.length, r.period, r.consensus);
}
h.write_json(std::io::stdout())?;   // the same JSON the command writes
```
`HybridOptions::default()` matches the command's defaults; its fields correspond to the command-line
options. `Hybrid` also holds the sequences read (`seqs`) and every FasTAN array before merging (`arrays`).
`hybrid` takes a file path because ULTRA reads the FASTA itself. The lower-level steps are public too:
`annotate` (FasTAN + consensus), `run_ultra`, `read_ultra_json` and `to_ultra_json` (merge and score).

## Performance
Benchmarked on the ULTRA paper's nine genomes and their shuffled copies, scored at periods ≤ 500 as in the
paper. Coverage is the fraction of the genome covered; false coverage is the same on the shuffled genome;
FDR is false coverage over coverage. "ultra" is this repository's `ultra -p 500` (0.2.x). All runs used one
thread per sequence; CPU is the total over all jobs.

| genome | ultrafast coverage (FDR) | ultra -p 500 coverage (FDR) |
|---|---|---|
| *C. trachomatis* | 0.66% (0.10) | 0.75% (0.09) |
| *P. torridus* | 2.68% (0.09) | 2.65% (0.09) |
| *T. aquaticus* | 3.79% (0.12) | 4.33% (0.10) |
| *H. marismortui* | 1.33% (0.12) | 1.42% (0.12) |
| *M. tuberculosis* | 4.15% (0.07) | 4.46% (0.07) |
| *S. coelicolor* | 6.37% (0.20) | 7.18% (0.19) |
| *P. falciparum* | 39.61% (0.29) | 40.95% (0.32) |
| *T. brucei* | 6.11% (0.01) | 6.25% (0.01) |
| *N. crassa* | 5.59% (0.03) | 5.49% (0.03) |
| total CPU | 1,749 s | 46,245 s |

- ultrafast covers 87-102% of what `ultra -p 500` covers. Its FDR is within 0.02 of ultra's: higher on six
  genomes, lower on three.
- It uses about 4% of the CPU time. On a synthetic file dense with long-unit arrays the gain is smaller
  (a reviewer measured 1 min 55 s against 4 min 15 s at `-t 8`): it depends on content.
- Against TRF and tantan on the same genomes, ultrafast covers 1.5-2.4 times what TRF covers and 1.6-6
  times what tantan covers. Its FDR is lower than both on every genome except *P. falciparum* (0.29,
  against 0.14 for TRF and 0.09 for tantan).
- On two simulated sets (640 arrays with units of 2-500 bp; 768 short arrays with units of 1-10 bp),
  ultrafast reports the exact unit for 54% and 68% of arrays. That is more than TRF, tantan, ULTRA 1.2.1
  or FasTAN alone. TRF reports the right period more often (85-88% of called arrays, against 78-84%),
  and TRF and FasTAN alone call fewer bases outside the arrays. Some "wrong" periods are multiples of the
  simulated unit in arrays where slippage duplicated several units at once, so the bases match better at
  the multiple (e.g. two near-identical 32 bp blocks of a diverged 8 bp unit); TRF, tantan and ULTRA
  without a period cap report those at the multiple too.

## Provenance and licence
- `src/fastan.rs` is a port of FasTAN 0.8 (`FasTAN.c`, and daligner's `align.c` wave aligner).
  Gene Myers' BSD-style notice is kept in the file.
- `src/consensus.rs` and `src/deltas.rs` are new code. The consensus vote follows the rule of TRF's
  `get_consensus` (Benson 1999; TRF 4.09 `tr30dat.c`): a column takes the plurality of A/C/G/T/gap, and
  an insertion is added where at least half the copies have one. That rule was reimplemented, not copied;
  the alignment and windowing around it differ from TRF's.
- The crate is BSD-3-Clause, like the rest of this repository.
