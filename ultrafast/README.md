# ultrafast

A hybrid tandem-repeat annotator built on this repository's ULTRA library. It writes ULTRA's JSON format.

- **FasTAN** (a Rust port of Gene Myers' tandem-array finder, `src/fastan.rs`) finds arrays of any period
  (up to about 30 kb) quickly. Each array gets a TRF-style iterative consensus unit (`src/consensus.rs`)
  and is split where the unit changes.
- **ULTRA** runs in-process at a short maximum period (`-p`, default 20), on its own threads alongside FasTAN.
  At short periods its sensitivity matters most and its cost is lowest.
- The two call sets are merged. Where they overlap, ULTRA's call is kept and FasTAN calls covered more than
  half by ULTRA calls are dropped (`--prefer-fastan` reverses this).

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
`ultrafast -h` lists every option. The main groups are:
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
`--delta-model te` adds `Matrix`.

With `-o FILE`, ULTRA's intermediate JSON goes to `FILE.ultra.partial`; without `-o`, it goes to the
system temp directory. Either way it is deleted after reading. ULTRA needs a FASTA file, so stdin works
only with `-p 0` or `--filler`.

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
Benchmarked against the ULTRA paper's nine genomes and their shuffled copies, scored at periods ≤ 500.
Full ULTRA here means ULTRA with `-p 500`.

| | Coverage | False coverage | FDR | Total CPU |
|---|---|---|---|---|
| **ultrafast** | 85–100% of full ULTRA | below full ULTRA | at or below full ULTRA | 1,671 s |
| **Full ULTRA** | (baseline) | (baseline) | (baseline) | 46,246 s |
| **TRF / tantan** | ultrafast has 1.5–3× theirs | n/a | ultrafast's is at or below theirs, except *P. falciparum* vs TRF | n/a |

On simulated arrays, ultrafast called the exact unit for 54% of arrays in the first simulation set and 68%
in the second. That is the most of the tools tested: TRF, tantan, ULTRA 1.2.1 and FasTAN alone.

## Provenance and licence
- `src/fastan.rs` is a port of FasTAN 0.8 (`FasTAN.c`, and daligner's `align.c` wave aligner).
  Gene Myers' BSD-style notice is kept in the file.
- `src/consensus.rs` and `src/deltas.rs` are new code. The consensus vote follows the rule of TRF's
  `get_consensus` (Benson 1999; TRF 4.09 `tr30dat.c`): a column takes the plurality of A/C/G/T/gap, and
  an insertion is added where at least half the copies have one. That rule was reimplemented, not copied;
  the alignment and windowing around it differ from TRF's.
- The crate is BSD-3-Clause, like the rest of this repository.
