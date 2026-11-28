//! ultrafast: a hybrid tandem-repeat annotator. FasTAN (a Rust port of Gene
//! Myers' tandem-array finder) finds arrays of any period quickly; ULTRA (this
//! repository's library) is run at a short maximum period alongside it; the two
//! call sets are merged and written as ULTRA's JSON, with per-base score deltas.
//! FasTAN arrays get a TRF-style iterative consensus unit (`consensus.rs`).
//!
//! ```no_run
//! let recs = ultrafast::read_fasta("genome.fa").unwrap();
//! let calls = ultrafast::annotate(&recs, &ultrafast::Options { consensus: true, ..Default::default() });
//! ```

pub mod consensus;
pub mod wrap;
pub mod fastan;

use rayon::prelude::*;
use std::io::{BufRead, BufReader};

pub use consensus::{refine, split, Consensus, Params, SplitParams};
pub use fastan::{tandem_arrays, TandemArray};

pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read a (possibly multi-record) FASTA file; sequences are upper-cased.
///
/// Positions are counted as ULTRA's reader counts them, so that FasTAN's and
/// ULTRA's coordinates agree: lines starting with ';' are comments, only
/// ASCII letters and digits are sequence (gaps, '*', spaces and the like are
/// skipped), and sequence before the first header is a record with an empty
/// name. The name is the first word of the header; any bytes are accepted
/// (a header that is not UTF-8 is converted leniently).
pub fn read_fasta(path: &str) -> std::io::Result<Vec<Record>> {
    let mut rd: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        Box::new(BufReader::new(std::fs::File::open(path)?))
    };
    let mut recs: Vec<Record> = Vec::new();
    let mut line: Vec<u8> = Vec::new();
    loop {
        line.clear();
        if rd.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        match line.first() {
            Some(b'>') => {
                let header = String::from_utf8_lossy(&line[1..]);
                recs.push(Record { name: header.split_whitespace().next().unwrap_or("").to_string(), seq: Vec::new() });
            }
            Some(b';') => {}
            _ => {
                let mut bases = line.iter().filter(|c| c.is_ascii_alphanumeric()).map(|c| c.to_ascii_uppercase()).peekable();
                if bases.peek().is_none() {
                    continue;
                }
                if recs.is_empty() {
                    recs.push(Record { name: String::new(), seq: Vec::new() });
                }
                recs.last_mut().expect("a record").seq.extend(bases);
            }
        }
    }
    Ok(recs)
}

/// Fails if two records holding sequence have the same name. ULTRA's calls
/// are matched to records by name, and the output names each call's record,
/// so calls on two such records could not be told apart. A record without
/// sequence gets no calls and may share its name.
fn check_names(recs: &[Record]) -> anyhow::Result<()> {
    let mut seen = std::collections::HashMap::new();
    for (i, r) in recs.iter().enumerate().filter(|(_, r)| !r.seq.is_empty()) {
        if let Some(j) = seen.insert(r.name.as_str(), i) {
            let empty = if r.name.is_empty() {
                "; sequence before the first header, and a header with nothing after '>', have the empty name"
            } else {
                ""
            };
            anyhow::bail!(
                "sequence name {:?} occurs more than once (records {} and {}), so calls on them could not be told apart; \
                 names are the first word of each header{empty}",
                r.name, j + 1, i + 1
            );
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
pub struct Options {
    /// Run the TRF-style consensus caller on each array.
    pub consensus: bool,
    pub cons: Params,
    /// Arrays with a period above this keep the representative copy as unit.
    pub max_cons_period: usize,
    /// Split arrays where the unit changes (needs `consensus`); None = don't.
    pub split: Option<SplitParams>,
}

impl Default for Options {
    fn default() -> Self {
        Options { consensus: false, cons: Params::default(), max_cons_period: 10_000, split: Some(SplitParams::default()) }
    }
}

#[derive(Clone, Debug)]
pub struct Call {
    pub seq: usize,
    /// Array [start, end), 0-based.
    pub start: usize,
    pub end: usize,
    /// FasTAN's period (mean alignment diagonal).
    pub period: usize,
    /// Differences in FasTAN's self-alignment.
    pub diffs: usize,
    /// Consensus unit if computed, else one copy from the middle of the array.
    pub unit: Vec<u8>,
    /// Result of the consensus caller, if run.
    pub cons: Option<Consensus>,
    /// Where the unit changes: (offset from start, unit) per part, as ULTRA's
    /// SubRepeats. Empty when the array wasn't split.
    pub parts: Vec<(usize, Vec<u8>)>,
}

impl Call {
    /// FasTAN identity of the self-alignment (as tanbed reports it, /1000).
    pub fn fastan_identity(&self) -> f64 {
        (1.0 - self.diffs as f64 / (self.end - self.start).max(1) as f64).max(0.0)
    }
}

/// The consensus seed for a re-vote of an ULTRA call: the copy of length
/// `period` of A, C, G and T only nearest the middle of `[start, end)`, or
/// `fallback` (ULTRA's own unit) when there is none. ULTRA calls can span N
/// runs, and a middle copy with N gave refine() a short seed (a wrong, shorter
/// period) or an N-only one (an "NNNNN" unit).
fn acgt_seed(seq: &[u8], start: usize, end: usize, period: usize, fallback: &[u8]) -> Vec<u8> {
    let p = period.clamp(1, end - start);
    // bad[i]: non-ACGT bases in seq[start..start + i]
    let mut bad = Vec::with_capacity(end - start + 1);
    bad.push(0usize);
    for &c in &seq[start..end] {
        let n = bad[bad.len() - 1] + !matches!(c, b'A' | b'C' | b'G' | b'T') as usize;
        bad.push(n);
    }
    let clean = |s: usize| bad[s - start + p] == bad[s - start];
    let mid = start + (end - start - p) / 2;
    let last = end - p;
    for d in 0..=last - start {
        for s in [mid.checked_add(d), mid.checked_sub(d)].into_iter().flatten() {
            if s >= start && s <= last && clean(s) {
                return seq[s..s + p].to_vec();
            }
        }
    }
    fallback.to_vec()
}

/// One copy of length `period` centred in the array (the consensus seed).
pub fn middle_copy(seq: &[u8], start: usize, end: usize, period: usize) -> &[u8] {
    let p = period.clamp(1, end - start);
    let s = start + (end - start - p) / 2;
    &seq[s..s + p]
}

/// `split` with each part's unit rotated to the phase of the part's first copy;
/// empty when there's a single part.
pub fn split_parts(array: &[u8], unit: &[u8], par: &Params, sp: &SplitParams) -> Vec<(usize, Vec<u8>)> {
    let p = split(array, unit, par, sp);
    if p.len() < 2 {
        return Vec::new();
    }
    let ends: Vec<usize> = p.iter().skip(1).map(|x| x.0).chain([array.len()]).collect();
    p.into_iter().zip(ends).map(|((o, u), e)| (o, phase_to_start(u, &array[o..e]))).collect()
}

/// Find tandem arrays in all records and optionally call consensus units.
pub fn annotate(recs: &[Record], opt: &Options) -> Vec<Call> {
    let mut text = Vec::with_capacity(recs.iter().map(|r| r.seq.len()).sum());
    let mut starts = Vec::with_capacity(recs.len());
    for r in recs {
        starts.push(text.len());
        text.extend_from_slice(&r.seq);
    }
    let arrays = tandem_arrays(&text, &starts);
    // longest jobs first for load balance
    let mut order: Vec<usize> = (0..arrays.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse((arrays[i].aend - arrays[i].bbeg) * arrays[i].period.max(1)));
    let mut calls: Vec<(usize, Call)> = order
        .par_iter()
        .with_max_len(1)
        .map(|&i| {
            let a = &arrays[i];
            let seq = &recs[a.seq].seq;
            let seed = middle_copy(seq, a.bbeg, a.aend, a.period);
            let cons = if opt.consensus && a.period <= opt.max_cons_period {
                Some(refine(&seq[a.bbeg..a.aend], seed, &opt.cons))
            } else {
                None
            };
            // rotated to the phase of the array's first copy, as ULTRA's units and the SubRepeats are
            let unit = phase_to_start(
                cons.as_ref().map(|c| c.unit.clone()).unwrap_or_else(|| seed.to_vec()),
                &seq[a.bbeg..a.aend],
            );
            let parts = match (&cons, &opt.split) {
                (Some(_), Some(sp)) => split_parts(&seq[a.bbeg..a.aend], &unit, &opt.cons, sp),
                _ => Vec::new(),
            };
            (i, Call { seq: a.seq, start: a.bbeg, end: a.aend, period: a.period, diffs: a.diffs, unit, cons, parts })
        })
        .collect();
    calls.sort_by_key(|c| c.0);
    calls.into_iter().map(|c| c.1).collect()
}

// ------------------------------------------------------------- ULTRA in-process

/// Read ULTRA's JSON output (`--json`, with or without `--show_delta`).
pub fn read_ultra_json(path: &str) -> anyhow::Result<Vec<UltraRecord>> {
    let j: UltraJson = serde_json::from_reader(BufReader::new(std::fs::File::open(path)?))?;
    Ok(j.repeats)
}

/// Collects ULTRA's repeats as records instead of writing them: the fields and
/// numbers of ULTRA's `--json --show_delta` writer (`writers/json.rs`), with
/// Score and deltas rounded to ULTRA's six decimals. A record whose
/// SubRepeats text does not parse is kept without SubRepeats, with a warning
/// on stderr: the field is optional (Aurora ignores it) and the record's other
/// fields do not depend on it.
struct Collect(std::sync::Arc<std::sync::Mutex<Vec<UltraRecord>>>);

impl ultra_repeats::writers::RepeatFileWriter for Collect {
    fn write_repeat(
        &mut self,
        _settings: &ultra_repeats::cli::Settings,
        r: &ultra_repeats::repeat::RepeatRegion,
        _out: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        use ultra_repeats::writers::{to_string_f32, JSONFileWriter};
        let num = |x: f32| to_string_f32(x).parse::<f64>().unwrap_or(f64::NAN);
        let upto_nul = |b: &[u8]| String::from_utf8_lossy(b.split(|&c| c == 0).next().unwrap_or(b)).into_owned();
        let len = r.repeat_length as usize;
        let sub_repeats = if r.splits.as_ref().is_some_and(|s| s.iter().any(|&x| x > 0)) {
            let txt = JSONFileWriter::sub_repeats_string(r);
            serde_json::from_slice(&txt).unwrap_or_else(|e| {
                eprintln!(
                    "ultrafast: warning: ULTRA's SubRepeats for {} at {} did not parse ({e}); the call is kept without them",
                    upto_nul(r.sequence_name.as_bytes()),
                    r.sequence_start
                );
                Vec::new()
            })
        } else {
            Vec::new()
        };
        let rec = UltraRecord {
            sequence_name: upto_nul(r.sequence_name.as_bytes()),
            start: r.sequence_start as usize,
            length: len,
            period: r.repeat_period as usize,
            score: num(r.region_score),
            consensus: upto_nul(&r.string_consensus),
            position_score_deltas: r.scores.as_ref().map(|v| v[..len.min(v.len())].iter().map(|&x| num(x)).collect()).unwrap_or_default(),
            matrix: None,
            sub_repeats,
            scored_alignment: None,
        };
        self.0.lock().expect("collector lock").push(rec);
        Ok(())
    }
}

/// ULTRA's calls at its defaults with maximum period `period` on the FASTA file
/// `input`: what `ultra -p period -t threads --json --show_delta` would write,
/// taken from this repository's library in memory (no files are written and
/// nothing is printed but the warning described at `Collect`).
pub fn run_ultra(input: &str, period: usize, threads: usize) -> anyhow::Result<Vec<UltraRecord>> {
    use ultra_repeats::{cli::Settings, ultra::Ultra};
    // ULTRA exits the process on an unreadable file; fail here with a message instead
    std::fs::File::open(input).map_err(|e| anyhow::anyhow!("{input}: {e}"))?;
    // The path is set directly, not passed through ULTRA's command-line rules
    // (which refuse a file name starting with '-'), and parsing never exits
    // the process.
    let argv: Vec<String> = ["ultra", "-p", &period.to_string(), "-t", &threads.to_string(), "--json", "--show_delta", "--hs"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut s = Settings::new();
    s.prepare_settings();
    s.in_file = input.to_string();
    let mut err = String::new();
    if s.try_parse_input(&argv, &mut err) != Ok(true) {
        anyhow::bail!("ULTRA rejected its settings ({}): {}", argv.join(" "), err.trim_end());
    }
    s.assign_settings();
    let got = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    {
        let mut u = Ultra::new(&s);
        // the JSON writer would print to stdout; collect instead
        u.writers = vec![Box::new(Collect(got.clone()))];
        u.analyze_file()
            .and_then(|()| u.output_repeats(true))
            .map_err(|e| anyhow::anyhow!("ULTRA on {input}: {e}"))?;
    }
    let out = std::mem::take(&mut *got.lock().expect("collector lock"));
    Ok(out)
}

// ------------------------------------------------------------- JSON output

pub mod deltas;

/// How per-base score deltas are computed (see `deltas.rs`).
pub enum DeltaModel {
    /// ULTRA 1.2.1's default model, bits.
    Ultra,
    /// Aurora's TE-alignment column scores, nats. `fixed` names one matrix to
    /// use for every array; otherwise each array gets the `<div>p<gc>g` matrix
    /// nearest its divergence and its sequence's GC (or `gc` if given).
    Te { matrices: Vec<deltas::Matrix>, fixed: Option<String>, gc: Option<f64>, window: usize },
}

/// One record of ULTRA's JSON (`--json --show_delta`), field names as ULTRA's.
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
#[serde(rename_all = "PascalCase")]
pub struct UltraRecord {
    pub sequence_name: String,
    /// 0-based, as ULTRA's.
    pub start: usize,
    pub length: usize,
    pub period: usize,
    pub score: f64,
    pub consensus: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub position_score_deltas: Vec<f64>,
    /// ultrafast only: matrix used by the TE delta model (ignored by Aurora).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub matrix: Option<String>,
    /// Parts where the unit changes, as ULTRA writes them (Start relative to the
    /// record's Start); absent for an unsplit array. Aurora 0.x ignores this field.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sub_repeats: Vec<SubRepeat>,
    /// ultrafast only, not written to the JSON: the alignment to `consensus` that the
    /// Score and deltas were computed along, for a call that ultrafast scored itself and
    /// did not split. BPAF output writes this alignment, so the record's Score describes
    /// the alignment stored with it. Its own `score` field is not set.
    #[serde(skip)]
    pub scored_alignment: Option<ultra_repeats::wrap_align::Alignment>,
}

/// A scoring path as the alignment type BPAF output takes.
fn to_alignment(ops: &[deltas::Op], period: usize) -> ultra_repeats::wrap_align::Alignment {
    use ultra_repeats::wrap_align::{Alignment, Op};
    let phase = ops.iter().find_map(|o| match o {
        deltas::Op::M(c) | deltas::Op::D(c) => Some(*c as usize),
        deltas::Op::I => None,
    });
    let mut runs: Vec<(Op, u64)> = Vec::new();
    for o in ops {
        let k = match o {
            deltas::Op::M(_) => Op::Aligned,
            deltas::Op::I => Op::Inserted,
            deltas::Op::D(_) => Op::Skipped,
        };
        match runs.last_mut() {
            Some(r) if r.0 == k => r.1 += 1,
            _ => runs.push((k, 1)),
        }
    }
    Alignment { period, phase: phase.unwrap_or(0), runs, score: 0 }
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
#[serde(rename_all = "PascalCase")]
pub struct SubRepeat {
    pub start: usize,
    pub consensus: String,
}

impl UltraRecord {
    /// The call aligned to its consensus, one part per sub-repeat (or one
    /// for the whole call): the alignment its Score was computed along when
    /// ultrafast scored the call (`scored_alignment`), else as `ultra` aligns
    /// its repeats, see [`ultra_repeats::wrap_align::align_parts`]. `seq` is the whole
    /// sequence the call is on; the call is clipped to it.
    pub fn consensus_alignment(&self, seq: &[u8]) -> Vec<ultra_repeats::wrap_align::Part> {
        let start = self.start.min(seq.len());
        let end = self.start.saturating_add(self.length).min(seq.len());
        if let Some(a) = self.scored_alignment.as_ref().filter(|a| self.sub_repeats.is_empty() && a.bases() == (end - start) as u64) {
            return vec![ultra_repeats::wrap_align::Part { start: 0, end: end - start, consensus: self.consensus.as_bytes().to_vec(), alignment: a.clone() }];
        }
        let parts: Vec<(usize, &[u8])> = if self.sub_repeats.is_empty() {
            vec![(0, self.consensus.as_bytes())]
        } else {
            self.sub_repeats.iter().map(|s| (s.start, s.consensus.as_bytes())).collect()
        };
        ultra_repeats::wrap_align::align_parts(&seq[start..end], &parts)
    }
}

fn subs_of(parts: &[(usize, Vec<u8>)]) -> Vec<SubRepeat> {
    parts.iter().map(|(o, u)| SubRepeat { start: *o, consensus: String::from_utf8_lossy(u).into_owned() }).collect()
}

#[derive(serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "PascalCase")]
pub struct UltraJson {
    pub repeats: Vec<UltraRecord>,
}

fn gc_pct(seq: &[u8]) -> f64 {
    let (mut gc, mut n) = (0usize, 0usize);
    for &c in seq {
        match consensus::code(c) {
            1 | 2 => {
                gc += 1;
                n += 1
            }
            0 | 3 => n += 1,
            _ => {}
        }
    }
    100.0 * gc as f64 / n.max(1) as f64
}

/// Deltas and their sum for seq[start..end] against `unit`.
pub fn score_region(seq: &[u8], seq_gc: f64, start: usize, end: usize, unit: &[u8], model: &DeltaModel, max_cells: usize)
    -> anyhow::Result<(Vec<f64>, Option<String>)> {
    score_region_path(seq, seq_gc, start, end, unit, model, max_cells).map(|(d, m, _)| (d, m))
}

/// [`score_region`], also returning the alignment path the deltas were computed along.
fn score_region_path(seq: &[u8], seq_gc: f64, start: usize, end: usize, unit: &[u8], model: &DeltaModel, max_cells: usize)
    -> anyhow::Result<(Vec<f64>, Option<String>, Vec<deltas::Op>)> {
    let arr = &seq[start..end];
    let costs = match model { DeltaModel::Ultra => deltas::ULTRA_COSTS, DeltaModel::Te { .. } => deltas::TE_COSTS };
    let ops = deltas::path(arr, unit, max_cells, costs);
    let (d, mname) = match model {
        DeltaModel::Ultra => (deltas::ultra_deltas(arr, unit.len(), &ops), None),
        DeltaModel::Te { matrices, fixed, gc, window } => {
            let m = match fixed {
                Some(n) => matrices.iter().find(|m| &m.name == n)
                    .ok_or_else(|| anyhow::anyhow!("matrix {n} not in the matrices file"))?,
                None => {
                    // substitutions + inserted bases + skipped columns
                    let diffs = ops.len() - count_matches(arr, unit, &ops);
                    let div = 100.0 * diffs as f64 / arr.len().max(1) as f64;
                    deltas::pick(matrices, div, gc.unwrap_or(seq_gc))
                        .ok_or_else(|| anyhow::anyhow!("no <div>p<gc>g matrix to pick from"))?
                }
            };
            (deltas::te_deltas(seq, start, unit, &ops, m, *window), Some(m.name.clone()))
        }
    };
    Ok((d, mname, ops))
}

fn count_matches(arr: &[u8], unit: &[u8], ops: &[deltas::Op]) -> usize {
    let u: Vec<u8> = unit.iter().map(|&c| consensus::code(c)).filter(|&c| c < 4).collect();
    let mut i = 0;
    let mut n = 0;
    for o in ops {
        match *o {
            deltas::Op::M(c) => {
                if consensus::code(arr[i]) == u[c as usize] {
                    n += 1;
                }
                i += 1;
            }
            deltas::Op::I => i += 1,
            deltas::Op::D(_) => {}
        }
    }
    n
}

/// Fraction of [s, e) covered by `ivs`: sorted, non-overlapping intervals.
fn covered(s: usize, e: usize, ivs: &[(usize, usize)]) -> f64 {
    let i0 = ivs.partition_point(|iv| iv.1 <= s);
    let mut cov = 0;
    for &(a, b) in &ivs[i0..] {
        if a >= e {
            break;
        }
        cov += b.min(e) - a.max(s);
    }
    cov as f64 / (e - s).max(1) as f64
}

fn merge(mut v: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    v.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(v.len());
    for (a, b) in v {
        match out.last_mut() {
            Some(l) if a <= l.1 => l.1 = l.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// Rotate `unit` to the phase of the array's first copy (fewest mismatches against
/// `array[..unit.len()]`, earliest rotation on ties), so units read as ULTRA's do.
pub fn phase_to_start(unit: Vec<u8>, array: &[u8]) -> Vec<u8> {
    let n = unit.len();
    if n < 2 || array.len() < n {
        return unit;
    }
    let up = |c: u8| c.to_ascii_uppercase();
    let best = (0..n)
        .min_by_key(|&r| (0..n).filter(|&i| up(unit[(r + i) % n]) != up(array[i])).count())
        .unwrap_or(0);
    let mut u = unit;
    u.rotate_left(best);
    u
}

/// Re-vote a filler call's units part by part: each of the filler's own parts (ULTRA's
/// SubRepeats; the whole array if it has none) gets a TRF-style consensus seeded with
/// its middle copy, and is then split further where its unit changes (`o.split`).
/// Returns the parts (empty when there's only one).
fn revote_parts(array: &[u8], parts: &[(usize, Vec<u8>)], period: usize, o: &Options) -> Vec<(usize, Vec<u8>)> {
    let bounds: Vec<(usize, usize, usize)> = if parts.len() < 2 {
        vec![(0, array.len(), period)]
    } else {
        let ends: Vec<usize> = parts.iter().skip(1).map(|x| x.0).chain([array.len()]).collect();
        parts.iter().zip(ends).filter(|(p, e)| *e > p.0).map(|(p, e)| (p.0, e, p.1.len().max(1))).collect()
    };
    let mut out: Vec<(usize, Vec<u8>)> = Vec::new();
    for (b, e, p) in bounds {
        let a = &array[b..e];
        let fallback = parts.iter().find(|x| x.0 == b).map(|x| x.1.clone()).unwrap_or_else(|| middle_copy(a, 0, a.len(), p).to_vec());
        let u = refine(a, &acgt_seed(a, 0, a.len(), p, &fallback), &o.cons).unit;
        match &o.split {
            Some(sp) => {
                let sub = split_parts(a, &u, &o.cons, sp);
                if sub.is_empty() {
                    out.push((b, phase_to_start(u, a)));
                } else {
                    out.extend(sub.into_iter().map(|(x, v)| (b + x, v)));
                }
            }
            None => out.push((b, phase_to_start(u, a))),
        }
    }
    // neighbouring parts that ended up with the same unit are one part
    let mut merged: Vec<(usize, Vec<u8>)> = Vec::new();
    for (o_, u) in out {
        if merged.last().is_some_and(|l| consensus_cyclic_eq(&l.1, &u)) {
            continue;
        }
        merged.push((o_, u));
    }
    if merged.len() < 2 { Vec::new() } else { merged }
}

/// Split a filler call where its unit changes, keeping the filler's units where nothing
/// changes: each of the filler's own parts (ULTRA's SubRepeats; the whole array if it
/// has none) is split with [`split_parts`], seeded with that part's unit. Parts whose
/// unit holds anything but ACGT are left as they are. Returns the parts (empty when
/// there's only one).
fn split_filler(array: &[u8], parts: &[(usize, Vec<u8>)], unit: &[u8], o: &Options) -> Vec<(usize, Vec<u8>)> {
    let Some(sp) = &o.split else { return parts.to_vec() };
    let bounds: Vec<(usize, usize, &[u8])> = if parts.len() < 2 {
        vec![(0, array.len(), unit)]
    } else {
        let ends: Vec<usize> = parts.iter().skip(1).map(|x| x.0).chain([array.len()]).collect();
        parts.iter().zip(ends).filter(|(p, e)| *e > p.0).map(|(p, e)| (p.0, e, p.1.as_slice())).collect()
    };
    let mut out: Vec<(usize, Vec<u8>)> = Vec::new();
    for (b, e, u) in bounds {
        let a = &array[b..e];
        let sub = if !u.is_empty() && u.iter().all(|&c| consensus::code(c) < 4) {
            split_parts(a, u, &o.cons, sp)
        } else {
            Vec::new()
        };
        if sub.is_empty() {
            out.push((b, u.to_vec()));
        } else {
            out.extend(sub.into_iter().map(|(x, v)| (b + x, v)));
        }
    }
    // neighbouring parts that ended up with the same unit are one part
    let mut merged: Vec<(usize, Vec<u8>)> = Vec::new();
    for (o_, u) in out {
        if merged.last().is_some_and(|l| consensus_cyclic_eq(&l.1, &u)) {
            continue;
        }
        merged.push((o_, u));
    }
    if merged.len() < 2 { Vec::new() } else { merged }
}

/// Fewest mismatches between `a` and any rotation of `b` (equal lengths).
fn cyclic_mismatches(a: &[u8], b: &[u8]) -> usize {
    let n = a.len();
    (0..n).map(|r| (0..n).filter(|&i| a[i] != b[(i + r) % n]).count()).min().unwrap_or(0)
}

/// A gap between two pieces may be this many units long; the bridging cost decides
/// within that.
const JOIN_MAX_GAP_UNITS: usize = 2;
/// Most a join may cost: see [`bridge_cost`]. On simulated arrays 90% of the breaks
/// ULTRA leaves inside one array cost 20 or less; on nine genomes most candidate joins
/// at periods above 30 cost over 60.
const JOIN_MAX_COST: i64 = 20;

/// What it costs to read the gap between two calls as part of one array: the
/// wraparound alignment score (matches +2, mismatches -2, indels -5 a base) of the end
/// of the first call and of the start of the second, each aligned to `unit` on its own,
/// minus the score of the two with the gap between them aligned as one stretch. Up to
/// five units (at least 60 bases) of each call are used. A gap that is a diverged copy
/// in phase costs little or nothing; unrelated sequence costs about a point a base, and
/// a second piece out of phase with the first costs its indels.
fn bridge_cost(seq: &[u8], a_start: usize, a_end: usize, b_start: usize, b_end: usize, unit: &[u8]) -> i64 {
    use ultra_repeats::wrap_align::{align_bytes, ULTRA_SCORES};
    let w = (5 * unit.len()).max(60);
    let t0 = a_end.saturating_sub(w).max(a_start);
    let h1 = (b_start + w).min(b_end);
    let s = |x: &[u8]| align_bytes(x, unit, ULTRA_SCORES).score;
    s(&seq[t0..a_end]) + s(&seq[b_start..h1]) - s(&seq[t0..h1])
}

/// Join filler calls that are pieces of one array. ULTRA's model sometimes leaves a
/// repeat at an indel or a cluster of mutations and re-enters it a few bases on, so a
/// long array is reported as several calls of the same period with gaps shorter than
/// a unit. Two calls on one sequence are joined when they have the same period, their
/// units differ in at most a tenth of their positions at the best rotation, the second
/// starts after the first ends and no more than two units on, and reading the gap as
/// part of the array costs no more than `JOIN_MAX_COST` (see [`bridge_cost`]).
/// Calls of other periods lying in between don't prevent a join. The joined call keeps
/// the longer piece's unit, has the pieces' units as SubRepeats where they differ, and
/// carries no deltas, so it is scored again over its whole span.
/// Returns the calls in the order given, joined ones at the place of their first piece.
pub fn join_filler<F: std::borrow::Borrow<UltraRecord>>(recs: &[Record], filler: &[F]) -> Vec<UltraRecord> {
    let filler: Vec<&UltraRecord> = filler.iter().map(|f| f.borrow()).collect();
    // a record without sequence may share its name with one that has
    let idx: std::collections::HashMap<&str, usize> =
        recs.iter().enumerate().filter(|(_, r)| !r.seq.is_empty()).map(|(i, r)| (r.name.as_str(), i)).collect();
    let name = |f: &UltraRecord| f.sequence_name.split_whitespace().next().unwrap_or("").to_string();
    let mut order: Vec<usize> = (0..filler.len()).collect();
    order.sort_by(|&a, &b| (name(filler[a]), filler[a].start, filler[a].length).cmp(&(name(filler[b]), filler[b].start, filler[b].length)));
    let mut out: Vec<Option<UltraRecord>> = filler.iter().map(|&f| Some(f.clone())).collect();
    // (sequence, period) -> index in `out` of the last call of that period
    let mut open: std::collections::HashMap<(String, usize), usize> = std::collections::HashMap::new();
    let debug = std::env::var_os("ULTRAFAST_JOIN_DEBUG").is_some();
    let max_cost = std::env::var("ULTRAFAST_JOIN_COST").ok().and_then(|v| v.parse().ok()).unwrap_or(JOIN_MAX_COST);
    for i in order {
        let f = filler[i].clone();
        let key = (name(&f), f.period);
        let join = open.get(&key).copied().filter(|&j| {
            let a = out[j].as_ref().expect("open call");
            let end = a.start + a.length;
            let (ua, ub) = (a.consensus.as_bytes(), f.consensus.as_bytes());
            if !(f.period >= 1
                && f.start >= end
                && f.start - end <= JOIN_MAX_GAP_UNITS * f.period
                && ua.len() == f.period
                && ub.len() == f.period
                && 10 * cyclic_mismatches(ua, ub) <= f.period)
            {
                return false;
            }
            let Some(seq) = idx.get(key.0.as_str()).map(|&si| recs[si].seq.as_slice()) else { return false };
            if f.start + f.length > seq.len() {
                return false;
            }
            let unit = if a.length >= f.length { ua } else { ub };
            let cost = bridge_cost(seq, a.start, end, f.start, f.start + f.length, unit);
            if debug {
                eprintln!("joincand\t{}\t{}\t{}\t{}\t{}\t{}", key.0, end, f.period, f.start - end, cost, cyclic_mismatches(ua, ub));
            }
            cost <= max_cost
        });
        let Some(j) = join else {
            open.insert(key, i);
            continue;
        };
        let a = out[j].take().expect("open call");
        out[i] = None;
        let off = f.start - a.start;
        let parts_of = |r: &UltraRecord, o: usize| -> Vec<SubRepeat> {
            if r.sub_repeats.is_empty() {
                vec![SubRepeat { start: o, consensus: r.consensus.clone() }]
            } else {
                r.sub_repeats.iter().map(|x| SubRepeat { start: o + x.start, consensus: x.consensus.clone() }).collect()
            }
        };
        let mut parts = parts_of(&a, 0);
        for x in parts_of(&f, off) {
            if !parts.last().is_some_and(|l| consensus_cyclic_eq(l.consensus.as_bytes(), x.consensus.as_bytes())) {
                parts.push(x);
            }
        }
        // the longer piece's unit, in the phase of the joined call's first copy
        let consensus = if a.length >= f.length {
            a.consensus.clone()
        } else {
            match idx.get(key.0.as_str()) {
                Some(&si) if a.start + f.period <= recs[si].seq.len() => {
                    String::from_utf8_lossy(&phase_to_start(f.consensus.as_bytes().to_vec(), &recs[si].seq[a.start..])).into_owned()
                }
                _ => a.consensus.clone(),
            }
        };
        out[j] = Some(UltraRecord {
            length: f.start + f.length - a.start,
            score: a.score + f.score,
            consensus,
            position_score_deltas: Vec::new(),
            sub_repeats: if parts.len() < 2 { Vec::new() } else { parts },
            scored_alignment: None,
            ..a
        });
    }
    out.into_iter().flatten().collect()
}

fn consensus_cyclic_eq(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && (a.is_empty() || (0..a.len()).any(|r| a[r..].iter().chain(&a[..r]).eq(b.iter())))
}

/// ULTRA-format records for FasTAN calls plus filler records (e.g. ULTRA -p 20
/// calls) with score >= `filler_min`, merged as `Merge` says. FasTAN records' deltas
/// (and filler records', unless kept) are computed with `model`; Score = their sum.
/// With `filler_cons`, each filler call's unit is re-voted with the TRF-style caller
/// (seed: the middle copy of the filler's period), as FasTAN calls are; its Period
/// becomes the new unit's length. Kept filler deltas are unaffected (ULTRA scores
/// each base against the previous copy, not the consensus).
/// Records are sorted by sequence (input order) and start.
pub fn to_ultra_json(recs: &[Record], calls: &[Call], filler: &[UltraRecord], merge_opt: &Merge,
                     model: &DeltaModel, max_cells: usize) -> anyhow::Result<Merged> {
    let Merge { filler_min, keep_filler_deltas, prefer_filler, filler_cap, filler_cons, filler_split, filler_join } = *merge_opt;
    let gcs: Vec<f64> = recs.par_iter().map(|r| gc_pct(&r.seq)).collect();
    // a record without sequence gets no calls, and may share its name with one that has
    let idx: std::collections::HashMap<&str, usize> =
        recs.iter().enumerate().filter(|(_, r)| !r.seq.is_empty()).map(|(i, r)| (r.name.as_str(), i)).collect();
    // filler calls passing the cutoff: (seq index, start, end, record)
    let mut n_unmatched = 0;
    let fl: Vec<(usize, usize, usize, &UltraRecord)> = filler
        .iter()
        .filter(|f| f.score >= filler_min)
        .filter_map(|f| {
            // ULTRA names a sequence by its whole header line; match on the first word
            let Some(&si) = idx.get(f.sequence_name.split_whitespace().next().unwrap_or("")) else {
                n_unmatched += 1;
                return None;
            };
            let (s, e) = (f.start, (f.start + f.length).min(recs[si].seq.len()));
            (e > s).then_some((si, s, e, f))
        })
        .collect();
    // filler calls (index into `fl`) by sequence, sorted by start, with the longest per
    // sequence (for overlap lookups)
    let mut by_seq: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); recs.len()];
    for (i, &(si, s, e, _)) in fl.iter().enumerate() {
        by_seq[si].push((s, e, i));
    }
    by_seq.iter_mut().for_each(|v| v.sort_unstable_by_key(|x| (x.0, x.1)));
    let longest: Vec<usize> = by_seq.iter().map(|v| v.iter().map(|x| x.1 - x.0).max().unwrap_or(0)).collect();
    let period_of = |c: &Call| if c.cons.is_some() { c.unit.len() } else { c.period };
    // A FasTAN call whose period is above the filler's cap is one the filler cannot
    // report, so it competes with the filler calls it overlaps, by score under ULTRA's
    // model. Each side is scored base by base: the FasTAN call with
    // `deltas::ultra_deltas`, a filler call with the larger of its own deltas and
    // `deltas::ultra_deltas` on its own unit (the two agree closely on arrays without
    // indels and differ around indels in either direction).
    //  * The FasTAN call wins when its whole score beats the filler calls' scores on the
    //    bases they share with it (the filler leaves the rest of the span as background,
    //    which scores nothing) by more than `UNIT_BASE_BITS` for every base its unit is
    //    longer than the unit of the filler call sharing the most bases. A winner is
    //    kept; the coverage rule below decides for the others.
    //  * A filler call is dropped when a winner covers more than half of it and
    //    outscores it on the bases they share. There the winner's first copy, which has
    //    no copy before it to be compared with, is scored against the copy after it.
    // A fragment of one larger filler call (it covers more than half of the FasTAN call,
    // which covers no more than half of it) is not compared.
    // The third field is the FasTAN call's deltas under ULTRA's model when the comparison computed
    // them: the final scoring below reuses them instead of aligning the array again.
    let verdicts: Vec<(bool, Vec<usize>, Option<(Vec<f64>, Vec<deltas::Op>)>)> = calls
        .par_iter()
        .with_max_len(1)
        .map(|c| {
            let none = (false, Vec::new(), None);
            let Some(cap) = filler_cap.filter(|_| prefer_filler) else { return none };
            let period = period_of(c);
            if period <= cap {
                return none;
            }
            let v = &by_seq[c.seq];
            let i0 = v.partition_point(|x| x.0 + longest[c.seq] <= c.start);
            let (cs, ce) = (c.start, c.end);
            let over: Vec<(usize, usize, usize)> = v[i0..].iter().take_while(|x| x.0 < ce).filter(|x| x.1 > cs).copied().collect();
            if over.is_empty() {
                return none;
            }
            let shared = |s: usize, e: usize| (s.max(cs), e.min(ce));
            if over.iter().any(|&(s, e, _)| {
                let (a, b) = shared(s, e);
                2 * (b - a) > ce - cs && 2 * (b - a) <= e - s
            }) {
                return none;
            }
            let seq = &recs[c.seq].seq;
            let score = |arr: &[u8], unit: &[u8]| deltas::ultra_deltas(arr, unit.len(), &deltas::path(arr, unit, max_cells, deltas::ULTRA_COSTS));
            let fastan_path = deltas::path(&seq[cs..ce], &c.unit, max_cells, deltas::ULTRA_COSTS);
            let fastan = deltas::ultra_deltas(&seq[cs..ce], c.unit.len(), &fastan_path);
            let mut filler_total = 0.0;
            let mut most_shared = (0, 0); // (shared bases, period) of the filler call sharing the most
            let theirs: Vec<f64> = over
                .iter()
                .map(|&(s, e, i)| {
                    let f = fl[i].3;
                    let (a, b) = shared(s, e);
                    let own: f64 = if f.position_score_deltas.len() == f.length {
                        f.position_score_deltas[a - f.start..b - f.start].iter().sum()
                    } else {
                        f.score * (b - a) as f64 / f.length.max(1) as f64
                    };
                    let same: f64 = score(&seq[s..e], f.consensus.as_bytes())[a - s..b - s].iter().sum();
                    filler_total += own.max(same);
                    most_shared = most_shared.max((b - a, f.period));
                    own.max(same)
                })
                .collect();
            let longer = period.saturating_sub(most_shared.1) as f64;
            let reuse = matches!(model, DeltaModel::Ultra).then(|| (fastan.clone(), fastan_path));
            if fastan.iter().sum::<f64>() <= filler_total + UNIT_BASE_BITS * longer {
                return (false, Vec::new(), reuse);
            }
            // the first copy against the copy after it: score the array read backwards
            let mut ours = fastan;
            let first = period.min(ours.len());
            if over.iter().any(|&(s, _, _)| s < cs + first) {
                let arr: Vec<u8> = seq[cs..ce].iter().rev().copied().collect();
                let unit: Vec<u8> = c.unit.iter().rev().copied().collect();
                let backward = score(&arr, &unit);
                for (x, y) in ours[..first].iter_mut().zip(backward.iter().rev()) {
                    *x = *y;
                }
            }
            let beaten = over
                .iter()
                .zip(&theirs)
                .filter(|&(&(s, e, _), &t)| {
                    let (a, b) = shared(s, e);
                    2 * (b - a) > e - s && ours[a - cs..b - cs].iter().sum::<f64>() > t
                })
                .map(|(x, _)| x.2)
                .collect();
            (true, beaten, reuse)
        })
        .collect();
    // with FasTAN preferred, all FasTAN calls together suppress filler calls
    let mut ivs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); recs.len()];
    if !prefer_filler {
        for c in calls {
            ivs[c.seq].push((c.start, c.end));
        }
    }
    let ivs: Vec<Vec<(usize, usize)>> = ivs.into_iter().map(merge).collect();
    let mut fivs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); recs.len()];
    if prefer_filler {
        for &(si, s, e, _) in &fl {
            fivs[si].push((s, e));
        }
    }
    let fivs: Vec<Vec<(usize, usize)>> = fivs.into_iter().map(merge).collect();
    // (seq index, start, end, unit, period, filler record whose deltas to keep, re-vote the unit, parts,
    //  is a filler call)
    let kept: Vec<bool> = calls
        .iter()
        .zip(&verdicts)
        .map(|(c, v)| !prefer_filler || v.0 || covered(c.start, c.end, &fivs[c.seq]) <= 0.5)
        .collect();
    let mut beaten = vec![false; fl.len()];
    for v in &verdicts {
        v.1.iter().for_each(|&i| beaten[i] = true);
    }
    type Job<'a> = (usize, usize, usize, Vec<u8>, usize, Option<&'a UltraRecord>, bool, Vec<(usize, Vec<u8>)>, bool, Option<&'a (Vec<f64>, Vec<deltas::Op>)>);
    let mut jobs: Vec<Job> = calls
        .iter()
        .enumerate()
        .filter(|&(i, _)| kept[i])
        .map(|(i, c)| (c.seq, c.start, c.end, c.unit.clone(), period_of(c), None, false, c.parts.clone(), false, verdicts[i].2.as_ref()))
        .collect();
    // The filler calls that remain. Pieces of one array are joined only now, so each
    // piece has met the score cutoff and the FasTAN calls on its own.
    let rest: Vec<&UltraRecord> = fl
        .iter()
        .enumerate()
        .filter(|&(i, &(si, s, e, _))| !if prefer_filler { beaten[i] } else { covered(s, e, &ivs[si]) > 0.5 })
        .map(|(_, x)| x.3)
        .collect();
    let joined: Vec<UltraRecord> = if filler_join { join_filler(recs, &rest) } else { Vec::new() };
    let rest: Vec<&UltraRecord> = if filler_join { joined.iter().collect() } else { rest };
    let mut nfill = 0;
    for f in rest {
        let si = idx[f.sequence_name.split_whitespace().next().unwrap_or("")];
        let (s, e) = (f.start, (f.start + f.length).min(recs[si].seq.len()));
        let keep = keep_filler_deltas && f.position_score_deltas.len() == f.length && e == f.start + f.length;
        let revote = filler_cons.is_some_and(|o| f.period >= 1 && f.period <= o.max_cons_period);
        let fparts = f.sub_repeats.iter().filter(|x| x.start < e - s).map(|x| (x.start, x.consensus.as_bytes().to_vec())).collect();
        jobs.push((si, s, e, f.consensus.as_bytes().to_vec(), f.period, keep.then_some(f), revote, fparts, true, None));
        nfill += 1;
    }
    jobs.sort_by_key(|j| (j.0, j.1, j.2));
    let out: anyhow::Result<Vec<UltraRecord>> = jobs
        .par_iter()
        .with_max_len(1)
        .map(|(si, s, e, unit, period, kept, revote, parts, is_filler, scored)| {
            let seq = &recs[*si].seq;
            let (unit, period, parts) = match (revote, filler_cons, filler_split) {
                (true, Some(o), _) => {
                    let seed = acgt_seed(seq, *s, *e, *period, unit);
                    let u = phase_to_start(refine(&seq[*s..*e], &seed, &o.cons).unit, &seq[*s..*e]);
                    if u.is_empty() || !u.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
                        // No usable vote: keep ULTRA's unit and period.
                        (unit.clone(), *period, parts.clone())
                    } else {
                        let n = u.len();
                        (u, n, revote_parts(&seq[*s..*e], parts, *period, o))
                    }
                }
                (_, _, Some(o)) if *is_filler => (unit.clone(), *period, split_filler(&seq[*s..*e], parts, unit, o)),
                _ => (unit.clone(), *period, parts.clone()),
            };
            let (unit, period) = (&unit, &period);
            let sub_repeats = subs_of(&parts);
            if let Some(f) = kept {
                return Ok(UltraRecord {
                    sequence_name: recs[*si].name.clone(),
                    matrix: None,
                    consensus: String::from_utf8_lossy(unit).into_owned(),
                    period: *period,
                    sub_repeats,
                    ..(*f).clone()
                });
            }
            // a FasTAN call keeps its unit and span from the comparison, so its deltas are the same
            let fresh;
            let (d, matrix, ops): (Vec<f64>, _, &[deltas::Op]) = match scored {
                Some((d, ops)) => (d.clone(), None, ops),
                None => {
                    let (d, m, ops) = score_region_path(seq, gcs[*si], *s, *e, unit, model, max_cells)?;
                    fresh = ops;
                    (d, m, &fresh)
                }
            };
            // the path's columns are the unit's own when it is plain ACGT
            let plain = !unit.is_empty() && unit.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T'));
            let scored_alignment = (sub_repeats.is_empty() && plain).then(|| to_alignment(ops, unit.len()));
            Ok(UltraRecord {
                sequence_name: recs[*si].name.clone(),
                start: *s,
                length: e - s,
                period: *period,
                score: round6(d.iter().sum()),
                consensus: String::from_utf8_lossy(unit).into_owned(),
                position_score_deltas: d.into_iter().map(round6).collect(),
                matrix,
                sub_repeats,
                scored_alignment,
            })
        })
        .collect();
    Ok(Merged { records: out?, n_filler: nfill, n_unmatched })
}

/// How `to_ultra_json` merges FasTAN's calls with the filler's (ULTRA's).
#[derive(Clone, Copy)]
pub struct Merge<'a> {
    /// Minimum filler Score.
    pub filler_min: f64,
    /// Keep the filler's own deltas and Score (else recompute them with the model).
    pub keep_filler_deltas: bool,
    /// Where they overlap, the filler call wins (else FasTAN's does).
    pub prefer_filler: bool,
    /// The filler's maximum period (ULTRA's -p). With `prefer_filler`, a FasTAN call
    /// whose period is above it is compared with the filler calls it overlaps by score
    /// under ULTRA's model, as `hybrid` describes. None: filler calls always win.
    pub filler_cap: Option<usize>,
    /// Re-vote filler units with the TRF-style caller.
    pub filler_cons: Option<&'a Options>,
    /// Split filler calls where their unit changes with these options' `split` (and
    /// `cons`), keeping the filler's units elsewhere; see `split_filler`. Not used for
    /// calls that `filler_cons` re-votes, which are split as they are re-voted.
    pub filler_split: Option<&'a Options>,
    /// Join the remaining filler calls that are pieces of one array; see [`join_filler`].
    pub filler_join: bool,
}

/// Result of `to_ultra_json`.
pub struct Merged {
    /// Output records, sorted by sequence and start.
    pub records: Vec<UltraRecord>,
    /// How many of them are filler calls.
    pub n_filler: usize,
    /// Filler calls dropped because their sequence is not in the FASTA.
    pub n_unmatched: usize,
}

// ------------------------------------------------------------- the hybrid

/// Settings for [`hybrid`]; `Default` is what the `ultrafast` command runs with no options.
pub struct HybridOptions {
    /// ULTRA's maximum period; 0 = don't run ULTRA (FasTAN only).
    pub ultra_period: usize,
    /// Threads for FasTAN (a rayon pool for this call) and ULTRA; 0 = all cores.
    pub threads: usize,
    /// FasTAN's consensus and splitting (`consensus` is always switched on).
    pub fastan: Options,
    /// Use these ULTRA-format calls instead of running ULTRA (e.g. read with
    /// [`read_ultra_json`] from a separate `ultra --json --show_delta` run).
    pub filler: Option<Vec<UltraRecord>>,
    /// Minimum ULTRA score for a short-period call.
    pub filler_min_score: f64,
    /// Join ULTRA calls of one period that lie less than a unit apart and have the
    /// same unit (pieces of one array); see [`join_filler`].
    pub join_ultra: bool,
    /// Where calls overlap, keep FasTAN's (default: keep ULTRA's).
    pub prefer_fastan: bool,
    /// Re-vote ULTRA's units with the TRF-style caller, as FasTAN's are.
    pub revote_filler: bool,
    /// Recompute ULTRA's deltas with `delta_model` (always done for `Te` and `rescore_only`).
    pub rescore_filler: bool,
    /// Don't run FasTAN: only rescore the `filler` calls.
    pub rescore_only: bool,
    /// Per-base score deltas for FasTAN calls.
    pub delta_model: DeltaModel,
    /// Multiply every delta and Score by this.
    pub delta_scale: f64,
}

impl Default for HybridOptions {
    fn default() -> Self {
        HybridOptions {
            ultra_period: 50,
            threads: 0,
            fastan: Options { consensus: true, ..Default::default() },
            filler: None,
            filler_min_score: f64::NEG_INFINITY,
            join_ultra: true,
            prefer_fastan: false,
            revote_filler: false,
            rescore_filler: false,
            rescore_only: false,
            delta_model: DeltaModel::Ultra,
            delta_scale: 1.0,
        }
    }
}

/// The BPAF scoring-system name of a Score on ULTRA's scale (bits); a call scored
/// with a TE matrix carries the matrix's name instead.
pub const BPAF_ULTRA_SCORING: &str = ultra_repeats::writers::bpaf::ULTRA_SCORING;

/// Result of [`hybrid`].
pub struct Hybrid {
    /// The input sequences (upper-cased).
    pub seqs: Vec<Record>,
    /// Every FasTAN array, before merging (`seq` indexes `seqs`).
    pub arrays: Vec<Call>,
    /// The merged calls, as ULTRA's JSON records, sorted by sequence and start.
    pub repeats: Vec<UltraRecord>,
    /// How many of `repeats` are ULTRA's calls (the rest are FasTAN's).
    pub n_ultra: usize,
    /// `filler` calls dropped because their sequence is not in the FASTA.
    pub n_unmatched: usize,
}

impl Hybrid {
    /// The merged calls as ULTRA's JSON (`{"Repeats": [...]}`).
    pub fn write_json<W: std::io::Write>(&self, w: W) -> anyhow::Result<()> {
        #[derive(serde::Serialize)]
        #[serde(rename_all = "PascalCase")]
        struct Out<'a> {
            repeats: &'a [UltraRecord],
        }
        serde_json::to_writer(w, &Out { repeats: &self.repeats })?;
        Ok(())
    }

    /// The merged calls as BPAF: each call aligned to its consensus unit, one
    /// record per call or per sub-repeat, with the call's Score (see
    /// [`ultra_repeats::bpaf`]). The units are stored in the file; the input
    /// sequences are not. The alignments are made here, on the current rayon
    /// pool.
    pub fn write_bpaf<W: std::io::Write>(&self, mut w: W) -> anyhow::Result<()> {
        use ultra_repeats::bpaf::{TandemRecord, TandemWriter};
        let seq_of: std::collections::HashMap<&str, &[u8]> =
            self.seqs.iter().filter(|r| !r.seq.is_empty()).map(|r| (r.name.as_str(), r.seq.as_slice())).collect();
        let mut bw = TandemWriter::new();
        bw.write_header(&mut w)?;
        // in batches, so the alignments held at once stay few
        for batch in self.repeats.chunks(1024) {
            let aligned: Vec<_> = batch
                .par_iter()
                .map(|r| seq_of.get(r.sequence_name.as_str()).map(|seq| (r.start.min(seq.len()), r.consensus_alignment(seq))))
                .collect();
            for (r, a) in batch.iter().zip(aligned) {
                let Some((start, parts)) = a else { continue };
                for part in &parts {
                    bw.write_record(
                        &mut w,
                        &TandemRecord {
                            query: &r.sequence_name,
                            q_start: (start + part.start + 1) as u64,
                            unit: &part.consensus,
                            alignment: &part.alignment,
                            score: Some(r.score as f32),
                            scoring: r.matrix.as_deref().unwrap_or(BPAF_ULTRA_SCORING),
                        },
                    )?;
                }
            }
        }
        bw.finish(&mut w)?;
        Ok(())
    }
}

/// Run the hybrid on a FASTA file: FasTAN on every sequence and, alongside it,
/// ULTRA at maximum period `opt.ultra_period`; then merge them. ULTRA calls keep
/// ULTRA's units, scores and deltas; FasTAN calls get a consensus unit and deltas
/// from `opt.delta_model`.
///
/// Merging (unless `opt.prefer_fastan`): a FasTAN call with a period up to ULTRA's
/// maximum is dropped when ULTRA calls cover more than half of it. A FasTAN call with
/// a longer period, which ULTRA cannot report, is kept when its score under ULTRA's
/// model is higher than the scores of the ULTRA calls it overlaps, on the bases they
/// share with it, by more than 2 bits for every base its unit is longer than the unit
/// of the ULTRA call sharing the most bases. It then replaces each ULTRA call that it
/// covers more than half of and outscores on the bases the two share (e.g. a short
/// microsatellite inside each copy of a satellite); other ULTRA calls inside it are
/// kept (e.g. a long microsatellite that interrupts a satellite). A fragment of one
/// larger ULTRA call is never kept this way. Calls can still overlap.
///
/// Threads: `opt.threads` (0 = all cores) in total, split between ULTRA and FasTAN,
/// which run at the same time; with 1 they run one after the other.
///
/// Sequence names (the first word of each header, which is the SequenceName
/// written) must be unique among records holding sequence. ULTRA reads the file
/// itself, so `input` must be a path; `"-"` (stdin) works only when ULTRA isn't
/// run (`ultra_period` 0 or `filler` given).
///
/// ```no_run
/// let h = ultrafast::hybrid("genome.fa", &ultrafast::HybridOptions::default()).unwrap();
/// for r in &h.repeats {
///     println!("{} {} {} {} {}", r.sequence_name, r.start, r.length, r.period, r.consensus);
/// }
/// ```
pub fn hybrid(input: &str, opt: &HybridOptions) -> anyhow::Result<Hybrid> {
    let n = if opt.threads > 0 { opt.threads } else { std::thread::available_parallelism().map_or(1, |n| n.get()) };
    let runs_ultra = opt.filler.is_none() && opt.ultra_period > 0;
    // ULTRA gets half the threads (FasTAN the rest) when both run at once
    let (ultra_threads, fastan_threads) = if n == 1 || !runs_ultra { (1, n) } else { (n / 2, n - n / 2) };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(fastan_threads).build()?;
    pool.install(|| hybrid_in_pool(input, opt, ultra_threads))
}

fn hybrid_in_pool(input: &str, opt: &HybridOptions, ultra_threads: usize) -> anyhow::Result<Hybrid> {
    let ultra_period = if opt.filler.is_some() { None } else { (opt.ultra_period > 0).then_some(opt.ultra_period) };
    if ultra_period.is_some() && input == "-" {
        anyhow::bail!("ULTRA needs a FASTA file, not stdin (or run FasTAN only: ultra_period 0)");
    }
    // ULTRA reads the input again after it has been read here, which works
    // only for a regular file: a pipe (<(...), /dev/stdin) would give ULTRA
    // nothing, silently, and a named FIFO would block.
    if ultra_period.is_some() && !std::fs::metadata(input).map(|m| m.is_file()).unwrap_or(true) {
        anyhow::bail!("ULTRA needs a regular FASTA file, not a pipe: {input} (or run FasTAN only: ultra_period 0)");
    }
    let fo = Options { consensus: true, ..opt.fastan.clone() };
    let seqs = read_fasta(input)?;
    check_names(&seqs)?;
    let fastan = || if opt.rescore_only { Vec::new() } else { annotate(&seqs, &fo) };
    let ultra = |p: usize| run_ultra(input, p, ultra_threads);
    let (arrays, ultra_calls) = match ultra_period {
        None => (fastan(), None),
        // one thread: one after the other
        Some(p) if opt.threads == 1 => {
            let u = ultra(p)?;
            (fastan(), Some(u))
        }
        // ULTRA (its own threads) runs alongside FasTAN (the rayon pool)
        Some(p) => std::thread::scope(|s| {
            let h = s.spawn(move || ultra(p));
            let arrays = fastan();
            let u = h.join().map_err(|_| anyhow::anyhow!("ULTRA thread panicked"))??;
            anyhow::Ok((arrays, Some(u)))
        })?,
    };
    let filler: &[UltraRecord] = match (&ultra_calls, &opt.filler) {
        (Some(u), _) => u,
        (None, Some(f)) => f,
        (None, None) => &[],
    };
    // the filler's maximum period: ULTRA's -p, or the longest period in a --filler file
    let filler_cap = match ultra_period {
        Some(p) => Some(p),
        None => filler.iter().map(|f| f.period).max(),
    };
    let keep = !(opt.rescore_filler || opt.rescore_only || matches!(opt.delta_model, DeltaModel::Te { .. }));
    let merge_opt = Merge {
        filler_min: opt.filler_min_score,
        keep_filler_deltas: keep,
        prefer_filler: !opt.prefer_fastan,
        filler_cap,
        filler_cons: opt.revote_filler.then_some(&fo),
        filler_split: fo.split.is_some().then_some(&fo),
        filler_join: opt.join_ultra,
    };
    let m = to_ultra_json(&seqs, &arrays, filler, &merge_opt, &opt.delta_model, fo.cons.max_cells)?;
    let mut repeats = m.records;
    if opt.delta_scale != 1.0 {
        for r in repeats.iter_mut() {
            r.position_score_deltas.iter_mut().for_each(|x| *x = round6(*x * opt.delta_scale));
            r.score = round6(r.score * opt.delta_scale);
        }
    }
    Ok(Hybrid { seqs, arrays, repeats, n_ultra: m.n_filler, n_unmatched: m.n_unmatched })
}

/// What a FasTAN call must gain over the filler calls it competes with for every
/// base its unit is longer than theirs: 2 bits, the information in one base. A longer
/// unit has more bases to fit to the same array, and without this charge a divergent
/// array is sometimes reported at two or three times its period.
const UNIT_BASE_BITS: f64 = 2.0;

/// Rounded to six decimals, as ULTRA writes its numbers.
fn round6(x: f64) -> f64 {
    (x * 1e6).round() / 1e6
}

#[cfg(test)]
mod tests {
    use super::{consensus_cyclic_eq, phase_to_start, split_filler, Options};

    #[test]
    fn phase_to_start_rotates_to_first_copy() {
        assert_eq!(phase_to_start(b"TGG".to_vec(), b"GGTGGTGGT"), b"GGT".to_vec());
        assert_eq!(phase_to_start(b"TG".to_vec(), b"gtgtgt"), b"GT".to_vec());
        assert_eq!(phase_to_start(b"ACGT".to_vec(), b"AC"), b"ACGT".to_vec());
    }

    #[test]
    fn split_filler_splits_and_keeps_units() {
        let o = Options { consensus: true, ..Default::default() };
        // (ACGTACGTAC)12 (ACGTTCGTAC)12: one change at a period of 10, which ULTRA's
        // splitter (5-copy windows, 85% join) leaves as one array
        let mut arr = b"ACGTACGTAC".repeat(12);
        arr.extend(b"ACGTTCGTAC".repeat(12));
        let parts = split_filler(&arr, &[], b"ACGTACGTAC", &o);
        assert_eq!(parts.len(), 2);
        assert!((115..=125).contains(&parts[1].0), "split at {}", parts[1].0);
        // units are phased to each part's first copy
        assert!(consensus_cyclic_eq(&parts[0].1, b"ACGTACGTAC"));
        assert!(consensus_cyclic_eq(&parts[1].1, b"ACGTTCGTAC"));
        // one unit throughout: no parts, the caller keeps the filler's unit
        let one = b"ACGTACGTAC".repeat(24);
        assert!(split_filler(&one, &[], b"ACGTACGTAC", &o).is_empty());
        // a unit with an ambiguity code is left alone
        assert!(split_filler(&arr, &[], b"ACGTNCGTAC", &o).is_empty());
        // no splitting asked for: the filler's own parts come back as they are
        let none = Options { split: None, ..o.clone() };
        let fp = vec![(0, b"ACGTACGTAC".to_vec()), (120, b"ACGTTCGTAC".to_vec())];
        assert_eq!(split_filler(&arr, &fp, b"ACGTACGTAC", &none), fp);
    }

    /// A call whose SubRepeats text is not valid JSON (here a consensus byte
    /// that is not UTF-8) is collected without SubRepeats. Collect used to
    /// return an error, on which ULTRA ended the whole process.
    #[test]
    fn collect_keeps_a_call_whose_sub_repeats_do_not_parse() {
        use ultra_repeats::{cli::Settings, repeat::RepeatRegion, ultra::Ultra};
        let mut s = Settings::new();
        s.run_without_reader = true;
        s.window_size = 1000;
        s.overlap = 0;
        s.assign_settings();
        let got = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut u = Ultra::new(&s);
        u.writers = vec![Box::new(super::Collect(got.clone()))];
        u.outs.truncate(1);
        let mut r = Box::new(RepeatRegion::new());
        r.sequence_name = "chr1".to_string();
        r.read_id = 0;
        r.repeat_period = 2;
        r.repeat_length = 10;
        r.region_score = 100.0;
        r.string_consensus = b"AC".to_vec();
        r.splits = Some(vec![4]);
        r.consensi = Some(vec![b"A\xff".to_vec(), b"AC".to_vec()]);
        u.out_repeats.push(r);
        u.output_repeats(true).unwrap();
        let got = got.lock().unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].sequence_name, "chr1");
        assert_eq!(got[0].length, 10);
        assert!(got[0].sub_repeats.is_empty());
    }

    #[test]
    fn join_filler_joins_pieces_of_one_array() {
        use super::{join_filler, Record, UltraRecord};
        let unit = b"ACGTTGCAAGCTTAGGCATC"; // period 20
        let call = |start: usize, length: usize| UltraRecord {
            sequence_name: "s".into(),
            start,
            length,
            period: 20,
            score: 100.0,
            consensus: String::from_utf8_lossy(unit).into_owned(),
            position_score_deltas: vec![1.0; length],
            matrix: None,
            sub_repeats: Vec::new(),
            scored_alignment: None,
        };
        // 10 copies, one copy with three substitutions that ULTRA left out, 10 copies
        let mut seq = unit.repeat(10);
        seq.extend(b"ACGATGCTAGCTTAGGGATC");
        seq.extend(unit.repeat(10));
        let recs = vec![Record { name: "s".into(), seq: seq.clone() }];
        let j = join_filler(&recs, &[call(0, 200), call(220, 200)]);
        assert_eq!(j.len(), 1);
        assert_eq!((j[0].start, j[0].length), (0, 420));
        assert!(j[0].position_score_deltas.is_empty() && j[0].sub_repeats.is_empty());
        // the same two calls with unrelated sequence between them stay apart
        let mut seq = unit.repeat(10);
        seq.extend(b"GGGGGGGGGGCCCCCCCCCCGGGGGGGGGGCCCCCCCC");
        seq.extend(unit.repeat(10));
        let recs = vec![Record { name: "s".into(), seq }];
        assert_eq!(join_filler(&recs, &[call(0, 200), call(238, 200)]).len(), 2);
        // and so do calls more than two units apart, or of different periods
        let recs = vec![Record { name: "s".into(), seq: unit.repeat(40) }];
        assert_eq!(join_filler(&recs, &[call(0, 200), call(260, 200)]).len(), 2);
        let mut other = call(200, 200);
        other.period = 10;
        assert_eq!(join_filler(&recs, &[call(0, 200), other]).len(), 2);
    }
}
