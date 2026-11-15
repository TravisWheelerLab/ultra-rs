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
pub fn read_fasta(path: &str) -> std::io::Result<Vec<Record>> {
    let rd: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        Box::new(BufReader::new(std::fs::File::open(path)?))
    };
    let mut recs: Vec<Record> = Vec::new();
    for line in rd.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            recs.push(Record { name: h.split_whitespace().next().unwrap_or("").to_string(), seq: Vec::new() });
        } else if let Some(r) = recs.last_mut() {
            r.seq.extend(line.trim_end().bytes().map(|c| c.to_ascii_uppercase()));
        }
    }
    Ok(recs)
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
/// Score and deltas rounded to ULTRA's six decimals.
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
            serde_json::from_slice(&txt).map_err(std::io::Error::other)?
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
        };
        self.0.lock().expect("collector lock").push(rec);
        Ok(())
    }
}

/// ULTRA's calls at its defaults with maximum period `period` on the FASTA file
/// `input`: what `ultra -p period -t threads --json --show_delta` would write,
/// taken from this repository's library in memory (no files are written and
/// nothing is printed). `upstream_compat` selects the C++-compatible scoring,
/// as `ULTRA_UPSTREAM_COMPAT` does for `ultra`.
pub fn run_ultra(input: &str, period: usize, threads: usize, upstream_compat: bool) -> anyhow::Result<Vec<UltraRecord>> {
    use ultra_repeats::{cli::Settings, ultra::Ultra};
    // ULTRA exits the process on an unreadable file; fail here with a message instead
    std::fs::File::open(input).map_err(|e| anyhow::anyhow!("{input}: {e}"))?;
    let argv: Vec<String> = ["ultra", "-p", &period.to_string(), "-t", &threads.to_string(), "--json", "--show_delta", "--hs", input]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut s = Settings::new();
    s.upstream_compat = upstream_compat;
    s.prepare_settings();
    if !s.parse_input(&argv) {
        anyhow::bail!("ULTRA rejected its settings: {}", argv.join(" "));
    }
    s.assign_settings();
    let got = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    {
        let mut u = Ultra::new(&s);
        // the JSON writer would print to stdout; collect instead
        u.writers = vec![Box::new(Collect(got.clone()))];
        u.analyze_file();
        u.output_repeats(true);
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
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
#[serde(rename_all = "PascalCase")]
pub struct SubRepeat {
    pub start: usize,
    pub consensus: String,
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
    Ok((d, mname))
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
        let u = refine(a, middle_copy(a, 0, a.len(), p), &o.cons).unit;
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
    let Merge { filler_min, keep_filler_deltas, prefer_filler, filler_cap, filler_cons } = *merge_opt;
    let gcs: Vec<f64> = recs.par_iter().map(|r| gc_pct(&r.seq)).collect();
    let idx: std::collections::HashMap<&str, usize> = recs.iter().enumerate().map(|(i, r)| (r.name.as_str(), i)).collect();
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
    // filler calls by sequence, sorted by start, with the longest per sequence (for overlap lookups)
    let mut by_seq: Vec<Vec<(usize, usize, &UltraRecord)>> = vec![Vec::new(); recs.len()];
    for &(si, s, e, f) in &fl {
        by_seq[si].push((s, e, f));
    }
    by_seq.iter_mut().for_each(|v| v.sort_unstable_by_key(|x| (x.0, x.1)));
    let longest: Vec<usize> = by_seq.iter().map(|v| v.iter().map(|x| x.1 - x.0).max().unwrap_or(0)).collect();
    let period_of = |c: &Call| if c.cons.is_some() { c.unit.len() } else { c.period };
    let ov = |a: (usize, usize), b: (usize, usize)| a.1.min(b.1).saturating_sub(a.0.max(b.0));
    // A FasTAN call whose period is above the filler's cap is one the filler cannot
    // report. It is protected (kept, and suppresses the filler calls it covers more than
    // half of) when it explains its bases better than the filler calls there: its score
    // under ULTRA's model beats the filler calls' deltas summed over the same bases. That
    // settles units just above the cap that ULTRA calls at the cap with an indel per copy
    // (FasTAN wins), a satellite with a microsatellite in each copy (FasTAN wins), and
    // short divergent arrays FasTAN called at a multiple of the period (ULTRA wins). A
    // fragment of one larger filler call (it covers more than half of the FasTAN call,
    // which covers no more than half of it) is never protected.
    let protected: Vec<bool> = calls
        .par_iter()
        .with_max_len(1)
        .map(|c| {
            let Some(cap) = filler_cap.filter(|_| prefer_filler) else { return false };
            if period_of(c) <= cap {
                return false;
            }
            let v = &by_seq[c.seq];
            let i0 = v.partition_point(|x| x.0 + longest[c.seq] <= c.start);
            let (cs, ce) = (c.start, c.end);
            let over: Vec<&(usize, usize, &UltraRecord)> = v[i0..].iter().take_while(|x| x.0 < ce).filter(|x| x.1 > cs).collect();
            if over.is_empty() {
                return true;
            }
            if over.iter().any(|&&(s, e, _)| {
                let o = ov((cs, ce), (s, e));
                2 * o > ce - cs && 2 * o <= e - s
            }) {
                return false;
            }
            let seq = &recs[c.seq].seq;
            let fastan: f64 = deltas::ultra_deltas(&seq[cs..ce], c.unit.len(), &deltas::path(&seq[cs..ce], &c.unit, max_cells, deltas::ULTRA_COSTS))
                .iter()
                .sum();
            let filler: f64 = over
                .iter()
                .map(|&&(s, e, f)| {
                    let (a, b) = (s.max(cs), e.min(ce));
                    if f.position_score_deltas.len() == f.length {
                        f.position_score_deltas[a - f.start..b - f.start].iter().sum::<f64>()
                    } else {
                        f.score * (b - a) as f64 / f.length.max(1) as f64
                    }
                })
                .sum();
            fastan > filler
        })
        .collect();
    // filler calls a single protected FasTAN call covers more than half of are dropped
    let mut prot: Vec<Vec<(usize, usize)>> = vec![Vec::new(); recs.len()];
    for (i, c) in calls.iter().enumerate() {
        if protected[i] {
            prot[c.seq].push((c.start, c.end));
        }
    }
    prot.iter_mut().for_each(|v| v.sort_unstable());
    let plongest: Vec<usize> = prot.iter().map(|v| v.iter().map(|x| x.1 - x.0).max().unwrap_or(0)).collect();
    let swallowed = |si: usize, s: usize, e: usize| {
        let v = &prot[si];
        let i0 = v.partition_point(|x| x.0 + plongest[si] <= s);
        v[i0..].iter().take_while(|x| x.0 < e).any(|&x| 2 * ov(x, (s, e)) > e - s)
    };
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
    // (seq index, start, end, unit, period, filler record whose deltas to keep, re-vote the unit, parts)
    type Job<'a> = (usize, usize, usize, Vec<u8>, usize, Option<&'a UltraRecord>, bool, Vec<(usize, Vec<u8>)>);
    let mut jobs: Vec<Job> = calls
        .iter()
        .enumerate()
        .filter(|&(i, c)| !prefer_filler || protected[i] || covered(c.start, c.end, &fivs[c.seq]) <= 0.5)
        .map(|(_, c)| (c.seq, c.start, c.end, c.unit.clone(), period_of(c), None, false, c.parts.clone()))
        .collect();
    let mut nfill = 0;
    for &(si, s, e, f) in &fl {
        if if prefer_filler { swallowed(si, s, e) } else { covered(s, e, &ivs[si]) > 0.5 } {
            continue;
        }
        let keep = keep_filler_deltas && f.position_score_deltas.len() == f.length && e == f.start + f.length;
        let revote = filler_cons.is_some_and(|o| f.period >= 1 && f.period <= o.max_cons_period);
        let fparts = f.sub_repeats.iter().filter(|x| x.start < e - s).map(|x| (x.start, x.consensus.as_bytes().to_vec())).collect();
        jobs.push((si, s, e, f.consensus.as_bytes().to_vec(), f.period, keep.then_some(f), revote, fparts));
        nfill += 1;
    }
    jobs.sort_by_key(|j| (j.0, j.1, j.2));
    let out: anyhow::Result<Vec<UltraRecord>> = jobs
        .par_iter()
        .with_max_len(1)
        .map(|(si, s, e, unit, period, kept, revote, parts)| {
            let seq = &recs[*si].seq;
            let (unit, period, parts) = match (revote, filler_cons) {
                (true, Some(o)) => {
                    let u = phase_to_start(refine(&seq[*s..*e], middle_copy(seq, *s, *e, *period), &o.cons).unit, &seq[*s..*e]);
                    let n = u.len();
                    (u, n, revote_parts(&seq[*s..*e], parts, *period, o))
                }
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
            let (d, matrix) = score_region(seq, gcs[*si], *s, *e, unit, model, max_cells)?;
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
    /// whose period is above it and that scores better under ULTRA's model than the
    /// filler calls on the same bases is kept, and filler calls it covers more than half
    /// of are dropped (a fragment of one larger filler call never is). None: filler calls
    /// always win.
    pub filler_cap: Option<usize>,
    /// Re-vote filler units with the TRF-style caller.
    pub filler_cons: Option<&'a Options>,
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
    /// Run ULTRA with the C++-compatible scoring (`ultra` does this when the
    /// environment variable `ULTRA_UPSTREAM_COMPAT` is set).
    pub ultra_upstream_compat: bool,
}

impl Default for HybridOptions {
    fn default() -> Self {
        HybridOptions {
            ultra_period: 20,
            threads: 0,
            fastan: Options { consensus: true, ..Default::default() },
            filler: None,
            filler_min_score: f64::NEG_INFINITY,
            prefer_fastan: false,
            revote_filler: false,
            rescore_filler: false,
            rescore_only: false,
            delta_model: DeltaModel::Ultra,
            delta_scale: 1.0,
            ultra_upstream_compat: false,
        }
    }
}

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
}

/// Run the hybrid on a FASTA file: FasTAN on every sequence and, alongside it,
/// ULTRA at maximum period `opt.ultra_period`; then merge them. ULTRA calls keep
/// ULTRA's units, scores and deltas; FasTAN calls get a consensus unit and deltas
/// from `opt.delta_model`.
///
/// Merging (unless `opt.prefer_fastan`): a FasTAN call with a period up to ULTRA's
/// maximum is dropped when ULTRA calls cover more than half of it. A FasTAN call with
/// a longer period, which ULTRA cannot report, is kept when its score under ULTRA's
/// model beats the ULTRA calls' deltas summed over the same bases, and then ULTRA
/// calls it covers more than half of are dropped (e.g. a microsatellite inside each
/// copy of a satellite). A fragment of one larger ULTRA call never wins this way.
/// Calls can still overlap in part.
///
/// Threads: `opt.threads` (0 = all cores) in total, split between ULTRA and FasTAN,
/// which run at the same time; with 1 they run one after the other.
///
/// Sequence names must be unique (the first word of each header, which is the
/// SequenceName written). ULTRA reads the file itself, so `input` must be a path;
/// `"-"` (stdin) works only when ULTRA isn't run (`ultra_period` 0 or `filler` given).
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
    let fo = Options { consensus: true, ..opt.fastan.clone() };
    let seqs = read_fasta(input)?;
    let mut seen = std::collections::HashSet::new();
    if let Some(r) = seqs.iter().find(|r| !seen.insert(r.name.as_str())) {
        anyhow::bail!("sequence name {:?} occurs more than once (names are the first word of each header)", r.name);
    }
    let fastan = || if opt.rescore_only { Vec::new() } else { annotate(&seqs, &fo) };
    let ultra = |p: usize| run_ultra(input, p, ultra_threads, opt.ultra_upstream_compat);
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

/// Rounded to six decimals, as ULTRA writes its numbers.
fn round6(x: f64) -> f64 {
    (x * 1e6).round() / 1e6
}

#[cfg(test)]
mod tests {
    use super::phase_to_start;

    #[test]
    fn phase_to_start_rotates_to_first_copy() {
        assert_eq!(phase_to_start(b"TGG".to_vec(), b"GGTGGTGGT"), b"GGT".to_vec());
        assert_eq!(phase_to_start(b"TG".to_vec(), b"gtgtgt"), b"GT".to_vec());
        assert_eq!(phase_to_start(b"ACGT".to_vec(), b"AC"), b"ACGT".to_vec());
    }
}
