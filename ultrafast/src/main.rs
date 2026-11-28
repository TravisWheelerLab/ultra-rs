use anyhow::Result;
use clap::Parser;
use std::io::Write;

/// ultrafast: hybrid tandem-repeat annotation with FasTAN and ULTRA.
///
/// FasTAN finds tandem arrays of any period; ULTRA runs in-process alongside it
/// at a short maximum period (-p, default 50). The result is written as ULTRA's
/// JSON (the fields of `ultra --json --show_delta`, Start 0-based), e.g. for
/// RepeatMasker's Aurora (-U). SequenceName is the first word of the header.
///
/// Merging: a FasTAN call with a period up to -p is dropped when ULTRA calls cover
/// more than half of it. A FasTAN call with a longer period, which ULTRA cannot
/// report, is kept when it scores better under ULTRA's model than the ULTRA calls it
/// overlaps, by 2 bits for every base its unit is longer than theirs. It then replaces
/// each ULTRA call that it covers more than half of and outscores on the bases they
/// share. Calls may still overlap. --prefer-fastan keeps every FasTAN call instead.
///
/// ULTRA calls keep ULTRA's own units, scores and per-base deltas. FasTAN calls
/// get a TRF-style iterative consensus unit, are split where the unit changes
/// (SubRepeats), and are scored per base with --delta-model:
///   ultra  ULTRA's default model, bits (the scale of ULTRA's own deltas)
///   te     Aurora's TE-alignment column scores, nats (needs --matrices)
#[derive(Parser)]
#[command(version)]
struct Args {
    /// FASTA input ('-' for stdin, only without ULTRA: -p 0 or --filler)
    input: String,
    /// Write the JSON here (default: stdout)
    #[arg(short = 'o', long, value_name = "FILE")]
    out: Option<String>,
    /// ULTRA's maximum period; 0 = don't run ULTRA (FasTAN only)
    #[arg(short = 'p', long, value_name = "PERIOD", default_value_t = 50, conflicts_with = "filler")]
    ultra_period: usize,
    /// Threads in total (0 = all cores), split between ULTRA and FasTAN, which run at
    /// the same time; with 1 they run one after the other
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
    /// Also write FasTAN's arrays as TSV (0-based half-open:
    /// seq start end period fastan_id cons_len copies cons_id rounds unit parts)
    #[arg(long, value_name = "FILE")]
    tsv: Option<String>,
    /// Also write the calls as BPAF, each aligned to its consensus unit (one record
    /// per call or sub-repeat; the units are stored in the file)
    #[arg(long, value_name = "FILE")]
    bpaf: Option<String>,

    /// Maximum consensus refinement rounds
    #[arg(long, default_value_t = 6, help_heading = "Consensus",
          value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..))]
    rounds: usize,
    /// Skip consensus for periods above this (the array keeps its middle copy)
    #[arg(long, default_value_t = 10_000, help_heading = "Consensus")]
    max_cons_period: usize,
    /// Don't split arrays where the unit changes (A...AB...B, same period): FasTAN's
    /// arrays, and ULTRA's calls beyond ULTRA's own splits
    #[arg(long, help_heading = "Consensus")]
    no_split: bool,
    /// Minimum edits a split must save
    #[arg(long, default_value_t = ultrafast::SplitParams::default().min_gain, help_heading = "Consensus",
          conflicts_with = "no_split")]
    split_min_gain: f64,
    /// ... and at least this fraction of the edits the two parts have against the parent unit
    #[arg(long, default_value_t = ultrafast::SplitParams::default().min_gain_frac, help_heading = "Consensus",
          conflicts_with = "no_split")]
    split_min_frac: f64,

    /// Score-delta model for FasTAN calls (ultra: the scale of ULTRA's own deltas)
    #[arg(long, value_enum, default_value_t = Model::Ultra, help_heading = "Score deltas")]
    delta_model: Model,
    /// Aurora-format matrices file (delta model te)
    #[arg(long, value_name = "FILE", help_heading = "Score deltas")]
    matrices: Option<String>,
    /// Use this matrix for every array (default: the <div>p<gc>g matrix nearest each
    /// array's divergence and its sequence's GC)
    #[arg(long, help_heading = "Score deltas")]
    matrix: Option<String>,
    /// GC percent for matrix choice, 0-100 (default: per sequence)
    #[arg(long, help_heading = "Score deltas")]
    gc: Option<f64>,
    /// Window for local base composition, bp (default 61, Aurora's)
    #[arg(long, help_heading = "Score deltas",
          value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..))]
    background_window: Option<usize>,
    /// Multiply every delta (and Score) by this: tunes how strongly tandem repeats
    /// compete with TE alignments in Aurora
    #[arg(long, default_value_t = 1.0, help_heading = "Score deltas")]
    delta_scale: f64,

    /// Take the short-period calls from this ULTRA-format JSON (e.g. a separate
    /// `ultra -p 20 --json --show_delta` run) instead of running ULTRA
    #[arg(long, value_name = "FILE", help_heading = "Merging")]
    filler: Option<String>,
    /// Minimum ULTRA score for a short-period call (default: keep all, as ULTRA does)
    #[arg(long, default_value_t = f64::NEG_INFINITY, allow_negative_numbers = true, help_heading = "Merging")]
    filler_min_score: f64,
    /// Don't join ULTRA calls of one period that lie less than a unit apart and have
    /// the same unit (ULTRA reports some long arrays in such pieces)
    #[arg(long, help_heading = "Merging")]
    no_join: bool,
    /// Where calls overlap, keep the FasTAN call instead: ULTRA calls covered more
    /// than half by FasTAN calls are dropped
    #[arg(long, help_heading = "Merging")]
    prefer_fastan: bool,
    /// Re-vote ULTRA's units with the TRF-style caller (seeded with the middle copy),
    /// as FasTAN's are (default: keep ULTRA's units)
    #[arg(long, help_heading = "Merging")]
    revote_filler: bool,
    /// Recompute the ULTRA calls' deltas with --delta-model instead of keeping ULTRA's
    /// own (always on for --delta-model te and --rescore-only)
    #[arg(long, help_heading = "Merging")]
    rescore_filler: bool,
    /// Don't run FasTAN: only rescore the --filler calls, e.g. to put stock ULTRA
    /// output on the TE scale
    #[arg(long, requires = "filler", help_heading = "Merging")]
    rescore_only: bool,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum Model {
    Ultra,
    Te,
}

fn main() -> Result<()> {
    let a = Args::parse();
    if matches!(a.delta_model, Model::Ultra) {
        for (set, name) in [
            (a.matrices.is_some(), "--matrices"),
            (a.matrix.is_some(), "--matrix"),
            (a.gc.is_some(), "--gc"),
            (a.background_window.is_some(), "--background-window"),
        ] {
            if set {
                anyhow::bail!("{name} applies only to --delta-model te");
            }
        }
    }
    // clap drops `requires = "filler"` when -p, which conflicts with --filler, is given
    if a.rescore_only && a.filler.is_none() {
        anyhow::bail!("--rescore-only needs --filler");
    }
    if !a.delta_scale.is_finite() {
        anyhow::bail!("--delta-scale must be a finite number, not {}", a.delta_scale);
    }
    if a.filler_min_score.is_nan() {
        anyhow::bail!("--filler-min-score must be a number");
    }
    if let Some(g) = a.gc {
        if !(0.0..=100.0).contains(&g) {
            anyhow::bail!("--gc must be a percentage from 0 to 100, not {g}");
        }
    }
    let delta_model = match a.delta_model {
        Model::Ultra => ultrafast::DeltaModel::Ultra,
        Model::Te => ultrafast::DeltaModel::Te {
            matrices: ultrafast::deltas::read_matrices(
                a.matrices.as_deref().ok_or_else(|| anyhow::anyhow!("--delta-model te needs --matrices"))?,
            )?,
            fixed: a.matrix.clone(),
            gc: a.gc,
            window: a.background_window.unwrap_or(61),
        },
    };
    let opt = ultrafast::HybridOptions {
        ultra_period: a.ultra_period,
        threads: a.threads,
        fastan: ultrafast::Options {
            consensus: true,
            cons: ultrafast::Params { rounds: a.rounds, ..Default::default() },
            max_cons_period: a.max_cons_period,
            split: (!a.no_split).then(|| ultrafast::SplitParams {
                min_gain: a.split_min_gain,
                min_gain_frac: a.split_min_frac,
                ..Default::default()
            }),
        },
        filler: a.filler.as_deref().map(ultrafast::read_ultra_json).transpose()?,
        filler_min_score: a.filler_min_score,
        prefer_fastan: a.prefer_fastan,
        join_ultra: !a.no_join,
        revote_filler: a.revote_filler,
        rescore_filler: a.rescore_filler,
        rescore_only: a.rescore_only,
        delta_model,
        delta_scale: a.delta_scale,
    };
    let h = ultrafast::hybrid(&a.input, &opt)?;
    if h.n_unmatched > 0 {
        let whose = if a.filler.is_some() { "--filler" } else { "ULTRA" };
        eprintln!("ultrafast: warning: {} {whose} calls name sequences not in {}; they were dropped", h.n_unmatched, a.input);
    }
    match a.out.as_deref() {
        // a regular file is written in full and then renamed into place; anything
        // else (e.g. /dev/null, a pipe) is written directly
        Some(path) if path != "-" && std::fs::metadata(path).map_or(true, |m| m.is_file()) => {
            let tmp = format!("{path}.partial");
            {
                let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
                h.write_json(&mut w)?;
                w.flush()?;
            }
            std::fs::rename(&tmp, path)?;
        }
        Some(path) if path != "-" => {
            let mut w = std::io::BufWriter::new(std::fs::OpenOptions::new().write(true).open(path)?);
            h.write_json(&mut w)?;
            w.flush()?;
        }
        _ => {
            let mut w = std::io::BufWriter::new(std::io::stdout().lock());
            h.write_json(&mut w)?;
            writeln!(w)?;
            w.flush()?;
        }
    }
    eprintln!(
        "ultrafast: {} FasTAN arrays ({} kept) + {} ULTRA calls kept",
        h.arrays.len(),
        h.repeats.len() - h.n_ultra,
        h.n_ultra
    );
    if let Some(path) = &a.tsv {
        write_tsv(path, &h.seqs, &h.arrays)?;
    }
    if let Some(path) = &a.bpaf {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
        // the alignments run on -t threads, as the calls did
        rayon::ThreadPoolBuilder::new().num_threads(a.threads).build()?.install(|| h.write_bpaf(&mut w))?;
        w.flush()?;
    }
    Ok(())
}

fn write_tsv(path: &str, recs: &[ultrafast::Record], calls: &[ultrafast::Call]) -> Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "#seq\tstart\tend\tperiod\tfastan_id\tcons_len\tcopies\tcons_id\trounds\tunit\tparts")?;
    for c in calls {
        let len = c.end - c.start;
        let (cid, rounds) = match &c.cons {
            Some(k) => (format!("{:.4}", 1.0 - k.ed as f64 / len as f64), k.rounds.to_string()),
            None => (".".into(), ".".into()),
        };
        let parts = if c.parts.is_empty() {
            ".".to_string()
        } else {
            c.parts.iter().map(|(o, u)| format!("{}:{}", o, String::from_utf8_lossy(u))).collect::<Vec<_>>().join(",")
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
            recs[c.seq].name,
            c.start,
            c.end,
            c.period,
            c.fastan_identity(),
            c.unit.len(),
            len as f64 / c.unit.len().max(1) as f64,
            cid,
            rounds,
            String::from_utf8_lossy(&c.unit),
            parts
        )?;
    }
    w.flush()?;
    Ok(())
}
