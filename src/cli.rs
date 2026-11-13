//! Command-line settings (`cli.hpp` / `cli.cpp`).
//!
//! The C++ code configures a vendored CLI11 parser inside `Settings`. Here the
//! parser is a small hand-written one that reproduces the CLI11 2.3.2 behaviour
//! ULTRA depends on (see `parse_args`), because the observable quirks (a
//! repeated option throwing "At Most 1 required", `--acgt` swallowing every
//! following non-option argument, float "flags", conversion through
//! `strtoll`/`strtold`, alias options being independent) are simpler to mirror
//! exactly than to coax out of clap. `Settings` stays plain data.

use std::io::Write;

/// The version this program reports in `--help` and in JSON output. It is
/// the version of the crate, set in `Cargo.toml`.
pub const ULTRA_VERSION_STRING: &str = env!("CARGO_PKG_VERSION");

/// Text CLI11 prints at the top of `--help` (the `CLI::App` description).
pub const ULTRA_BANNER: &str = concat!(
    "\n",
    "=================================================\n",
    "(U)ltra (L)ocates (T)andemly (R)epetitive (A)reas\n",
    "     Daniel R. Olson and Travis J. Wheeler\n",
    "                 Version ",
    env!("CARGO_PKG_VERSION"),
    "\n",
    "     Use '--cite' for citation instructions\n",
    "=================================================\n",
);

const CITE_TEXT: &str = "BibTeX: \n\
@article{10.1093/bioadv/vbae149,\n\
\x20   author = {Olson, Daniel R and Wheeler, Travis J},\n\
\x20   title = {ULTRA-effective labeling of tandem repeats in genomic sequence},\n\
\x20   journal = {Bioinformatics Advances},\n\
\x20   volume = {4},\n\
\x20   number = {1},\n\
\x20   pages = {vbae149},\n\
\x20   year = {2024},\n\
\x20   month = {10},\n\
\x20   issn = {2635-0041},\n\
\x20   doi = {10.1093/bioadv/vbae149},\n\
\x20   url = {https://doi.org/10.1093/bioadv/vbae149},\n\
\x20   eprint = {https://academic.oup.com/bioinformaticsadvances/article-pdf/4/1/vbae149/60779841/vbae149.pdf},\n\
}\n";

const HELP_BODY: &str = r####"Usage: {PROG} [OPTIONS] [input_file]

Positionals:
  input_file TEXT             DNA FASTA input file

Options:
  -h,--help                   Print this help message and exit


Input:
  -r,--read_all               Read entire input file into memory (disables streaming input)


Output:
  -o,--out TEXT               Output file path
  --disable_streaming_out     Disables streaming output; no output will be created until all analysis has been completed
  -c,--show_counts            Output #copies, #substitutions, #insertions, #deletions
  --pval                      Use p-values instead of scores in BED output
  --pval_loc                  The exponential location used for converting scores to p-values.
  --pval_scale                The exponential scale used for converting scores to p-values
  --tsv                       Use TSV output format
  --json                      Use JSON output format
  --bed                       Use BED output format
  --max_consensus INT [1000000] 
                              The maximum length of consensus pattern to include in output
  --show_seq                  Output repetitive region
  --show_delta                Show positional score deltas in JSON output
  --show_trace                Show Viterbi trace in JSON output
  --show_wid                  Show sequence window IDs in JSON output
  --show_logo                 Show logo numbers in JSON output
  --hs,--hide_settings        Do not output settings
  --suppress                  Do not output BED or JSON annotation
  --disable_summary           Disable summary statistics
  --fdr                       Estimate the False Discovery rate (runtime will be twice as long)


Masking:
  --mask TEXT                 File path to save a masked FASTA
  --nmask                     Use n masking instead of case masking


System:
  -t,--threads INT [1]        Number of threads to use
  --win_size INT              Manually set sequence window size
  --overlap INT               Manually set sequence window overlap size
  --windows INT [1024]        Number of sequence windows to store at once
  --mem                       Display memory requirements for current settings and exit


Filter:
  -s,--score FLOAT [-100.0]   Minimum reportable repeat score
  --min_length UINT [10]      Minimum reportable repeat length
  --min_unit UINT [2]         Minimum reportable number of repeat units


Parameter Tuning:
  --tune                      Tune parameters using a small search grid before running (see README)
  --tune_medium               Tune parameters before running (see README)
  --tune_large                Tune parameters using a larger search grid before running (see README)
  --tune_file TEXT            Use custom parameter search during tuning (see README)
  --tune_indel                Enable indels while tuning
  --tune_fdr FLOAT [0.05]     FDR to be tuned against (see README)
  --tune_only                 Tune parameters and don't run (see README)


Model:
  -p,--period UINT [100]      Maximum detectable repeat period
  -i,--inserts UINT [10]      Maximum number of insertion states
  -d,--deletes UINT [10]      Maximum number of deletion states


Probabilities:
  --at FLOAT [0.5]            Expected AT content
  --acgt FLOAT x 4            Expected frequency of A C G T
  -m,--match FLOAT [0.7]      Expected conservation in repeats
  --decay FLOAT [0.85]        Decay penalty applied to repetitive states
  --nr FLOAT [0.01]           Probability of a repeat starting
  --rn FLOAT [0.05]           Probability of a repeat stopping
  --ri FLOAT [0.02]           Probability of an insertion occurring
  --rd FLOAT [0.02]           Probability of a deletion occurring
  --ii FLOAT [0.02]           Probability of consecutive insertions occurring
  --dd FLOAT [0.02]           Probability of consecutive deletions occurring


Splitting and Naming:
  --no_split                  Do not perform repeat splitting
  --max_split INT [10]        The maximum repeat period to perform repeat splitting
  --split_threshold FLOAT [3.5] 
                              Split threshold value
  --split_depth UINT [5]      Number of repeat units to use in repeat splitting
  --min_split_window UINT [15] 
                              Minimum repeat split window size

For additional information see README
Supported by: NIH NIGMS P20GM103546 and NIH NHGRI U24HG010136

"####;

/// All ULTRA settings. Field names, types and defaults follow the C++ struct.
#[derive(Clone, Debug)]
pub struct Settings {
    pub args: String,

    pub cite: bool,

    pub in_file: String,
    pub read_all: bool,

    pub out_file: String,
    pub pval: bool,
    pub disable_summary: bool,
    pub disable_streaming_out: bool,
    pub p_value_loc: f32,
    pub p_value_scale: f32,
    pub p_value_freq: f32,
    pub estimate_fdr: bool,

    pub max_consensus_period: i32,
    pub ultra_out: bool,
    pub bed_out: bool,
    pub json_out: bool,

    pub show_counts: bool,
    pub show_seq: bool,
    pub show_deltas: bool,
    pub show_trace: bool,
    pub show_wid: bool,
    pub show_logo_nums: bool,
    pub hide_settings: bool,
    pub suppress_out: bool,

    pub produce_mask: bool,
    pub mask_file: String,
    pub mask_with_n: bool,

    pub threads: i32,
    pub window_size: i64,
    pub overlap: i64,
    pub windows: i64,
    pub show_memory: bool,

    pub min_score: f32,
    pub min_length: u64,
    pub min_units: u64,

    pub tune_fdr: f64,
    pub tune: bool,
    pub tune_medium: bool,
    pub tune_large: bool,
    pub tune_only: bool,
    pub tune_indels: bool,
    pub tune_param_path: String,

    pub max_period: u64,
    pub max_insert: u64,
    pub max_delete: u64,

    pub a_freq: f32,
    pub c_freq: f32,
    pub g_freq: f32,
    pub t_freq: f32,
    pub acgt: Vec<f32>,
    pub at: f32,
    pub match_probability: f32,

    pub period_decay: f32,
    pub transition_nr: f32,
    pub transition_rn: f32,
    pub transition_ri: f32,
    pub transition_rd: f32,
    pub transition_ii: f32,
    pub transition_dd: f32,

    pub no_split: bool,
    pub max_split: i32,
    pub split_threshold: f32,
    pub split_depth: u64,
    pub min_split_window: u64,
    pub max_namable_period: u64,
    pub max_highfi_naming: u64,

    pub run_without_reader: bool,

    /// Set by `set_multi_option`: when an option is given more than once the
    /// last occurrence wins instead of raising an error (CLI11 `TakeLast`).
    /// Stands in for the parser state the C++ struct keeps inside `CLI::App`.
    pub multi_option: bool,

    /// Build each consensus as the C++ does (ULTRA 1.2.2 with the fixes of
    /// pull request 95) instead of with the indel test and realignment of
    /// `RepeatRegion::create_logo`. No command-line option sets this. The
    /// `ultra` binary turns it on when the environment variable
    /// `ULTRA_UPSTREAM_CONSENSUS` is set to anything but `0` or the empty
    /// string, which is how the differential tests compare the rest of the
    /// output with the C++ byte for byte.
    pub upstream_consensus: bool,
}

impl Default for Settings {
    fn default() -> Self {
        // Unsuffixed C++ literals assigned to float are doubles narrowed to float.
        Settings {
            args: String::new(),
            cite: false,
            in_file: String::new(),
            read_all: false,
            out_file: String::new(),
            pval: false,
            disable_summary: false,
            disable_streaming_out: false,
            p_value_loc: 4.27294921875_f64 as f32,
            p_value_scale: 1.8913602828979492_f64 as f32,
            p_value_freq: 0.000330256_f64 as f32,
            estimate_fdr: false,
            max_consensus_period: 1000000,
            ultra_out: false,
            bed_out: false,
            json_out: false,
            show_counts: false,
            show_seq: false,
            show_deltas: false,
            show_trace: false,
            show_wid: false,
            show_logo_nums: false,
            hide_settings: false,
            suppress_out: false,
            produce_mask: false,
            mask_file: String::new(),
            mask_with_n: false,
            threads: 1,
            window_size: -1,
            overlap: -1,
            windows: 1024,
            show_memory: false,
            min_score: -100.0_f64 as f32,
            min_length: 10,
            min_units: 2,
            tune_fdr: 0.05,
            tune: false,
            tune_medium: false,
            tune_large: false,
            tune_only: false,
            tune_indels: false,
            tune_param_path: String::new(),
            max_period: 100,
            max_insert: 10,
            max_delete: 10,
            a_freq: 0.25_f64 as f32,
            c_freq: 0.25_f64 as f32,
            g_freq: 0.25_f64 as f32,
            t_freq: 0.25_f64 as f32,
            acgt: Vec::new(),
            at: 0.5_f64 as f32,
            match_probability: 0.7_f64 as f32,
            period_decay: 0.85_f64 as f32,
            transition_nr: 0.01_f64 as f32,
            transition_rn: 0.05_f64 as f32,
            transition_ri: 0.02_f64 as f32,
            transition_rd: 0.02_f64 as f32,
            transition_ii: 0.02_f64 as f32,
            transition_dd: 0.02_f64 as f32,
            no_split: false,
            max_split: 10,
            split_threshold: 3.5_f64 as f32,
            split_depth: 5,
            min_split_window: 15,
            max_namable_period: 50,
            max_highfi_naming: 20,
            run_without_reader: false,
            multi_option: false,
            upstream_consensus: false,
        }
    }
}

/// Why `try_parse_input` stopped before returning a pass/fail verdict. Each
/// variant is a place where the C++ either calls `exit` or lets a CLI11
/// exception escape to `main`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CliExit {
    /// `--help`: print the text to stdout, exit 0.
    Help(String),
    /// Unrecognised argument: print "Unrecognized flag or argument: <msg>" to stderr, exit 0.
    Extras(String),
    /// `--cite`: print the BibTeX entry to stdout, exit 0.
    Cite,
    /// A CLI11 conversion/mismatch exception: `main` prints
    /// "Standard exception caught: <msg>" to stderr and returns -1 (exit 255).
    Exception(String),
}

// ---------------------------------------------------------------------------
// Value conversion (CLI11 `lexical_cast`)
// ---------------------------------------------------------------------------

/// C `isspace` in the "C" locale.
fn c_isspace(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Result of `strtoull`/`strtoll` with base 0: value, bytes consumed (0 if no
/// digits), and whether the magnitude overflowed 64 bits.
fn strto_base0(s: &[u8]) -> (u64, bool, usize, bool) {
    let mut i = 0;
    while i < s.len() && c_isspace(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut radix = 10u32;
    if i + 2 < s.len()
        && s[i] == b'0'
        && (s[i + 1] == b'x' || s[i + 1] == b'X')
        && (s[i + 2] as char).is_ascii_hexdigit()
    {
        radix = 16;
        i += 2;
    } else if i < s.len() && s[i] == b'0' {
        radix = 8;
    }
    let start = i;
    let mut mag: u128 = 0;
    let mut overflow = false;
    while i < s.len() {
        let d = match (s[i] as char).to_digit(radix) {
            Some(d) => d,
            None => break,
        };
        if !overflow {
            mag = mag * radix as u128 + d as u128;
            if mag > u64::MAX as u128 {
                overflow = true;
            }
        }
        i += 1;
    }
    if i == start {
        // No digits: glibc consumes nothing (a lone "0x" is handled above as "0").
        return (0, neg, 0, false);
    }
    (if overflow { u64::MAX } else { mag as u64 }, neg, i, overflow)
}

/// CLI11 integral conversion for unsigned 64-bit targets.
fn lex_u64(s: &str) -> Option<u64> {
    let b = s.as_bytes();
    if b.is_empty() || b[0] == b'-' {
        return None;
    }
    let (mag, neg, used, overflow) = strto_base0(b);
    if overflow {
        return None; // ERANGE
    }
    if used == b.len() {
        // strtoull negates modulo 2^64 (reachable via leading whitespace, e.g. " -5").
        return Some(if neg { mag.wrapping_neg() } else { mag });
    }
    None
}

/// `strtoll` result, or None on ERANGE / no conversion.
fn strtoll_full(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    let (mag, neg, used, overflow) = strto_base0(b);
    if overflow || used != b.len() {
        return None;
    }
    if neg {
        if mag > (i64::MAX as u64) + 1 {
            None
        } else {
            Some((mag as i64).wrapping_neg())
        }
    } else if mag > i64::MAX as u64 {
        None
    } else {
        Some(mag as i64)
    }
}

/// CLI11 integral conversion for signed targets (`"true"` becomes 1).
fn lex_i64(s: &str) -> Option<i64> {
    if s.is_empty() {
        return None;
    }
    if let Some(v) = strtoll_full(s) {
        return Some(v);
    }
    if s == "true" {
        return Some(1);
    }
    None
}

fn lex_i32(s: &str) -> Option<i32> {
    if s.is_empty() {
        return None;
    }
    match strtoll_full(s) {
        Some(v) => {
            if v == v as i32 as i64 {
                Some(v as i32)
            } else if s == "true" {
                Some(1)
            } else {
                None
            }
        }
        None => {
            if s == "true" {
                Some(1)
            } else {
                None
            }
        }
    }
}

/// Parse a hexadecimal floating literal body (after the sign), as accepted by `strtod`.
fn parse_hex_float(s: &str) -> Option<f64> {
    let b = s.as_bytes();
    if b.len() < 3 || b[0] != b'0' || (b[1] != b'x' && b[1] != b'X') {
        return None;
    }
    let mut i = 2;
    let mut mant: f64 = 0.0;
    let mut digits = 0;
    let mut exp: i32 = 0;
    while i < b.len() && (b[i] as char).is_ascii_hexdigit() {
        mant = mant * 16.0 + (b[i] as char).to_digit(16).unwrap() as f64;
        digits += 1;
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && (b[i] as char).is_ascii_hexdigit() {
            mant = mant * 16.0 + (b[i] as char).to_digit(16).unwrap() as f64;
            exp -= 4;
            digits += 1;
            i += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if i < b.len() && (b[i] == b'p' || b[i] == b'P') {
        let mut j = i + 1;
        let mut eneg = false;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            eneg = b[j] == b'-';
            j += 1;
        }
        let st = j;
        let mut e: i32 = 0;
        while j < b.len() && b[j].is_ascii_digit() {
            e = e.saturating_mul(10).saturating_add((b[j] - b'0') as i32);
            j += 1;
        }
        if j == st {
            return None;
        }
        exp = exp.saturating_add(if eneg { -e } else { e });
        i = j;
    }
    if i != b.len() {
        return None;
    }
    Some(mant * 2f64.powi(exp))
}

/// Strip what `strtod` tolerates before the number (leading whitespace) and split the sign.
fn float_prefix(s: &str) -> Option<(bool, &str)> {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && c_isspace(b[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        neg = b[i] == b'-';
        i += 1;
    }
    let rest = &s[i..];
    if rest.is_empty() {
        None
    } else {
        Some((neg, rest))
    }
}

/// CLI11 float conversion (`strtold`, then narrowing to `f32`).
fn lex_f32(s: &str) -> Option<f32> {
    if s.is_empty() {
        return None;
    }
    let (neg, rest) = float_prefix(s)?;
    let v = if rest.len() > 1 && rest.as_bytes()[0] == b'0' && (rest.as_bytes()[1] | 0x20) == b'x' {
        parse_hex_float(rest)? as f32
    } else if rest.as_bytes()[0] == b'+' || rest.as_bytes()[0] == b'-' {
        return None;
    } else {
        rest.parse::<f32>().ok()?
    };
    Some(if neg { -v } else { v })
}

/// CLI11 float conversion for `double` targets.
fn lex_f64(s: &str) -> Option<f64> {
    if s.is_empty() {
        return None;
    }
    let (neg, rest) = float_prefix(s)?;
    let v = if rest.len() > 1 && rest.as_bytes()[0] == b'0' && (rest.as_bytes()[1] | 0x20) == b'x' {
        parse_hex_float(rest)?
    } else if rest.as_bytes()[0] == b'+' || rest.as_bytes()[0] == b'-' {
        return None;
    } else {
        rest.parse::<f64>().ok()?
    };
    Some(if neg { -v } else { v })
}

/// CLI11 `to_flag_value`; `Err(false)` is std::invalid_argument, `Err(true)` is std::out_of_range.
fn to_flag_value(val: &str) -> Result<i64, bool> {
    if val == "true" {
        return Ok(1);
    }
    if val == "false" {
        return Ok(-1);
    }
    let low = val.to_ascii_lowercase();
    let lb = low.as_bytes();
    if lb.len() == 1 {
        return match lb[0] {
            c @ b'1'..=b'9' => Ok((c - b'0') as i64),
            b'0' | b'f' | b'n' | b'-' => Ok(-1),
            b't' | b'y' | b'+' => Ok(1),
            _ => Err(false),
        };
    }
    match low.as_str() {
        "true" | "on" | "yes" | "enable" => Ok(1),
        "false" | "off" | "no" | "disable" => Ok(-1),
        _ => {
            // std::stoll: skip whitespace, sign, decimal digits, ignore the rest.
            let mut i = 0;
            while i < lb.len() && c_isspace(lb[i]) {
                i += 1;
            }
            let mut neg = false;
            if i < lb.len() && (lb[i] == b'+' || lb[i] == b'-') {
                neg = lb[i] == b'-';
                i += 1;
            }
            let st = i;
            let mut mag: u128 = 0;
            while i < lb.len() && lb[i].is_ascii_digit() {
                mag = (mag * 10 + (lb[i] - b'0') as u128).min(u128::from(u64::MAX) * 4);
                i += 1;
            }
            if i == st {
                return Err(false);
            }
            if neg {
                if mag > (i64::MAX as u128) + 1 {
                    Err(true)
                } else {
                    Ok(-(mag as i128) as i64)
                }
            } else if mag > i64::MAX as u128 {
                Err(true)
            } else {
                Ok(mag as i64)
            }
        }
    }
}

/// CLI11 boolean conversion.
fn lex_bool(s: &str) -> Option<bool> {
    match to_flag_value(s) {
        Ok(v) => Some(v > 0),
        Err(false) => None,
        Err(true) => Some(!s.starts_with('-')),
    }
}

/// Types an option can be converted to, with CLI11's empty-string rule
/// (an empty argument assigns the type's zero value).
trait Lex: Default {
    fn lex(s: &str) -> Option<Self>;
}
impl Lex for i32 {
    fn lex(s: &str) -> Option<Self> {
        lex_i32(s)
    }
}
impl Lex for i64 {
    fn lex(s: &str) -> Option<Self> {
        lex_i64(s)
    }
}
impl Lex for u64 {
    fn lex(s: &str) -> Option<Self> {
        lex_u64(s)
    }
}
impl Lex for f32 {
    fn lex(s: &str) -> Option<Self> {
        lex_f32(s)
    }
}
impl Lex for f64 {
    fn lex(s: &str) -> Option<Self> {
        lex_f64(s)
    }
}
impl Lex for String {
    fn lex(s: &str) -> Option<Self> {
        Some(s.to_string())
    }
}

fn lex_assign<T: Lex>(s: &str) -> Option<T> {
    if s.is_empty() {
        Some(T::default())
    } else {
        T::lex(s)
    }
}

// ---------------------------------------------------------------------------
// Option table (registration order of `prepare_settings`)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Help,
    Positional,
    /// Takes no value; boolean target.
    Flag,
    /// A flag whose target is a float (`--pval_loc`, `--pval_scale`): only usable as `--name=number`.
    FlagF32,
    /// Takes exactly one value.
    Value,
    /// `--acgt`: one or more values, exactly 4 expected.
    Multi4,
}

struct OptDef {
    shorts: &'static [u8],
    longs: &'static [&'static str],
    kind: Kind,
    /// CLI11 group "" (alias options): their name is empty in error messages.
    hidden: bool,
    type_name: &'static str,
    set: fn(&mut Settings, &[String]) -> bool,
}

impl OptDef {
    /// CLI11 `Option::get_name()` as used in error messages.
    fn name(&self) -> String {
        if self.hidden {
            String::new()
        } else if let Some(l) = self.longs.first() {
            format!("--{}", l)
        } else if let Some(s) = self.shorts.first() {
            format!("-{}", *s as char)
        } else {
            "input_file".to_string()
        }
    }
}

macro_rules! flag {
    ($s:expr, $l:expr, $hid:expr, $f:ident) => {
        OptDef {
            shorts: $s,
            longs: $l,
            kind: Kind::Flag,
            hidden: $hid,
            type_name: "",
            set: |st, r| match lex_bool(&r[0]) {
                Some(v) => {
                    st.$f = v;
                    true
                }
                None => false,
            },
        }
    };
}
macro_rules! flagf {
    ($s:expr, $l:expr, $f:ident) => {
        OptDef {
            shorts: $s,
            longs: $l,
            kind: Kind::FlagF32,
            hidden: false,
            type_name: "",
            set: |st, r| match lex_f32(&r[0]) {
                Some(v) => {
                    st.$f = v;
                    true
                }
                None => false,
            },
        }
    };
}
macro_rules! val {
    ($s:expr, $l:expr, $hid:expr, $f:ident, $ty:ty, $tn:expr) => {
        OptDef {
            shorts: $s,
            longs: $l,
            kind: Kind::Value,
            hidden: $hid,
            type_name: $tn,
            set: |st, r| match lex_assign::<$ty>(&r[0]) {
                Some(v) => {
                    st.$f = v;
                    true
                }
                None => false,
            },
        }
    };
}

fn option_table() -> Vec<OptDef> {
    vec![
        OptDef {
            shorts: b"h",
            longs: &["help"],
            kind: Kind::Help,
            hidden: false,
            type_name: "",
            set: |_, _| true,
        },
        // Input
        OptDef {
            shorts: b"",
            longs: &[],
            kind: Kind::Positional,
            hidden: false,
            type_name: "TEXT",
            set: |st, r| {
                st.in_file = r[0].clone();
                true
            },
        },
        flag!(b"r", &["read_all"], false, read_all),
        flag!(b"", &["readall"], true, read_all),
        // Output
        val!(b"o", &["out"], false, out_file, String, "TEXT"),
        flag!(b"", &["disable_streaming_out"], false, disable_streaming_out),
        flag!(b"c", &["show_counts"], false, show_counts),
        flag!(b"", &["pval"], false, pval),
        flagf!(b"", &["pval_loc"], p_value_loc),
        flagf!(b"", &["pval_scale"], p_value_scale),
        flag!(b"", &["tsv"], false, ultra_out),
        flag!(b"", &["json"], false, json_out),
        flag!(b"", &["bed"], false, bed_out),
        val!(b"", &["max_consensus"], false, max_consensus_period, i32, "INT"),
        flag!(b"", &["show_seq"], false, show_seq),
        flag!(b"", &["show_delta"], false, show_deltas),
        flag!(b"", &["showdelta"], true, show_deltas),
        flag!(b"", &["show_trace"], false, show_trace),
        flag!(b"", &["showtrace"], true, show_trace),
        flag!(b"", &["show_wid"], false, show_wid),
        flag!(b"", &["showwid"], true, show_wid),
        flag!(b"", &["show_logo"], false, show_logo_nums),
        flag!(b"", &["showlogo"], true, show_logo_nums),
        flag!(b"", &["hs", "hide_settings"], false, hide_settings),
        flag!(b"", &["hidesettings"], true, hide_settings),
        flag!(b"", &["suppress"], false, suppress_out),
        flag!(b"", &["disable_summary"], false, disable_summary),
        flag!(b"", &["fdr"], false, estimate_fdr),
        // Masking
        val!(b"", &["mask"], false, mask_file, String, "TEXT"),
        flag!(b"", &["nmask"], false, mask_with_n),
        // System
        val!(b"t", &["threads"], false, threads, i32, "INT"),
        val!(b"", &["win_size"], false, window_size, i64, "INT"),
        val!(b"", &["winsize"], true, window_size, i64, "INT"),
        val!(b"", &["overlap"], false, overlap, i64, "INT"),
        val!(b"", &["windows"], false, windows, i64, "INT"),
        flag!(b"", &["mem"], false, show_memory),
        // Filter
        val!(b"s", &["score"], false, min_score, f32, "FLOAT"),
        val!(b"", &["min_length"], false, min_length, u64, "UINT"),
        val!(b"", &["minlen"], true, min_length, u64, "UINT"),
        val!(b"", &["min_unit"], false, min_units, u64, "UINT"),
        val!(b"", &["minunit"], true, min_units, u64, "UINT"),
        // Tune
        flag!(b"", &["tune"], false, tune),
        flag!(b"", &["tune_medium"], false, tune_medium),
        flag!(b"", &["tune_large"], false, tune_large),
        val!(b"", &["tune_file"], false, tune_param_path, String, "TEXT"),
        flag!(b"", &["tune_indel"], false, tune_indels),
        flag!(b"", &["tune_indels"], true, tune_indels),
        val!(b"", &["tune_fdr"], false, tune_fdr, f64, "FLOAT"),
        flag!(b"", &["tune_only"], false, tune_only),
        // Model
        val!(b"p", &["period"], false, max_period, u64, "UINT"),
        val!(b"i", &["inserts"], false, max_insert, u64, "UINT"),
        val!(b"d", &["deletes"], false, max_delete, u64, "UINT"),
        // Probabilities
        val!(b"", &["at"], false, at, f32, "FLOAT"),
        OptDef {
            shorts: b"",
            longs: &["acgt"],
            kind: Kind::Multi4,
            hidden: false,
            type_name: "FLOAT",
            set: |st, r| {
                let mut v = Vec::new();
                if !(r.len() == 1 && r[0] == "{}") {
                    for e in r {
                        match lex_assign::<f32>(e) {
                            Some(x) => v.push(x),
                            None => return false,
                        }
                    }
                }
                st.acgt = v;
                true
            },
        },
        val!(b"m", &["match"], false, match_probability, f32, "FLOAT"),
        val!(b"", &["decay"], false, period_decay, f32, "FLOAT"),
        val!(b"", &["nr"], false, transition_nr, f32, "FLOAT"),
        val!(b"", &["rn"], false, transition_rn, f32, "FLOAT"),
        val!(b"", &["ri"], false, transition_ri, f32, "FLOAT"),
        val!(b"", &["rd"], false, transition_rd, f32, "FLOAT"),
        val!(b"", &["ii"], false, transition_ii, f32, "FLOAT"),
        val!(b"", &["dd"], false, transition_dd, f32, "FLOAT"),
        // Splitting
        flag!(b"", &["no_split"], false, no_split),
        flag!(b"", &["nosplit"], true, no_split),
        val!(b"", &["max_split"], false, max_split, i32, "INT"),
        val!(b"", &["maxsplit"], true, max_split, i32, "INT"),
        val!(b"", &["split_threshold"], false, split_threshold, f32, "FLOAT"),
        val!(b"", &["splitval"], true, split_threshold, f32, "FLOAT"),
        val!(b"", &["split_depth"], false, split_depth, u64, "UINT"),
        val!(b"", &["splitdepth"], true, split_depth, u64, "UINT"),
        val!(b"", &["min_split_window"], false, min_split_window, u64, "UINT"),
        val!(b"", &["minsplitwin"], true, min_split_window, u64, "UINT"),
        flag!(b"", &["cite"], true, cite),
    ]
}

/// CLI11 `valid_first_char`.
fn valid_first_char(c: u8) -> bool {
    c != b'-' && c != b'!' && c != b' ' && c != b'\n'
}

enum Cls<'a> {
    Mark,
    /// `--name` or `--name=value`
    Long(&'a str, &'a str),
    /// `-n` plus the rest of the argument
    Short(u8, &'a str),
    None,
}

/// CLI11 `_recognize`. Digits after a single dash never name an option (no
/// option here has a digit short name), so `-5` is a positional value.
fn classify(s: &str) -> Cls<'_> {
    let b = s.as_bytes();
    if s == "--" {
        return Cls::Mark;
    }
    if b.len() > 2 && &b[..2] == b"--" && valid_first_char(b[2]) {
        return match s.find('=') {
            Some(loc) => Cls::Long(&s[2..loc], &s[loc + 1..]),
            None => Cls::Long(&s[2..], ""),
        };
    }
    if b.len() > 1 && b[0] == b'-' && valid_first_char(b[1]) {
        if b[1].is_ascii_digit() {
            return Cls::None;
        }
        // Non-ASCII second byte cannot match any option; keep the raw remainder lossless.
        let rest = s.get(2..).unwrap_or("");
        return Cls::Short(b[1], rest);
    }
    Cls::None
}

/// Parse `argv` (program name first) into `settings`, mirroring CLI11's
/// `App::parse` for the option set above. Errors are returned as `CliExit`.
fn parse_args(settings: &mut Settings, argv: &[String]) -> Result<(), CliExit> {
    let defs = option_table();
    let multi = settings.multi_option;
    let mut results: Vec<Vec<String>> = vec![Vec::new(); defs.len()];
    // Arguments CLI11 could not place; the flag marks a "--" separator.
    let mut missing: Vec<(bool, String)> = Vec::new();

    let find = |long: bool, name: &str| -> Option<usize> {
        defs.iter().position(|o| {
            if long {
                o.longs.iter().any(|l| *l == name)
            } else {
                name.len() == 1 && o.shorts.contains(&name.as_bytes()[0])
            }
        })
    };

    let args: &[String] = if argv.is_empty() { &[] } else { &argv[1..] };
    let mut i = 0;
    let mut positional_only = false;
    while i < args.len() {
        let current = &args[i];
        i += 1;
        let cls = if positional_only { Cls::None } else { classify(current) };
        let (idx, value) = match cls {
            Cls::Mark => {
                positional_only = true;
                missing.push((true, "--".to_string()));
                continue;
            }
            Cls::None => {
                // Positional: the first free positional option takes it, otherwise it is extra.
                match defs.iter().position(|o| o.kind == Kind::Positional) {
                    Some(p) if results[p].is_empty() => results[p].push(current.clone()),
                    _ => missing.push((false, current.clone())),
                }
                continue;
            }
            Cls::Long(name, value) => match find(true, name) {
                Some(p) => (p, value),
                None => {
                    missing.push((false, current.clone()));
                    continue;
                }
            },
            Cls::Short(name, _rest) => {
                let key = [name];
                match std::str::from_utf8(&key).ok().and_then(|n| find(false, n)) {
                    // `rest` is never non-empty for a matching short option: the
                    // argument preprocessor rejects `-xyz` before parsing.
                    Some(p) => (p, ""),
                    None => {
                        missing.push((false, current.clone()));
                        continue;
                    }
                }
            }
        };
        let def = &defs[idx];
        match def.kind {
            Kind::Help | Kind::Flag | Kind::FlagF32 => {
                // get_flag_value: an empty value means "true".
                results[idx].push(if value.is_empty() {
                    "true".to_string()
                } else {
                    value.to_string()
                });
            }
            Kind::Value | Kind::Multi4 => {
                let mut collected = 0;
                if !value.is_empty() {
                    results[idx].push(value.to_string());
                    collected += 1;
                }
                // The minimum (one value) is taken from the next argument whatever it looks like.
                if collected < 1 {
                    if i < args.len() {
                        results[idx].push(args[i].clone());
                        i += 1;
                    } else {
                        return Err(CliExit::Exception(format!(
                            "{}: 1 required {} missing",
                            def.name(),
                            def.type_name
                        )));
                    }
                }
                if def.kind == Kind::Multi4 {
                    // allow_extra_args: swallow every following argument that is not an option.
                    while i < args.len() && matches!(classify(&args[i]), Cls::None) {
                        results[idx].push(args[i].clone());
                        i += 1;
                    }
                    // "--" ends the list and is eaten.
                    if i < args.len() && args[i] == "--" {
                        i += 1;
                    }
                }
            }
            Kind::Positional => unreachable!("positionals are not named"),
        }
    }

    // _process_callbacks: options in registration order; errors surface in that order.
    for (idx, def) in defs.iter().enumerate() {
        let all = &results[idx];
        if all.is_empty() {
            continue;
        }
        let (items_min, items_max) = match def.kind {
            Kind::Multi4 => (4usize, 4usize),
            Kind::Value | Kind::Positional => (1, 1),
            _ => (0, 0),
        };
        let mut send: &[String] = all;
        let throw_policy = !multi && matches!(def.kind, Kind::Value | Kind::Multi4 | Kind::Positional);
        if throw_policy {
            if all.len() < items_min {
                return Err(CliExit::Exception(format!(
                    "{}: At least {} required but received {}",
                    def.name(),
                    items_min,
                    all.len()
                )));
            }
            if all.len() > items_max {
                return Err(CliExit::Exception(format!(
                    "{}: At Most {} required but received {}",
                    def.name(),
                    items_max,
                    all.len()
                )));
            }
        } else {
            // TakeLast: keep the final max(items_max, 1) results.
            let trim = items_max.max(1).min(all.len());
            send = &all[all.len() - trim..];
        }
        if !(def.set)(settings, send) {
            return Err(CliExit::Exception(format!(
                "Could not convert: {} = {}",
                def.name(),
                all.join(",")
            )));
        }
    }

    // _process_help_flags
    if !results[0].is_empty() {
        let prog = argv.first().map(String::as_str).unwrap_or("");
        return Err(CliExit::Help(help_text(prog)));
    }

    // _process_extras
    if missing.iter().any(|(mark, _)| !*mark) {
        let list: Vec<&str> = missing.iter().rev().map(|(_, s)| s.as_str()).collect();
        let lead = if missing.len() > 1 {
            "The following arguments were not expected: "
        } else {
            "The following argument was not expected: "
        };
        return Err(CliExit::Extras(format!("{}{}", lead, list.join(" "))));
    }
    Ok(())
}

/// `app.help()`: banner, usage line for `prog`, then the option listing.
pub fn help_text(prog: &str) -> String {
    format!("{}\n{}", ULTRA_BANNER, HELP_BODY.replace("{PROG}", prog))
}

/// `std::to_string(float)`: `%f` of the value promoted to double.
fn to_string_f(v: f32) -> String {
    fmt_f(v as f64)
}

/// printf `%f`, including how glibc spells inf/nan.
fn fmt_f(v: f64) -> String {
    if v.is_nan() {
        if v.is_sign_negative() { "-nan" } else { "nan" }.to_string()
    } else if v.is_infinite() {
        if v < 0.0 { "-inf" } else { "inf" }.to_string()
    } else {
        format!("{:.6}", v)
    }
}

impl Settings {
    pub fn new() -> Settings {
        Settings::default()
    }

    /// In C++ this registers every option on the embedded `CLI::App`. The Rust
    /// parser is built on demand inside `parse_input`, so nothing is needed here.
    pub fn prepare_settings(&mut self) {}

    /// Make a repeated option take its last value instead of raising an error,
    /// so tuning runs can append overriding arguments after the user's own.
    pub fn set_multi_option(&mut self) {
        self.multi_option = true;
    }

    /// Parse arguments (`argv[0]` is the program name) and validate them.
    /// Returns false (after writing a message to stderr) when validation fails.
    /// Like the C++, it exits the process for `--help`, `--cite`, unknown
    /// arguments (code 0) and CLI11 exceptions (code 255, printed the way
    /// `main` does).
    pub fn parse_input(&mut self, argv: &[String]) -> bool {
        let mut err = String::new();
        let r = self.try_parse_input(argv, &mut err);
        eprint!("{}", err);
        match r {
            Ok(passed) => passed,
            Err(CliExit::Help(text)) => {
                print!("{}", text);
                let _ = std::io::stdout().flush();
                std::process::exit(0)
            }
            Err(CliExit::Extras(msg)) => {
                eprintln!("Unrecognized flag or argument: {}", msg);
                std::process::exit(0)
            }
            Err(CliExit::Cite) => {
                print!("{}", CITE_TEXT);
                let _ = std::io::stdout().flush();
                std::process::exit(0)
            }
            Err(CliExit::Exception(msg)) => {
                let _ = std::io::stdout().flush();
                eprintln!("Standard exception caught: {}", msg);
                std::process::exit(255)
            }
        }
    }

    /// `parse_input` without the exits: everything the C++ writes to stderr
    /// before returning is appended to `err`, and the exit paths come back as
    /// `Err(CliExit)`.
    pub fn try_parse_input(&mut self, argv: &[String], err: &mut String) -> Result<bool, CliExit> {
        for (i, a) in argv.iter().enumerate() {
            if i > 0 {
                self.args.push(' ');
            }
            self.args.push_str(a);

            // Argument preprocessor
            let b = a.as_bytes();
            if b.len() >= 3 && b[0] == b'-' && b[1].is_ascii_alphabetic() {
                err.push_str(&format!(
                    "Argument '{}' is not allowed (long arguments begin with --, filenames may not begin with -)\n",
                    a
                ));
                return Ok(false);
            }
        }

        parse_args(self, argv)?;

        if self.cite {
            return Err(CliExit::Cite);
        }

        let mut passed = true;
        let mut fail = |msg: &str| {
            err.push_str(msg);
            passed = false;
        };

        if self.in_file.is_empty() && !self.show_memory && !self.run_without_reader {
            fail("Input file required.\n");
        }

        if !self.json_out {
            if self.show_deltas {
                fail("--show_delta is only available with --json\n");
            }
            if self.show_trace {
                fail("--show_trace is only available with --json\n");
            }
            if self.show_logo_nums {
                fail("--show_logo is only available with --json\n");
            }
        }

        if (self.ultra_out || self.json_out || self.bed_out) && self.suppress_out {
            fail("--suppress is incompatible with --tsv, --json, and --bed\n");
        }

        if self.mask_file.is_empty() && self.mask_with_n {
            fail("--maskn requires an output file path provided to --mask\n");
        }

        if self.threads < 1 {
            fail("--threads must be at least 1.\n");
        }

        // long long vs unsigned long long: the signed side converts to unsigned.
        if self.window_size != -1 && (self.window_size as u64) < self.max_period {
            fail("--winsize cannot be smaller than --period\n");
        }

        if self.windows < self.threads as i64 {
            fail("--windows must be at least as large as --threads\n");
        }

        if self.max_period < 1 {
            fail("--period must be at least 1\n");
        }

        if (self.at as f64) < 0.0 || (self.at as f64) > 1.0 {
            fail("--at must be >= 0.0 and <= 1.0\n");
        }

        if self.acgt.len() == 4 {
            let mut sum: f32 = 0.0;
            for i in 0..4 {
                sum += self.acgt[i];
                if self.acgt[i] < 0.0 || (self.acgt[i] as f64) > 1.0 {
                    fail(&format!("--acgt argument [{}] must be >= 0.0 and <= 1.0\n", i));
                }
            }

            // Upstream 1.2.2 calls the integer abs() here, which truncates the
            // difference first and lets any sum strictly between 0 and 2 pass.
            if (1.0_f64 - sum as f64).abs() > 0.001 {
                fail("--acgt arguments must sum to 1.0 \u{b1} 0.001\n");
            }
        }

        let frac_ok = |v: f32, lo_inclusive: bool| -> bool {
            let v = v as f64;
            if lo_inclusive {
                !(v < 0.0 || v > 1.0)
            } else {
                !(v <= 0.0 || v > 1.0)
            }
        };
        if !frac_ok(self.match_probability, false) {
            fail("--match must be > 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.period_decay, false) {
            fail("--decay must be > 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_nr, true) {
            fail("--nr must be >= 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_rn, true) {
            fail("--rn must be >= 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_ri, true) {
            fail("--ri must be >= 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_rd, true) {
            fail("--rd must be >= 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_ii, true) {
            fail("--ii must be >= 0.0 and <= 1.0\n");
        }
        if !frac_ok(self.transition_dd, true) {
            fail("--dd must be >= 0.0 and <= 1.0\n");
        }

        if self.split_threshold <= 0.0 {
            fail("--splitval must be > 0.0\n");
        }

        if self.tune_fdr < 0.0 || self.tune_fdr > 1.0 {
            fail("--tune_fdr must be >= 0 and <= 1.0\n");
        }

        if self.tune_only || self.tune_medium || self.tune_large || self.tune_indels {
            self.tune = true;
        }

        if self.tune_medium && self.tune_large {
            fail("Cannot use both --tune_medium and --tune_large\n");
        }

        if !self.tune_param_path.is_empty() && (self.tune_medium || self.tune_large) {
            fail("Cannot use both --tune_file and (--tune_small or --tune_large).\n");
        }

        let mut c = 0;
        if self.ultra_out {
            c += 1;
        }
        if self.json_out {
            c += 1;
        }
        if self.bed_out {
            c += 1;
        }

        if c > 1 && self.out_file.is_empty() {
            fail("Output file path must be provided when using multiple output formats .\n");
        }

        if !(self.ultra_out || self.json_out || self.bed_out) {
            self.ultra_out = true;
        }

        Ok(passed)
    }

    /// Parse `argv` followed by the whitespace-split tokens of `arg_str`
    /// (used by tuning; pair with `set_multi_option` so later options win).
    pub fn parse_multi_input(&mut self, argv: &[String], arg_str: &str) -> bool {
        let extra = string_to_args(arg_str);
        let combined = combine_args(argv, &extra);
        self.parse_input(&combined)
    }

    /// Number of DP states per column. Also sets `windows = -1` when `read_all` is set.
    pub fn calculate_num_states(&mut self) -> i32 {
        let mut num_of_states: i32 = 0;

        // `int i` compared with `unsigned long long` limits: i converts to unsigned.
        let mut i: i32 = 0;
        while (i as u64) < self.max_period {
            num_of_states += 1;

            if (i as u64) <= self.max_insert {
                if i > 1 {
                    num_of_states += i - 1;
                }
            } else {
                // int += unsigned long long: computed in 64 bits, truncated back to int.
                num_of_states = num_of_states.wrapping_add(self.max_insert as i32);
            }

            if (i as u64) <= self.max_delete {
                if i > 1 {
                    num_of_states += i - 1;
                }
            } else {
                num_of_states = num_of_states.wrapping_add(self.max_delete as i32);
            }
            i += 1;
        }

        if self.read_all {
            self.windows = -1;
        }

        num_of_states
    }

    /// Derive window size, overlap and base frequencies from the parsed settings.
    pub fn assign_settings(&mut self) {
        if !self.mask_file.is_empty() {
            self.produce_mask = true;
        }

        if self.overlap == -1 {
            self.overlap = self.max_period as i64;
        }

        if self.window_size == -1 {
            let num_states = self.calculate_num_states();
            // int * unsigned long long: the int sign-extends to 64 bits.
            let period_memory: u64 = (num_states as i64 as u64).wrapping_mul(self.max_period);
            let factor: u64 = if period_memory <= 1000 {
                10000 // itty bitty models use ~40 mb per thread
            } else if period_memory <= 10000 {
                2000 // tiny models use 80 mb per thread
            } else if period_memory <= 200000 {
                1000 // small models use less than 160 mb per thread
            } else if period_memory <= 2000000 {
                500 // medium models use less than 1 GB per thread
            } else if period_memory <= 20000000 {
                100
            } else if period_memory <= 50000000 {
                50
            } else {
                25 // large models use less than 4 GB per thread
            };
            self.window_size = factor.wrapping_mul(self.max_period) as i64;
        }

        self.a_freq = (self.at as f64 / 2.0) as f32;
        self.t_freq = (self.at as f64 / 2.0) as f32;
        self.c_freq = ((1.0 - self.at as f64) / 2.0) as f32;
        self.g_freq = ((1.0 - self.at as f64) / 2.0) as f32;

        if self.acgt.len() == 4 {
            let mut sum: f32 = 0.0;
            for i in 0..4 {
                sum += self.acgt[i];
            }
            self.a_freq = self.acgt[0] / sum;
            self.c_freq = self.acgt[1] / sum;
            self.g_freq = self.acgt[2] / sum;
            self.t_freq = self.acgt[3] / sum;
        }

        if self.no_split {
            self.max_split = 0;
        }
    }

    /// The text `print_memory_usage` writes to stdout.
    pub fn memory_usage_string(&mut self) -> String {
        let num_states = self.calculate_num_states();
        let mut dp_size = num_states as u64;
        // All the C++ arithmetic here is on 64-bit integers that may wrap.
        let total_window: i64 = self.window_size.wrapping_add(self.overlap.wrapping_mul(2));
        dp_size = dp_size.wrapping_mul(total_window as u64);
        let mut s = String::new();
        s.push_str("----------------------------\n");
        s.push_str(&format!("Maximum repeat period: {}\n", self.max_period));
        s.push_str(&format!("Number of model states: {}\n", num_states));
        s.push_str(&format!("Total sequence window size: {}\n", total_window));
        s.push_str(&format!("DP matrix cells: {}\n", dp_size));
        s.push_str(&format!("Threads: {}\n", self.threads));
        dp_size = dp_size.wrapping_mul(self.threads as u64);
        dp_size = dp_size.wrapping_mul(4);
        s.push_str(&format!("Sequence window queue length: {}\n", self.windows as u64));
        dp_size = dp_size.wrapping_add(total_window.wrapping_mul(self.windows) as u64);
        let gb_size: f32 = ((dp_size as f32) as f64 / (1024.0 * 1024.0 * 1024.0)) as f32;
        if gb_size as f64 > 0.09 {
            s.push_str(&format!("Total size: {:.2} GB ({} bytes)\n", gb_size as f64, dp_size));
        } else {
            s.push_str(&format!("Total size: ({} bytes)\n", dp_size));
        }
        s.push_str("----------------------------\n");
        s.push_str(
            "*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices\n",
        );
        s.push_str("============================\n");
        s
    }

    pub fn print_memory_usage(&mut self) {
        print!("{}", self.memory_usage_string());
    }

    /// Settings as the JSON fragment written into output files.
    pub fn json_string(&self) -> String {
        let mut j = String::new();
        fn s(j: &mut String, name: &str, v: &str) {
            j.push_str(&format!("\"{}\": \"{}\",\n", name, v));
        }
        fn b(j: &mut String, name: &str, v: bool) {
            j.push_str(&format!("\"{}\": {},\n", name, if v { "true" } else { "false" }));
        }
        fn n<T: std::fmt::Display>(j: &mut String, name: &str, v: T) {
            j.push_str(&format!("\"{}\": {},\n", name, v));
        }
        fn f(j: &mut String, name: &str, v: f32) {
            j.push_str(&format!("\"{}\": {},\n", name, to_string_f(v)));
        }
        s(&mut j, "args", &self.args);

        s(&mut j, "in_file", &self.in_file);
        b(&mut j, "read_all", self.read_all);

        s(&mut j, "out_file", &self.out_file);
        b(&mut j, "pval", self.pval);
        f(&mut j, "p_value_loc", self.p_value_loc);
        f(&mut j, "p_value_scale", self.p_value_scale);
        f(&mut j, "p_value_freq", self.p_value_freq);

        b(&mut j, "ultra_out", self.ultra_out);
        b(&mut j, "json_out", self.json_out);
        b(&mut j, "bed_out", self.bed_out);
        n(&mut j, "max_consensus_period", self.max_consensus_period);
        b(&mut j, "show_seq", self.show_seq);
        b(&mut j, "show_deltas", self.show_deltas);
        b(&mut j, "show_trace", self.show_trace);
        b(&mut j, "show_wid", self.show_wid);
        b(&mut j, "show_logo_nums", self.show_logo_nums);
        b(&mut j, "suppress_out", self.suppress_out);

        b(&mut j, "produce_mask", self.produce_mask);
        s(&mut j, "mask_file", &self.mask_file);
        b(&mut j, "mask_with_n", self.mask_with_n);

        n(&mut j, "threads", self.threads);
        n(&mut j, "window_size", self.window_size);
        n(&mut j, "overlap", self.overlap);
        n(&mut j, "windows", self.windows);

        f(&mut j, "min_score", self.min_score);
        n(&mut j, "min_length", self.min_length);
        n(&mut j, "min_units", self.min_units);

        n(&mut j, "max_period", self.max_period);
        n(&mut j, "max_insert", self.max_insert);
        n(&mut j, "max_delete", self.max_delete);

        f(&mut j, "a_freq", self.a_freq);
        f(&mut j, "c_freq", self.c_freq);
        f(&mut j, "g_freq", self.g_freq);
        f(&mut j, "t_freq", self.t_freq);
        f(&mut j, "match_probability", self.match_probability);

        f(&mut j, "period_decay", self.period_decay);
        f(&mut j, "transition_nr", self.transition_nr);
        f(&mut j, "transition_rn", self.transition_rn);
        f(&mut j, "transition_ri", self.transition_ri);
        f(&mut j, "transition_rd", self.transition_rd);
        f(&mut j, "transition_ii", self.transition_ii);
        f(&mut j, "transition_dd", self.transition_dd);

        b(&mut j, "no_split", self.no_split);
        n(&mut j, "max_split", self.max_split);
        f(&mut j, "split_threshold", self.split_threshold);
        n(&mut j, "split_depth", self.split_depth);
        n(&mut j, "min_split_window", self.min_split_window);
        j.pop();
        j.pop();
        j.push('\n');
        j
    }
}

fn tune_line(m: f32, at: f32, start: f32, stop: f32) -> String {
    let mut p = format!("-m {} ", to_string_f(m));
    p += &format!("--at {} ", to_string_f(at));
    p += &format!("--rn {} ", to_string_f(start));
    p += &format!("--nr {}", to_string_f(stop));
    p
}

fn f32s(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

/// Shared by the small and medium grids: start/stop values are paired by index.
fn paired_grid(matches: &[f64], ats: &[f64]) -> Vec<String> {
    let (matches, ats) = (f32s(matches), f32s(ats));
    let repeat_start = f32s(&[0.001, 0.01]);
    let repeat_stop = f32s(&[0.005, 0.05]);
    let mut settings = Vec::new();
    for &m in &matches {
        for &at in &ats {
            for i in 0..2 {
                settings.push(tune_line(m, at, repeat_start[i], repeat_stop[i]));
            }
        }
    }
    settings
}

/// Small tuning grid (18 parameter strings).
pub fn small_tune_settings() -> Vec<String> {
    paired_grid(&[0.6, 0.75, 0.9], &[0.4, 0.5, 0.6])
}

/// Medium tuning grid (40 parameter strings).
pub fn medium_tune_settings() -> Vec<String> {
    paired_grid(&[0.6, 0.7, 0.8, 0.9], &[0.3, 0.4, 0.5, 0.6, 0.7])
}

/// Large tuning grid (252 parameter strings).
pub fn large_tune_settings() -> Vec<String> {
    let matches = f32s(&[0.6, 0.7, 0.8, 0.9]);
    let ats = f32s(&[0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7]);
    let repeat_start = f32s(&[0.001, 0.005, 0.01]);
    let repeat_stop = f32s(&[0.005, 0.01, 0.05]);
    let mut settings = Vec::new();
    for &m in &matches {
        for &at in &ats {
            for &rs in &repeat_start {
                for &rt in &repeat_stop {
                    settings.push(tune_line(m, at, rs, rt));
                }
            }
        }
    }
    settings
}

/// Read parameter strings, one per non-empty line, from a tune file. Each line
/// is validated by parsing it as a full command line (so, as in the C++, its
/// first token is consumed as the program name). Exits with code 0 on a bad
/// file or line.
pub fn tune_settings_for_path(path: &str) -> Vec<String> {
    let mut settings = Vec::new();
    let mut bytes = Vec::new();
    match std::fs::File::open(path) {
        Ok(mut f) => {
            // A directory opens but reads nothing, like ifstream.
            let _ = std::io::Read::read_to_end(&mut f, &mut bytes);
        }
        Err(_) => {
            eprintln!("Failed to open settings file: {}", path);
            std::process::exit(0);
        }
    }

    let text = String::from_utf8_lossy(&bytes).into_owned();
    let mut lines: Vec<&str> = text.split('\n').collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    for (i, line) in lines.iter().enumerate() {
        let line_num = i + 1; // we 1 index line numbers
        if !line.is_empty() {
            let argv = string_to_args(line);
            let mut test_settings = Settings::new();
            test_settings.prepare_settings();
            if !test_settings.parse_input(&argv) {
                eprintln!("Invalid arguments on line {} in tune file. \"{}\"", line_num, line);
                std::process::exit(0);
            }
            settings.push(line.to_string());
        }
    }
    settings
}

/// Split on whitespace (C `isspace`: space, \t, \n, \v, \f, \r), like `istringstream >> string`.
pub fn string_to_args(s: &str) -> Vec<String> {
    s.split(|c: char| c.is_ascii() && c_isspace(c as u8))
        .filter(|t| !t.is_empty())
        .map(String::from)
        .collect()
}

/// Concatenate two argument lists.
pub fn combine_args(a: &[String], b: &[String]) -> Vec<String> {
    let mut v = Vec::with_capacity(a.len() + b.len());
    v.extend_from_slice(a);
    v.extend_from_slice(b);
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    /// Parse like the C++ binary up to `assign_settings`; report the outcome class
    /// the way the reference run was classified.
    fn classify_parse(a: &[&str]) -> (String, String) {
        let mut s = Settings::new();
        let mut err = String::new();
        match s.try_parse_input(&argv(a), &mut err) {
            Ok(true) => ("Ok".into(), String::new()),
            Ok(false) => ("Fail".into(), err),
            Err(CliExit::Help(_)) => ("Help".into(), String::new()),
            Err(CliExit::Cite) => ("Cite".into(), String::new()),
            Err(CliExit::Extras(m)) => ("Extras".into(), m),
            Err(CliExit::Exception(m)) => ("Exc".into(), m),
        }
    }

    #[test]
    fn parse_cases_match_reference() {
        for (a, kind, text) in PARSE_CASES {
            let (k, t) = classify_parse(a);
            assert_eq!(&k, kind, "outcome for {:?}", a);
            let expect = if *kind == "Fail" {
                text.to_string()
            } else {
                text.trim_end_matches('\n').to_string()
            };
            let got = if *kind == "Fail" {
                t
            } else {
                t.trim_end_matches('\n').to_string()
            };
            if *kind == "Fail" || *kind == "Exc" || *kind == "Extras" {
                assert_eq!(got, expect, "message for {:?}", a);
            }
        }
    }

    #[test]
    fn json_cases_match_reference() {
        for (a, expected) in JSON_CASES {
            let mut s = Settings::new();
            let mut err = String::new();
            assert_eq!(s.try_parse_input(&argv(a), &mut err), Ok(true), "{:?}: {}", a, err);
            s.assign_settings();
            let full = format!(
                "{{\"Version\": \"{}\", \n\"Parameters\": {{\n{}}}}}\n",
                ULTRA_VERSION_STRING,
                s.json_string()
            );
            let expected = expected.replace("@VERSION@", ULTRA_VERSION_STRING);
            assert_eq!(full, expected, "json for {:?}", a);
        }
    }

    #[test]
    fn memory_cases_match_reference() {
        for (a, expected) in MEM_CASES {
            let mut s = Settings::new();
            let mut err = String::new();
            assert_eq!(s.try_parse_input(&argv(a), &mut err), Ok(true), "{:?}: {}", a, err);
            s.assign_settings();
            assert_eq!(&s.memory_usage_string(), expected, "memory for {:?}", a);
        }
    }

    #[test]
    fn num_states_known_values() {
        let mut s = Settings::new();
        assert_eq!(s.calculate_num_states(), 1970);
        s.max_period = 5;
        assert_eq!(s.calculate_num_states(), 17);
        s.read_all = true;
        assert_eq!(s.windows, 1024);
        s.calculate_num_states();
        assert_eq!(s.windows, -1);
    }

    #[test]
    fn string_to_args_and_combine() {
        assert_eq!(
            string_to_args("  -m 0.6\t--at\n0.4\x0b--rn  1 \r\x0c"),
            argv(&["-m", "0.6", "--at", "0.4", "--rn", "1"])
        );
        assert!(string_to_args("").is_empty());
        assert!(string_to_args(" \t ").is_empty());
        assert_eq!(combine_args(&argv(&["a", "b"]), &argv(&["c"])), argv(&["a", "b", "c"]));
        assert_eq!(combine_args(&[], &[]), Vec::<String>::new());
    }

    #[test]
    fn multi_option_last_wins() {
        let a = argv(&["./ultra", "ex.fa", "--at", "0.3", "-p", "5"]);
        let mut s = Settings::new();
        s.set_multi_option();
        let mut err = String::new();
        let combined = combine_args(&a, &string_to_args("-m 0.6 --at 0.4 --rn 0.001 --nr 0.005 --period 7"));
        assert_eq!(s.try_parse_input(&combined, &mut err), Ok(true), "{}", err);
        assert_eq!(s.at, 0.4_f32);
        assert_eq!(s.max_period, 7);
        assert_eq!(s.match_probability, 0.6_f32);
        assert_eq!(s.transition_rn, 0.001_f32);
        assert_eq!(s.transition_nr, 0.005_f32);
        assert_eq!(s.in_file, "ex.fa");

        // acgt keeps its last four values; a second positional is still an error.
        let mut s = Settings::new();
        s.set_multi_option();
        let c = argv(&[
            "u", "f", "--acgt", "0.1", "0.2", "0.3", "0.4", "--acgt", "0.4", "0.3", "0.2", "0.1",
        ]);
        assert_eq!(s.try_parse_input(&c, &mut err), Ok(true));
        assert_eq!(s.acgt, vec![0.4_f32, 0.3, 0.2, 0.1]);
        let mut s = Settings::new();
        s.set_multi_option();
        assert_eq!(
            s.try_parse_input(&argv(&["u", "f", "g"]), &mut err),
            Err(CliExit::Extras("The following argument was not expected: g".to_string()))
        );
    }

    #[test]
    fn without_multi_option_duplicate_is_error() {
        let mut s = Settings::new();
        let mut err = String::new();
        assert_eq!(
            s.try_parse_input(&argv(&["u", "f", "--at", "0.3", "--at", "0.4"]), &mut err),
            Err(CliExit::Exception(
                "--at: At Most 1 required but received 2".to_string()
            ))
        );
    }

    #[test]
    fn tune_grids() {
        let s = small_tune_settings();
        assert_eq!(s.len(), 18);
        assert_eq!(s[0], "-m 0.600000 --at 0.400000 --rn 0.001000 --nr 0.005000");
        assert_eq!(s[1], "-m 0.600000 --at 0.400000 --rn 0.010000 --nr 0.050000");
        assert_eq!(s[17], "-m 0.900000 --at 0.600000 --rn 0.010000 --nr 0.050000");
        assert_eq!(s[6], "-m 0.750000 --at 0.400000 --rn 0.001000 --nr 0.005000");
        assert_eq!(medium_tune_settings().len(), 40);
        let l = large_tune_settings();
        assert_eq!(l.len(), 252);
        assert_eq!(l[0], "-m 0.600000 --at 0.300000 --rn 0.001000 --nr 0.005000");
        assert_eq!(l[1], "-m 0.600000 --at 0.300000 --rn 0.001000 --nr 0.010000");
        assert_eq!(l[9], "-m 0.600000 --at 0.350000 --rn 0.001000 --nr 0.005000");
        assert_eq!(l[251], "-m 0.900000 --at 0.700000 --rn 0.010000 --nr 0.050000");
    }

    #[test]
    fn help_text_shape() {
        let h = help_text("./ultra");
        assert!(
            h.starts_with(
                "
=====
"
            ) || h.starts_with(
                "
====="
            )
        );
        assert!(h.contains("Usage: ./ultra [OPTIONS] [input_file]\n"));
        assert!(h.contains("  --acgt FLOAT x 4            Expected frequency of A C G T\n"));
        assert!(h.ends_with("NIH NHGRI U24HG010136\n\n"));
    }

    #[test]
    fn numeric_conversions() {
        assert_eq!(lex_u64("0x10"), Some(16));
        assert_eq!(lex_u64("010"), Some(8));
        assert_eq!(lex_u64("08"), None);
        assert_eq!(lex_u64("-1"), None);
        assert_eq!(lex_u64("99999999999999999999"), None);
        assert_eq!(lex_i32("true"), Some(1));
        assert_eq!(lex_i32("2147483648"), None);
        assert_eq!(lex_i32(" 3"), Some(3));
        assert_eq!(lex_i64("-9223372036854775808"), Some(i64::MIN));
        assert_eq!(lex_f32("0x1"), Some(1.0));
        assert_eq!(lex_f32(" .5"), Some(0.5));
        assert_eq!(lex_f32("1e"), None);
        assert_eq!(lex_f32("abc"), None);
        assert!(lex_f32("nan").unwrap().is_nan());
        assert_eq!(lex_bool("T"), Some(true));
        assert_eq!(lex_bool("0"), Some(false));
        assert_eq!(lex_bool("maybe"), None);
        assert_eq!(lex_bool("12abc"), Some(true));
        assert_eq!(lex_bool("-1"), Some(false));
        assert_eq!(fmt_f(-0.0), "-0.000000");
    }

    const JSON_CASES: &[(&[&str], &str)] = &[
        (
            &["./ultra", "ex.fa"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-p", "20", "-i", "3", "-d", "4"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -p 20 -i 3 -d 4",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 40000,
"overlap": 20,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 20,
"max_insert": 3,
"max_delete": 4,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--at", "0.3"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --at 0.3",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.150000,
"c_freq": 0.350000,
"g_freq": 0.350000,
"t_freq": 0.150000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--acgt", "0.1", "0.2", "0.3", "0.4"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --acgt 0.1 0.2 0.3 0.4",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.100000,
"c_freq": 0.200000,
"g_freq": 0.300000,
"t_freq": 0.400000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--acgt", "0.3", "0.2", "0.2", "0.3", "--at", "0.9"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --acgt 0.3 0.2 0.2 0.3 --at 0.9",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.300000,
"c_freq": 0.200000,
"g_freq": 0.200000,
"t_freq": 0.300000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--win_size", "5000", "--overlap", "50"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --win_size 5000 --overlap 50",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 5000,
"overlap": 50,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-t", "4", "--windows", "100"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -t 4 --windows 100",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 4,
"window_size": 100000,
"overlap": 100,
"windows": 100,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--json", "-o", "out1"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --json -o out1",
"in_file": "ex.fa",
"read_all": false,
"out_file": "out1",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": false,
"json_out": true,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--bed", "--pval", "-o", "out2"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --bed --pval -o out2",
"in_file": "ex.fa",
"read_all": false,
"out_file": "out2",
"pval": true,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": false,
"json_out": false,
"bed_out": true,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--tsv", "--json", "--bed", "-o", "out3"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --tsv --json --bed -o out3",
"in_file": "ex.fa",
"read_all": false,
"out_file": "out3",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": true,
"bed_out": true,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--json",
                "--show_seq",
                "--show_delta",
                "--show_trace",
                "--show_wid",
                "--show_logo",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --json --show_seq --show_delta --show_trace --show_wid --show_logo",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": false,
"json_out": true,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": true,
"show_deltas": true,
"show_trace": true,
"show_wid": true,
"show_logo_nums": true,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--suppress"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --suppress",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": true,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mask", "m.fa", "--nmask"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --mask m.fa --nmask",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": true,
"mask_file": "m.fa",
"mask_with_n": true,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "-s",
                "-20.5",
                "--min_length",
                "15",
                "--min_unit",
                "3",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -s -20.5 --min_length 15 --min_unit 3",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -20.500000,
"min_length": 15,
"min_units": 3,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra", "ex.fa", "-m", "0.8", "--decay", "0.9", "--nr", "0.02", "--rn", "0.1", "--ri", "0.03",
                "--rd", "0.04", "--ii", "0.05", "--dd", "0.06",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -m 0.8 --decay 0.9 --nr 0.02 --rn 0.1 --ri 0.03 --rd 0.04 --ii 0.05 --dd 0.06",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.800000,
"period_decay": 0.900000,
"transition_nr": 0.020000,
"transition_rn": 0.100000,
"transition_ri": 0.030000,
"transition_rd": 0.040000,
"transition_ii": 0.050000,
"transition_dd": 0.060000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--no_split"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --no_split",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": true,
"max_split": 0,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--max_split",
                "5",
                "--split_threshold",
                "2.25",
                "--split_depth",
                "7",
                "--min_split_window",
                "20",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --max_split 5 --split_threshold 2.25 --split_depth 7 --min_split_window 20",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 5,
"split_threshold": 2.250000,
"split_depth": 7,
"min_split_window": 20
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-r"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -r",
"in_file": "ex.fa",
"read_all": true,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": -1,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--max_consensus", "50"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --max_consensus 50",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 50,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--pval_loc=2.5", "--pval_scale=1.5"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --pval_loc=2.5 --pval_scale=1.5",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 2.500000,
"p_value_scale": 1.500000,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-p", "5"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -p 5",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 50000,
"overlap": 5,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 5,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-p", "12", "-i", "2", "-d", "2"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -p 12 -i 2 -d 2",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 120000,
"overlap": 12,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 12,
"max_insert": 2,
"max_delete": 2,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--minlen",
                "7",
                "--minunit",
                "4",
                "--winsize",
                "3000",
                "--maxsplit",
                "3",
                "--splitval",
                "1.5",
                "--splitdepth",
                "9",
                "--minsplitwin",
                "30",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --minlen 7 --minunit 4 --winsize 3000 --maxsplit 3 --splitval 1.5 --splitdepth 9 --minsplitwin 30",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 3000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 7,
"min_units": 4,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 3,
"split_threshold": 1.500000,
"split_depth": 9,
"min_split_window": 30
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "--at", "0.333333"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa --at 0.333333",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 100,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.166666,
"c_freq": 0.333333,
"g_freq": 0.333333,
"t_freq": 0.166666,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "-s",
                "1e2",
                "-t",
                "0x2",
                "--windows",
                "010",
                "--overlap",
                "33",
            ],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -s 1e2 -t 0x2 --windows 010 --overlap 33",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 2,
"window_size": 100000,
"overlap": 33,
"windows": 8,
"min_score": 100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 100,
"max_insert": 10,
"max_delete": 10,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": false,
"max_split": 10,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
        (
            &["./ultra", "ex.fa", "-p", "50", "-i", "0", "-d", "0", "--nosplit"],
            r####"{"Version": "@VERSION@", 
"Parameters": {
"args": "./ultra ex.fa -p 50 -i 0 -d 0 --nosplit",
"in_file": "ex.fa",
"read_all": false,
"out_file": "",
"pval": false,
"p_value_loc": 4.272949,
"p_value_scale": 1.891360,
"p_value_freq": 0.000330,
"ultra_out": true,
"json_out": false,
"bed_out": false,
"max_consensus_period": 1000000,
"show_seq": false,
"show_deltas": false,
"show_trace": false,
"show_wid": false,
"show_logo_nums": false,
"suppress_out": false,
"produce_mask": false,
"mask_file": "",
"mask_with_n": false,
"threads": 1,
"window_size": 100000,
"overlap": 50,
"windows": 1024,
"min_score": -100.000000,
"min_length": 10,
"min_units": 2,
"max_period": 50,
"max_insert": 0,
"max_delete": 0,
"a_freq": 0.250000,
"c_freq": 0.250000,
"g_freq": 0.250000,
"t_freq": 0.250000,
"match_probability": 0.700000,
"period_decay": 0.850000,
"transition_nr": 0.010000,
"transition_rn": 0.050000,
"transition_ri": 0.020000,
"transition_rd": 0.020000,
"transition_ii": 0.020000,
"transition_dd": 0.020000,
"no_split": true,
"max_split": 0,
"split_threshold": 3.500000,
"split_depth": 5,
"min_split_window": 15
}}
"####,
        ),
    ];

    const MEM_CASES: &[(&[&str], &str)] = &[
        (
            &["./ultra", "ex.fa", "--mem"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 1970
Total sequence window size: 100200
DP matrix cells: 197394000
Threads: 1
Sequence window queue length: 1024
Total size: 0.83 GB (892180800 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5"],
            r####"----------------------------
Maximum repeat period: 5
Number of model states: 17
Total sequence window size: 50010
DP matrix cells: 850170
Threads: 1
Sequence window queue length: 1024
Total size: (54610920 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "10", "-i", "2", "-d", "3"],
            r####"----------------------------
Maximum repeat period: 10
Number of model states: 46
Total sequence window size: 100020
DP matrix cells: 4600920
Threads: 1
Sequence window queue length: 1024
Total size: 0.11 GB (120824160 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "30"],
            r####"----------------------------
Maximum repeat period: 30
Number of model states: 500
Total sequence window size: 30060
DP matrix cells: 15030000
Threads: 1
Sequence window queue length: 1024
Total size: (90901440 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "150"],
            r####"----------------------------
Maximum repeat period: 150
Number of model states: 3020
Total sequence window size: 75300
DP matrix cells: 227406000
Threads: 1
Sequence window queue length: 1024
Total size: 0.92 GB (986731200 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "300", "-i", "0", "-d", "0"],
            r####"----------------------------
Maximum repeat period: 300
Number of model states: 300
Total sequence window size: 300600
DP matrix cells: 90180000
Threads: 1
Sequence window queue length: 1024
Total size: 0.62 GB (668534400 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "3000", "-i", "0", "-d", "0"],
            r####"----------------------------
Maximum repeat period: 3000
Number of model states: 3000
Total sequence window size: 306000
DP matrix cells: 918000000
Threads: 1
Sequence window queue length: 1024
Total size: 3.71 GB (3985344000 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5000", "-i", "0", "-d", "0"],
            r####"----------------------------
Maximum repeat period: 5000
Number of model states: 5000
Total sequence window size: 260000
DP matrix cells: 1300000000
Threads: 1
Sequence window queue length: 1024
Total size: 5.09 GB (5466240000 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "10000", "-i", "0", "-d", "0"],
            r####"----------------------------
Maximum repeat period: 10000
Number of model states: 10000
Total sequence window size: 270000
DP matrix cells: 2700000000
Threads: 1
Sequence window queue length: 1024
Total size: 10.32 GB (11076480000 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-t", "4"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 1970
Total sequence window size: 100200
DP matrix cells: 197394000
Threads: 4
Sequence window queue length: 1024
Total size: 3.04 GB (3260908800 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--windows", "10"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 1970
Total sequence window size: 100200
DP matrix cells: 197394000
Threads: 1
Sequence window queue length: 10
Total size: 0.74 GB (790578000 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-r"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 1970
Total sequence window size: 100200
DP matrix cells: 197394000
Threads: 1
Sequence window queue length: 18446744073709551615
Total size: 0.74 GB (789475800 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--win_size", "777", "--overlap", "7"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 1970
Total sequence window size: 791
DP matrix cells: 1558270
Threads: 1
Sequence window queue length: 1024
Total size: (7043064 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "2"],
            r####"----------------------------
Maximum repeat period: 2
Number of model states: 2
Total sequence window size: 20004
DP matrix cells: 40008
Threads: 1
Sequence window queue length: 1024
Total size: (20644128 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "1", "-i", "0", "-d", "0"],
            r####"----------------------------
Maximum repeat period: 1
Number of model states: 1
Total sequence window size: 10002
DP matrix cells: 10002
Threads: 1
Sequence window queue length: 1024
Total size: (10282056 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--mem",
                "-t",
                "3",
                "--windows",
                "5",
                "-p",
                "7",
                "-i",
                "1",
                "-d",
                "20",
            ],
            r####"----------------------------
Maximum repeat period: 7
Number of model states: 27
Total sequence window size: 70014
DP matrix cells: 1890378
Threads: 3
Sequence window queue length: 5
Total size: (23034606 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "100", "-i", "100", "-d", "100"],
            r####"----------------------------
Maximum repeat period: 100
Number of model states: 9802
Total sequence window size: 50200
DP matrix cells: 492060400
Threads: 1
Sequence window queue length: 1024
Total size: 1.88 GB (2019646400 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--mem",
                "-p",
                "40",
                "-i",
                "50",
                "-d",
                "3",
                "--win_size",
                "99999",
                "-t",
                "8",
            ],
            r####"----------------------------
Maximum repeat period: 40
Number of model states: 892
Total sequence window size: 100079
DP matrix cells: 89270468
Threads: 8
Sequence window queue length: 1024
Total size: 2.76 GB (2959135872 bytes)
----------------------------
*Actual memory usage will be slightly greater due to output repeat queue and repeat split matrices
============================
"####,
        ),
    ];

    const PARSE_CASES: &[(&[&str], &str, &str)] = &[
        (&["./ultra", "--help"], "Help", r####""####),
        (&["./ultra", "-h"], "Help", r####""####),
        (&["./ultra", "ex.fa", "--mem", "-h"], "Help", r####""####),
        (
            &["./ultra", "--version"],
            "Extras",
            r####"The following argument was not expected: --version"####,
        ),
        (
            &["./ultra", "--bogus"],
            "Extras",
            r####"The following argument was not expected: --bogus"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--bogus"],
            "Extras",
            r####"The following argument was not expected: --bogus"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-x"],
            "Extras",
            r####"The following argument was not expected: -x"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-pt", "5"],
            "Fail",
            r####"Argument '-pt' is not allowed (long arguments begin with --, filenames may not begin with -)
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p5"],
            "Fail",
            r####"Argument '-p5' is not allowed (long arguments begin with --, filenames may not begin with -)
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-o=x"],
            "Fail",
            r####"Argument '-o=x' is not allowed (long arguments begin with --, filenames may not begin with -)
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-rc"],
            "Fail",
            r####"Argument '-rc' is not allowed (long arguments begin with --, filenames may not begin with -)
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--Mem"],
            "Extras",
            r####"The following argument was not expected: --Mem"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--me"],
            "Extras",
            r####"The following argument was not expected: --me"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "extra"],
            "Extras",
            r####"The following argument was not expected: extra"####,
        ),
        (
            &["./ultra", "ex.fa", "extra", "--mem"],
            "Extras",
            r####"The following argument was not expected: extra"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-"],
            "Extras",
            r####"The following argument was not expected: -"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--", "-x"],
            "Extras",
            r####"The following arguments were not expected: -x --"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--", "extra"],
            "Extras",
            r####"The following arguments were not expected: extra --"####,
        ),
        (&["./ultra", "--mem", "--", "ex.fa"], "Ok", r####""####),
        (&["./ultra", "--mem", "-5"], "Ok", r####""####),
        (
            &["./ultra", "--mem", "-.5"],
            "Extras",
            r####"The following argument was not expected: -.5"####,
        ),
        (
            &["./ultra", "--mem", "-_a"],
            "Extras",
            r####"The following argument was not expected: -_a"####,
        ),
        (
            &["./ultra", "--mem", "--5"],
            "Extras",
            r####"The following argument was not expected: --5"####,
        ),
        (
            &["./ultra", "--mem", "--=x"],
            "Extras",
            r####"The following argument was not expected: --=x"####,
        ),
        (&["./ultra", "--mem"], "Ok", r####""####),
        (&["./ultra", "ex.fa"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=false"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem=maybe"],
            "Exc",
            r####"Could not convert: --mem = maybe"####,
        ),
        (&["./ultra", "ex.fa", "--mem=yes"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem="], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=0"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=12abc"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=-1"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=T"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=n"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--mem=false"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem=false", "--mem"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "abc"],
            "Exc",
            r####"Could not convert: --period = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "-5"],
            "Exc",
            r####"Could not convert: --period = -5"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5.0"],
            "Exc",
            r####"Could not convert: --period = 5.0"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "1e2"],
            "Exc",
            r####"Could not convert: --period = 1e2"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "-p", "0x10"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "-p", "010"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "08"],
            "Exc",
            r####"Could not convert: --period = 08"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", ""],
            "Fail",
            r####"--period must be at least 1
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "0"],
            "Fail",
            r####"--period must be at least 1
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "true"],
            "Exc",
            r####"Could not convert: --period = true"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "99999999999999999999"],
            "Exc",
            r####"Could not convert: --period = 99999999999999999999"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "-t", "true"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-t", "99999999999"],
            "Exc",
            r####"Could not convert: --threads = 99999999999"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-t", "-1"],
            "Fail",
            r####"--threads must be at least 1.
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-t", "0"],
            "Fail",
            r####"--threads must be at least 1.
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--windows", "-1"],
            "Fail",
            r####"--windows must be at least as large as --threads
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--windows", "0", "-t", "1"],
            "Fail",
            r####"--windows must be at least as large as --threads
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-t", "2", "--windows", "1"],
            "Fail",
            r####"--windows must be at least as large as --threads
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5", "-p", "6"],
            "Exc",
            r####"--period: At Most 1 required but received 2"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5", "--period", "6"],
            "Exc",
            r####"--period: At Most 1 required but received 2"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--period=5", "--period=6"],
            "Exc",
            r####"--period: At Most 1 required but received 2"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--mem",
                "--winsize",
                "50",
                "--win_size",
                "60",
                "-p",
                "3",
            ],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--winsize", "abc"],
            "Exc",
            r####"Could not convert:  = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--minlen", "abc"],
            "Exc",
            r####"Could not convert:  = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--at", "0.3", "--at", "0.4"],
            "Exc",
            r####"--at: At Most 1 required but received 2"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-o"],
            "Exc",
            r####"--out: 1 required TEXT missing"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "-o", "--json"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-p"],
            "Exc",
            r####"--period: 1 required UINT missing"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--period="],
            "Exc",
            r####"--period: 1 required UINT missing"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--period=", "7"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--at", "1.5"],
            "Fail",
            r####"--at must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--at", "-0.1"],
            "Fail",
            r####"--at must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--at", "abc"],
            "Exc",
            r####"Could not convert: --at = abc"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--at", "nan"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--at", "inf"],
            "Fail",
            r####"--at must be >= 0.0 and <= 1.0
"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--at", "0x1"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--at", " 0.5"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--at", ""], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "-s", "inf"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "-s", "-50"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-s", "-50", "--score=-3"],
            "Exc",
            r####"--score: At Most 1 required but received 2"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-s=-3"],
            "Fail",
            r####"Argument '-s=-3' is not allowed (long arguments begin with --, filenames may not begin with -)
"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--score=-3"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--pval_loc", "3.0"],
            "Exc",
            r####"Could not convert: --pval_loc = true"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--pval_loc"],
            "Exc",
            r####"Could not convert: --pval_loc = true"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--pval_loc", "--pval_loc"],
            "Exc",
            r####"Could not convert: --pval_loc = true,true"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--pval_loc=1"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--pval_loc=x"],
            "Exc",
            r####"Could not convert: --pval_loc = x"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.1", "0.2", "0.3", "0.4"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.1", "0.2", "0.3"],
            "Exc",
            r####"--acgt: At least 4 required but received 3"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.1", "0.2", "0.3", "0.4", "0.5"],
            "Exc",
            r####"--acgt: At Most 4 required but received 5"####,
        ),
        (
            &[
                "./ultra", "ex.fa", "--mem", "--acgt", "0.1", "0.2", "0.3", "0.4", "extra",
            ],
            "Exc",
            r####"--acgt: At Most 4 required but received 5"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt=0.1", "0.2", "0.3", "0.4"],
            "Ok",
            r####""####,
        ),
        (
            &[
                "./ultra", "ex.fa", "--mem", "--acgt", "0.25", "0.25", "0.25", "0.25", "--acgt", "0.1", "0.2", "0.3",
                "0.4",
            ],
            "Exc",
            r####"--acgt: At Most 4 required but received 8"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.25", "-0.25", "0.5", "0.5"],
            "Fail",
            r####"--acgt argument [1] must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &[
                "./ultra", "ex.fa", "--mem", "--acgt", "0.25", "0.25", "0.25", "0.25", "--", "extra",
            ],
            "Extras",
            r####"The following argument was not expected: extra"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt"],
            "Exc",
            r####"--acgt: 1 required FLOAT missing"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.5"],
            "Exc",
            r####"--acgt: At least 4 required but received 1"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.25", "0.25", "0.25", "abc"],
            "Exc",
            r####"Could not convert: --acgt = 0.25,0.25,0.25,abc"####,
        ),
        (
            // Upstream 1.2.2 accepts this sum of 0.4 (integer abs()).
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.1", "0.1", "0.1", "0.1"],
            "Fail",
            r####"--acgt arguments must sum to 1.0 ± 0.001
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.9", "0.9", "0.9", "0.9"],
            "Fail",
            r####"--acgt arguments must sum to 1.0 ± 0.001
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "1", "1", "1", "1"],
            "Fail",
            r####"--acgt arguments must sum to 1.0 ± 0.001
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0", "0", "0", "0"],
            "Fail",
            r####"--acgt arguments must sum to 1.0 ± 0.001
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "2", "-1", "0", "0"],
            "Fail",
            r####"--acgt argument [0] must be >= 0.0 and <= 1.0
--acgt argument [1] must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--acgt", "0.5", "0.5", "0.5", "-0.5"],
            "Fail",
            r####"--acgt argument [3] must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "--mem", "--acgt", "0.25", "0.25", "0.25", "0.25", "ex.fa"],
            "Exc",
            r####"--acgt: At Most 4 required but received 5"####,
        ),
        (
            &[
                "./ultra", "--mem", "--acgt", "0.25", "0.25", "0.25", "0.25", "--", "ex.fa",
            ],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--json", "--suppress"],
            "Fail",
            r####"--suppress is incompatible with --tsv, --json, and --bed
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tsv", "--suppress"],
            "Fail",
            r####"--suppress is incompatible with --tsv, --json, and --bed
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--bed", "--suppress"],
            "Fail",
            r####"--suppress is incompatible with --tsv, --json, and --bed
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tsv", "--json"],
            "Fail",
            r####"Output file path must be provided when using multiple output formats .
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tsv", "--json", "-o", "x"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--show_delta"],
            "Fail",
            r####"--show_delta is only available with --json
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--show_trace"],
            "Fail",
            r####"--show_trace is only available with --json
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--show_logo"],
            "Fail",
            r####"--show_logo is only available with --json
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--showdelta", "--showtrace", "--showlogo"],
            "Fail",
            r####"--show_delta is only available with --json
--show_trace is only available with --json
--show_logo is only available with --json
"####,
        ),
        (
            &[
                "./ultra",
                "ex.fa",
                "--mem",
                "--json",
                "--showdelta",
                "--showtrace",
                "--showlogo",
                "--showwid",
            ],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--nmask"],
            "Fail",
            r####"--maskn requires an output file path provided to --mask
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--nmask", "--mask", "m.fa"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--win_size", "5", "-p", "10"],
            "Fail",
            r####"--winsize cannot be smaller than --period
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--win_size", "-5", "-p", "10"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--win_size", "-1", "-p", "10"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--win_size", "10", "-p", "10"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-m", "0"],
            "Fail",
            r####"--match must be > 0.0 and <= 1.0
"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "-m", "1"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-m", "1.1"],
            "Fail",
            r####"--match must be > 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--decay", "0"],
            "Fail",
            r####"--decay must be > 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--decay", "1.01"],
            "Fail",
            r####"--decay must be > 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--nr", "-0.5"],
            "Fail",
            r####"--nr must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--rn", "2"],
            "Fail",
            r####"--rn must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--ri", "2"],
            "Fail",
            r####"--ri must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--rd", "2"],
            "Fail",
            r####"--rd must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--ii", "2"],
            "Fail",
            r####"--ii must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--dd", "2"],
            "Fail",
            r####"--dd must be >= 0.0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--split_threshold", "0"],
            "Fail",
            r####"--splitval must be > 0.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--splitval", "-1"],
            "Fail",
            r####"--splitval must be > 0.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_fdr", "2"],
            "Fail",
            r####"--tune_fdr must be >= 0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_fdr", "-0.1"],
            "Fail",
            r####"--tune_fdr must be >= 0 and <= 1.0
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_fdr", "abc"],
            "Exc",
            r####"Could not convert: --tune_fdr = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_medium", "--tune_large"],
            "Fail",
            r####"Cannot use both --tune_medium and --tune_large
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_medium", "--tune_file", "x"],
            "Fail",
            r####"Cannot use both --tune_file and (--tune_small or --tune_large).
"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--tune_large", "--tune_file", "x"],
            "Fail",
            r####"Cannot use both --tune_file and (--tune_small or --tune_large).
"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--tune_file", "x"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--tune_file", ""], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--tune_only"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--tune_indel"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--tune_indels"], "Ok", r####""####),
        (&["./ultra", "ex.fa", "--mem", "--cite"], "Cite", r####""####),
        (&["./ultra", "--cite"], "Cite", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "-h", "-p", "abc"],
            "Exc",
            r####"Could not convert: --period = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "-p", "abc", "-h"],
            "Exc",
            r####"Could not convert: --period = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p", "5", "-t", "abc", "-h"],
            "Exc",
            r####"Could not convert: --threads = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-p"],
            "Exc",
            r####"--period: 1 required UINT missing"####,
        ),
        (
            &["./ultra", "-p", "abc", "ex.fa", "--mem"],
            "Exc",
            r####"Could not convert: --period = abc"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--hs"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--hidesettings", "--hide_settings"],
            "Ok",
            r####""####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--readall", "-r"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--read_all=false", "--readall"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--max_consensus", "-1"],
            "Ok",
            r####""####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--max_consensus", "99999999999"],
            "Exc",
            r####"Could not convert: --max_consensus = 99999999999"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--max_split", "-3"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--min_length", "-1"],
            "Exc",
            r####"Could not convert: --min_length = -1"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-d", "5.5"],
            "Exc",
            r####"Could not convert: --deletes = 5.5"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "-i", "abc"],
            "Exc",
            r####"Could not convert: --inserts = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--overlap", "abc"],
            "Exc",
            r####"Could not convert: --overlap = abc"####,
        ),
        (&["./ultra", "ex.fa", "--mem", "--overlap", "-7"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "--mem", "--min_split_window", "abc"],
            "Exc",
            r####"Could not convert: --min_split_window = abc"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "--split_depth", "-2"],
            "Exc",
            r####"Could not convert: --split_depth = -2"####,
        ),
        (
            &["./ultra", "ex.fa", "--mem", "/dev/null", "/dev/null"],
            "Extras",
            r####"The following arguments were not expected: /dev/null /dev/null"####,
        ),
        (&["./ultra", "-p", "5", "ex.fa", "--mem"], "Ok", r####""####),
        (
            &["./ultra", "ex.fa", "-p", "5", "-p", "6"],
            "Exc",
            r####"--period: At Most 1 required but received 2"####,
        ),
    ];
}
