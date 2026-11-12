//! Repeat output writers (`RepeatFileWriter.hpp` and the TSV, JSON and BED
//! writers), plus the output sink they write to and the printf-style number
//! formatting shared with `main`.
//!
//! The C++ writers keep an `Ultra *owner` and a `FILE *out`. Here the writer
//! gets the settings (everything it reads through `owner`, including the
//! p-value parameters) and the sink as arguments on every call.

pub mod bed;
pub mod json;
pub mod tab;

pub use bed::BEDFileWriter;
pub use json::JSONFileWriter;
pub use tab::TabFileWriter;

use std::fs::File;
use std::io::{self, BufWriter, Stdout, Write};

use crate::cli::Settings;
use crate::repeat::RepeatRegion;

/// C++ `RepeatFileWriter`: the interface every output format implements.
/// The default methods do nothing, as in the C++ base class.
pub trait RepeatFileWriter: Send {
    /// Writes the format's header, if any.
    fn initialize_writer(&mut self, _settings: &Settings, _out: &mut dyn Write) -> io::Result<()> {
        Ok(())
    }
    /// Writes one repeat.
    fn write_repeat(
        &mut self,
        _settings: &Settings,
        _repeat: &RepeatRegion,
        _out: &mut dyn Write,
    ) -> io::Result<()> {
        Ok(())
    }
    /// Writes the format's trailer, if any.
    fn end_writer(&mut self, _settings: &Settings, _out: &mut dyn Write) -> io::Result<()> {
        Ok(())
    }
}

/// An output stream: the C++ `FILE *` that is either `stdout` or an
/// `fopen`ed file. Both are buffered; `flush` is called where the C++
/// `fflush`es, and dropping the sink flushes (and closes) it as `fclose`
/// would.
///
/// Stdout output is buffered here and pushed into `std::io::stdout()` on
/// every flush. Ultra flushes at the end of every `OutputRepeats`, before any
/// other component prints, so stdout text keeps the C++ order.
pub enum OutSink {
    Stdout(BufWriter<Stdout>),
    File(BufWriter<File>),
}

impl OutSink {
    pub fn stdout() -> Self {
        OutSink::Stdout(BufWriter::with_capacity(1 << 16, io::stdout()))
    }

    pub fn file(f: File) -> Self {
        OutSink::File(BufWriter::with_capacity(1 << 16, f))
    }

    pub fn is_stdout(&self) -> bool {
        matches!(self, OutSink::Stdout(_))
    }
}

impl Write for OutSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutSink::Stdout(w) => w.write(buf),
            OutSink::File(w) => w.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            OutSink::Stdout(w) => w.write_all(buf),
            OutSink::File(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutSink::Stdout(w) => w.flush(),
            OutSink::File(w) => w.flush(),
        }
    }
}

/// C++ `PvalForScore`: `exp(-(score - loc) / scale) * freq`, in double.
#[allow(clippy::neg_multiply)] // written as upstream
pub fn pval_for_score(settings: &Settings, score: f32) -> f64 {
    let loc = settings.p_value_loc as f64;
    let scale = settings.p_value_scale as f64;
    let freq = settings.p_value_freq as f64;
    (-1.0 * (score as f64 - loc) / scale).exp() * freq
}

/// C++ `Log2PvalForScore` (`period` is unused, as upstream).
#[allow(clippy::neg_multiply)] // written as upstream
pub fn log2_pval_for_score(settings: &Settings, score: f32, _period: f32) -> f64 {
    let loc = settings.p_value_loc as f64;
    let scale = settings.p_value_scale as f64;
    let freq = settings.p_value_freq as f64;
    ((-1.0 * (score as f64 - loc) / scale).exp() * freq).log2()
}

/// glibc spelling of a non-finite value for `%f`/`%g`, or `None`.
fn non_finite(v: f64) -> Option<String> {
    if v.is_nan() {
        Some(if v.is_sign_negative() { "-nan" } else { "nan" }.to_string())
    } else if v.is_infinite() {
        Some(if v < 0.0 { "-inf" } else { "inf" }.to_string())
    } else {
        None
    }
}

/// printf `%.<prec>f` (glibc), including `nan`/`-nan`/`inf`.
pub fn fmt_f(v: f64, prec: usize) -> String {
    match non_finite(v) {
        Some(s) => s,
        None => format!("{:.*}", prec, v),
    }
}

/// `std::to_string(float)`: `%f` of the value widened to double.
pub fn to_string_f32(v: f32) -> String {
    fmt_f(v as f64, 6)
}

/// `std::to_string(double)`: `%f`.
pub fn to_string_f64(v: f64) -> String {
    fmt_f(v, 6)
}

/// Removes trailing zeros after a decimal point, and the point itself if
/// nothing follows it (the `%g` rule without `#`).
fn strip_fraction_zeros(s: &mut String) {
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
}

/// printf `%g` (glibc, precision 6).
pub fn fmt_g(v: f64) -> String {
    const P: i32 = 6;
    if let Some(s) = non_finite(v) {
        return s;
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    // The exponent after rounding to P significant digits decides the style.
    let e_form = format!("{:.*e}", (P - 1) as usize, v);
    let (mantissa, exp) = e_form.split_once('e').expect("exponent form");
    let x: i32 = exp.parse().expect("exponent");
    if (-4..P).contains(&x) {
        let mut s = format!("{:.*}", (P - 1 - x) as usize, v);
        strip_fraction_zeros(&mut s);
        s
    } else {
        let mut m = mantissa.to_string();
        strip_fraction_zeros(&mut m);
        format!("{}e{}{:02}", m, if x < 0 { '-' } else { '+' }, x.abs())
    }
}

/// The TSV/BED writers' name cleanup: the sequence name up to (not
/// including) the first byte outside `[A-Za-z0-9-_.:*#]`.
pub(crate) fn clean_sequence_name(name: &str) -> &str {
    let b = name.as_bytes();
    for (i, &c) in b.iter().enumerate() {
        if c.is_ascii_alphanumeric()
            || c == b'-'
            || c == b'_'
            || c == b'.'
            || c == b':'
            || c == b'*'
            || c == b'#'
        {
            continue;
        }
        // Every earlier byte is ASCII, so `i` is a char boundary.
        return &name[..i];
    }
    name
}

/// `%s` of a C++ `std::string` via `c_str()`: stops at the first NUL.
pub(crate) fn c_str(s: &[u8]) -> &[u8] {
    match s.iter().position(|&c| c == 0) {
        Some(i) => &s[..i],
        None => s,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected strings come from glibc printf (gcc, x86-64).

    #[test]
    fn fmt_g_matches_glibc() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0"),
            (-0.0, "-0"),
            (1.0, "1"),
            (100000.0, "100000"),
            (1000000.0, "1e+06"),
            (123456789.0, "1.23457e+08"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (0.000123456789, "0.000123457"),
            (1.5e-300, "1.5e-300"),
            (2.5, "2.5"),
            (0.1, "0.1"),
            (999999.5, "1e+06"),
            (-3.25e-7, "-3.25e-07"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (0.330256, "0.330256"),
            (1.234565e-10, "1.23457e-10"),
        ];
        for &(v, want) in cases {
            assert_eq!(fmt_g(v), want, "{v:e}");
        }
    }

    #[test]
    fn fmt_f_matches_glibc() {
        assert_eq!(fmt_f(1.0 / 3.0, 6), "0.333333");
        assert_eq!(fmt_f(-0.0000001, 6), "-0.000000");
        assert_eq!(fmt_f(2.5, 0), "2");
        assert_eq!(fmt_f(0.12345, 4), "0.1235");
        assert_eq!(fmt_f(f64::NAN.copysign(-1.0), 4), "-nan");
        assert_eq!(fmt_f(f64::NAN.copysign(1.0), 4), "nan");
        assert_eq!(to_string_f32(42.2_f32), "42.200001");
    }

    #[test]
    fn name_cleanup() {
        assert_eq!(clean_sequence_name("seq1 synthetic"), "seq1");
        assert_eq!(clean_sequence_name("a-b_c.d:e*f#g"), "a-b_c.d:e*f#g");
        assert_eq!(clean_sequence_name("crlf\r"), "crlf");
        assert_eq!(clean_sequence_name("caf\u{e9}x"), "caf");
        assert_eq!(clean_sequence_name(""), "");
    }

    #[test]
    fn pval_matches_cpp() {
        // Bit patterns and %g text from the same expression compiled with gcc.
        let s = Settings::default();
        let cases: &[(f32, u64, &str)] = &[
            (42.2, 0x3D66B7FF6AF05773, "6.45705e-13"),
            (3.0, 0x3F45367204DAFAAA, "0.00064736"),
            (-5.5, 0x3FADA93CC41E2FB3, "0.0579318"),
            (1208.5, 0x05CD504667E8840D, "1.00931e-280"),
        ];
        for &(score, bits, text) in cases {
            let p = pval_for_score(&s, score);
            assert_eq!(p.to_bits(), bits, "{score}");
            assert_eq!(fmt_g(p), text);
        }
    }
}
