//! JSON output (`JSONFileWriter.hpp` / `JSONFileWriter.cpp`).

use std::io::{self, Write};

use super::{c_str, json_escape, pval_for_score, to_string_f32, to_string_f64, RepeatFileWriter};
use crate::cli::Settings;
use crate::repeat::RepeatRegion;

/// C++ `JSONFileWriter`.
#[derive(Default)]
pub struct JSONFileWriter {
    repeats_output: i32,
}

impl JSONFileWriter {
    pub fn new() -> Self {
        JSONFileWriter { repeats_output: 0 }
    }

    fn output_json_key_value(
        out: &mut dyn Write,
        key: &str,
        value: &[u8],
        quotes: bool,
    ) -> io::Result<()> {
        // `%s` of `value.c_str()`.
        let value = c_str(value);
        if !quotes {
            write!(out, ",\n\"{}\": ", key)?;
            out.write_all(value)
        } else {
            write!(out, ",\n\"{}\": \"", key)?;
            out.write_all(value)?;
            out.write_all(b"\"")
        }
    }

    /// C++ `StringForSubRepeat`. The unused `frac_rep`/`subScore` are dropped.
    pub fn string_for_sub_repeat(r: &RepeatRegion, split_index: i32, start_pos: i32) -> Vec<u8> {
        Self::string_for_sub_repeat_hiding(r, split_index, start_pos, None)
    }

    /// [`string_for_sub_repeat`](Self::string_for_sub_repeat), with the
    /// consensus hidden as `--max_consensus` hides it when `hide` is
    /// `Some(max_consensus_period)`: "." for a period above it, no consensus
    /// below 0.
    fn string_for_sub_repeat_hiding(r: &RepeatRegion, split_index: i32, start_pos: i32, hide: Option<i32>) -> Vec<u8> {
        let splits = r.splits.as_ref().expect("splits is null");
        let consensi = r.consensi.as_ref().expect("consensi is null");
        let mut s: Vec<u8> = b"{".to_vec();

        let start = start_pos;
        let mut consensus_position = splits.len();
        if split_index >= 0 {
            consensus_position = split_index as usize;
        }

        s.extend_from_slice(b"\"Start\": ");
        s.extend_from_slice(start.to_string().as_bytes());
        match hide {
            None => {
                s.extend_from_slice(b",\n\"Consensus\": \"");
                // `+= consensi->at(i).c_str()`: stops at a NUL.
                s.extend_from_slice(c_str(&consensi[consensus_position]));
                s.extend_from_slice(b"\"");
            }
            Some(max) if max >= 0 => s.extend_from_slice(b",\n\"Consensus\": \".\""),
            Some(_) => {}
        }
        s.push(b'}');
        s
    }

    /// C++ `SubRepeatsString`.
    pub fn sub_repeats_string(r: &RepeatRegion) -> Vec<u8> {
        Self::sub_repeats_string_hiding(r, None)
    }

    /// [`sub_repeats_string`](Self::sub_repeats_string) with the consensi
    /// hidden as in [`string_for_sub_repeat_hiding`].
    ///
    /// [`string_for_sub_repeat_hiding`]: Self::string_for_sub_repeat_hiding
    fn sub_repeats_string_hiding(r: &RepeatRegion, hide: Option<i32>) -> Vec<u8> {
        let splits = r.splits.as_ref().expect("splits is null");
        let mut s: Vec<u8> = b"[".to_vec();

        let mut start: i32 = 0;
        for (i, &split_i) in splits.iter().enumerate() {
            if split_i >= 0 {
                if start > 0 {
                    s.extend_from_slice(b",\n");
                }
                s.extend_from_slice(&Self::string_for_sub_repeat_hiding(r, i as i32, start, hide));
                start = split_i;
            }
        }

        if start > 0 {
            s.extend_from_slice(b",\n");
        }
        s.extend_from_slice(&Self::string_for_sub_repeat_hiding(r, -1, start, hide));
        s.push(b']');
        s
    }
}

impl RepeatFileWriter for JSONFileWriter {
    fn initialize_writer(&mut self, _settings: &Settings, out: &mut dyn Write) -> io::Result<()> {
        out.write_all(b"{\"Repeats\": [\n")
    }

    // Indexed loops: entries up to `repeat_length` must exist, as the C++ assumes.
    #[allow(clippy::needless_range_loop)]
    fn write_repeat(
        &mut self,
        settings: &Settings,
        repeat: &RepeatRegion,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        if self.repeats_output > 0 {
            out.write_all(b",\n\n")?;
        }

        self.repeats_output += 1;

        // The header line up to a NUL (`%s` of `c_str()`), without the `\r`
        // that a file with CRLF line endings leaves at its end.
        let name = repeat.sequence_name.split('\0').next().unwrap_or("");
        let name = name.strip_suffix('\r').unwrap_or(name);
        out.write_all(b"{\"SequenceName\": \"")?;
        out.write_all(json_escape(name).as_bytes())?;
        out.write_all(b"\"")?;

        let kv = |out: &mut dyn Write, key: &str, value: String| {
            Self::output_json_key_value(out, key, value.as_bytes(), false)
        };

        kv(out, "Start", repeat.sequence_start.to_string())?;
        kv(out, "Length", repeat.repeat_length.to_string())?;
        kv(out, "Period", repeat.repeat_period.to_string())?;
        kv(out, "Score", to_string_f32(repeat.region_score))?;
        if settings.pval {
            let pval = pval_for_score(settings, repeat.region_score);
            // %g as in the TSV: the C++ prints %f, so every p-value below
            // 5e-7 (the best-scoring repeats) came out as 0.000000.
            let text = if settings.upstream_compat { to_string_f64(pval) } else { super::fmt_g(pval) };
            kv(out, "PVal", text)?;
        }

        if settings.show_counts {
            // unsigned long - int + int, then / int, all in unsigned long.
            let copies = repeat
                .repeat_length
                .wrapping_sub(repeat.insertions as i64 as u64)
                .wrapping_add(repeat.deletions as i64 as u64)
                / (repeat.repeat_period as i64 as u64);
            kv(out, "Copies", copies.to_string())?;
            kv(out, "Substitutions", repeat.mismatches.to_string())?;
            kv(out, "Insertions", repeat.insertions.to_string())?;
            kv(out, "Deletions", repeat.deletions.to_string())?;
        }

        // --max_consensus applies to JSON as to TSV and BED (the C++ ignores
        // it here): "." for a period above it, no consensus below 0.
        let hide_consensus = !settings.upstream_compat && settings.max_consensus_period < repeat.repeat_period;
        if !hide_consensus {
            Self::output_json_key_value(out, "Consensus", &repeat.string_consensus, true)?;
        } else if settings.max_consensus_period >= 0 {
            Self::output_json_key_value(out, "Consensus", b".", true)?;
        }

        // --show_wid: the C++ accepts it and prints nothing.
        if settings.show_wid && !settings.upstream_compat {
            kv(out, "WindowID", repeat.read_id.to_string())?;
        }

        if settings.show_seq {
            Self::output_json_key_value(out, "Sequence", &repeat.sequence, true)?;
        }

        if settings.show_trace {
            Self::output_json_key_value(out, "Traceback", &repeat.traceback, true)?;
        }

        if settings.show_logo_nums {
            if let Some(ln) = &repeat.logo_numbers {
                // Upstream 1.2.2 writes entry 0 twice.
                let mut logo_numbers = format!("\"{}", ln[0]);
                for i in 1..repeat.repeat_length as usize {
                    logo_numbers.push(',');
                    logo_numbers.push_str(&ln[i].to_string());
                }
                logo_numbers.push('"');
                kv(out, "LogoNumbers", logo_numbers)?;
            }
        }

        if settings.show_deltas {
            if let Some(scores) = &repeat.scores {
                let mut deltas = String::from("[");
                for i in 0..repeat.repeat_length as usize {
                    if i > 0 {
                        deltas.push(',');
                    }
                    deltas.push_str(&to_string_f32(scores[i]));
                }
                deltas.push(']');
                kv(out, "PositionScoreDeltas", deltas)?;
            }
        }

        let mut number_of_valid_splits = 0;
        if let Some(splits) = &repeat.splits {
            for &s in splits {
                if s > 0 {
                    number_of_valid_splits += 1;
                }
            }
        }

        if number_of_valid_splits > 0 {
            let sub_repeats = Self::sub_repeats_string_hiding(repeat, hide_consensus.then_some(settings.max_consensus_period));
            Self::output_json_key_value(out, "SubRepeats", &sub_repeats, false)?;
        }
        out.write_all(b"}")
    }

    fn end_writer(&mut self, _settings: &Settings, out: &mut dyn Write) -> io::Result<()> {
        out.write_all(b"]\n}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> RepeatRegion {
        let mut r = RepeatRegion::new();
        r.sequence_name = "chr1 description".into();
        r.sequence_start = 100;
        r.repeat_length = 4;
        r.repeat_period = 2;
        r.region_score = 7.5;
        r.string_consensus = b"CA".to_vec();
        r
    }

    #[test]
    fn two_repeats_and_trailer() {
        let s = Settings::default();
        let mut w = JSONFileWriter::new();
        let mut out = Vec::new();
        w.initialize_writer(&s, &mut out).unwrap();
        w.write_repeat(&s, &region(), &mut out).unwrap();
        w.write_repeat(&s, &region(), &mut out).unwrap();
        w.end_writer(&s, &mut out).unwrap();
        let one = "{\"SequenceName\": \"chr1 description\",\n\"Start\": 100,\n\"Length\": 4,\n\
                   \"Period\": 2,\n\"Score\": 7.500000,\n\"Consensus\": \"CA\"}";
        assert_eq!(
            String::from_utf8(out).unwrap(),
            format!("{{\"Repeats\": [\n{one},\n\n{one}]\n}}")
        );
    }

    /// --show_wid prints the window (read) id; --max_consensus hides the
    /// consensus as in TSV and BED: "." above it, nothing below 0. The C++
    /// does neither.
    #[test]
    fn window_id_and_hidden_consensus() {
        let mut r = region();
        r.read_id = 7;
        r.splits = Some(vec![2]);
        r.consensi = Some(vec![b"CA".to_vec(), b"TA".to_vec()]);
        let text = |s: &Settings| {
            let mut out = Vec::new();
            JSONFileWriter::new().write_repeat(s, &r, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        };
        let mut s = Settings::default();
        s.json_out = true;
        s.show_wid = true;
        assert!(text(&s).contains("\"WindowID\": 7"), "{}", text(&s));
        s.max_consensus_period = 1;
        let t = text(&s);
        assert!(t.contains("\"Consensus\": \".\",") && !t.contains("\"CA\"") && !t.contains("\"TA\""), "{t}");
        s.max_consensus_period = -1;
        assert!(!text(&s).contains("\"Consensus\""), "{}", text(&s));
        s.upstream_compat = true;
        let t = text(&s);
        assert!(t.contains("\"Consensus\": \"CA\"") && !t.contains("WindowID"), "{t}");
    }

    /// Small p-values keep their digits (%g, as in the TSV); the C++ prints
    /// %f, which turns them into 0.000000.
    #[test]
    fn pval_is_printed_with_g() {
        let mut s = Settings::default();
        s.pval = true;
        let mut r = region();
        r.region_score = 42.2;
        let text = |s: &Settings| {
            let mut out = Vec::new();
            JSONFileWriter::new().write_repeat(s, &r, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        };
        assert!(text(&s).contains("\"PVal\": 6.45705e-13,"), "{}", text(&s));
        s.upstream_compat = true;
        assert!(text(&s).contains("\"PVal\": 0.000000,"), "{}", text(&s));
    }

    #[test]
    fn sequence_name_is_escaped() {
        let s = Settings::default();
        let mut out = Vec::new();
        let mut r = region();
        r.sequence_name = "id \"quoted\" back\\slash\ttab\r".into();
        JSONFileWriter::new().write_repeat(&s, &r, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.starts_with(
                "{\"SequenceName\": \"id \\\"quoted\\\" back\\\\slash\\ttab\",\n\"Start\": 100,"
            ),
            "{text}"
        );
    }

    #[test]
    fn empty_output() {
        let s = Settings::default();
        let mut w = JSONFileWriter::new();
        let mut out = Vec::new();
        w.initialize_writer(&s, &mut out).unwrap();
        w.end_writer(&s, &mut out).unwrap();
        assert_eq!(out, b"{\"Repeats\": [\n]\n}");
    }

    #[test]
    fn details_and_subrepeats() {
        let mut s = Settings::default();
        s.show_logo_nums = true;
        s.show_deltas = true;
        s.show_counts = true;
        s.pval = true;
        let mut r = region();
        r.logo_numbers = Some(vec![0, 1, -1, 0]);
        r.scores = Some(vec![1.0, 0.5, -0.25, 2.0]);
        r.splits = Some(vec![2, -1]);
        r.consensi = Some(vec![b"CA".to_vec(), b"CA".to_vec(), b"TA".to_vec()]);
        let mut w = JSONFileWriter::new();
        let mut out = Vec::new();
        w.write_repeat(&s, &r, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let pval = crate::writers::fmt_g(pval_for_score(&s, r.region_score));
        assert!(text.contains(&format!(",\n\"PVal\": {pval},")), "{text}");
        assert_eq!(pval, "5.99591e-05");
        assert!(
            text.contains(",\n\"Copies\": 2,\n\"Substitutions\": 0"),
            "{text}"
        );
        assert!(
            text.contains(",\n\"LogoNumbers\": \"0,1,-1,0\""),
            "{text}"
        );
        assert!(
            text.contains(",\n\"PositionScoreDeltas\": [1.000000,0.500000,-0.250000,2.000000]"),
            "{text}"
        );
        assert!(
            text.ends_with(
                ",\n\"SubRepeats\": [{\"Start\": 0,\n\"Consensus\": \"CA\"},\n\
                 {\"Start\": 2,\n\"Consensus\": \"TA\"}]}"
            ),
            "{text}"
        );
    }
}
