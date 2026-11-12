//! JSON output (`JSONFileWriter.hpp` / `JSONFileWriter.cpp`).

use std::io::{self, Write};

use super::{c_str, pval_for_score, to_string_f32, to_string_f64, RepeatFileWriter};
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
        s.extend_from_slice(b",\n\"Consensus\": \"");
        // `+= consensi->at(i).c_str()`: stops at a NUL.
        s.extend_from_slice(c_str(&consensi[consensus_position]));
        s.extend_from_slice(b"\"}");
        s
    }

    /// C++ `SubRepeatsString`.
    pub fn sub_repeats_string(r: &RepeatRegion) -> Vec<u8> {
        let splits = r.splits.as_ref().expect("splits is null");
        let mut s: Vec<u8> = b"[".to_vec();

        let mut start: i32 = 0;
        for (i, &split_i) in splits.iter().enumerate() {
            if split_i >= 0 {
                if start > 0 {
                    s.extend_from_slice(b",\n");
                }
                s.extend_from_slice(&Self::string_for_sub_repeat(r, i as i32, start));
                start = split_i;
            }
        }

        if start > 0 {
            s.extend_from_slice(b",\n");
        }
        s.extend_from_slice(&Self::string_for_sub_repeat(r, -1, start));
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

        out.write_all(b"{\"SequenceName\": \"")?;
        out.write_all(c_str(repeat.sequence_name.as_bytes()))?;
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
            kv(out, "PVal", to_string_f64(pval))?;
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

        Self::output_json_key_value(out, "Consensus", &repeat.string_consensus, true)?;

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
            let sub_repeats = Self::sub_repeats_string(repeat);
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
        assert!(text.contains(",\n\"PVal\": 0.000060,"), "{text}");
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
