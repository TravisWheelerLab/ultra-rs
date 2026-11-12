//! TSV output (`TabFileWriter.hpp` / `TabFileWriter.cpp`).

use std::io::{self, Write};

use super::{c_str, clean_sequence_name, fmt_g, to_string_f32, RepeatFileWriter};
use crate::cli::Settings;
use crate::repeat::RepeatRegion;

/// C++ `TabFileWriter`.
#[derive(Default)]
pub struct TabFileWriter;

impl TabFileWriter {
    pub fn new() -> Self {
        TabFileWriter
    }
}

impl RepeatFileWriter for TabFileWriter {
    fn initialize_writer(&mut self, settings: &Settings, out: &mut dyn Write) -> io::Result<()> {
        out.write_all(b"SeqID")?;
        out.write_all(b"\tStart")?;
        out.write_all(b"\tEnd")?;
        out.write_all(b"\tPeriod")?;
        out.write_all(b"\tScore")?;
        // Upstream quirk: a comma, not a tab.
        if settings.pval {
            out.write_all(b",PValue")?;
        }
        // Upstream quirk: the header tests `>= 0` but rows test `!= 0`, so
        // `--max_consensus 0` gets a Consensus header and no column.
        if settings.max_consensus_period >= 0 {
            out.write_all(b"\tConsensus")?;
        }

        if settings.show_counts {
            out.write_all(b"\t#copies\t#substitutions\t#insertions\t#deletions")?;
        }

        if settings.max_split >= 0 {
            out.write_all(b"\t#Subrepeats")?;
            out.write_all(b"\tSubrepeatStarts")?;
            if settings.max_consensus_period >= 0 {
                out.write_all(b"\tSubrepeatConsensi")?;
            }
        }

        if settings.show_seq {
            out.write_all(b"\tSequence")?;
        }

        out.write_all(b"\n")
    }

    fn write_repeat(
        &mut self,
        settings: &Settings,
        repeat: &RepeatRegion,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        let name = clean_sequence_name(&repeat.sequence_name);

        // Columns 1 (name) 2 (start) 3 (end) 4 (period) 5 (score)
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            name,
            repeat.sequence_start,
            repeat.sequence_start.wrapping_add(repeat.repeat_length),
            repeat.repeat_period,
            to_string_f32(repeat.region_score)
        )?;
        if settings.pval {
            write!(
                out,
                ",{}",
                fmt_g(super::pval_for_score(settings, repeat.region_score))
            )?;
        }

        let mut rep_con: &[u8] = b".";
        // Same test as the header. Upstream 1.2.2 tests `!= 0` here, so with
        // `--max_consensus 0` the rows have one column fewer than the header.
        if settings.max_consensus_period >= 0 {
            if settings.max_consensus_period >= repeat.repeat_period
                && !repeat.string_consensus.is_empty()
            {
                rep_con = &repeat.string_consensus;
            }

            out.write_all(b"\t")?;
            out.write_all(c_str(rep_con))?;
        }

        if settings.show_counts {
            // unsigned long - int + int, then / int, all in unsigned long.
            let copies = repeat
                .repeat_length
                .wrapping_sub(repeat.insertions as i64 as u64)
                .wrapping_add(repeat.deletions as i64 as u64)
                / (repeat.repeat_period as i64 as u64);

            write!(
                out,
                "\t{}\t{}\t{}\t{}",
                copies, repeat.mismatches, repeat.insertions, repeat.deletions
            )?;
        }

        if settings.max_split >= 0 {
            let mut starts: Vec<u8> = b"0".to_vec();
            let mut consensi: Vec<u8> = Vec::new();
            let mut number_of_valid_splits: i32 = 0;

            match &repeat.splits {
                Some(splits) if !splits.is_empty() => {
                    // `sizes` is built in the C++ but never printed.
                    for &split_i in splits.iter() {
                        if split_i > 0 {
                            starts.push(b',');
                            starts.extend_from_slice(split_i.to_string().as_bytes());
                            number_of_valid_splits += 1;
                        }
                    }

                    if settings.max_consensus_period >= 0 {
                        let rc = repeat.consensi.as_ref().expect("consensi is null");
                        for i in 0..rc.len() {
                            let mut con: &[u8] = b".";
                            if settings.max_consensus_period >= repeat.repeat_period && rc.len() > i
                            {
                                if i > 0 && rc[i] == rc[i - 1] {
                                    continue;
                                }
                                con = &rc[i];
                            }
                            if i > 0 {
                                consensi.push(b',');
                            }
                            consensi.extend_from_slice(c_str(con));
                        }
                    }

                    write!(out, "\t{}\t", number_of_valid_splits + 1)?;
                    out.write_all(&starts)?;
                    if settings.max_consensus_period >= 0 {
                        out.write_all(b"\t")?;
                        out.write_all(&consensi)?;
                    }
                }
                _ => {
                    out.write_all(b"\t1\t0")?;
                    if settings.max_consensus_period >= 0 {
                        out.write_all(b"\t")?;
                        out.write_all(c_str(rep_con))?;
                    }
                }
            }
        }

        if settings.show_seq {
            out.write_all(b"\t")?;
            out.write_all(c_str(&repeat.sequence))?;
        }

        out.write_all(b"\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> RepeatRegion {
        let mut r = RepeatRegion::new();
        r.sequence_name = "chr1 description".into();
        r.sequence_start = 100;
        r.repeat_length = 30;
        r.repeat_period = 3;
        r.region_score = 42.2;
        r.string_consensus = b"CAG".to_vec();
        r.sequence = b"CAGCAGCAG".to_vec();
        r.insertions = 1;
        r.deletions = 2;
        r.mismatches = 4;
        r
    }

    fn render(settings: &Settings, r: &RepeatRegion) -> String {
        let mut out = Vec::new();
        let mut w = TabFileWriter::new();
        w.initialize_writer(settings, &mut out).unwrap();
        w.write_repeat(settings, r, &mut out).unwrap();
        w.end_writer(settings, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn default_columns() {
        let s = Settings::default();
        assert_eq!(
            render(&s, &region()),
            "SeqID\tStart\tEnd\tPeriod\tScore\tConsensus\t#Subrepeats\tSubrepeatStarts\tSubrepeatConsensi\n\
             chr1\t100\t130\t3\t42.200001\tCAG\t1\t0\tCAG\n"
        );
    }

    #[test]
    fn splits_counts_and_sequence() {
        let mut s = Settings::default();
        s.show_counts = true;
        s.show_seq = true;
        s.pval = true;
        let mut r = region();
        r.splits = Some(vec![9, -1, 21]);
        r.consensi = Some(vec![
            b"CAG".to_vec(),
            b"CAG".to_vec(),
            b"CTG".to_vec(),
            b"CTG".to_vec(),
        ]);
        let text = render(&s, &r);
        assert_eq!(
            text.lines().nth(1).unwrap(),
            "chr1\t100\t130\t3\t42.200001,6.45705e-13\tCAG\t10\t4\t1\t2\t3\t0,9,21\tCAG,CTG\tCAGCAGCAG"
        );
    }

    #[test]
    fn max_consensus_zero_keeps_the_consensus_column() {
        let mut s = Settings::default();
        s.max_consensus_period = 0;
        s.max_split = 0;
        let text = render(&s, &region());
        assert_eq!(
            text,
            "SeqID\tStart\tEnd\tPeriod\tScore\tConsensus\t#Subrepeats\tSubrepeatStarts\tSubrepeatConsensi\n\
             chr1\t100\t130\t3\t42.200001\t.\t1\t0\t.\n"
        );
    }
}
