//! BED output (`BEDFileWriter.hpp` / `BEDFileWriter.cpp`).

use std::io::{self, Write};

use super::{c_str, clean_sequence_name, RepeatFileWriter};
use crate::cli::Settings;
use crate::repeat::RepeatRegion;

/// C++ `BEDFileWriter`. No header or trailer.
#[derive(Default)]
pub struct BEDFileWriter {
    /// Unused upstream.
    pub output_multiline_split_repeats: bool,
}

impl BEDFileWriter {
    pub fn new() -> Self {
        BEDFileWriter {
            output_multiline_split_repeats: false,
        }
    }
}

impl RepeatFileWriter for BEDFileWriter {
    fn write_repeat(
        &mut self,
        settings: &Settings,
        repeat: &RepeatRegion,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        let name = clean_sequence_name(&repeat.sequence_name);

        // Columns 1 (name) 2 (start) 3 (end)
        write!(
            out,
            "{}\t{}\t{}",
            name,
            repeat.sequence_start,
            repeat.sequence_start.wrapping_add(repeat.repeat_length)
        )?;

        let period = repeat.repeat_period.to_string();
        let mut rep_con: &[u8] = period.as_bytes();
        if settings.max_consensus_period >= repeat.repeat_period
            && !repeat.string_consensus.is_empty()
        {
            rep_con = &repeat.string_consensus;
        }

        out.write_all(b"\t")?;
        out.write_all(c_str(rep_con))?;
        out.write_all(b"\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_format() {
        let mut s = Settings::default();
        let mut r = RepeatRegion::new();
        r.sequence_name = "chr2\tx".into();
        r.sequence_start = 5;
        r.repeat_length = 12;
        r.repeat_period = 4;
        r.string_consensus = b"ACGT".to_vec();
        let mut out = Vec::new();
        let mut w = BEDFileWriter::new();
        w.initialize_writer(&s, &mut out).unwrap();
        w.write_repeat(&s, &r, &mut out).unwrap();
        s.max_consensus_period = 3;
        w.write_repeat(&s, &r, &mut out).unwrap();
        w.end_writer(&s, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr2\t5\t17\tACGT\nchr2\t5\t17\t4\n"
        );
    }
}
