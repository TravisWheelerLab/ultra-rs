//! BPAF output: each repeat aligned to its consensus (see [`crate::bpaf`]).
//! Not in upstream ULTRA.

use std::io::{self, Write};

use super::RepeatFileWriter;
use crate::bpaf::{TandemRecord, TandemWriter};
use crate::cli::Settings;
use crate::repeat::RepeatRegion;

/// The `mats` entry for ULTRA's repeat score.
pub const ULTRA_SCORING: &str = "ultra";

/// Writes one record per repeat, or per sub-repeat of a split repeat, each
/// with the repeat's score. The alignment is made here, when the repeat is
/// written, so runs without BPAF output do no extra work.
#[derive(Default)]
pub struct BPAFFileWriter {
    writer: TandemWriter,
}

impl BPAFFileWriter {
    pub fn new() -> Self {
        Self::default()
    }
}

/// The first word of a FASTA header line: the sequence's id.
pub fn sequence_id(name: &str) -> &str {
    name.split(|c: char| c.is_whitespace() || c == '\0')
        .next()
        .unwrap_or("")
}

impl RepeatFileWriter for BPAFFileWriter {
    fn initialize_writer(&mut self, _settings: &Settings, out: &mut dyn Write) -> io::Result<()> {
        self.writer.write_header(out)
    }

    fn write_repeat(
        &mut self,
        _settings: &Settings,
        repeat: &RepeatRegion,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        let Some(parts) = repeat.consensus_alignment() else {
            return Ok(());
        };
        let query = sequence_id(&repeat.sequence_name);
        for part in &parts {
            self.writer.write_record(
                out,
                &TandemRecord {
                    query,
                    q_start: repeat.sequence_start + part.start as u64 + 1,
                    unit: &part.consensus,
                    alignment: &part.alignment,
                    score: Some(repeat.region_score),
                    scoring: ULTRA_SCORING,
                },
            )?;
        }
        Ok(())
    }

    fn end_writer(&mut self, _settings: &Settings, out: &mut dyn Write) -> io::Result<()> {
        self.writer.finish(out)
    }
}
