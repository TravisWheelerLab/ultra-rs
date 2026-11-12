//! A windowed driver for [`UModel`]. Port of `umodeldriver.hpp` /
//! `umodeldriver.cpp`.
//!
//! Upstream compiles this file but nothing calls it: `Ultra` drives the model
//! itself (`Ultra::AnalyzeSequenceWindow`). It is kept for parity.

use crate::sequence_window::SequenceWindow;
use crate::umodel::UModel;

/// A candidate repeat region (C++ `str_region`).
#[derive(Debug, Clone)]
pub struct StrRegion {
    // Where in the sequence this occurs
    pub sequence_start: i32,
    pub sequence_end: i32,

    pub starting_period: i32,

    pub start_score: f64,
    pub max_score: f64,
    pub d_score: f64,

    // Traceback/score info
    pub ar_start: i32,
    pub ar_end: i32,
    pub ar_length: i32,

    /// Never set upstream (C++ `char *traceback = NULL`).
    pub traceback: Option<Vec<u8>>,
    /// Never set upstream (C++ `double *scores = NULL`).
    pub scores: Option<Vec<f64>>,
    pub tracer: String,
    pub window_index: i32,
    pub window_start: i32,
    pub window_end: i32,

    /// Uninitialised in C++ until assigned; 0 here.
    pub last_adjustment: f64,
    /// Uninitialised in C++ until assigned; 0 here.
    pub total_adjustment: f64,

    pub continues: bool,
}

impl Default for StrRegion {
    fn default() -> Self {
        StrRegion {
            sequence_start: 0,
            sequence_end: 0,
            starting_period: 0,
            start_score: 0.0,
            max_score: -10000.0,
            d_score: 0.0,
            ar_start: 0,
            ar_end: 0,
            ar_length: 0,
            traceback: None,
            scores: None,
            tracer: String::new(),
            window_index: 0,
            window_start: 0,
            window_end: 0,
            last_adjustment: 0.0,
            total_adjustment: 0.0,
            continues: false,
        }
    }
}

/// C++ `UModelDriver`. Borrows the model; the C++ `matrix` member is
/// `model.matrix`.
#[derive(Debug)]
pub struct UModelDriver<'a> {
    pub model: &'a mut UModel,

    /// `None` until the first [`analyze_window`](Self::analyze_window).
    pub repeats: Option<Vec<StrRegion>>,

    pub window_count: i32,

    pub window_size: i32,
    pub threshold: f64,

    pub window_start: i32,
    pub window_end: i32,

    /// `matrix.length + 2` entries (uninitialised in C++, zeroed here).
    pub current_traceback: Vec<i32>,

    pub adjustment_threshold: f64,
}

impl<'a> UModelDriver<'a> {
    pub fn new(m: &'a mut UModel, ws: i32) -> Self {
        let n = m.matrix.length + 2;
        UModelDriver {
            model: m,
            repeats: None,
            window_count: 0,
            window_size: ws,
            threshold: 10.0,
            window_start: 0,
            window_end: 0,
            current_traceback: vec![0; n],
            adjustment_threshold: 10000.0,
        }
    }

    pub fn analyze_sequence(&mut self, sequence: &SequenceWindow, start: i32, end: i32) {
        let mut tstart = start;
        let mut tend = end;

        if start >= end {
            return;
        }
        let mut number_of_passes = (end - start) / self.window_size;

        if number_of_passes * self.window_size < (end - start) {
            number_of_passes += 1;
        }

        let mlen = self.model.matrix.length as u64;

        for pass in 0..number_of_passes {
            self.window_count += 1;
            if pass == 0 {
                // int + unsigned long long, narrowed to int
                tend = (tstart as i64 as u64).wrapping_add(mlen) as i32;

                if tend > end {
                    tend = end;
                }

                for i in tstart..tend {
                    self.model.calculate_current_column(sequence, i, None);
                }
            } else {
                if tstart >= end {
                    return;
                }

                for i in tstart..tend {
                    self.model.calculate_current_column(sequence, i, None);
                }
            }

            let should_flush = pass == number_of_passes - 1;

            let mut length = mlen as i32;
            if tend - start < length {
                length = tend - start;
            }

            let mut nuc_start = (tend as i64 as u64).wrapping_sub(mlen) as i32;
            if nuc_start < 0 {
                nuc_start = 0;
            }

            self.analyze_window(length, nuc_start, should_flush);

            // abs(float) resolves to the float overload; compared in double.
            let v = self.model.matrix.previous_score_column()[0].abs();
            if v as f64 > self.adjustment_threshold {
                let n = (mlen.wrapping_sub(self.window_size as i64 as u64)).wrapping_add(2) as i32;
                self.model.matrix.adjust_score_to_zero(n);
            }

            tstart = tend;
            tend = tstart + self.window_size;

            if tend > end {
                tend = end;
            }
        }
    }

    pub fn analyze_window(&mut self, length: i32, start: i32, flush: bool) {
        let repeats = self.repeats.get_or_insert_with(Vec::new);
        let matrix = &mut self.model.matrix;

        let length = length - 1;

        matrix.current_traceback_column_mut()[0] = 0;
        matrix.forward_traceback(&mut self.current_traceback, length, 0);

        let mut end = length;
        if !flush {
            // int - (unsigned long long - int), narrowed to int
            let k = (matrix.length as u64).wrapping_sub(self.window_size as i64 as u64);
            end = (end as i64 as u64).wrapping_sub(k) as i32;
        }

        let mut rep = StrRegion::default();

        let mut p = length - 1;
        let mut i: i32 = 0;
        while i < end {
            let mut c = self.current_traceback[i as usize];

            if c != 0 {
                rep.start_score = matrix.previous_score(c, p) as f64;
                rep.max_score = crate::umatrix::NEG_INF;
                rep.window_start = i;
                rep.starting_period = matrix.cell_descriptions[c as usize].order;
                rep.sequence_start = i + start - rep.starting_period;

                i += 1;
                p -= 1;
                while i < end {
                    c = self.current_traceback[i as usize];

                    let s = matrix.previous_score(c, p) as f64;
                    if s > rep.max_score {
                        rep.max_score = s;
                    }

                    if c == 0 {
                        break;
                    }
                    i += 1;
                    p -= 1;
                }

                rep.window_end = i;
                rep.ar_end = i;
                rep.sequence_end = start + i;
                rep.d_score = rep.max_score - rep.start_score;
                rep.window_index = self.window_count;
                rep.last_adjustment = matrix.last_score_adjustment as f64;
                rep.total_adjustment = matrix.total_score_adjustment as f64;

                let mut prev = StrRegion::default();
                if !repeats.is_empty() {
                    prev = repeats.last().unwrap().clone();
                }

                if rep.d_score >= self.threshold && rep.sequence_end > prev.sequence_end {
                    if !repeats.is_empty() && prev.sequence_end >= rep.sequence_start {
                        prev.continues = true;
                        repeats.pop();
                        repeats.push(prev);
                    }

                    repeats.push(rep.clone());
                }
            }
            i += 1;
            p -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbol::symbol_for_char;

    /// Expected values were dumped from the upstream C++ driven the same way.
    #[test]
    fn analyze_window_matches_cpp() {
        let s = "ACGTTGCANNNNGCAGCAGCAGCAGCATCAGCAGCAGCGCAGCAGCAGGTTACGATATATATANATATATATCGATGCTAGTACCAGTCAGTTAGCAGTCAGTTAGCAGTCATCAGTTAGCAGGTCAGTTAGCACCGTA";
        let mlen = s.len() as i32 + 3;
        let mut m = UModel::new(7, 0, 0, mlen);
        m.period_decay = 0.85;
        m.set_match_probabilities(0.7);
        m.set_atcg_probabilities(0.25, 0.25, 0.25, 0.25);
        m.tp_match_to_zero = 0.05;
        m.tp_match_to_insertion = 0.02;
        m.tp_match_to_deletion = 0.02;
        m.tp_consecutive_insertion = 0.02;
        m.tp_consecutive_deletion = 0.02;
        m.calculate_scores();

        let mut w = SequenceWindow::new(mlen as usize, 0);
        w.length = s.len();
        for (i, c) in s.bytes().enumerate() {
            w.seq_mem[i] = symbol_for_char(c);
        }
        m.matrix.restart_matrix();
        for i in 0..s.len() as i32 {
            m.calculate_current_column(&w, i, None);
        }
        m.calculate_current_column_without_emission(None);
        m.matrix.adjust_score_to_zero(5);
        m.matrix.adjust_score_to_zero_from(2, 3);

        let mut d = UModelDriver::new(&mut m, mlen - 4);
        d.threshold = 1.0;
        d.analyze_window(s.len() as i32, 3, true);
        d.analyze_window(s.len() as i32, 5, false);

        let got: Vec<_> = d
            .repeats
            .as_ref()
            .unwrap()
            .iter()
            .map(|r| {
                (
                    r.sequence_start,
                    r.sequence_end,
                    r.starting_period,
                    r.start_score.to_bits(),
                    r.max_score.to_bits(),
                    r.d_score.to_bits(),
                    r.window_start,
                    r.total_adjustment.to_bits(),
                    r.continues,
                )
            })
            .collect();
        let adj = 0xc034928480000000u64;
        assert_eq!(
            got,
            vec![
                (
                    14,
                    50,
                    3,
                    0xc037474d80000000,
                    0x4009534400000000,
                    0x403a71b600000000,
                    14,
                    adj,
                    false
                ),
                (
                    56,
                    74,
                    2,
                    0xc01a844fc0000000,
                    0x401fc44900000000,
                    0x402d244c60000000,
                    55,
                    adj,
                    true
                ),
                (
                    58,
                    76,
                    2,
                    0xc01a844fc0000000,
                    0x401fc44900000000,
                    0x402d244c60000000,
                    55,
                    adj,
                    false
                ),
            ]
        );
    }
}
