//! Cell descriptions and the rotating DP matrix used by ULTRA. Port of
//! `umatrix.hpp` / `umatrix.cpp`.
//!
//! `UMatrix` holds two flat matrices:
//! - a score matrix of `f32`, `cells_per_column` rows per column;
//! - a traceback matrix of `i32`, `max_period + 1` rows per column. Row `k`
//!   of a traceback column records which cell order-`k` was reached from, so
//!   it does not describe every cell, only the background row (row 0) and one
//!   row per period.
//!
//! The C++ keeps per-column pointer tables (`scoreColumns`,
//! `tracebackColumns`) and `current*`/`previous*` column pointers. Here a
//! column is addressed by index and the slices come from accessor methods
//! ([`UMatrix::score_column`], [`UMatrix::current_score_column`], ...).
//!
//! Both matrices have `length + 1` columns of storage. The C++ only fills
//! pointer-table entries `0..length`, so column `length` is backed by memory
//! but has no valid `scoreColumns` entry there; no code path reaches it when
//! the matrix is sized as `Ultra` sizes it.

/// Score used for impossible transitions. A double in C++; it is converted to
/// `f32` wherever it is stored into a float.
pub const NEG_INF: f64 = -10000000000000.0;

/// `NEG_INF` as stored into a C++ `float`.
pub const NEG_INF_F32: f32 = NEG_INF as f32;

/// Cell type (C++ `cell_type`, an `unsigned char`).
pub type CellType = u8;

pub const CT_BACKGROUND: CellType = 0;
pub const CT_MATCH: CellType = 1;
pub const CT_INSERTION: CellType = 2;
pub const CT_DELETION: CellType = 3;

/// Description of one row of a score column (C++ `cell`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Cell {
    /// What type of cell this is (`CT_*`).
    pub cell_type: CellType,
    /// The row of the parent match cell (0 for background and match cells).
    pub parent_index: i32,
    /// Period (or the parent period, for indel cells).
    pub order: i32,
    /// Which indel this row represents (1-based; 0 for non-indel cells).
    pub indel_number: i32,
}

/// C++ `PreviousColumn`: `currentColumnIndex - d`, wrapped by `length` when
/// negative. The arithmetic goes through `unsigned long long` and is narrowed
/// to `int` as in the C++.
#[inline(always)]
pub(crate) fn previous_column_of(current_column_index: usize, length: usize, d: i32) -> i32 {
    // d = currentColumnIndex - d  (computed in unsigned long long, stored in int)
    let mut d = (current_column_index as u64).wrapping_sub(d as i64 as u64) as i32;
    if d < 0 {
        // d += length  (int + unsigned long long, narrowed back to int)
        d = (d as i64 as u64).wrapping_add(length as u64) as i32;
    }
    d
}

/// The DP matrix (C++ `UMatrix`).
#[derive(Debug, Clone)]
pub struct UMatrix {
    /// Max period is included.
    pub max_period: usize,
    /// Max insertion is included.
    pub max_insertions: usize,
    /// Max deletion is included.
    pub max_deletions: usize,

    /// Number of columns in the matrix.
    pub length: usize,
    /// Number of rows in a score column.
    pub cells_per_column: usize,

    /// Flat score storage, `cells_per_column * (length + 1)` floats. Column
    /// `c` occupies `[c * cells_per_column, (c + 1) * cells_per_column)`.
    pub score_matrix: Vec<f32>,
    /// Flat traceback storage, `(max_period + 1) * (length + 1)` ints. Column
    /// `c` occupies `[c * (max_period + 1), (c + 1) * (max_period + 1))`.
    pub traceback_matrix: Vec<i32>,

    pub current_column_index: usize,
    pub previous_column_index: usize,

    /// One entry per row, `cells_per_column` entries.
    pub cell_descriptions: Vec<Cell>,

    pub last_score_adjustment: f32,
    pub total_score_adjustment: f32,

    /// Per-column best row, filled by [`UMatrix::calculate_traceback`].
    /// `length + 2` entries, zeroed at construction and never reset.
    pub traceback: Vec<i32>,
}

impl UMatrix {
    /// C++ `UMatrix(maximumPeriod, maximumInsertions, maximumDeletions,
    /// matrixLength)`.
    pub fn new(
        maximum_period: i32,
        maximum_insertions: i32,
        maximum_deletions: i32,
        matrix_length: i32,
    ) -> Self {
        // int -> unsigned long long, as the C++ assignments do.
        let mut m = UMatrix {
            max_period: maximum_period as i64 as u64 as usize,
            max_insertions: maximum_insertions as i64 as u64 as usize,
            max_deletions: maximum_deletions as i64 as u64 as usize,
            length: matrix_length as i64 as u64 as usize,
            cells_per_column: 0,
            score_matrix: Vec::new(),
            traceback_matrix: Vec::new(),
            current_column_index: 0,
            previous_column_index: 0,
            cell_descriptions: Vec::new(),
            last_score_adjustment: 0.0,
            total_score_adjustment: 0.0,
            traceback: Vec::new(),
        };
        m.create_matrix();
        m
    }

    /// Rows per traceback column (`max_period + 1`).
    #[inline(always)]
    pub fn traceback_rows(&self) -> usize {
        self.max_period + 1
    }

    /// Score column `col` (C++ `scoreColumns[col]`).
    #[inline]
    pub fn score_column(&self, col: usize) -> &[f32] {
        let n = self.cells_per_column;
        &self.score_matrix[col * n..(col + 1) * n]
    }

    #[inline]
    pub fn score_column_mut(&mut self, col: usize) -> &mut [f32] {
        let n = self.cells_per_column;
        &mut self.score_matrix[col * n..(col + 1) * n]
    }

    /// Traceback column `col` (C++ `tracebackColumns[col]`), `max_period + 1`
    /// entries.
    #[inline]
    pub fn traceback_column(&self, col: usize) -> &[i32] {
        let n = self.traceback_rows();
        &self.traceback_matrix[col * n..(col + 1) * n]
    }

    #[inline]
    pub fn traceback_column_mut(&mut self, col: usize) -> &mut [i32] {
        let n = self.traceback_rows();
        &mut self.traceback_matrix[col * n..(col + 1) * n]
    }

    /// C++ `scoreColumns[col][row]`.
    #[inline]
    pub fn score_at(&self, col: usize, row: usize) -> f32 {
        self.score_matrix[col * self.cells_per_column + row]
    }

    /// C++ `tracebackColumns[col][row]`.
    #[inline]
    pub fn traceback_at(&self, col: usize, row: usize) -> i32 {
        self.traceback_matrix[col * self.traceback_rows() + row]
    }

    /// C++ `currentScoreColumn`.
    #[inline]
    pub fn current_score_column(&self) -> &[f32] {
        self.score_column(self.current_column_index)
    }

    #[inline]
    pub fn current_score_column_mut(&mut self) -> &mut [f32] {
        self.score_column_mut(self.current_column_index)
    }

    /// C++ `previousScoreColumn`.
    #[inline]
    pub fn previous_score_column(&self) -> &[f32] {
        self.score_column(self.previous_column_index)
    }

    #[inline]
    pub fn previous_score_column_mut(&mut self) -> &mut [f32] {
        self.score_column_mut(self.previous_column_index)
    }

    /// C++ `currentTracebackColumn`.
    #[inline]
    pub fn current_traceback_column(&self) -> &[i32] {
        self.traceback_column(self.current_column_index)
    }

    #[inline]
    pub fn current_traceback_column_mut(&mut self) -> &mut [i32] {
        self.traceback_column_mut(self.current_column_index)
    }

    /// C++ `previousTracebackColumn`.
    #[inline]
    pub fn previous_traceback_column(&self) -> &[i32] {
        self.traceback_column(self.previous_column_index)
    }

    #[inline]
    pub fn previous_traceback_column_mut(&mut self) -> &mut [i32] {
        self.traceback_column_mut(self.previous_column_index)
    }

    /// Advances the current/previous column indices by one. The C++ wrap-around
    /// check is commented out upstream, so this never wraps and always returns
    /// `false`.
    pub fn move_matrix_forward(&mut self) -> bool {
        let wrap = false;

        self.previous_column_index = self.current_column_index;
        self.current_column_index += 1;

        wrap
    }

    /// C++ `AdjustScoreToZero(int numberOfColumns)`. Shifts every column by
    /// minus the previous column's background score. `number_of_columns` is
    /// unused, as upstream (the loop runs over all `length` columns).
    pub fn adjust_score_to_zero(&mut self, _number_of_columns: i32) {
        // previousScoreColumn[0] * -1 is float * int, computed in float.
        self.last_score_adjustment = self.previous_score_column()[0] * -1.0f32;
        self.total_score_adjustment += self.last_score_adjustment;

        // for (int i = 0; i < length; ++i): int compared with unsigned long long.
        let mut i: i32 = 0;
        while (i as i64 as u64) < self.length as u64 {
            let c = self.previous_column(i + 1);
            self.adjust_column(c, self.last_score_adjustment);
            i += 1;
        }
    }

    /// C++ `AdjustScoreToZero(int startingAt, int numberOfColumns)`. Resets
    /// `number_of_columns` columns from `starting_at` to row 0 = 0 and every
    /// other row = `NEG_INF`.
    pub fn adjust_score_to_zero_from(&mut self, starting_at: i32, number_of_columns: i32) {
        for i in 0..number_of_columns {
            let p = (i + starting_at) as usize;
            let col = self.score_column_mut(p);
            for v in col.iter_mut().skip(1) {
                *v = NEG_INF_F32;
            }
            col[0] = 0.0;
        }
    }

    /// Adds `adjustment` to every row of `column`.
    pub fn adjust_column(&mut self, column: i32, adjustment: f32) {
        for v in self.score_column_mut(column as usize) {
            *v += adjustment;
        }
    }

    /// Score at `row` of the column `d` before the current one, or `NEG_INF`
    /// if `d >= length`. `d = 0` is the current column.
    #[inline]
    pub fn previous_score(&self, row: i32, d: i32) -> f32 {
        let mut s = NEG_INF_F32;
        // int d compared with unsigned long long length.
        if (d as i64 as u64) < self.length as u64 {
            let c = self.previous_column(d);
            s = self.score_at(c as usize, row as usize);
        }
        s
    }

    /// Traceback at `row` of the column `d` before the current one, or 0 if
    /// `d >= length`.
    #[inline]
    pub fn previous_traceback(&self, row: i32, d: i32) -> i32 {
        let mut r = 0;
        if (d as i64 as u64) < self.length as u64 {
            let c = self.previous_column(d);
            r = self.traceback_at(c as usize, row as usize);
        }
        r
    }

    /// `current_column_index - d`, wrapped by `length` when negative.
    #[inline]
    pub fn previous_column(&self, d: i32) -> i32 {
        previous_column_of(self.current_column_index, self.length, d)
    }

    /// C++ `ForwardTraceback`. Walks the traceback matrix backwards from
    /// `previous_column_index` starting at `row`, writing a forward-oriented
    /// traceback into `trace_array[0..=window_length]` with indels moved to
    /// where they occurred.
    ///
    /// The C++ allocates `trace_array` when passed NULL; here the caller must
    /// pass a slice of at least `window_length + 1` entries. Returns `false`
    /// (after printing the C++ message to stdout) where the C++ returns NULL.
    pub fn forward_traceback(&self, trace_array: &mut [i32], window_length: i32, row: i32) -> bool {
        // int windowLength compared with unsigned long long length.
        if (window_length as i64 as u64) > self.length as u64 {
            println!("(ForwardTraceback) Length greater than matrix length, ret NULL");
            return false;
        }

        let mut row = row;
        let length = self.length as i32;
        let mut pos = self.previous_column_index as i32;

        let mut i: i32 = 0;
        while i <= window_length {
            let desc = self.cell_descriptions[row as usize];
            let row_order = desc.order;

            if desc.cell_type == CT_MATCH || desc.cell_type == CT_BACKGROUND {
                trace_array[(window_length - i) as usize] = row;
                row = self.traceback_at(pos as usize, row_order as usize);
            } else if desc.cell_type == CT_INSERTION {
                let characters_back = desc.order + desc.indel_number;
                for _ in 0..characters_back {
                    trace_array[(window_length - i) as usize] = desc.parent_index;
                    i += 1;

                    if i > window_length {
                        return true;
                    }

                    pos -= 1;
                    if pos < 0 {
                        pos = length - 1;
                    }
                }

                trace_array[(window_length - i) as usize] = row;
                row = desc.parent_index;
            } else if desc.cell_type == CT_DELETION {
                let characters_back = desc.order;
                for _ in 0..characters_back {
                    trace_array[(window_length - i) as usize] = desc.parent_index;
                    i += 1;

                    if i > window_length {
                        return true;
                    }

                    pos -= 1;
                    if pos < 0 {
                        pos = length - 1;
                    }
                }

                trace_array[(window_length - i) as usize] = row;
                row = desc.parent_index;
            }

            pos -= 1;
            if pos < 0 {
                pos = length - 1;
            }
            i += 1;
        }

        true
    }

    /// C++ `RestartMatrix`: makes column 0 the previous column and column 1
    /// the current one, and zeroes column 0 of both matrices. Other columns
    /// and `traceback` are left as they are.
    pub fn restart_matrix(&mut self) {
        self.current_column_index = 1;
        self.previous_column_index = 0;

        let cpc = self.cells_per_column;
        self.score_matrix[..cpc].fill(0.0);
        let tbr = self.traceback_rows();
        self.traceback_matrix[..tbr].fill(0);
    }

    /// C++ `CalculateTraceback`. Picks the best row of column `start_column`
    /// (row 0, or the first match row with a strictly greater score) and walks
    /// back to column 0, writing the row of every column into `traceback`.
    ///
    /// Unlike [`forward_traceback`](Self::forward_traceback), indel cells do
    /// not consult the traceback matrix: the walk continues at the parent
    /// match row. An indel that reaches column 1 or earlier returns early,
    /// leaving the low entries of `traceback` as they were from the previous
    /// call.
    pub fn calculate_traceback(&mut self, start_column: usize) {
        let mut row: usize = 0;
        let start_col = self.score_column(start_column);
        let mut best_value = start_col[0];

        for i in 1..self.cells_per_column {
            if self.cell_descriptions[i].cell_type == CT_MATCH && start_col[i] > best_value {
                row = i;
                best_value = start_col[i];
            }
        }

        let mut i: usize = 0;
        while i <= start_column {
            let desc = self.cell_descriptions[row];
            let row_order = desc.order;

            if desc.cell_type == CT_MATCH || desc.cell_type == CT_BACKGROUND {
                self.traceback[start_column - i] = row as i32;
                // int -> unsigned long long
                row =
                    self.traceback_at(start_column - i, row_order as usize) as i64 as u64 as usize;
            } else if desc.cell_type == CT_INSERTION {
                let characters_back = desc.order + desc.indel_number;
                for _ in 0..characters_back {
                    self.traceback[start_column - i] = desc.parent_index;
                    i += 1;

                    if i >= start_column {
                        return;
                    }
                }

                self.traceback[start_column - i] = row as i32;
                row = desc.parent_index as usize;
            } else if desc.cell_type == CT_DELETION {
                let characters_back = desc.order;
                for _ in 0..characters_back {
                    self.traceback[start_column - i] = desc.parent_index;
                    i += 1;

                    if i >= start_column {
                        return;
                    }
                }

                self.traceback[start_column - i] = row as i32;
                row = desc.parent_index as usize;
            }
            i += 1;
        }
    }

    /// C++ `CreateMatrix`: builds the cell descriptions and allocates and
    /// initialises the matrices.
    ///
    /// Row layout: row 0 is the background cell. Then for each period
    /// `i` in `1..max_period`: one `CT_MATCH` row (order `i`), then
    /// `min(max_insertions, i - 1)` `CT_INSERTION` rows (indel 1, 2, ...),
    /// then `min(max_deletions, i - 1)` `CT_DELETION` rows, all with
    /// `parent_index` = that match row.
    pub fn create_matrix(&mut self) {
        let max_period = self.max_period as u64;
        let max_ins = self.max_insertions as u64;
        let max_del = self.max_deletions as u64;

        // Count the cells. Comparisons are int vs unsigned long long.
        let mut number_of_cells: i32 = 0;
        let mut i: i32 = 0;
        while (i as i64 as u64) < max_period {
            number_of_cells += 1;

            if (i as i64 as u64) <= max_ins {
                if i > 1 {
                    number_of_cells += i - 1;
                }
            } else {
                number_of_cells = (number_of_cells as i64 as u64).wrapping_add(max_ins) as i32;
            }

            if (i as i64 as u64) <= max_del {
                if i > 1 {
                    number_of_cells += i - 1;
                }
            } else {
                number_of_cells = (number_of_cells as i64 as u64).wrapping_add(max_del) as i32;
            }
            i += 1;
        }

        // The C++ allocates number_of_cells + 1 descriptions; the last is
        // never written or read. Here the vector is truncated to the cells
        // actually built.
        let mut desc: Vec<Cell> = Vec::with_capacity(number_of_cells as usize + 1);

        desc.push(Cell {
            cell_type: CT_BACKGROUND,
            parent_index: 0,
            order: 0,
            indel_number: 0,
        });

        let mut i: i32 = 1;
        while (i as i64 as u64) < max_period {
            let match_index = desc.len() as i32;

            desc.push(Cell {
                cell_type: CT_MATCH,
                parent_index: 0,
                order: i,
                indel_number: 0,
            });

            // unsigned long long -> int
            let mut num_ins = max_ins as i32;
            let mut num_del = max_del as i32;

            if num_ins >= i {
                num_ins = i - 1;
            }

            if num_del >= i {
                num_del = i - 1;
            }

            for j in 0..num_ins {
                desc.push(Cell {
                    cell_type: CT_INSERTION,
                    parent_index: match_index,
                    order: i,
                    indel_number: j + 1,
                });
            }

            for j in 0..num_del {
                desc.push(Cell {
                    cell_type: CT_DELETION,
                    parent_index: match_index,
                    order: i,
                    indel_number: j + 1,
                });
            }
            i += 1;
        }

        self.cells_per_column = desc.len();
        self.cell_descriptions = desc;

        // Rust zero-fills where the C++ mallocs uninitialised memory. Every
        // cell is written before it is read on the paths ULTRA uses.
        self.score_matrix = vec![0.0; self.cells_per_column * (self.length + 1)];
        self.traceback_matrix = vec![0; (self.max_period + 1) * (self.length + 1)];

        self.traceback = vec![0; self.length + 2];

        self.current_column_index = 1;
        self.previous_column_index = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_layout_matches_cpp() {
        // maxPeriod 5, maxInsertions 2, maxDeletions 1
        let m = UMatrix::new(5, 2, 1, 8);
        let got: Vec<(u8, i32, i32, i32)> = m
            .cell_descriptions
            .iter()
            .map(|c| (c.cell_type, c.parent_index, c.order, c.indel_number))
            .collect();
        let want = vec![
            (0, 0, 0, 0),
            (1, 0, 1, 0),
            (1, 0, 2, 0),
            (2, 2, 2, 1),
            (3, 2, 2, 1),
            (1, 0, 3, 0),
            (2, 5, 3, 1),
            (2, 5, 3, 2),
            (3, 5, 3, 1),
            (1, 0, 4, 0),
            (2, 9, 4, 1),
            (2, 9, 4, 2),
            (3, 9, 4, 1),
        ];
        assert_eq!(got, want);
        assert_eq!(m.cells_per_column, 13);
        assert_eq!(m.score_matrix.len(), 13 * 9);
        assert_eq!(m.traceback_matrix.len(), 6 * 9);
        assert_eq!(m.traceback.len(), 10);
        assert_eq!((m.current_column_index, m.previous_column_index), (1, 0));
    }

    #[test]
    fn previous_column_wraps() {
        let mut m = UMatrix::new(3, 0, 0, 10);
        m.current_column_index = 2;
        assert_eq!(m.previous_column(1), 1);
        assert_eq!(m.previous_column(3), 9);
        assert_eq!(m.previous_score(0, 10), NEG_INF_F32);
        assert_eq!(m.previous_score(0, -1), NEG_INF_F32);
        assert_eq!(m.previous_traceback(0, 10), 0);
        assert_eq!(NEG_INF_F32.to_bits(), 0xd51184e7);
    }

    #[test]
    fn forward_traceback_rejects_long_windows() {
        let m = UMatrix::new(3, 0, 0, 10);
        let mut t = vec![0; 20];
        assert!(!m.forward_traceback(&mut t, 11, 0));
        assert!(!m.forward_traceback(&mut t, -1, 0));
    }
}
