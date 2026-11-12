//! Rust port of ULTRA (ULTRA Locates Tandemly Repetitive Areas).
//!
//! Upstream: <https://github.com/TravisWheelerLab/ULTRA>. Each module
//! corresponds to one C++ translation unit of the same name.

pub mod cli;
pub mod fasta_reader;
pub mod file_reader;
pub mod mask;
pub mod repeat;
pub mod repeat_splitter;
pub mod sequence_window;
pub mod std_sort;
pub mod symbol;
pub mod ultra;
pub mod umatrix;
pub mod umodel;
pub mod umodel_driver;
pub mod writers;
