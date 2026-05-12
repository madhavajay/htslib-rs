//! Pure Rust HTSlib-compatible building blocks.
//!
//! The initial target is a Rust-only API with test parity against the selected
//! HTSlib suite. Format-specific parsing and writing should use noodles by
//! default, with this crate filling the HTSlib-compatible behavior around it.

pub mod alignment_compat;
pub mod bgzf_compat;
pub mod endian;
pub mod error;
pub mod expr;
pub mod faidx_compat;
pub mod fastq_compat;
pub mod format;
pub mod hfile_compat;
pub mod hts_defs;
pub mod hts_os;
pub mod index_compat;
pub mod io;
pub mod kbitset;
pub mod khash;
pub mod klist;
pub mod kroundup;
pub mod ksort;
pub mod kstring;
pub mod log_compat;
pub mod math;
pub mod probaln;
pub mod regidx;
pub mod region;
pub mod sequence;
pub mod tabix_compat;
pub mod text;
pub mod time;
pub mod variant;
pub mod variant_io_compat;

pub use noodles::{bam, bcf, bgzf, core, cram, csi, fasta, fastq, refget, sam, tabix, vcf};
