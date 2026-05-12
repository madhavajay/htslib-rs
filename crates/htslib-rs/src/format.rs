//! Format classification and local file detection.

use std::{
    fmt, fs,
    io::{self, Read},
    path::Path,
};

use flate2::read::MultiGzDecoder;

use crate::error::{ContextError, ResultExt as _};

/// High-level HTSlib format categories.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Category {
    /// The category could not be determined.
    Unknown,
    /// Sequence alignment data, e.g. SAM, BAM, or CRAM.
    SequenceData,
    /// Variant data, e.g. VCF or BCF.
    VariantData,
    /// An index associated with a data file.
    IndexFile,
    /// Coordinate intervals or regions.
    RegionList,
}

/// Exact file formats recognized by the initial detector.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Exact {
    /// The exact format could not be determined.
    Unknown,
    /// An empty file.
    Empty,
    /// SAM.
    Sam,
    /// BAM.
    Bam,
    /// BAI.
    Bai,
    /// CRAM.
    Cram,
    /// CRAI.
    Crai,
    /// VCF.
    Vcf,
    /// BCF.
    Bcf,
    /// CSI.
    Csi,
    /// GZI.
    Gzi,
    /// TBI.
    Tbi,
    /// BED.
    Bed,
    /// FASTA.
    Fasta,
    /// FASTQ.
    Fastq,
    /// FAI.
    Fai,
    /// FQI.
    Fqi,
}

/// Compression recognized by the initial detector.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    /// The stream is not compressed.
    None,
    /// Gzip-compatible compression.
    Gzip,
    /// Blocked gzip compression.
    Bgzf,
}

/// A detected format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Format {
    /// The high-level category.
    pub category: Category,
    /// The exact format.
    pub exact: Exact,
    /// The detected compression.
    pub compression: Compression,
}

impl Format {
    /// Creates a detected format.
    pub const fn new(category: Category, exact: Exact, compression: Compression) -> Self {
        Self {
            category,
            exact,
            compression,
        }
    }
}

impl Default for Format {
    fn default() -> Self {
        Self::new(Category::Unknown, Exact::Unknown, Compression::None)
    }
}

/// Errors returned while detecting a file format.
#[derive(Debug)]
pub enum DetectError {
    /// An I/O error occurred.
    Io(io::Error),
    /// An error with operation or path context occurred.
    Context(ContextError),
}

impl fmt::Display for DetectError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "failed to read input: {e}"),
            Self::Context(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for DetectError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Context(e) => Some(e),
        }
    }
}

impl From<io::Error> for DetectError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<ContextError> for DetectError {
    fn from(e: ContextError) -> Self {
        Self::Context(e)
    }
}

/// Detects the format of a local file.
pub fn detect_path<P>(src: P) -> Result<Format, DetectError>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let mut file =
        fs::File::open(path).with_path_context("open input for format detection", path)?;
    let mut buf = [0; 65536];
    let n = file
        .read(&mut buf)
        .with_path_context("read input for format detection", path)?;

    let mut format = if n == 0 {
        Format::new(Category::Unknown, Exact::Empty, Compression::None)
    } else {
        detect_bytes(&buf[..n])
    };

    if matches!(format.compression, Compression::Bgzf | Compression::Gzip) {
        let inner_format = detect_compressed_path(path, format.compression)
            .with_path_context("inspect compressed input for format detection", path)?;

        if inner_format.exact != Exact::Unknown {
            format = inner_format;
        }
    }

    Ok(format)
}

/// Detects a format from the start of a byte stream.
pub fn detect_bytes(buf: &[u8]) -> Format {
    if buf.is_empty() {
        return Format::new(Category::Unknown, Exact::Empty, Compression::None);
    }

    if starts_with(buf, b"BAM\x01") {
        return Format::new(Category::SequenceData, Exact::Bam, Compression::None);
    }

    if starts_with(buf, b"CRAM") {
        return Format::new(Category::SequenceData, Exact::Cram, Compression::None);
    }

    if starts_with(buf, b"BCF\x02") {
        return Format::new(Category::VariantData, Exact::Bcf, Compression::None);
    }

    if starts_with(buf, b"BAI\x01") {
        return Format::new(Category::IndexFile, Exact::Bai, Compression::None);
    }

    if starts_with(buf, b"CSI\x01") {
        return Format::new(Category::IndexFile, Exact::Csi, Compression::None);
    }

    if starts_with(buf, b"TBI\x01") {
        return Format::new(Category::IndexFile, Exact::Tbi, Compression::None);
    }

    if is_bgzf(buf) {
        return Format::new(Category::Unknown, Exact::Unknown, Compression::Bgzf);
    }

    if starts_with(buf, &[0x1f, 0x8b]) {
        return Format::new(Category::Unknown, Exact::Unknown, Compression::Gzip);
    }

    detect_text(buf)
}

fn detect_text(buf: &[u8]) -> Format {
    let first = first_non_empty_line(buf);

    match first {
        Some(line) if line.starts_with(b"@HD") || line.starts_with(b"@SQ") => {
            Format::new(Category::SequenceData, Exact::Sam, Compression::None)
        }
        Some(line) if line.starts_with(b"##fileformat=VCF") || line.starts_with(b"#CHROM") => {
            Format::new(Category::VariantData, Exact::Vcf, Compression::None)
        }
        Some(line) if line.starts_with(b">") => {
            Format::new(Category::SequenceData, Exact::Fasta, Compression::None)
        }
        Some(line) if line.starts_with(b"@") && looks_like_fastq(buf) => {
            Format::new(Category::SequenceData, Exact::Fastq, Compression::None)
        }
        _ if has_bed_record(buf) => {
            Format::new(Category::RegionList, Exact::Bed, Compression::None)
        }
        _ => Format::default(),
    }
}

fn starts_with(buf: &[u8], magic: &[u8]) -> bool {
    buf.len() >= magic.len() && &buf[..magic.len()] == magic
}

fn first_non_empty_line(buf: &[u8]) -> Option<&[u8]> {
    buf.split(|&b| b == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .find(|line| !line.is_empty())
}

fn is_bgzf(buf: &[u8]) -> bool {
    if buf.len() < 18 || !starts_with(buf, &[0x1f, 0x8b]) {
        return false;
    }

    let flags = buf[3];

    if flags & 0x04 == 0 {
        return false;
    }

    let xlen = u16::from_le_bytes([buf[10], buf[11]]) as usize;
    let mut i = 12;
    let end = 12 + xlen;

    if end > buf.len() {
        return false;
    }

    while i + 4 <= end {
        let si1 = buf[i];
        let si2 = buf[i + 1];
        let slen = u16::from_le_bytes([buf[i + 2], buf[i + 3]]) as usize;
        i += 4;

        if i + slen > end {
            return false;
        }

        if si1 == b'B' && si2 == b'C' && slen == 2 {
            return true;
        }

        i += slen;
    }

    false
}

fn looks_like_fastq(buf: &[u8]) -> bool {
    let mut lines = buf
        .split(|&b| b == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line));

    matches!(
        (lines.next(), lines.next(), lines.next(), lines.next()),
        (Some(name), Some(seq), Some(plus), Some(qual))
            if name.starts_with(b"@")
                && !seq.is_empty()
                && plus.starts_with(b"+")
                && qual.len() >= seq.len()
    )
}

fn looks_like_bed_line(line: &[u8]) -> bool {
    if line.starts_with(b"#")
        || line.starts_with(b"track")
        || line.starts_with(b"browser")
        || line.starts_with(b"@")
    {
        return false;
    }

    let mut fields = line.split(|&b| b == b'\t');
    matches!(
        (fields.next(), fields.next(), fields.next()),
        (Some(name), Some(start), Some(end))
            if !name.is_empty() && is_ascii_u64(start) && is_ascii_u64(end)
    )
}

fn is_ascii_u64(s: &[u8]) -> bool {
    !s.is_empty() && s.iter().all(u8::is_ascii_digit)
}

fn has_bed_record(buf: &[u8]) -> bool {
    buf.split(|&b| b == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .any(looks_like_bed_line)
}

fn detect_compressed_path(path: &Path, compression: Compression) -> Result<Format, DetectError> {
    let file = fs::File::open(path)?;
    let mut reader = MultiGzDecoder::new(file);
    let mut buf = [0; 65536];
    let n = reader.read(&mut buf)?;

    let mut format = if n == 0 {
        Format::new(Category::Unknown, Exact::Empty, compression)
    } else {
        detect_bytes(&buf[..n])
    };

    format.compression = compression;

    Ok(format)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_empty() {
        assert_eq!(
            detect_bytes(b""),
            Format::new(Category::Unknown, Exact::Empty, Compression::None)
        );
    }

    #[test]
    fn detects_text_formats() {
        assert_eq!(
            detect_bytes(b"@HD\tVN:1.6\n"),
            Format::new(Category::SequenceData, Exact::Sam, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"##fileformat=VCFv4.3\n"),
            Format::new(Category::VariantData, Exact::Vcf, Compression::None)
        );
        assert_eq!(
            detect_bytes(b">sq0\nACGT\n"),
            Format::new(Category::SequenceData, Exact::Fasta, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"@r0\nACGT\n+\n!!!!\n"),
            Format::new(Category::SequenceData, Exact::Fastq, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"sq0\t0\t4\n"),
            Format::new(Category::RegionList, Exact::Bed, Compression::None)
        );
    }

    #[test]
    fn detects_binary_magic_numbers() {
        assert_eq!(
            detect_bytes(b"BAM\x01"),
            Format::new(Category::SequenceData, Exact::Bam, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"CRAM"),
            Format::new(Category::SequenceData, Exact::Cram, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"BCF\x02\x02"),
            Format::new(Category::VariantData, Exact::Bcf, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"BAI\x01"),
            Format::new(Category::IndexFile, Exact::Bai, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"CSI\x01"),
            Format::new(Category::IndexFile, Exact::Csi, Compression::None)
        );
        assert_eq!(
            detect_bytes(b"TBI\x01"),
            Format::new(Category::IndexFile, Exact::Tbi, Compression::None)
        );
    }
}
