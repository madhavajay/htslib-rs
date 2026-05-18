//! HTSlib-compatible BGZF helpers built on noodles BGZF primitives.

use std::{
    fs::File,
    io::{self, BufReader, Read, Seek, Write},
    num::NonZero,
    path::Path,
};

use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use noodles::{bam, bgzf};

/// A BGZF/GZIP uncompressed offset index.
pub type GziIndex = bgzf::gzi::Index;

/// Number of bits used for the uncompressed block offset in a BGZF virtual offset.
pub const UNCOMPRESSED_OFFSET_BITS: u64 = 16;

/// Mask for the uncompressed block offset part of a BGZF virtual offset.
pub const UNCOMPRESSED_OFFSET_MASK: u64 = 0xffff;

/// Maximum compressed file offset representable in a BGZF virtual offset.
pub const MAX_COMPRESSED_OFFSET: u64 = (1 << 48) - 1;

/// The canonical empty BGZF block used as an end-of-file marker.
pub const EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Compression kind used by HTSlib BGZF open modes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CompressionKind {
    /// Plain uncompressed bytes (`wu`).
    Uncompressed,
    /// Standard gzip stream (`wg`).
    Gzip,
    /// BGZF stream (`w`/`w0`..`w9`).
    Bgzf,
}

/// Error returned when creating a BGZF virtual offset fails.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VirtualOffsetError {
    /// The compressed file offset is larger than 48 bits.
    CompressedOffsetOverflow,
}

impl std::fmt::Display for VirtualOffsetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CompressedOffsetOverflow => {
                f.write_str("the compressed offset is larger than 2^48 - 1")
            }
        }
    }
}

impl std::error::Error for VirtualOffsetError {}

/// Creates an HTSlib-style BGZF virtual offset.
///
/// This matches the `bgzf_tell` layout: the compressed block address occupies
/// the high 48 bits, and the uncompressed in-block offset occupies the low 16
/// bits.
pub const fn virtual_offset(
    compressed_offset: u64,
    uncompressed_offset: u16,
) -> Result<u64, VirtualOffsetError> {
    if compressed_offset > MAX_COMPRESSED_OFFSET {
        return Err(VirtualOffsetError::CompressedOffsetOverflow);
    }

    Ok((compressed_offset << UNCOMPRESSED_OFFSET_BITS) | uncompressed_offset as u64)
}

/// Splits an HTSlib-style BGZF virtual offset into `(compressed, uncompressed)` parts.
pub const fn virtual_offset_parts(offset: u64) -> (u64, u16) {
    (
        offset >> UNCOMPRESSED_OFFSET_BITS,
        (offset & UNCOMPRESSED_OFFSET_MASK) as u16,
    )
}

/// Converts an HTSlib-style BGZF virtual offset to a noodles virtual position.
pub fn virtual_position_from_offset(offset: u64) -> bgzf::VirtualPosition {
    bgzf::VirtualPosition::from(offset)
}

/// Converts a noodles virtual position to an HTSlib-style BGZF virtual offset.
pub fn virtual_offset_from_position(position: bgzf::VirtualPosition) -> u64 {
    position.into()
}

/// Returns whether a byte stream ends with the canonical BGZF EOF marker.
pub fn has_eof_marker(data: &[u8]) -> bool {
    data.ends_with(&EOF_MARKER)
}

/// Reads all uncompressed bytes from a BGZF stream.
pub fn read_all<R>(reader: R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut reader = bgzf::io::Reader::new(reader);
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    Ok(buf)
}

/// Reads all uncompressed bytes from a BGZF stream using noodles' worker-count API.
///
/// This is the Rust equivalent of HTSlib's BGZF thread-pool read path for the
/// current Rust-only target: callers select a nonzero worker count, and the
/// implementation delegates decompression scheduling to `noodles-bgzf`.
pub fn read_all_with_worker_count<R>(reader: R, worker_count: NonZero<usize>) -> io::Result<Vec<u8>>
where
    R: Read + Send + 'static,
{
    let mut reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, reader);
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    Ok(buf)
}

/// Detects whether a local byte stream is plain, gzip, or BGZF-compressed.
pub fn detect_compression_kind(data: &[u8]) -> CompressionKind {
    const GZIP_ID1: u8 = 0x1f;
    const GZIP_ID2: u8 = 0x8b;
    const FLG_FEXTRA: u8 = 4;
    const BGZF_EXTRA_ID: &[u8] = b"BC";

    if data.len() >= 18
        && data[0] == GZIP_ID1
        && data[1] == GZIP_ID2
        && data[3] & FLG_FEXTRA != 0
        && &data[12..14] == BGZF_EXTRA_ID
    {
        CompressionKind::Bgzf
    } else if data.len() >= 2 && data[0] == GZIP_ID1 && data[1] == GZIP_ID2 {
        CompressionKind::Gzip
    } else {
        CompressionKind::Uncompressed
    }
}

/// Reads all bytes from a plain, gzip, or BGZF local byte stream.
pub fn read_auto(data: &[u8]) -> io::Result<Vec<u8>> {
    read_auto_with_worker_count(data, None)
}

/// Reads all bytes from a plain, gzip, or BGZF local byte stream.
///
/// When `worker_count` is set and the stream is BGZF-compressed, decompression
/// is delegated to noodles' multithreaded BGZF reader. Plain and gzip streams
/// are unaffected by the worker count.
pub fn read_auto_with_worker_count(
    data: &[u8],
    worker_count: Option<NonZero<usize>>,
) -> io::Result<Vec<u8>> {
    match detect_compression_kind(data) {
        CompressionKind::Uncompressed => Ok(data.to_vec()),
        CompressionKind::Gzip => {
            let mut reader = GzDecoder::new(data);
            let mut buf = Vec::new();
            reader.read_to_end(&mut buf)?;
            Ok(buf)
        }
        CompressionKind::Bgzf => match worker_count {
            Some(worker_count) => {
                read_all_with_worker_count(io::Cursor::new(data.to_vec()), worker_count)
            }
            None => read_all(data),
        },
    }
}

/// Writes all bytes to a BGZF stream and returns the finished inner writer.
pub fn write_all<W>(writer: W, data: &[u8]) -> io::Result<W>
where
    W: Write,
{
    let mut writer = bgzf::io::Writer::new(writer);
    writer.write_all(data)?;
    writer.finish()
}

/// Writes all bytes to a BGZF stream using noodles' worker-count API.
///
/// This is the Rust equivalent of HTSlib's BGZF thread-pool write path for the
/// current Rust-only target. Passing `NonZero<usize>` makes the unsupported
/// zero-worker case unrepresentable.
pub fn write_all_with_worker_count<W>(
    writer: W,
    data: &[u8],
    worker_count: NonZero<usize>,
) -> io::Result<W>
where
    W: Write + Send + 'static,
{
    let mut writer = bgzf::io::MultithreadedWriter::with_worker_count(worker_count, writer);
    writer.write_all(data)?;
    writer.finish()
}

/// Writes bytes using an HTSlib BGZF open-mode compression kind.
pub fn write_all_with_kind(data: &[u8], kind: CompressionKind) -> io::Result<Vec<u8>> {
    write_all_with_kind_and_worker_count(data, kind, None)
}

/// Writes bytes using an HTSlib BGZF open-mode compression kind.
///
/// When `worker_count` is set and `kind` is BGZF, compression is delegated to
/// noodles' multithreaded BGZF writer. Plain and gzip modes are unaffected.
pub fn write_all_with_kind_and_worker_count(
    data: &[u8],
    kind: CompressionKind,
    worker_count: Option<NonZero<usize>>,
) -> io::Result<Vec<u8>> {
    match kind {
        CompressionKind::Uncompressed => Ok(data.to_vec()),
        CompressionKind::Gzip => {
            let mut writer = GzEncoder::new(Vec::new(), Compression::default());
            writer.write_all(data)?;
            writer.finish()
        }
        CompressionKind::Bgzf => match worker_count {
            Some(worker_count) => write_all_with_worker_count(Vec::new(), data, worker_count),
            None => write_all(Vec::new(), data),
        },
    }
}

/// Reads a GZI index from a reader.
pub fn read_gzi<R>(reader: R) -> io::Result<GziIndex>
where
    R: Read,
{
    let mut reader = bgzf::gzi::io::Reader::new(reader);

    reader.read_index()
}

/// Reads a GZI index from a local file.
pub fn read_gzi_from_path<P>(src: P) -> io::Result<GziIndex>
where
    P: AsRef<Path>,
{
    bgzf::gzi::fs::read(src)
}

/// Writes a GZI index to a writer.
pub fn write_gzi<W>(writer: W, index: &GziIndex) -> io::Result<W>
where
    W: Write,
{
    let mut writer = bgzf::gzi::io::Writer::new(writer);

    writer.write_index(index)?;

    Ok(writer.into_inner())
}

/// Writes a GZI index to a local file.
pub fn write_gzi_to_path<P>(dst: P, index: &GziIndex) -> io::Result<()>
where
    P: AsRef<Path>,
{
    bgzf::gzi::fs::write(dst, index)
}

/// Builds a GZI index from a BGZF stream.
pub fn build_gzi<R>(reader: &mut R) -> io::Result<GziIndex>
where
    R: Read + Seek,
{
    let mut compressed_offset = reader.stream_position()?;
    let mut uncompressed_offset = 0;
    let mut offsets = Vec::new();

    loop {
        let mut header = [0; 18];

        match reader.read_exact(&mut header) {
            Ok(()) => {}
            Err(e)
                if e.kind() == io::ErrorKind::UnexpectedEof && header.iter().all(|b| *b == 0) =>
            {
                break;
            }
            Err(e) => return Err(e),
        }

        validate_bgzf_header(&header)?;

        let block_size = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;

        if block_size < 26 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BGZF block size",
            ));
        }

        let tail_len = usize::try_from(block_size - 18)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut tail = vec![0; tail_len];
        reader.read_exact(&mut tail)?;

        let isize_offset = tail
            .len()
            .checked_sub(4)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing BGZF ISIZE"))?;
        let isize = u32::from_le_bytes([
            tail[isize_offset],
            tail[isize_offset + 1],
            tail[isize_offset + 2],
            tail[isize_offset + 3],
        ]);

        if isize > 0 {
            if compressed_offset != 0 {
                offsets.push((compressed_offset, uncompressed_offset));
            }

            uncompressed_offset += u64::from(isize);
        }

        compressed_offset = reader.stream_position()?;
    }

    Ok(GziIndex::from(offsets))
}

/// Resolves an uncompressed offset to a BGZF virtual position using a GZI index.
pub fn query_gzi(index: &GziIndex, uncompressed_offset: u64) -> io::Result<bgzf::VirtualPosition> {
    index.query(uncompressed_offset)
}

fn validate_bgzf_header(header: &[u8; 18]) -> io::Result<()> {
    const GZIP_ID1: u8 = 0x1f;
    const GZIP_ID2: u8 = 0x8b;
    const CM_DEFLATE: u8 = 8;
    const FLG_FEXTRA: u8 = 4;
    const SI1: u8 = b'B';
    const SI2: u8 = b'C';
    const SLEN: [u8; 2] = 2u16.to_le_bytes();

    let is_bgzf = header[0] == GZIP_ID1
        && header[1] == GZIP_ID2
        && header[2] == CM_DEFLATE
        && header[3] & FLG_FEXTRA != 0
        && header[12] == SI1
        && header[13] == SI2
        && header[14..16] == SLEN;

    if is_bgzf {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ))
    }
}

/// Result of attempting the block-level BAM concatenation fast path
/// for one input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BamFrameCopy {
    /// The input's alignment-section BGZF frames were copied verbatim.
    Copied,
    /// The input's BAM header does not end on a BGZF block boundary, so
    /// nothing was written; the caller must fall back to a
    /// record-level re-encode for this input.
    NotBlockAligned,
}

/// Appends the alignment-section BGZF frames of the BAM file at `src`
/// to `out` **verbatim** (no re-deflate), iff the input's BAM header
/// ends exactly on a BGZF block boundary.
///
/// This is `samtools cat`'s BAM fast path: the bulk of every input is
/// copied as already-compressed BGZF frames instead of being decoded
/// and re-encoded. The input's own trailing BGZF EOF block is dropped;
/// the caller is responsible for writing a single
/// [`noodles::bgzf::io::EOF`] block after the last input and for
/// emitting the (shared) BAM header before the first input.
///
/// Returns [`BamFrameCopy::NotBlockAligned`] **without writing any
/// bytes** when the header shares its final BGZF block with alignment
/// data (so a verbatim copy would corrupt the stream). The caller must
/// then fall back to a correct record-level copy for that input. This
/// keeps the optimization strictly safe: it only ever fast-paths
/// inputs where the boundary is provable.
pub fn append_bam_alignment_frames<W>(src: &Path, out: &mut W) -> io::Result<BamFrameCopy>
where
    W: Write,
{
    use noodles::bgzf::io::{EOF, read_raw_frame};

    // Decode just the BAM header to learn where the alignment section
    // begins. `virtual_position()` after `read_header` is
    // (compressed = file offset of the current BGZF block,
    //  uncompressed = byte offset within that decoded block).
    let mut header_reader = bam::io::Reader::new(File::open(src)?);
    header_reader.read_header()?;
    let position = header_reader.get_ref().virtual_position();
    if position.uncompressed() != 0 {
        // Header ends mid-block (shares a block with alignments); a
        // raw frame copy is impossible without re-deflating that
        // block. Signal the caller to fall back.
        return Ok(BamFrameCopy::NotBlockAligned);
    }
    let alignment_start = position.compressed();

    // Re-open and copy raw frames starting at the alignment block
    // boundary. A one-frame lookahead lets us drop the input's own
    // trailing EOF block while preserving any (pathological) empty
    // block that is not the terminator.
    let mut raw = BufReader::new(File::open(src)?);
    raw.seek(io::SeekFrom::Start(alignment_start))?;

    let mut pending: Option<Vec<u8>> = None;
    let mut buf = Vec::new();
    while read_raw_frame(&mut raw, &mut buf)? {
        if let Some(frame) = pending.take() {
            out.write_all(&frame)?;
        }
        pending = Some(std::mem::take(&mut buf));
    }
    if let Some(last) = pending
        && last.as_slice() != EOF
    {
        out.write_all(&last)?;
    }

    Ok(BamFrameCopy::Copied)
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::num::NonZero;

    use super::{
        CompressionKind, MAX_COMPRESSED_OFFSET, VirtualOffsetError, build_gzi,
        detect_compression_kind, query_gzi, read_all, read_auto, read_auto_with_worker_count,
        read_gzi, virtual_offset, virtual_offset_from_position, virtual_offset_parts,
        virtual_position_from_offset, write_all, write_all_with_kind,
        write_all_with_kind_and_worker_count, write_gzi,
    };

    #[test]
    fn test_virtual_offset_round_trip() {
        let offset = virtual_offset(57, 6086).unwrap();

        assert_eq!(offset, 3741638);
        assert_eq!(virtual_offset_parts(offset), (57, 6086));
    }

    #[test]
    fn test_virtual_offset_rejects_compressed_overflow() {
        assert_eq!(
            virtual_offset(MAX_COMPRESSED_OFFSET + 1, 0),
            Err(VirtualOffsetError::CompressedOffsetOverflow)
        );
    }

    #[test]
    fn test_noodles_virtual_position_conversion() {
        let offset = virtual_offset(399103671, 321).unwrap();
        let position = virtual_position_from_offset(offset);

        assert_eq!(virtual_offset_from_position(position), offset);
        assert_eq!(position.compressed(), 399103671);
        assert_eq!(position.uncompressed(), 321);
    }

    #[test]
    fn test_read_and_write_all() {
        const PLAIN: &[u8] = include_bytes!("../../../htslib/test/bgziptest.txt");
        const BGZF: &[u8] = include_bytes!("../../../htslib/test/bgziptest.txt.gz");

        assert_eq!(read_all(Cursor::new(BGZF)).unwrap(), PLAIN);

        let encoded = write_all(Vec::new(), PLAIN).unwrap();
        assert_eq!(read_all(Cursor::new(encoded)).unwrap(), PLAIN);
    }

    #[test]
    fn test_htslib_compression_modes() {
        const PLAIN: &[u8] = include_bytes!("../../../htslib/test/bgziptest.txt");

        for kind in [
            CompressionKind::Uncompressed,
            CompressionKind::Gzip,
            CompressionKind::Bgzf,
        ] {
            let encoded = write_all_with_kind(PLAIN, kind).unwrap();

            assert_eq!(detect_compression_kind(&encoded), kind);
            assert_eq!(read_auto(&encoded).unwrap(), PLAIN);
        }
    }

    #[test]
    fn test_htslib_bgzf_worker_count_paths() {
        const PLAIN: &[u8] = include_bytes!("../../../htslib/test/bgziptest.txt");
        let worker_count = NonZero::new(2).unwrap();

        let encoded =
            write_all_with_kind_and_worker_count(PLAIN, CompressionKind::Bgzf, Some(worker_count))
                .unwrap();

        assert_eq!(detect_compression_kind(&encoded), CompressionKind::Bgzf);
        assert_eq!(
            read_auto_with_worker_count(&encoded, Some(worker_count)).unwrap(),
            PLAIN
        );
    }

    #[test]
    fn test_read_write_and_query_gzi() {
        const GZI: &[u8] = include_bytes!("../../../htslib/test/bgziptest.txt.gz.gzi");

        let index = read_gzi(Cursor::new(GZI)).unwrap();
        assert_eq!(index.as_ref().len(), 5);

        let position = query_gzi(&index, 1000).unwrap();
        assert!(position.compressed() > 0);

        let encoded = write_gzi(Vec::new(), &index).unwrap();
        assert_eq!(encoded, GZI);
    }

    #[test]
    fn test_build_gzi() {
        let plain = vec![b'N'; 200_000];

        let encoded = write_all(Vec::new(), &plain).unwrap();
        let index = build_gzi(&mut Cursor::new(encoded)).unwrap();
        let position = query_gzi(&index, 70_000).unwrap();

        assert!(position.compressed() > 0);
    }

    #[test]
    fn append_bam_alignment_frames_reconstructs_records() {
        use std::path::PathBuf;

        use noodles::bam;
        use noodles::bgzf::io::EOF;

        use super::{BamFrameCopy, append_bam_alignment_frames};

        let src = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../htslib/test/range.bam");

        // Reference: header + all records decoded normally.
        fn names<R: std::io::Read>(
            reader: &mut bam::io::Reader<noodles::bgzf::io::Reader<R>>,
        ) -> Vec<Option<Vec<u8>>> {
            reader
                .records()
                .map(|r| r.unwrap().name().map(|n| n.to_vec()))
                .collect()
        }

        let mut reader = bam::io::Reader::new(std::fs::File::open(&src).unwrap());
        let header = reader.read_header().unwrap();
        let original = names(&mut reader);
        assert!(!original.is_empty());

        // The header's own (block-aligned) BGZF frames, taken verbatim
        // from the source file up to the alignment boundary.
        let mut probe = bam::io::Reader::new(std::fs::File::open(&src).unwrap());
        probe.read_header().unwrap();
        let boundary = probe.get_ref().virtual_position().compressed() as usize;
        let src_bytes = std::fs::read(&src).unwrap();
        let mut rebuilt = src_bytes[..boundary].to_vec();

        // Append the alignment frames twice (cat of the file with
        // itself) then a single shared EOF.
        for _ in 0..2 {
            let mut sink = Vec::new();
            assert_eq!(
                append_bam_alignment_frames(&src, &mut sink).unwrap(),
                BamFrameCopy::Copied,
                "range.bam is samtools-written and must be block-aligned"
            );
            rebuilt.extend_from_slice(&sink);
        }
        rebuilt.extend_from_slice(&EOF);

        // The reconstructed stream must decode to exactly two copies of
        // the original records, with the same header.
        let mut out = bam::io::Reader::new(Cursor::new(rebuilt));
        let rebuilt_header = out.read_header().unwrap();
        assert_eq!(
            rebuilt_header.reference_sequences().len(),
            header.reference_sequences().len()
        );
        let rebuilt_records = names(&mut out);

        let mut expected = original.clone();
        expected.extend(original.clone());
        assert_eq!(rebuilt_records, expected);
    }
}
