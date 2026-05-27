//! HTSlib-compatible index helpers backed by noodles index readers and builders.

use std::{
    ffi::OsString,
    fs::File,
    io::{self, BufRead, Read, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use crate::{bam, bcf, bgzf, core::Position, cram, csi, sam, tabix, tabix_compat, vcf};

/// A BAM BAI index.
pub type BaiIndex = bam::bai::Index;

/// A coordinate-sorted CSI index.
pub type CsiIndex = csi::Index;

/// A tabix TBI index.
pub type TbiIndex = tabix::Index;

/// A CRAM CRAI index.
pub type CraiIndex = cram::crai::Index;

/// The HTSlib explicit-index delimiter.
pub const INDEX_DELIMITER: &str = "##idx##";

/// HTSlib index formats used for associated-index lookup.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndexFormat {
    /// CSI index.
    Csi,
    /// BAI index.
    Bai,
    /// TBI index.
    Tbi,
    /// CRAI index.
    Crai,
}

impl IndexFormat {
    fn extension(self) -> &'static str {
        match self {
            Self::Csi => ".csi",
            Self::Bai => ".bai",
            Self::Tbi => ".tbi",
            Self::Crai => ".crai",
        }
    }
}

/// A located associated index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LocatedIndex {
    /// Index path.
    pub path: PathBuf,
    /// Index format inferred from the selected path.
    pub format: IndexFormat,
}

/// Reads a BAI index from a local file.
pub fn read_bai<P>(src: P) -> io::Result<BaiIndex>
where
    P: AsRef<Path>,
{
    bam::bai::fs::read(src)
}

/// Reads a CSI index from a local file.
pub fn read_csi<P>(src: P) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    csi::fs::read(src)
}

/// Reads a TBI index from a local file.
pub fn read_tbi<P>(src: P) -> io::Result<TbiIndex>
where
    P: AsRef<Path>,
{
    tabix::fs::read(src)
}

/// Reads a CRAI index from a local file.
pub fn read_cram_crai<P>(src: P) -> io::Result<CraiIndex>
where
    P: AsRef<Path>,
{
    cram::crai::fs::read(src)
}

/// Returns the HTSlib local associated-index candidate paths in lookup order.
pub fn associated_index_candidates<P>(src: P, format: IndexFormat) -> Vec<LocatedIndex>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();

    if let Some(index_path) = explicit_index_path(src) {
        return vec![LocatedIndex {
            format: infer_index_format(&index_path).unwrap_or(format),
            path: index_path,
        }];
    }

    let mut candidates = Vec::new();
    push_index_candidates(&mut candidates, src, IndexFormat::Csi);

    match format {
        IndexFormat::Csi => {}
        IndexFormat::Bai => push_index_candidates(&mut candidates, src, IndexFormat::Bai),
        IndexFormat::Tbi => push_index_candidates(&mut candidates, src, IndexFormat::Tbi),
        IndexFormat::Crai => push_index_candidates(&mut candidates, src, IndexFormat::Crai),
    }

    candidates
}

/// Locates the first existing local associated index using HTSlib lookup order.
pub fn locate_associated_index<P>(src: P, format: IndexFormat) -> Option<LocatedIndex>
where
    P: AsRef<Path>,
{
    associated_index_candidates(src, format)
        .into_iter()
        .find(|candidate| candidate.path.exists())
}

/// Returns the data-file path with any HTSlib explicit-index suffix removed.
pub fn associated_data_path<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    let src = src.as_ref();
    let src = src.to_string_lossy();
    let (data, _) = src.split_once(INDEX_DELIMITER).unwrap_or((&src, ""));

    PathBuf::from(data)
}

/// Reads a BAM-associated BAI or CSI index using HTSlib lookup order.
pub fn read_associated_bam_index<P>(src: P) -> io::Result<Box<dyn csi::BinningIndex>>
where
    P: AsRef<Path>,
{
    let located = locate_associated_index(&src, IndexFormat::Bai)
        .ok_or_else(|| missing_associated_index_error(src.as_ref(), IndexFormat::Bai))?;

    match located.format {
        IndexFormat::Csi => read_csi(located.path).map(|index| Box::new(index) as _),
        IndexFormat::Bai => read_bai(located.path).map(|index| Box::new(index) as _),
        IndexFormat::Tbi | IndexFormat::Crai => Err(unsupported_index_format_error(located)),
    }
}

/// Reads a VCF-associated TBI or CSI index using HTSlib lookup order.
pub fn read_associated_vcf_index<P>(src: P) -> io::Result<Box<dyn csi::BinningIndex>>
where
    P: AsRef<Path>,
{
    let located = locate_associated_index(&src, IndexFormat::Tbi)
        .ok_or_else(|| missing_associated_index_error(src.as_ref(), IndexFormat::Tbi))?;

    match located.format {
        IndexFormat::Csi => read_csi(located.path).map(|index| Box::new(index) as _),
        IndexFormat::Tbi => read_tbi(located.path).map(|index| Box::new(index) as _),
        IndexFormat::Bai | IndexFormat::Crai => Err(unsupported_index_format_error(located)),
    }
}

/// Reads a BCF-associated CSI index using HTSlib lookup order.
pub fn read_associated_bcf_index<P>(src: P) -> io::Result<Box<dyn csi::BinningIndex>>
where
    P: AsRef<Path>,
{
    let located = locate_associated_index(&src, IndexFormat::Csi)
        .ok_or_else(|| missing_associated_index_error(src.as_ref(), IndexFormat::Csi))?;

    match located.format {
        IndexFormat::Csi => read_csi(located.path).map(|index| Box::new(index) as _),
        IndexFormat::Bai | IndexFormat::Tbi | IndexFormat::Crai => {
            Err(unsupported_index_format_error(located))
        }
    }
}

/// Reads a CRAM-associated CRAI index using HTSlib lookup order.
pub fn read_associated_cram_index<P>(src: P) -> io::Result<CraiIndex>
where
    P: AsRef<Path>,
{
    let located = locate_associated_index(&src, IndexFormat::Crai)
        .ok_or_else(|| missing_associated_index_error(src.as_ref(), IndexFormat::Crai))?;

    match located.format {
        IndexFormat::Crai => read_cram_crai(located.path),
        IndexFormat::Csi | IndexFormat::Bai | IndexFormat::Tbi => {
            Err(unsupported_index_format_error(located))
        }
    }
}

/// Builds a BAI index for a coordinate-sorted BAM file.
///
/// Unlike `noodles_bam::fs::index`, this does **not** require the SAM header
/// to carry `@HD SO:coordinate`: upstream `samtools index` indexes
/// coordinate-ordered BAMs whose header omits the sort-order tag (e.g.
/// `test/dat/test_input_1_{a,b}.bam`), so the data order — not the header
/// annotation — is authoritative. The record loop mirrors the CSI/SAM-BAI
/// builders, which already index without the header check.
pub fn build_bai<P>(src: P) -> io::Result<BaiIndex>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    build_bai_from_bam_reader(&mut reader)
}

/// Builds a BAI index using a multithreaded BGZF reader.
pub fn build_bai_with_worker_count<P>(src: P, worker_count: NonZero<usize>) -> io::Result<BaiIndex>
where
    P: AsRef<Path>,
{
    let reader = File::open(src)?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, reader);
    let mut reader = bam::io::Reader::from(decoder);
    build_bai_from_bam_reader(&mut reader)
}

fn build_bai_from_bam_reader<R>(reader: &mut bam::io::Reader<R>) -> io::Result<BaiIndex>
where
    R: Read + bgzf::io::Read,
{
    let header = reader.read_header()?;
    let mut indexer = alignment_bai_indexer();
    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_position,
            end_position,
        );
        let alignment_context = bam_alignment_context(&record)?;

        indexer.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

/// Builds a CRAI index for a coordinate-sorted CRAM file.
pub fn build_cram_crai<P>(src: P) -> io::Result<CraiIndex>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();
    match cram_reference_repository_from_header_uri(src)? {
        Some(repository) => cram::fs::index_with_reference_sequence_repository(src, repository),
        None => cram::fs::index(src),
    }
}

fn cram_reference_repository_from_header_uri(
    src: &Path,
) -> io::Result<Option<crate::fasta::Repository>> {
    let mut reader = File::open(src).map(cram::io::Reader::new)?;
    let header = reader.read_header()?;

    for (_, reference_sequence) in header.reference_sequences() {
        let Some(uri) = reference_sequence
            .other_fields()
            .get(&sam::header::record::value::map::reference_sequence::tag::URI)
        else {
            continue;
        };
        let uri = String::from_utf8_lossy(uri);
        let path = uri.strip_prefix("file://").unwrap_or(&uri);
        let mut path = PathBuf::from(path);
        if !path.is_absolute()
            && let Some(parent) = src.parent()
        {
            path = parent.join(path);
        }
        if path.is_file() {
            return crate::alignment_compat::cram_reference_repository_from_fasta_path(path)
                .map(Some);
        }
    }

    Ok(None)
}

/// Builds a CSI index for a coordinate-sorted BAM file.
pub fn build_bam_csi<P>(src: P) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    build_bam_csi_with_min_shift(src, 14)
}

/// Builds a CSI index for a coordinate-sorted BAM file with a custom min_shift.
pub fn build_bam_csi_with_min_shift<P>(src: P, min_shift: u8) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    build_bam_csi_from_reader(&mut reader, min_shift)
}

/// Builds a CSI index using a multithreaded BGZF reader.
pub fn build_bam_csi_with_worker_count<P>(
    src: P,
    worker_count: NonZero<usize>,
) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    build_bam_csi_with_min_shift_and_worker_count(src, 14, worker_count)
}

/// Builds a CSI index with a custom min_shift using a multithreaded BGZF reader.
pub fn build_bam_csi_with_min_shift_and_worker_count<P>(
    src: P,
    min_shift: u8,
    worker_count: NonZero<usize>,
) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    let reader = File::open(src)?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, reader);
    let mut reader = bam::io::Reader::from(decoder);
    build_bam_csi_from_reader(&mut reader, min_shift)
}

fn build_bam_csi_from_reader<R>(
    reader: &mut bam::io::Reader<R>,
    min_shift: u8,
) -> io::Result<CsiIndex>
where
    R: Read + bgzf::io::Read,
{
    let header = reader.read_header()?;
    // Size the CSI depth from the largest reference so very large
    // references (e.g. > 2^29, which BAI cannot address) get enough bin
    // levels — matching upstream CSI auto-sizing. The SAM-CSI builder
    // already does this; the BAM path previously used a fixed depth of 5.
    let depth = alignment_csi_depth_for_header(&header, min_shift);
    let mut indexer = alignment_csi_indexer_with_depth(min_shift, depth);
    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_position,
            end_position,
        );
        let alignment_context = bam_alignment_context(&record)?;

        indexer.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

/// Builds a CSI index for a coordinate-sorted BGZF-compressed SAM file.
pub fn build_sam_csi<P>(src: P) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    build_sam_csi_with_min_shift(src, 14)
}

/// Builds a CSI index for a coordinate-sorted BGZF-compressed SAM file with a custom min_shift.
pub fn build_sam_csi_with_min_shift<P>(src: P, min_shift: u8) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(sam::io::Reader::new)?;
    build_sam_csi_from_reader(&mut reader, min_shift)
}

/// Builds a CSI index for a BGZF-compressed SAM file using a multithreaded BGZF reader.
pub fn build_sam_csi_with_worker_count<P>(
    src: P,
    worker_count: NonZero<usize>,
) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    build_sam_csi_with_min_shift_and_worker_count(src, 14, worker_count)
}

/// Builds a CSI index for a BGZF-compressed SAM file with a custom min_shift using a multithreaded BGZF reader.
pub fn build_sam_csi_with_min_shift_and_worker_count<P>(
    src: P,
    min_shift: u8,
    worker_count: NonZero<usize>,
) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    let reader = File::open(src)?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, reader);
    let mut reader = sam::io::Reader::new(decoder);
    build_sam_csi_from_reader(&mut reader, min_shift)
}

fn build_sam_csi_from_reader<R>(
    reader: &mut sam::io::Reader<R>,
    min_shift: u8,
) -> io::Result<CsiIndex>
where
    R: BufRead + bgzf::io::Read,
{
    let header = reader.read_header()?;
    let depth = alignment_csi_depth_for_header(&header, min_shift);
    let mut indexer = alignment_csi_indexer_with_depth(min_shift, depth);
    let mut record = sam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_position,
            end_position,
        );
        let alignment_context = sam_alignment_context(&header, &record)?;

        indexer.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

/// Builds a BAI index for a coordinate-sorted BGZF-compressed SAM file.
pub fn build_sam_bai<P>(src: P) -> io::Result<BaiIndex>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(sam::io::Reader::new)?;
    build_sam_bai_from_reader(&mut reader)
}

/// Builds a BAI index for a BGZF-compressed SAM file using a multithreaded BGZF reader.
pub fn build_sam_bai_with_worker_count<P>(
    src: P,
    worker_count: NonZero<usize>,
) -> io::Result<BaiIndex>
where
    P: AsRef<Path>,
{
    let reader = File::open(src)?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, reader);
    let mut reader = sam::io::Reader::new(decoder);
    build_sam_bai_from_reader(&mut reader)
}

fn build_sam_bai_from_reader<R>(reader: &mut sam::io::Reader<R>) -> io::Result<BaiIndex>
where
    R: BufRead + bgzf::io::Read,
{
    let header = reader.read_header()?;
    let mut indexer = alignment_bai_indexer();
    let mut record = sam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_position,
            end_position,
        );
        let alignment_context = sam_alignment_context(&header, &record)?;

        indexer.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

/// Writes a BAI index to a local file.
pub fn write_bai<P>(dst: P, index: &BaiIndex) -> io::Result<()>
where
    P: AsRef<Path>,
{
    bam::bai::fs::write(dst, index)
}

/// Builds a BGZF-compressed VCF stream and a matching TBI index.
pub fn build_vcf_tbi<R, W>(reader: R, writer: W) -> io::Result<(W, TbiIndex)>
where
    R: BufRead,
    W: Write,
{
    tabix_compat::build_bgzf_and_index(reader, writer, tabix_compat::TextFormat::Vcf)
}

/// Builds a BGZF-compressed VCF stream and a matching CSI index.
pub fn build_vcf_csi<R, W>(reader: R, writer: W) -> io::Result<(W, CsiIndex)>
where
    R: BufRead,
    W: Write,
{
    tabix_compat::build_bgzf_and_csi(reader, writer, tabix_compat::TextFormat::Vcf)
}

/// Builds a BGZF-compressed VCF stream and a matching CSI index with a custom min_shift.
pub fn build_vcf_csi_with_min_shift<R, W>(
    reader: R,
    writer: W,
    min_shift: u8,
) -> io::Result<(W, CsiIndex)>
where
    R: BufRead,
    W: Write,
{
    tabix_compat::build_bgzf_and_csi_with_min_shift(
        reader,
        writer,
        tabix_compat::TextFormat::Vcf,
        min_shift,
    )
}

/// Builds a CSI index for an existing BGZF-compressed VCF file using the
/// default CSI `min_shift` of 14.
///
/// Walks the existing file's BGZF virtual offsets — does not rewrite the
/// data. Equivalent to `bcftools index -c file.vcf.gz`.
pub fn build_vcf_csi_from_path<P>(src: P) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    tabix_compat::build_csi_from_bgzf_path(src, tabix_compat::TextFormat::Vcf)
}

/// Builds a CSI index for an existing BGZF-compressed VCF file with a custom
/// `min_shift`.
pub fn build_vcf_csi_from_path_with_min_shift<P>(src: P, min_shift: u8) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    tabix_compat::build_csi_from_bgzf_path_with_min_shift(
        src,
        tabix_compat::TextFormat::Vcf,
        min_shift,
    )
}

/// Builds a TBI index for an existing BGZF-compressed VCF file. Equivalent to
/// `bcftools index -t file.vcf.gz` (or `tabix -p vcf file.vcf.gz`).
pub fn build_vcf_tbi_from_path<P>(src: P) -> io::Result<TbiIndex>
where
    P: AsRef<Path>,
{
    tabix_compat::build_tbi_from_bgzf_path(src, tabix_compat::TextFormat::Vcf)
}

/// Builds a CSI index for a local BCF file.
pub fn build_bcf_csi<P>(src: P) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    build_bcf_csi_with_min_shift(src, 14)
}

/// Builds a CSI index for a local BCF file with a custom min_shift.
pub fn build_bcf_csi_with_min_shift<P>(src: P, min_shift: u8) -> io::Result<CsiIndex>
where
    P: AsRef<Path>,
{
    use vcf::variant::Record as _;

    let mut reader = File::open(src).map(bcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut indexer = alignment_csi_indexer(min_shift);
    let mut record = bcf::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_position,
            end_position,
        );
        let reference_sequence_id = record.reference_sequence_id()?;
        let start = record
            .variant_start()
            .transpose()?
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing position"))?;
        let end = record.variant_end(&header)?;

        indexer.add_record(Some((reference_sequence_id, start, end, true)), chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.contigs().len()))
}

/// Writes a TBI index to a local file.
pub fn write_tbi<P>(dst: P, index: &TbiIndex) -> io::Result<()>
where
    P: AsRef<Path>,
{
    tabix::fs::write(dst, index)
}

/// Writes a CSI index to a local file.
pub fn write_csi<P>(dst: P, index: &CsiIndex) -> io::Result<()>
where
    P: AsRef<Path>,
{
    csi::fs::write(dst, index)
}

/// Writes a CRAI index to a local file.
pub fn write_cram_crai<P>(dst: P, index: &CraiIndex) -> io::Result<()>
where
    P: AsRef<Path>,
{
    cram::crai::fs::write(dst, index)
}

/// Returns the number of indexed reference sequences in a BAI index.
pub fn bai_reference_sequence_count(index: &BaiIndex) -> usize {
    index.reference_sequences().len()
}

/// Returns the number of indexed reference sequences in a CSI index.
pub fn csi_reference_sequence_count(index: &CsiIndex) -> usize {
    index.reference_sequences().len()
}

/// Returns the number of indexed reference sequences in a TBI index.
pub fn tbi_reference_sequence_count(index: &TbiIndex) -> usize {
    index.reference_sequences().len()
}

/// Returns the number of records in a CRAI index.
pub fn crai_record_count(index: &CraiIndex) -> usize {
    index.len()
}

fn push_index_candidates(candidates: &mut Vec<LocatedIndex>, src: &Path, format: IndexFormat) {
    candidates.push(LocatedIndex {
        path: append_extension(src, format.extension()),
        format,
    });

    let replaced = replace_extension(src, format.extension());

    if replaced != candidates.last().expect("candidate").path {
        candidates.push(LocatedIndex {
            path: replaced,
            format,
        });
    }
}

fn append_extension(src: &Path, extension: &str) -> PathBuf {
    let mut path = OsString::from(src);
    path.push(extension);
    PathBuf::from(path)
}

fn replace_extension(src: &Path, extension: &str) -> PathBuf {
    let mut path = src.to_path_buf();
    path.set_extension(extension.trim_start_matches('.'));
    path
}

fn explicit_index_path(src: &Path) -> Option<PathBuf> {
    let src = src.to_string_lossy();
    let (_, index) = src.split_once(INDEX_DELIMITER)?;
    Some(PathBuf::from(index))
}

fn infer_index_format(path: &Path) -> Option<IndexFormat> {
    let path = path.to_string_lossy();

    if path.ends_with(IndexFormat::Csi.extension()) {
        Some(IndexFormat::Csi)
    } else if path.ends_with(IndexFormat::Bai.extension()) {
        Some(IndexFormat::Bai)
    } else if path.ends_with(IndexFormat::Tbi.extension()) {
        Some(IndexFormat::Tbi)
    } else if path.ends_with(IndexFormat::Crai.extension()) {
        Some(IndexFormat::Crai)
    } else {
        None
    }
}

fn alignment_csi_indexer(
    min_shift: u8,
) -> csi::binning_index::Indexer<csi::binning_index::index::reference_sequence::index::BinnedIndex>
{
    alignment_csi_indexer_with_depth(min_shift, 5)
}

fn alignment_csi_indexer_with_depth(
    min_shift: u8,
    depth: u8,
) -> csi::binning_index::Indexer<csi::binning_index::index::reference_sequence::index::BinnedIndex>
{
    csi::binning_index::Indexer::new(min_shift, depth)
}

fn alignment_csi_depth_for_header(header: &sam::Header, min_shift: u8) -> u8 {
    const DEFAULT_DEPTH: u8 = 5;

    let max_reference_len = header
        .reference_sequences()
        .values()
        .map(|reference_sequence| usize::from(reference_sequence.length()))
        .max()
        .unwrap_or_default() as u128;

    let mut depth = DEFAULT_DEPTH;

    while max_reference_len > csi_max_position(min_shift, depth) {
        depth += 1;
    }

    depth
}

fn csi_max_position(min_shift: u8, depth: u8) -> u128 {
    let bit_count = u32::from(min_shift) + 3 * u32::from(depth);
    (1u128 << bit_count) - 1
}

fn alignment_bai_indexer()
-> csi::binning_index::Indexer<csi::binning_index::index::reference_sequence::index::LinearIndex> {
    csi::binning_index::Indexer::default()
}

fn bam_alignment_context(
    record: &bam::Record,
) -> io::Result<Option<(usize, Position, Position, bool)>> {
    use sam::alignment::Record as _;

    let context = match (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => {
            let is_mapped = !record.flags().is_unmapped();
            Some((id, start, end, is_mapped))
        }
        _ => None,
    };

    Ok(context)
}

fn sam_alignment_context(
    header: &sam::Header,
    record: &sam::Record,
) -> io::Result<Option<(usize, Position, Position, bool)>> {
    use sam::alignment::Record as _;

    let context = match (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => {
            let is_mapped = !record.flags()?.is_unmapped();
            Some((id, start, end, is_mapped))
        }
        _ => None,
    };

    Ok(context)
}

fn missing_associated_index_error(src: &Path, format: IndexFormat) -> io::Error {
    let candidates: Vec<_> = associated_index_candidates(src, format)
        .into_iter()
        .map(|candidate| candidate.path)
        .collect();

    io::Error::new(
        io::ErrorKind::NotFound,
        format!("missing associated index; tried {candidates:?}"),
    )
}

fn unsupported_index_format_error(located: LocatedIndex) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "unsupported associated index format {:?} for {}",
            located.format,
            located.path.display()
        ),
    )
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;
    use std::path::PathBuf;

    use crate::csi::BinningIndex;

    use super::{
        IndexFormat, associated_data_path, bai_reference_sequence_count, build_bai,
        build_bai_with_worker_count, build_bam_csi_with_worker_count, build_bcf_csi, build_vcf_csi,
        build_vcf_csi_with_min_shift, build_vcf_tbi, csi_reference_sequence_count,
        locate_associated_index, read_bai, read_csi, read_tbi, tbi_reference_sequence_count,
    };

    fn fixture(path: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(path)
    }

    #[test]
    fn test_read_bai() {
        let index = read_bai(fixture("repos/htslib/test/index.bam.bai")).unwrap();

        assert!(bai_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_read_bam_csi() {
        let index = read_csi(fixture("repos/htslib/test/index.bam.csi")).unwrap();

        assert!(csi_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_read_vcf_csi() {
        let index = read_csi(fixture("repos/htslib/test/index.vcf.gz.csi")).unwrap();

        assert!(csi_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_read_tbi() {
        let index = read_tbi(fixture("repos/htslib/test/index.vcf.gz.tbi")).unwrap();

        assert!(tbi_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_build_bai() {
        let index = build_bai(fixture("repos/htslib/test/range.bam")).unwrap();

        assert!(bai_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_build_bam_indexes_with_worker_count() {
        let worker_count = NonZero::new(2).unwrap();

        let bai =
            build_bai_with_worker_count(fixture("repos/htslib/test/range.bam"), worker_count).unwrap();
        assert!(bai_reference_sequence_count(&bai) > 0);

        let csi = build_bam_csi_with_worker_count(fixture("repos/htslib/test/range.bam"), worker_count)
            .unwrap();
        assert!(csi_reference_sequence_count(&csi) > 0);
    }

    #[test]
    fn test_build_vcf_tbi_and_csi() {
        let vcf = std::fs::File::open(fixture("repos/htslib/test/index.vcf")).unwrap();
        let (_, tbi) = build_vcf_tbi(std::io::BufReader::new(vcf), Vec::new()).unwrap();
        assert!(tbi_reference_sequence_count(&tbi) > 0);

        let vcf = std::fs::File::open(fixture("repos/htslib/test/index.vcf")).unwrap();
        let (_, csi) = build_vcf_csi(std::io::BufReader::new(vcf), Vec::new()).unwrap();
        assert!(csi_reference_sequence_count(&csi) > 0);

        let vcf = std::fs::File::open(fixture("repos/htslib/test/index.vcf")).unwrap();
        let (_, csi) =
            build_vcf_csi_with_min_shift(std::io::BufReader::new(vcf), Vec::new(), 10).unwrap();
        assert_eq!(csi.min_shift(), 10);
    }

    #[test]
    fn test_build_bcf_csi() {
        let index = build_bcf_csi(fixture("repos/htslib/test/tabix/vcf_file.bcf")).unwrap();

        assert!(csi_reference_sequence_count(&index) > 0);
    }

    #[test]
    fn test_associated_index_candidates() {
        let candidates = super::associated_index_candidates("sample.bam", IndexFormat::Bai);
        let paths: Vec<_> = candidates
            .iter()
            .map(|candidate| (candidate.path.clone(), candidate.format))
            .collect();

        assert_eq!(
            paths,
            vec![
                (PathBuf::from("sample.bam.csi"), IndexFormat::Csi),
                (PathBuf::from("sample.csi"), IndexFormat::Csi),
                (PathBuf::from("sample.bam.bai"), IndexFormat::Bai),
                (PathBuf::from("sample.bai"), IndexFormat::Bai),
            ]
        );

        let candidates = super::associated_index_candidates("sample.vcf.gz", IndexFormat::Tbi);
        let paths: Vec<_> = candidates
            .iter()
            .map(|candidate| (candidate.path.clone(), candidate.format))
            .collect();

        assert_eq!(
            paths,
            vec![
                (PathBuf::from("sample.vcf.gz.csi"), IndexFormat::Csi),
                (PathBuf::from("sample.vcf.csi"), IndexFormat::Csi),
                (PathBuf::from("sample.vcf.gz.tbi"), IndexFormat::Tbi),
                (PathBuf::from("sample.vcf.tbi"), IndexFormat::Tbi),
            ]
        );
    }

    #[test]
    fn test_locate_associated_index_prefers_csi() {
        let dir =
            std::env::temp_dir().join(format!("htslib-rs-index-locate-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir(&dir).unwrap();

        let bam = dir.join("sample.bam");
        let csi = dir.join("sample.bam.csi");
        let bai = dir.join("sample.bam.bai");
        std::fs::write(&bam, []).unwrap();
        std::fs::write(&csi, []).unwrap();
        std::fs::write(&bai, []).unwrap();

        let located = locate_associated_index(&bam, IndexFormat::Bai).unwrap();

        std::fs::remove_dir_all(&dir).unwrap();

        assert_eq!(located.path, csi);
        assert_eq!(located.format, IndexFormat::Csi);
    }

    #[test]
    fn test_explicit_index_delimiter() {
        let candidates =
            super::associated_index_candidates("sample.bam##idx##custom.bai", IndexFormat::Bai);

        assert_eq!(
            candidates,
            [super::LocatedIndex {
                path: PathBuf::from("custom.bai"),
                format: IndexFormat::Bai
            }]
        );

        assert_eq!(
            associated_data_path("sample.bam##idx##custom.bai"),
            PathBuf::from("sample.bam")
        );
    }
}
