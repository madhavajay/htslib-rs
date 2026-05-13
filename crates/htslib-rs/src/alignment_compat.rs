//! HTSlib-compatible alignment I/O helpers backed by noodles readers.

use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Read, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use crate::{
    bam, bgzf,
    core::Region,
    cram,
    error::{IoResultExt as _, io_error_with_path},
    expr::{Filter, Value as ExprValue},
    fasta,
    index_compat::{
        associated_data_path, build_bai, build_cram_crai, read_associated_bam_index,
        read_associated_cram_index, write_bai, write_cram_crai,
    },
    probaln::{ProbalnParams, probaln_glocal},
    sam,
};

/// A SAM header shared by SAM, BAM, and CRAM streams.
pub type Header = sam::Header;

/// A synchronized multi-input pileup column.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SynchronizedPileupColumn {
    /// Reference sequence name.
    pub reference_name: String,
    /// 1-based reference position.
    pub position: usize,
    /// Total depth across all inputs.
    pub total_depth: usize,
    /// Per-input depth at this position.
    pub depths_by_input: Vec<usize>,
    /// Per-input pileup base strings.
    pub bases_by_input: Vec<String>,
    /// Per-input pileup quality strings.
    pub qualities_by_input: Vec<String>,
}

type SynchronizedPileupSite = (String, usize);
type SynchronizedPileupEntry = (usize, usize, usize);

/// A mutable adapter for HTSlib-style SAM header operations.
pub struct SamHeaderAdapter<'a> {
    header: &'a mut Header,
}

impl<'a> SamHeaderAdapter<'a> {
    /// Creates a mutable header adapter.
    pub fn new(header: &'a mut Header) -> Self {
        Self { header }
    }

    /// Returns the number of reference sequences in the header.
    pub fn reference_sequence_count(&self) -> usize {
        self.header.reference_sequences().len()
    }

    /// Returns the length of a reference sequence by name.
    pub fn reference_sequence_len(&self, name: &str) -> Option<usize> {
        self.header
            .reference_sequences()
            .get(name.as_bytes())
            .map(|reference_sequence| usize::from(reference_sequence.length()))
    }

    /// Inserts or replaces an `@SQ` reference sequence and returns the replaced length.
    pub fn insert_reference_sequence(
        &mut self,
        name: &str,
        len: usize,
    ) -> io::Result<Option<usize>> {
        use sam::header::record::value::{Map, map::ReferenceSequence};

        let len = NonZero::new(len).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference sequence length must be nonzero",
            )
        })?;

        Ok(self
            .header
            .reference_sequences_mut()
            .insert(name.as_bytes().into(), Map::<ReferenceSequence>::new(len))
            .map(|reference_sequence| usize::from(reference_sequence.length())))
    }

    /// Removes an `@SQ` reference sequence and returns its length.
    pub fn remove_reference_sequence(&mut self, name: &str) -> Option<usize> {
        self.header
            .reference_sequences_mut()
            .shift_remove(name.as_bytes())
            .map(|reference_sequence| usize::from(reference_sequence.length()))
    }
}

/// A mutable adapter for HTSlib-style SAM record core field operations.
pub struct SamRecordAdapter<'a> {
    record: &'a mut sam::alignment::RecordBuf,
}

impl<'a> SamRecordAdapter<'a> {
    /// Creates a mutable record adapter.
    pub fn new(record: &'a mut sam::alignment::RecordBuf) -> Self {
        Self { record }
    }

    /// Returns the query template name as UTF-8 lossily decoded text.
    pub fn name(&self) -> Option<String> {
        self.record
            .name()
            .map(|name| String::from_utf8_lossy(name).into_owned())
    }

    /// Sets or clears the query template name.
    pub fn set_name(&mut self, name: Option<&str>) {
        *self.record.name_mut() = name.map(Into::into);
    }

    /// Returns the raw SAM flag bits.
    pub fn flags(&self) -> u16 {
        self.record.flags().bits()
    }

    /// Sets the raw SAM flag bits.
    pub fn set_flags(&mut self, flags: u16) {
        *self.record.flags_mut() = flags.into();
    }

    /// Returns the reference sequence ID.
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.record.reference_sequence_id()
    }

    /// Sets or clears the reference sequence ID.
    pub fn set_reference_sequence_id(&mut self, reference_sequence_id: Option<usize>) {
        *self.record.reference_sequence_id_mut() = reference_sequence_id;
    }

    /// Returns the 1-based alignment start.
    pub fn alignment_start(&self) -> Option<usize> {
        self.record.alignment_start().map(usize::from)
    }

    /// Sets or clears the 1-based alignment start.
    pub fn set_alignment_start(&mut self, alignment_start: Option<usize>) -> io::Result<()> {
        *self.record.alignment_start_mut() =
            position_from_one_based("alignment start", alignment_start)?;

        Ok(())
    }

    /// Returns the raw mapping quality, using 255 for a missing value.
    pub fn mapping_quality(&self) -> u8 {
        self.record
            .mapping_quality()
            .map(|mapping_quality| mapping_quality.get())
            .unwrap_or(255)
    }

    /// Sets the raw mapping quality, using 255 to clear the value.
    pub fn set_mapping_quality(&mut self, mapping_quality: u8) {
        use sam::alignment::record::MappingQuality;

        *self.record.mapping_quality_mut() = MappingQuality::new(mapping_quality);
    }

    /// Returns the mate reference sequence ID.
    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.record.mate_reference_sequence_id()
    }

    /// Sets or clears the mate reference sequence ID.
    pub fn set_mate_reference_sequence_id(&mut self, mate_reference_sequence_id: Option<usize>) {
        *self.record.mate_reference_sequence_id_mut() = mate_reference_sequence_id;
    }

    /// Returns the 1-based mate alignment start.
    pub fn mate_alignment_start(&self) -> Option<usize> {
        self.record.mate_alignment_start().map(usize::from)
    }

    /// Sets or clears the 1-based mate alignment start.
    pub fn set_mate_alignment_start(
        &mut self,
        mate_alignment_start: Option<usize>,
    ) -> io::Result<()> {
        *self.record.mate_alignment_start_mut() =
            position_from_one_based("mate alignment start", mate_alignment_start)?;

        Ok(())
    }

    /// Returns the template length.
    pub fn template_length(&self) -> i32 {
        self.record.template_length()
    }

    /// Sets the template length.
    pub fn set_template_length(&mut self, template_length: i32) {
        *self.record.template_length_mut() = template_length;
    }
}

fn position_from_one_based(
    field: &'static str,
    position: Option<usize>,
) -> io::Result<Option<crate::core::Position>> {
    position
        .map(|n| {
            crate::core::Position::new(n).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("{field} must be nonzero"),
                )
            })
        })
        .transpose()
}

/// A typed SAM auxiliary field value for HTSlib-style record mutation.
#[derive(Clone, Debug, PartialEq)]
pub enum SamAuxValue {
    /// A character (`A`).
    Character(u8),
    /// An 8-bit integer (`c`).
    Int8(i8),
    /// An 8-bit unsigned integer (`C`).
    UInt8(u8),
    /// A 16-bit integer (`s`).
    Int16(i16),
    /// A 16-bit unsigned integer (`S`).
    UInt16(u16),
    /// A 32-bit integer (`i`).
    Int32(i32),
    /// A 32-bit unsigned integer (`I`).
    UInt32(u32),
    /// A single-precision floating-point value (`f`).
    Float(f32),
    /// A string (`Z`).
    String(String),
    /// A hex string (`H`).
    Hex(String),
    /// An 8-bit integer array (`B:c`).
    Int8Array(Vec<i8>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8Array(Vec<u8>),
    /// A 16-bit integer array (`B:s`).
    Int16Array(Vec<i16>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16Array(Vec<u16>),
    /// A 32-bit integer array (`B:i`).
    Int32Array(Vec<i32>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32Array(Vec<u32>),
    /// A single-precision floating-point array (`B:f`).
    FloatArray(Vec<f32>),
}

/// Returns a typed auxiliary field from a mutable SAM-compatible record buffer.
pub fn sam_aux_get(record: &sam::alignment::RecordBuf, tag: [u8; 2]) -> Option<SamAuxValue> {
    use sam::alignment::record::data::field::Tag;

    record.data().get(&Tag::from(tag)).map(SamAuxValue::from)
}

/// Inserts or replaces a typed auxiliary field in a mutable SAM-compatible record buffer.
pub fn sam_aux_insert(
    record: &mut sam::alignment::RecordBuf,
    tag: [u8; 2],
    value: SamAuxValue,
) -> Option<SamAuxValue> {
    use sam::alignment::record::data::field::Tag;

    record
        .data_mut()
        .insert(Tag::from(tag), value.into())
        .map(|(_, value)| SamAuxValue::from(&value))
}

/// Removes a typed auxiliary field from a mutable SAM-compatible record buffer.
pub fn sam_aux_remove(record: &mut sam::alignment::RecordBuf, tag: [u8; 2]) -> Option<SamAuxValue> {
    use sam::alignment::record::data::field::Tag;

    record
        .data_mut()
        .remove(&Tag::from(tag))
        .map(|(_, value)| SamAuxValue::from(&value))
}

impl From<&sam::alignment::record_buf::data::field::Value> for SamAuxValue {
    fn from(value: &sam::alignment::record_buf::data::field::Value) -> Self {
        use sam::alignment::record_buf::data::field::{Value, value::Array};

        match value {
            Value::Character(n) => Self::Character(*n),
            Value::Int8(n) => Self::Int8(*n),
            Value::UInt8(n) => Self::UInt8(*n),
            Value::Int16(n) => Self::Int16(*n),
            Value::UInt16(n) => Self::UInt16(*n),
            Value::Int32(n) => Self::Int32(*n),
            Value::UInt32(n) => Self::UInt32(*n),
            Value::Float(n) => Self::Float(*n),
            Value::String(s) => Self::String(s.to_string()),
            Value::Hex(s) => Self::Hex(s.to_string()),
            Value::Array(Array::Int8(values)) => Self::Int8Array(values.clone()),
            Value::Array(Array::UInt8(values)) => Self::UInt8Array(values.clone()),
            Value::Array(Array::Int16(values)) => Self::Int16Array(values.clone()),
            Value::Array(Array::UInt16(values)) => Self::UInt16Array(values.clone()),
            Value::Array(Array::Int32(values)) => Self::Int32Array(values.clone()),
            Value::Array(Array::UInt32(values)) => Self::UInt32Array(values.clone()),
            Value::Array(Array::Float(values)) => Self::FloatArray(values.clone()),
        }
    }
}

impl From<SamAuxValue> for sam::alignment::record_buf::data::field::Value {
    fn from(value: SamAuxValue) -> Self {
        use sam::alignment::record_buf::data::field::{Value, value::Array};

        match value {
            SamAuxValue::Character(n) => Value::Character(n),
            SamAuxValue::Int8(n) => Value::Int8(n),
            SamAuxValue::UInt8(n) => Value::UInt8(n),
            SamAuxValue::Int16(n) => Value::Int16(n),
            SamAuxValue::UInt16(n) => Value::UInt16(n),
            SamAuxValue::Int32(n) => Value::Int32(n),
            SamAuxValue::UInt32(n) => Value::UInt32(n),
            SamAuxValue::Float(n) => Value::Float(n),
            SamAuxValue::String(s) => Value::String(s.into()),
            SamAuxValue::Hex(s) => Value::Hex(s.into()),
            SamAuxValue::Int8Array(values) => Value::Array(Array::Int8(values)),
            SamAuxValue::UInt8Array(values) => Value::Array(Array::UInt8(values)),
            SamAuxValue::Int16Array(values) => Value::Array(Array::Int16(values)),
            SamAuxValue::UInt16Array(values) => Value::Array(Array::UInt16(values)),
            SamAuxValue::Int32Array(values) => Value::Array(Array::Int32(values)),
            SamAuxValue::UInt32Array(values) => Value::Array(Array::UInt32(values)),
            SamAuxValue::FloatArray(values) => Value::Array(Array::Float(values)),
        }
    }
}

/// A format-neutral alignment record summary for parity tests.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentRecordSummary {
    name: Option<Vec<u8>>,
    flags: sam::alignment::record::Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<crate::core::Position>,
    mapping_quality: Option<u8>,
    cigar: Vec<sam::alignment::record::cigar::Op>,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<crate::core::Position>,
    template_length: i32,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

/// Split FASTA/FASTQ text outputs for paired-read extraction.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FastxSplitText {
    pub read1: String,
    pub read2: String,
    pub singleton: String,
}

impl AlignmentRecordSummary {
    /// Returns the template length.
    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    /// Returns the alignment flags as a `u16`.
    pub fn flags_u16(&self) -> u16 {
        self.flags.bits()
    }

    /// Returns the typed alignment flags.
    pub fn flags(&self) -> sam::alignment::record::Flags {
        self.flags
    }

    /// Returns the 0-based reference sequence id, if any.
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    /// Returns the 0-based mate reference sequence id, if any.
    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id
    }

    /// Returns the mapping quality, if any.
    pub fn mapping_quality(&self) -> Option<u8> {
        self.mapping_quality
    }
}

/// Reads a SAM header from a buffered reader.
pub fn read_sam_header<R>(reader: R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut reader = sam::io::Reader::new(reader);

    reader.read_header()
}

/// Reads a SAM header from a local file.
pub fn read_sam_header_from_path<P>(src: P) -> io::Result<Header>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(read_sam_header)
}

/// Counts SAM records from a buffered reader.
pub fn count_sam_records<R>(reader: R) -> io::Result<usize>
where
    R: BufRead,
{
    let mut reader = sam::io::Reader::new(reader);
    reader.read_header()?;

    reader
        .records()
        .try_fold(0, |n, result| result.map(|_| n + 1))
}

/// Counts SAM records from a local file.
pub fn count_sam_records_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(count_sam_records)
}

/// Reads SAM input without producing output and returns the number of records seen.
pub fn benchmark_sam_view_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    count_sam_records_from_path(src)
}

/// Reads SAM records from a local file into format-neutral summaries.
pub fn summarize_sam_records_from_path<P>(src: P) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;

    reader
        .records()
        .map(|result| result.and_then(|record| summarize_alignment_record(&header, &record)))
        .collect()
}

/// Writes a SAM view of the first `limit` records, including the header.
pub fn view_sam_text_from_path_with_limit<P>(src: P, limit: Option<usize>) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_sam_text_from_path_with_limit_and_parse_errors(src, limit, false)
}

/// Writes a SAM view while optionally ignoring malformed SAM records.
pub fn view_sam_text_from_path_with_limit_and_parse_errors<P>(
    src: P,
    limit: Option<usize>,
    ignore_parse_errors: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(Vec::new());
    let mut remaining = limit.unwrap_or(usize::MAX);

    writer.write_header(&header)?;

    for result in reader.records() {
        if remaining == 0 {
            break;
        }

        let record = match result {
            Ok(record) => record,
            Err(e) if ignore_parse_errors => {
                drop(e);
                continue;
            }
            Err(e) => return Err(e),
        };

        writer.write_record(&header, &record)?;
        remaining -= 1;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a SAM view of the first `limit` records from a BGZF-compressed SAM file.
pub fn view_bgzf_sam_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(bgzf::io::Reader::new)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes SAM input as BGZF-compressed SAM, including the header and first `limit` records.
pub fn write_bgzf_sam_from_path_with_limit<P, W>(
    src: P,
    dst: W,
    limit: Option<usize>,
) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let bgzf_writer = bgzf::io::Writer::new(dst);
    let mut writer = sam::io::Writer::new(bgzf_writer);

    writer.write_header(&header)?;

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    writer.into_inner().finish()
}

/// Writes a FASTQ view of the first `limit` SAM records.
pub fn view_sam_as_fastq_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fastq_record(&mut writer, &record)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of the first `limit` SAM records, optionally appending read numbers.
pub fn view_sam_as_fastq_text_from_path_with_limit_and_suffix<P>(
    src: P,
    limit: Option<usize>,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fastq_record_with_suffix(&mut writer, &record, append_read_number)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of SAM records passing flag filters.
pub fn view_sam_as_fastq_text_from_path_with_flag_filter<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fastq_record(&mut writer, &record)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of filtered SAM records, optionally appending read numbers.
pub fn view_sam_as_fastq_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    view_sam_as_fastq_text_from_reader_with_flag_filter_and_suffix(
        &mut reader,
        require_flags,
        exclude_flags,
        exclude_all_flags,
        append_read_number,
    )
}

/// Writes a FASTQ view of filtered SAM records from any buffered reader.
pub fn view_sam_as_fastq_text_from_reader_with_flag_filter_and_suffix<R>(
    reader: &mut sam::io::Reader<R>,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    R: BufRead,
{
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fastq_record_with_suffix(&mut writer, &record, append_read_number)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of the first `limit` SAM records.
pub fn view_sam_as_fasta_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fasta_record(&mut writer, &record)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of the first `limit` SAM records, optionally appending read numbers.
pub fn view_sam_as_fasta_text_from_path_with_limit_and_suffix<P>(
    src: P,
    limit: Option<usize>,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fasta_record_with_suffix(&mut writer, &record, append_read_number)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of SAM records passing flag filters.
pub fn view_sam_as_fasta_text_from_path_with_flag_filter<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fasta_record(&mut writer, &record)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of filtered SAM records, optionally appending read numbers.
pub fn view_sam_as_fasta_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    view_sam_as_fasta_text_from_reader_with_flag_filter_and_suffix(
        &mut reader,
        require_flags,
        exclude_flags,
        exclude_all_flags,
        append_read_number,
    )
}

/// Writes a FASTA view of filtered SAM records from any buffered reader.
pub fn view_sam_as_fasta_text_from_reader_with_flag_filter_and_suffix<R>(
    reader: &mut sam::io::Reader<R>,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    R: BufRead,
{
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fasta_record_with_suffix(&mut writer, &record, append_read_number)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes split FASTQ views of SAM records passing flag filters.
pub fn view_sam_as_fastq_split_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    view_sam_as_fastq_split_text_from_reader_with_flag_filter_and_suffix(
        &mut reader,
        require_flags,
        exclude_flags,
        exclude_all_flags,
        append_read_number,
    )
}

/// Writes split FASTQ views of filtered SAM records from any buffered reader.
pub fn view_sam_as_fastq_split_text_from_reader_with_flag_filter_and_suffix<R>(
    reader: &mut sam::io::Reader<R>,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    R: BufRead,
{
    view_sam_as_fastq_split_text_from_reader_with_flag_filter_suffix_and_aux(
        reader,
        require_flags,
        exclude_flags,
        exclude_all_flags,
        append_read_number,
        None,
    )
}

/// Writes split FASTQ views of filtered SAM records from any buffered reader, preserving selected aux tags.
pub fn view_sam_as_fastq_split_text_from_reader_with_flag_filter_suffix_and_aux<R>(
    reader: &mut sam::io::Reader<R>,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
    aux_tags: Option<&[[u8; 2]]>,
) -> io::Result<FastxSplitText>
where
    R: BufRead,
{
    let _header = reader.read_header()?;
    let mut split = FastxSplitBuffers::default();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_split_fastq_record_with_aux(&mut split, &record, append_read_number, aux_tags)?;
        }
    }

    split.into_text()
}

/// Writes split FASTA views of SAM records passing flag filters.
pub fn view_sam_as_fasta_split_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    view_sam_as_fasta_split_text_from_reader_with_flag_filter_and_suffix(
        &mut reader,
        require_flags,
        exclude_flags,
        exclude_all_flags,
        append_read_number,
    )
}

/// Writes split FASTA views of filtered SAM records from any buffered reader.
pub fn view_sam_as_fasta_split_text_from_reader_with_flag_filter_and_suffix<R>(
    reader: &mut sam::io::Reader<R>,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    R: BufRead,
{
    let _header = reader.read_header()?;
    let mut split = FastxSplitBuffers::default();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_split_fasta_record(&mut split, &record, append_read_number)?;
        }
    }

    split.into_text()
}

/// Applies existing `BQ:Z` BAQ tags to SAM quality strings and renames them to `ZQ:Z`.
pub fn apply_existing_baq_from_sam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    convert_existing_baq_from_sam_path(src, true)
}

/// Reverts existing `ZQ:Z` BAQ tags from SAM quality strings and renames them to `BQ:Z`.
pub fn revert_existing_baq_from_sam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    convert_existing_baq_from_sam_path(src, false)
}

/// Recalculates non-extended BAQ tags for a local SAM file.
///
/// This mirrors the default `sam_prob_realn` path used by HTSlib's `calmd` tests for local SAM
/// input: existing BAQ tags are left unchanged, reference skips and records without match CIGAR
/// operations are no-ops, and recalculated BAQ is emitted as `BQ:Z`.
pub fn recalculate_baq_from_sam_path<P, Q>(sam_src: P, reference_src: Q) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    recalculate_baq_from_sam_path_with_options(sam_src, reference_src, false, false, false)
}

/// Recalculates non-extended BAQ tags, applies them to qualities, and emits `ZQ:Z`.
pub fn recalculate_and_apply_baq_from_sam_path<P, Q>(
    sam_src: P,
    reference_src: Q,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    recalculate_baq_from_sam_path_with_options(sam_src, reference_src, true, false, false)
}

/// Forces non-extended BAQ recalculation, replacing any existing `BQ:Z` tag.
pub fn force_recalculate_baq_from_sam_path<P, Q>(sam_src: P, reference_src: Q) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    recalculate_baq_from_sam_path_with_options(sam_src, reference_src, false, true, false)
}

/// Recalculates extended BAQ tags for a local SAM file.
pub fn recalculate_extended_baq_from_sam_path<P, Q>(
    sam_src: P,
    reference_src: Q,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    recalculate_baq_from_sam_path_with_options(sam_src, reference_src, false, false, true)
}

/// Calculates a BAQ tag for a single mpileup-style read/reference alignment.
///
/// `alignment_start` is 0-based. The returned string is the raw `BQ:Z`/`ZQ:Z`
/// payload, or `None` when HTSlib would skip BAQ for the CIGAR shape.
pub fn mpileup_baq_from_alignment(
    sequence: &[u8],
    qualities: &[u8],
    reference: &[u8],
    alignment_start: usize,
    cigar: &str,
    extended: bool,
) -> io::Result<Option<String>> {
    if sequence.len() != qualities.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "sequence and quality lengths differ",
        ));
    }

    let cigar = parse_sam_cigar(cigar)?;
    calculate_baq(
        sequence,
        qualities,
        reference,
        alignment_start,
        &cigar,
        extended,
    )
}

/// Scores a read against a candidate haplotype using the `bam2bcf_*` indel
/// realignment `probaln_glocal` parameters.
pub fn mpileup_indel_alignment_score(
    reference: &[u8],
    query: &[u8],
    qualities: Option<&[u8]>,
) -> io::Result<i32> {
    if qualities.is_some_and(|qualities| qualities.len() != query.len()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "quality length differs from query length",
        ));
    }

    let translated_ref = reference
        .iter()
        .map(|base| probaln_base_code(*base))
        .collect::<Vec<_>>();
    let translated_query = query
        .iter()
        .map(|base| probaln_base_code(*base))
        .collect::<Vec<_>>();
    let clipped_qualities = qualities.map(|qualities| {
        qualities
            .iter()
            .map(|quality| (*quality).clamp(7, 30))
            .collect::<Vec<_>>()
    });
    let result = probaln_glocal(
        &translated_ref,
        &translated_query,
        clipped_qualities.as_deref(),
        ProbalnParams {
            d: 1e-4,
            e: 1e-2,
            bw: 10,
        },
        false,
    )?;

    Ok(result.likelihood)
}

fn recalculate_baq_from_sam_path_with_options<P, Q>(
    sam_src: P,
    reference_src: Q,
    apply: bool,
    force: bool,
    extended: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequences = read_fasta_sequences(reference_src)?;
    let text = std::fs::read_to_string(sam_src)?;
    let mut out = String::new();

    for line in text.lines() {
        if line.starts_with('@') {
            out.push_str(line);
        } else {
            out.push_str(&recalculate_baq_sam_line(
                line,
                &reference_sequences,
                apply,
                force,
                extended,
            )?);
        }

        out.push('\n');
    }

    Ok(out)
}

/// Computes the HTSlib `sam_cap_mapq` mapping-quality cap for a SAM record.
pub fn cap_mapping_quality_for_record<R>(
    header: &Header,
    record: &R,
    reference_sequence: &[u8],
    threshold: i32,
) -> io::Result<i32>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::cigar::op::Kind;

    let threshold = if threshold < 0 { 40 } else { threshold };
    let Some(alignment_start) = record.alignment_start().transpose()? else {
        return Ok(-1);
    };

    record.reference_sequence_id(header).transpose()?;

    let sequence = record.sequence().iter().collect::<Vec<_>>();
    let mut quality_scores = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<_>>>()?;

    if quality_scores.is_empty() {
        quality_scores.resize(sequence.len(), u8::MAX);
    }

    let mut reference_index = usize::from(alignment_start) - 1;
    let mut query_index = 0usize;
    let mut mismatch_count = 0usize;
    let mut mismatch_quality_sum = 0usize;
    let mut len = 0usize;
    let mut clipped_quality_sum = 0usize;

    'ops: for result in record.cigar().iter() {
        let op = result?;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let mut consumed = 0usize;

                for offset in 0..op.len() {
                    let ref_base = match reference_sequence.get(reference_index + offset) {
                        Some(base) if *base != 0 => *base,
                        _ => break,
                    };
                    let read_base =
                        sequence.get(query_index + offset).copied().ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "CIGAR qpos is out of range")
                        })?;
                    let quality = quality_scores
                        .get(query_index + offset)
                        .copied()
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "quality score qpos is out of range",
                            )
                        })?;
                    let read_code = nt16_code(read_base);
                    let ref_code = nt16_code(ref_base);

                    if ref_code != 15 && read_code != 15 && quality >= 13 {
                        len += 1;

                        if read_code != 0 && read_code != ref_code {
                            mismatch_count += 1;
                            mismatch_quality_sum += usize::from(quality.min(33));
                        }
                    }

                    consumed += 1;
                }

                if consumed < op.len() {
                    break 'ops;
                }

                reference_index += op.len();
                query_index += op.len();
                len += op.len();
            }
            Kind::Deletion => {
                for offset in 0..op.len() {
                    match reference_sequence.get(reference_index + offset) {
                        Some(base) if *base != 0 => {}
                        _ => break 'ops,
                    }
                }

                reference_index += op.len();
            }
            Kind::SoftClip => {
                for offset in 0..op.len() {
                    clipped_quality_sum += usize::from(
                        quality_scores
                            .get(query_index + offset)
                            .copied()
                            .ok_or_else(|| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    "soft clip qpos is out of range",
                                )
                            })?,
                    );
                }

                query_index += op.len();
            }
            Kind::HardClip => {
                clipped_quality_sum += 13 * op.len();
            }
            Kind::Insertion => {
                query_index += op.len();
            }
            Kind::Skip => {
                reference_index += op.len();
            }
            Kind::Pad => {}
        }
    }

    let mut t = 1.0;
    let len = len as f64;

    for i in 0..mismatch_count {
        t *= len / (i + 1) as f64;
    }

    let t = mismatch_quality_sum as f64 - 4.343 * t.ln() + clipped_quality_sum as f64 / 5.0;

    if t > f64::from(threshold) {
        return Ok(-1);
    }

    let t = t.max(0.0);
    let cap = ((f64::from(threshold) - t) / f64::from(threshold)).sqrt() * f64::from(threshold);

    Ok((cap + 0.499) as i32)
}

/// Computes `sam_cap_mapq` caps for all records in a local SAM file.
pub fn cap_mapping_qualities_from_sam_path<P, Q>(
    sam_src: P,
    reference_src: Q,
    threshold: i32,
) -> io::Result<Vec<i32>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequences = read_fasta_sequences(reference_src)?;
    let mut reader = File::open(sam_src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut caps = Vec::new();

    for result in reader.records() {
        let record = result?;
        let Some(reference_sequence_id) = record.reference_sequence_id(&header).transpose()? else {
            caps.push(-1);
            continue;
        };
        let (reference_sequence_name, _) = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "reference sequence ID is out of range",
                )
            })?;
        let reference_sequence = reference_sequences
            .get(reference_sequence_name.as_ref() as &[u8])
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "reference sequence {} not found",
                        String::from_utf8_lossy(reference_sequence_name)
                    ),
                )
            })?;

        caps.push(cap_mapping_quality_for_record(
            &header,
            &record,
            reference_sequence,
            threshold,
        )?);
    }

    Ok(caps)
}

fn nt16_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

fn convert_existing_baq_from_sam_path<P>(src: P, apply: bool) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let text = std::fs::read_to_string(src)?;
    let mut out = String::new();

    for line in text.lines() {
        if line.starts_with('@') {
            out.push_str(line);
        } else {
            out.push_str(&convert_existing_baq_sam_line(line, apply)?);
        }

        out.push('\n');
    }

    Ok(out)
}

fn convert_existing_baq_sam_line(line: &str, apply: bool) -> io::Result<String> {
    let tag_prefix = if apply { "BQ:Z:" } else { "ZQ:Z:" };
    let replacement_prefix = if apply { "ZQ:Z:" } else { "BQ:Z:" };
    let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();

    if fields.len() < 11 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM record has fewer than 11 fields",
        ));
    }

    let Some(tag_index) = fields
        .iter()
        .position(|field| field.starts_with(tag_prefix))
    else {
        return Ok(line.to_owned());
    };

    let tag = fields[tag_index]
        .strip_prefix(tag_prefix)
        .expect("tag prefix was matched")
        .to_owned();

    if fields[10].len() != tag.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BAQ tag length does not match SAM quality length",
        ));
    }

    let mut quality = Vec::with_capacity(fields[10].len());

    for (qual, baq) in fields[10].bytes().zip(tag.bytes()) {
        if !(33..=126).contains(&qual) || baq < 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid SAM quality or BAQ character",
            ));
        }

        let qual = qual - 33;
        let baq = baq - 64;
        let adjusted = if apply {
            qual.saturating_sub(baq)
        } else {
            qual.checked_add(baq).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "quality score overflow")
            })?
        };

        quality.push(
            adjusted.checked_add(33).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "quality score overflow")
            })?,
        );
    }

    fields[10] =
        String::from_utf8(quality).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    fields[tag_index] = format!("{replacement_prefix}{tag}");

    Ok(fields.join("\t"))
}

#[derive(Clone, Copy, Debug)]
struct SamCigarOp {
    len: usize,
    op: u8,
}

fn recalculate_baq_sam_line(
    line: &str,
    reference_sequences: &HashMap<Vec<u8>, Vec<u8>>,
    apply: bool,
    force: bool,
    extended: bool,
) -> io::Result<String> {
    let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();

    if fields.len() < 11 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM record has fewer than 11 fields",
        ));
    }

    let has_bq = fields
        .iter()
        .skip(11)
        .any(|field| field.starts_with("BQ:Z:") || field.starts_with("BQ:B:"));
    let has_zq = fields
        .iter()
        .skip(11)
        .any(|field| field.starts_with("ZQ:Z:"));

    if !force && (has_bq || has_zq) {
        return Ok(line.to_owned());
    }

    if force && has_bq {
        let mut kept_fields = fields.drain(..11).collect::<Vec<_>>();
        kept_fields.extend(
            fields
                .into_iter()
                .filter(|field| !field.starts_with("BQ:Z:") && !field.starts_with("BQ:B:")),
        );
        fields = kept_fields;
    }

    if has_zq {
        return Ok(fields.join("\t"));
    }

    let flags = fields[1].parse::<u16>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid SAM flags: {e}"),
        )
    })?;

    if flags & 0x04 != 0 || fields[9] == "*" || fields[10] == "*" {
        return Ok(line.to_owned());
    }

    let pos = fields[3].parse::<usize>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid SAM position: {e}"),
        )
    })?;

    if pos == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM position must be nonzero",
        ));
    }

    let cigar = parse_sam_cigar(&fields[5])?;
    let sequence = fields[9].as_bytes();
    let qualities = sam_quality_scores(fields[10].as_bytes())?;

    if sequence.len() != qualities.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM sequence and quality lengths differ",
        ));
    }

    let Some(reference) = reference_sequences.get(fields[2].as_bytes()) else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("reference sequence {} not found", fields[2]),
        ));
    };

    let Some(baq) = calculate_baq(sequence, &qualities, reference, pos - 1, &cigar, extended)?
    else {
        return Ok(line.to_owned());
    };

    if apply {
        fields[10] = apply_baq_to_quality_string(fields[10].as_bytes(), baq.as_bytes())?;
        fields.push(format!("ZQ:Z:{baq}"));
    } else {
        fields.push(format!("BQ:Z:{baq}"));
    }

    Ok(fields.join("\t"))
}

fn parse_sam_cigar(s: &str) -> io::Result<Vec<SamCigarOp>> {
    if s == "*" {
        return Ok(Vec::new());
    }

    let mut ops = Vec::new();
    let mut len = 0usize;
    let mut have_len = false;

    for b in s.bytes() {
        if b.is_ascii_digit() {
            have_len = true;
            len = len
                .checked_mul(10)
                .and_then(|n| n.checked_add(usize::from(b - b'0')))
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CIGAR overflow"))?;
        } else {
            if !have_len {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CIGAR op is missing a length",
                ));
            }

            match b {
                b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X' => {
                    ops.push(SamCigarOp { len, op: b });
                    len = 0;
                    have_len = false;
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid CIGAR op {}", char::from(b)),
                    ));
                }
            }
        }
    }

    if have_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CIGAR has trailing length without op",
        ));
    }

    Ok(ops)
}

fn sam_quality_scores(s: &[u8]) -> io::Result<Vec<u8>> {
    s.iter()
        .map(|b| {
            if !(33..=126).contains(b) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid SAM quality character",
                ));
            }

            Ok(b - 33)
        })
        .collect()
}

fn apply_baq_to_quality_string(qualities: &[u8], baq: &[u8]) -> io::Result<String> {
    if qualities.len() != baq.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BAQ tag length does not match SAM quality length",
        ));
    }

    let adjusted = qualities
        .iter()
        .zip(baq)
        .map(|(quality, baq)| {
            if !(33..=126).contains(quality) || *baq < 64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid SAM quality or BAQ character",
                ));
            }

            let quality = quality - 33;
            let baq_adjustment = baq - 64;

            Ok(quality.saturating_sub(baq_adjustment) + 33)
        })
        .collect::<io::Result<Vec<_>>>()?;

    String::from_utf8(adjusted).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn calculate_baq(
    sequence: &[u8],
    qualities: &[u8],
    reference: &[u8],
    alignment_start: usize,
    cigar: &[SamCigarOp],
    extended: bool,
) -> io::Result<Option<String>> {
    let mut reference_index = alignment_start;
    let mut query_index = 0usize;
    let mut query_begin = None;
    let mut query_end = 0usize;
    let mut reference_begin = None;
    let mut reference_end = 0usize;

    for op in cigar {
        match op.op {
            b'M' | b'=' | b'X' => {
                query_begin.get_or_insert(query_index);
                reference_begin.get_or_insert(reference_index);
                query_end = query_index.saturating_add(op.len);
                reference_end = reference_index.saturating_add(op.len);
                reference_index = reference_index.saturating_add(op.len);
                query_index = query_index.saturating_add(op.len);
            }
            b'S' | b'I' => query_index = query_index.saturating_add(op.len),
            b'D' => reference_index = reference_index.saturating_add(op.len),
            b'N' => return Ok(None),
            b'H' | b'P' => {}
            _ => unreachable!("CIGAR op validated by parse_sam_cigar"),
        }
    }

    let (query_begin, reference_begin) = match (query_begin, reference_begin) {
        (Some(query_begin), Some(reference_begin)) => (query_begin, reference_begin),
        _ => return Ok(None),
    };

    let mut bandwidth = 7usize;
    let reference_span = reference_end.saturating_sub(reference_begin);
    let query_span = query_end.saturating_sub(query_begin);

    if reference_span.abs_diff(query_span) > bandwidth {
        bandwidth = reference_span.abs_diff(query_span) + 3;
    }

    let mut window_begin = reference_begin.saturating_sub(query_begin + bandwidth / 2);
    let mut window_end = reference_end.saturating_add(sequence.len() - query_end + bandwidth / 2);

    if window_end
        .saturating_sub(window_begin)
        .saturating_sub(sequence.len())
        > bandwidth
    {
        let delta = (window_end - window_begin - sequence.len() - bandwidth) / 2;
        window_begin += delta;
        window_end -= delta;
    }

    window_end = window_end.min(reference.len());

    let translated_ref = reference[window_begin..window_end]
        .iter()
        .map(|base| probaln_base_code(*base))
        .collect::<Vec<_>>();
    let translated_seq = sequence
        .iter()
        .map(|base| probaln_base_code(*base))
        .collect::<Vec<_>>();
    let params = ProbalnParams {
        d: if sequence.len() > 1000 { 1e-7 } else { 0.001 },
        e: 0.1,
        bw: bandwidth,
    };
    let result = probaln_glocal(
        &translated_ref,
        &translated_seq,
        Some(qualities),
        params,
        true,
    )?;
    let state = result.state.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "probaln did not return MAP states",
        )
    })?;
    let posterior_quality = result.posterior_quality.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "probaln did not return posterior qualities",
        )
    })?;
    let mut baq = qualities.to_vec();

    if extended {
        calculate_extended_baq_scores(
            cigar,
            sequence,
            window_begin,
            &state,
            &posterior_quality,
            &mut baq,
            alignment_start,
        );
    } else {
        calculate_default_baq_scores(
            cigar,
            sequence,
            window_begin,
            &state,
            &posterior_quality,
            &mut baq,
            alignment_start,
        );
    }

    let baq = if extended {
        qualities
            .iter()
            .zip(baq)
            .map(|(quality, baq)| 64 + quality.saturating_sub(baq))
            .collect::<Vec<_>>()
    } else {
        qualities
            .iter()
            .zip(baq)
            .map(|(quality, baq)| quality - baq + 64)
            .collect::<Vec<_>>()
    };

    String::from_utf8(baq)
        .map(Some)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn calculate_default_baq_scores(
    cigar: &[SamCigarOp],
    sequence: &[u8],
    window_begin: usize,
    state: &[i32],
    posterior_quality: &[u8],
    baq: &mut [u8],
    alignment_start: usize,
) {
    let mut reference_index = alignment_start;
    let mut query_index = 0usize;

    for op in cigar {
        match op.op {
            b'M' | b'=' | b'X' => {
                let len = op.len.min(sequence.len().saturating_sub(query_index));

                for i in query_index..query_index + len {
                    let expected_reference_state =
                        reference_index - window_begin + (i - query_index);

                    if (state[i] & 3) != 0 || (state[i] >> 2) as usize != expected_reference_state {
                        baq[i] = 0;
                    } else {
                        baq[i] = baq[i].min(posterior_quality[i]);
                    }
                }

                reference_index = reference_index.saturating_add(op.len);
                query_index = query_index.saturating_add(op.len);
            }
            b'S' | b'I' => query_index = query_index.saturating_add(op.len),
            b'D' => reference_index = reference_index.saturating_add(op.len),
            b'H' | b'P' => {}
            b'N' => unreachable!("reference skips returned before BAQ calculation"),
            _ => unreachable!("CIGAR op validated by parse_sam_cigar"),
        }
    }
}

fn calculate_extended_baq_scores(
    cigar: &[SamCigarOp],
    sequence: &[u8],
    window_begin: usize,
    state: &[i32],
    posterior_quality: &[u8],
    baq: &mut [u8],
    alignment_start: usize,
) {
    let mut reference_index = alignment_start;
    let mut query_index = 0usize;
    let mut previous_match_len = 0usize;
    let mut left = vec![0; sequence.len()];
    let mut right = vec![0; sequence.len()];

    for (op_index, op) in cigar.iter().enumerate() {
        let is_match = matches!(op.op, b'M' | b'=' | b'X');
        let next_is_match = cigar
            .get(op_index + 1)
            .is_some_and(|next| matches!(next.op, b'M' | b'=' | b'X'));

        if is_match && next_is_match {
            previous_match_len = previous_match_len.saturating_add(op.len);
            continue;
        }

        match op.op {
            b'M' | b'=' | b'X' => {
                let len = op
                    .len
                    .saturating_add(previous_match_len)
                    .min(sequence.len().saturating_sub(query_index));
                previous_match_len = 0;

                if len > 0 {
                    for i in query_index..query_index + len {
                        let expected_reference_state =
                            reference_index - window_begin + (i - query_index);

                        baq[i] = if (state[i] & 3) != 0
                            || (state[i] >> 2) as usize != expected_reference_state
                        {
                            0
                        } else {
                            posterior_quality[i]
                        };
                    }

                    left[query_index] = baq[query_index];
                    for i in query_index + 1..query_index + len {
                        left[i] = baq[i].max(left[i - 1]);
                    }

                    let last = query_index + len - 1;
                    right[last] = baq[last];
                    for i in (query_index..last).rev() {
                        right[i] = baq[i].max(right[i + 1]);
                    }

                    for i in query_index..query_index + len {
                        baq[i] = left[i].min(right[i]);
                    }
                }

                reference_index = reference_index.saturating_add(op.len);
                query_index = query_index.saturating_add(op.len);
            }
            b'S' | b'I' => {
                let len = op.len.min(sequence.len().saturating_sub(query_index));
                query_index = query_index.saturating_add(len);
            }
            b'D' => reference_index = reference_index.saturating_add(op.len),
            b'H' | b'P' => {}
            b'N' => unreachable!("reference skips returned before BAQ calculation"),
            _ => unreachable!("CIGAR op validated by parse_sam_cigar"),
        }
    }
}

fn probaln_base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

fn write_fastq_record<W, R>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    write_fastq_record_with_suffix(writer, record, false)
}

fn write_fastq_record_with_suffix<W, R>(
    writer: &mut W,
    record: &R,
    append_read_number: bool,
) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    let name = fastx_record_name(record)?;
    let name = append_fastx_read_number(name, record, append_read_number)?;
    write_fastq_record_with_name_and_aux(writer, record, &name, None)
}

fn write_fastq_record_with_suffix_and_aux<W, R>(
    writer: &mut W,
    record: &R,
    append_read_number: bool,
    aux_tags: Option<&[[u8; 2]]>,
) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    let name = fastx_record_name(record)?;
    let name = append_fastx_read_number(name, record, append_read_number)?;
    write_fastq_record_with_name_and_aux(writer, record, &name, aux_tags)
}

fn write_fastq_record_with_name_and_aux<W, R>(
    writer: &mut W,
    record: &R,
    name: &str,
    aux_tags: Option<&[[u8; 2]]>,
) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    let sequence = sequence_string(record).unwrap_or_default();
    let quality = fastq_quality_scores_string(record)?;

    if sequence.len() != quality.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTQ quality length differs from sequence length for {name}"),
        ));
    }

    write!(writer, "@{name}")?;
    for field in fastq_aux_fields(record, aux_tags)? {
        write!(writer, "\t{field}")?;
    }
    writeln!(writer)?;
    writeln!(writer, "{sequence}")?;
    writeln!(writer, "+")?;
    writeln!(writer, "{quality}")?;

    Ok(())
}

fn fastq_aux_fields<R>(record: &R, aux_tags: Option<&[[u8; 2]]>) -> io::Result<Vec<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    let Some(aux_tags) = aux_tags else {
        return Ok(Vec::new());
    };

    record
        .data()
        .iter()
        .filter_map(|result| match result {
            Ok((tag, value)) => {
                let tag_bytes = <[u8; 2]>::from(tag);
                aux_tags
                    .iter()
                    .any(|wanted| wanted == &tag_bytes)
                    .then_some((tag_bytes, value))
                    .map(Ok)
            }
            Err(e) => Some(Err(e)),
        })
        .filter_map(|result| match result {
            Ok((tag, value)) => format_fastq_aux_field(&tag, value).map(Ok),
            Err(e) => Some(Err(e)),
        })
        .collect()
}

fn format_fastq_aux_field(
    tag: &[u8; 2],
    value: sam::alignment::record::data::field::Value<'_>,
) -> Option<String> {
    use sam::alignment::record::data::field::Value;

    let tag = std::str::from_utf8(tag).ok()?;
    match value {
        Value::Character(n) => Some(format!("{tag}:A:{}", char::from(n))),
        Value::Int8(n) => Some(format!("{tag}:i:{n}")),
        Value::UInt8(n) => Some(format!("{tag}:i:{n}")),
        Value::Int16(n) => Some(format!("{tag}:i:{n}")),
        Value::UInt16(n) => Some(format!("{tag}:i:{n}")),
        Value::Int32(n) => Some(format!("{tag}:i:{n}")),
        Value::UInt32(n) => Some(format!("{tag}:i:{n}")),
        Value::Float(n) => Some(format!("{tag}:f:{n:e}")),
        Value::String(s) => Some(format!("{tag}:Z:{}", String::from_utf8_lossy(s))),
        Value::Hex(s) => Some(format!("{tag}:H:{}", String::from_utf8_lossy(s))),
        Value::Array(_) => None,
    }
}

#[derive(Default)]
struct FastxSplitBuffers {
    read1: Vec<u8>,
    read2: Vec<u8>,
    singleton: Vec<u8>,
}

impl FastxSplitBuffers {
    fn into_text(self) -> io::Result<FastxSplitText> {
        Ok(FastxSplitText {
            read1: String::from_utf8(self.read1)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            read2: String::from_utf8(self.read2)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            singleton: String::from_utf8(self.singleton)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        })
    }
}

enum FastxSplitTarget {
    Read1,
    Read2,
    Singleton,
}

fn fastx_split_target<R>(record: &R) -> io::Result<FastxSplitTarget>
where
    R: sam::alignment::Record + ?Sized,
{
    let flags = record.flags()?;
    if flags.is_first_segment() {
        Ok(FastxSplitTarget::Read1)
    } else if flags.is_last_segment() {
        Ok(FastxSplitTarget::Read2)
    } else {
        Ok(FastxSplitTarget::Singleton)
    }
}

fn write_split_fastq_record<R>(
    split: &mut FastxSplitBuffers,
    record: &R,
    append_read_number: bool,
) -> io::Result<()>
where
    R: sam::alignment::Record + ?Sized,
{
    match fastx_split_target(record)? {
        FastxSplitTarget::Read1 => {
            write_fastq_record_with_suffix(&mut split.read1, record, append_read_number)
        }
        FastxSplitTarget::Read2 => {
            write_fastq_record_with_suffix(&mut split.read2, record, append_read_number)
        }
        FastxSplitTarget::Singleton => {
            write_fastq_record_with_suffix(&mut split.singleton, record, append_read_number)
        }
    }
}

fn write_split_fastq_record_with_aux<R>(
    split: &mut FastxSplitBuffers,
    record: &R,
    append_read_number: bool,
    aux_tags: Option<&[[u8; 2]]>,
) -> io::Result<()>
where
    R: sam::alignment::Record + ?Sized,
{
    match fastx_split_target(record)? {
        FastxSplitTarget::Read1 => write_fastq_record_with_suffix_and_aux(
            &mut split.read1,
            record,
            append_read_number,
            aux_tags,
        ),
        FastxSplitTarget::Read2 => write_fastq_record_with_suffix_and_aux(
            &mut split.read2,
            record,
            append_read_number,
            aux_tags,
        ),
        FastxSplitTarget::Singleton => write_fastq_record_with_suffix_and_aux(
            &mut split.singleton,
            record,
            append_read_number,
            aux_tags,
        ),
    }
}

fn write_split_fasta_record<R>(
    split: &mut FastxSplitBuffers,
    record: &R,
    append_read_number: bool,
) -> io::Result<()>
where
    R: sam::alignment::Record + ?Sized,
{
    match fastx_split_target(record)? {
        FastxSplitTarget::Read1 => {
            write_fasta_record_with_suffix(&mut split.read1, record, append_read_number)
        }
        FastxSplitTarget::Read2 => {
            write_fasta_record_with_suffix(&mut split.read2, record, append_read_number)
        }
        FastxSplitTarget::Singleton => {
            write_fasta_record_with_suffix(&mut split.singleton, record, append_read_number)
        }
    }
}

fn record_passes_flag_filter<R>(
    record: &R,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
) -> io::Result<bool>
where
    R: sam::alignment::Record + ?Sized,
{
    let flag = record.flags()?.bits();
    Ok(
        (require_flags == 0 || (flag & require_flags) == require_flags)
            && (exclude_flags == 0 || (flag & exclude_flags) == 0)
            && (exclude_all_flags == 0 || (flag & exclude_all_flags) != exclude_all_flags),
    )
}

fn write_fasta_record<W, R>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    write_fasta_record_with_suffix(writer, record, false)
}

fn write_fasta_record_with_suffix<W, R>(
    writer: &mut W,
    record: &R,
    append_read_number: bool,
) -> io::Result<()>
where
    W: Write,
    R: sam::alignment::Record + ?Sized,
{
    let name = fastx_record_name(record)?;
    let name = append_fastx_read_number(name, record, append_read_number)?;
    let sequence = sequence_string(record).unwrap_or_default();

    writeln!(writer, ">{name}")?;
    writeln!(writer, "{sequence}")?;

    Ok(())
}

fn fastx_record_name<R>(record: &R) -> io::Result<String>
where
    R: sam::alignment::Record + ?Sized,
{
    record
        .name()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing FASTX record name"))
}

fn append_fastx_read_number<R>(
    mut name: String,
    record: &R,
    append_read_number: bool,
) -> io::Result<String>
where
    R: sam::alignment::Record + ?Sized,
{
    if append_read_number {
        let flags = record.flags()?;
        if flags.is_first_segment() {
            name.push_str("/1");
        } else if flags.is_last_segment() {
            name.push_str("/2");
        }
    }

    Ok(name)
}

/// Counts SAM records matching an HTSlib-style filter expression.
pub fn count_sam_records_matching_filter_from_path<P>(src: P, filter: &str) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(|reader| count_sam_records_matching_filter(reader, filter))
}

/// Counts SAM records from a buffered reader matching an HTSlib-style filter expression.
pub fn count_sam_records_matching_filter<R>(reader: R, filter: &str) -> io::Result<usize>
where
    R: BufRead,
{
    let mut reader = sam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut count = 0;

    for result in reader.records() {
        let record = result?;
        let context = SamFilterContext::new(&header, &record)?;
        let value = filter
            .eval_with(|symbol| context.lookup(symbol))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        if value.truth() {
            count += 1;
        }
    }

    Ok(count)
}

/// Counts BAM records matching an HTSlib-style filter expression.
pub fn count_bam_records_matching_filter_from_path<P>(src: P, filter: &str) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| count_bam_records_matching_filter(reader, filter))
}

/// Counts BAM records from a reader matching an HTSlib-style filter expression.
pub fn count_bam_records_matching_filter<R>(reader: R, filter: &str) -> io::Result<usize>
where
    R: Read,
{
    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut count = 0;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            count += 1;
        }
    }

    Ok(count)
}

/// Counts CRAM records matching an HTSlib-style filter expression using a FASTA reference.
pub fn count_cram_records_matching_filter_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
    filter: &str,
) -> io::Result<usize>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        count_cram_records_matching_filter_with_reference_repository(
            reader,
            reference_sequence_repository,
            filter,
        )
    })
}

/// Counts CRAM records from a reader matching an HTSlib-style filter expression.
pub fn count_cram_records_matching_filter_with_reference<R, Q>(
    reader: R,
    reference_src: Q,
    filter: &str,
) -> io::Result<usize>
where
    R: Read,
    Q: AsRef<Path>,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    count_cram_records_matching_filter_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
    )
}

fn count_cram_records_matching_filter_with_reference_repository<R>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    filter: &str,
) -> io::Result<usize>
where
    R: Read,
{
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut count = 0;

    for result in reader.records(&header) {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            count += 1;
        }
    }

    Ok(count)
}

/// Writes SAM records matching an HTSlib-style filter expression, including the header.
pub fn view_sam_text_matching_filter_from_path<P>(src: P, filter: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(|reader| view_sam_text_matching_filter(reader, filter))
}

/// Writes SAM records from a buffered reader matching an HTSlib-style filter expression.
pub fn view_sam_text_matching_filter<R>(reader: R, filter: &str) -> io::Result<String>
where
    R: BufRead,
{
    let mut reader = sam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_record(&header, &record)?;
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes BAM records matching an HTSlib-style filter expression as SAM text, including the header.
pub fn view_bam_as_sam_text_matching_filter_from_path<P>(src: P, filter: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| view_bam_as_sam_text_matching_filter(reader, filter))
}

/// Writes BAM records from a reader matching an HTSlib-style filter expression as SAM text.
pub fn view_bam_as_sam_text_matching_filter<R>(reader: R, filter: &str) -> io::Result<String>
where
    R: Read,
{
    use sam::alignment::io::Write as _;

    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes CRAM records matching an HTSlib-style filter expression as SAM text, including the header.
pub fn view_cram_as_sam_text_matching_filter_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
    filter: &str,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        view_cram_as_sam_text_matching_filter_with_reference_repository(
            reader,
            reference_sequence_repository,
            filter,
        )
    })
}

/// Writes CRAM records from a reader matching a filter expression as SAM text.
pub fn view_cram_as_sam_text_matching_filter_with_reference<R, Q>(
    reader: R,
    reference_src: Q,
    filter: &str,
) -> io::Result<String>
where
    R: Read,
    Q: AsRef<Path>,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    view_cram_as_sam_text_matching_filter_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
    )
}

fn view_cram_as_sam_text_matching_filter_with_reference_repository<R>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    filter: &str,
) -> io::Result<String>
where
    R: Read,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn record_matches_filter<R>(header: &Header, record: &R, filter: &Filter) -> io::Result<bool>
where
    R: sam::alignment::Record + ?Sized,
{
    let context = SamFilterContext::new(header, record)?;
    let value = filter
        .eval_with(|symbol| context.lookup(symbol))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    Ok(value.truth())
}

struct SamFilterContext {
    qname: Option<String>,
    rname: Option<String>,
    pos: Option<usize>,
    flag: u16,
    mapq: Option<u8>,
    qlen: usize,
    cigar: Option<String>,
    seq: Option<String>,
    qual: Option<String>,
    library: Option<String>,
    rlen: usize,
    sclen: usize,
    hclen: usize,
    data: Vec<([u8; 2], SamFilterFieldValue)>,
}

enum SamFilterFieldValue {
    Number(f64),
    String(String),
}

impl SamFilterContext {
    fn new<R>(header: &Header, record: &R) -> io::Result<Self>
    where
        R: sam::alignment::Record + ?Sized,
    {
        let qname = record
            .name()
            .map(|name| String::from_utf8_lossy(name).into_owned());
        let rname = record
            .reference_sequence(header)
            .transpose()?
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned());
        let pos = record.alignment_start().transpose()?.map(usize::from);
        let flag = record.flags()?.bits();
        let mapq = record.mapping_quality().transpose()?.map(|mapq| mapq.get());
        let sequence_len = record.sequence().len();
        let qlen = if sequence_len == 0 {
            record.cigar().read_length()?
        } else {
            sequence_len
        };
        let cigar = cigar_string(record)?;
        let seq = sequence_string(record);
        let qual = quality_scores_string(record)?;
        let library = read_group_library(header, record)?;
        let cigar_ops = cigar_ops(record)?;
        let rlen = cigar_ops
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        let sclen = soft_clip_len(&cigar_ops);
        let hclen = hard_clip_len(&cigar_ops);
        let data = data_values(record)?;

        Ok(Self {
            qname,
            rname,
            pos,
            flag,
            mapq,
            qlen,
            cigar,
            seq,
            qual,
            library,
            rlen,
            sclen,
            hclen,
            data,
        })
    }

    fn lookup(&self, symbol: &str) -> Option<(ExprValue, usize)> {
        let (value, len) = if symbol.starts_with("qname") {
            (string_or_undefined(self.qname.as_deref()), 5)
        } else if symbol.starts_with("rname") {
            (string_or_undefined(self.rname.as_deref()), 5)
        } else if symbol.starts_with("cigar") {
            (string_or_undefined(self.cigar.as_deref()), 5)
        } else if symbol.starts_with("seq") {
            (string_or_undefined(self.seq.as_deref()), 3)
        } else if symbol.starts_with("qual") {
            (string_or_undefined(self.qual.as_deref()), 4)
        } else if symbol.starts_with("library") {
            (
                ExprValue::string(self.library.as_deref().unwrap_or_default()),
                7,
            )
        } else if symbol.starts_with("pos") {
            (number_or_undefined(self.pos.map(|n| n as f64)), 3)
        } else if symbol.starts_with("flag") {
            (ExprValue::number(f64::from(self.flag)), 4)
        } else if symbol.starts_with("mapq") {
            (number_or_undefined(self.mapq.map(f64::from)), 4)
        } else if symbol.starts_with("qlen") {
            (ExprValue::number(self.qlen as f64), 4)
        } else if symbol.starts_with("rlen") {
            (ExprValue::number(self.rlen as f64), 4)
        } else if symbol.starts_with("sclen") {
            (ExprValue::number(self.sclen as f64), 5)
        } else if symbol.starts_with("hclen") {
            (ExprValue::number(self.hclen as f64), 5)
        } else if let Some((tag, len)) = parse_bracketed_tag(symbol) {
            (self.aux_value(&tag), len)
        } else {
            return None;
        };
        Some((value, len))
    }

    fn aux_value(&self, tag: &[u8; 2]) -> ExprValue {
        self.data
            .iter()
            .find_map(|(candidate, value)| (candidate == tag).then_some(value))
            .map(|value| match value {
                SamFilterFieldValue::Number(n) => ExprValue::number(*n),
                SamFilterFieldValue::String(s) => ExprValue::string(s),
            })
            .unwrap_or_else(ExprValue::undefined)
    }
}

fn string_or_undefined(value: Option<&str>) -> ExprValue {
    value.map_or_else(ExprValue::undefined, ExprValue::string)
}

fn number_or_undefined(value: Option<f64>) -> ExprValue {
    value.map_or_else(ExprValue::undefined, ExprValue::number)
}

fn cigar_ops<R>(record: &R) -> io::Result<Vec<sam::alignment::record::cigar::Op>>
where
    R: sam::alignment::Record + ?Sized,
{
    record.cigar().iter().collect()
}

fn cigar_string<R>(record: &R) -> io::Result<Option<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    let cigar = record.cigar();

    if cigar.is_empty() {
        return Ok(None);
    }

    let mut s = String::new();

    for result in cigar.iter() {
        let op = result?;
        s.push_str(&op.len().to_string());
        s.push(match op.kind() {
            sam::alignment::record::cigar::op::Kind::Match => 'M',
            sam::alignment::record::cigar::op::Kind::Insertion => 'I',
            sam::alignment::record::cigar::op::Kind::Deletion => 'D',
            sam::alignment::record::cigar::op::Kind::Skip => 'N',
            sam::alignment::record::cigar::op::Kind::SoftClip => 'S',
            sam::alignment::record::cigar::op::Kind::HardClip => 'H',
            sam::alignment::record::cigar::op::Kind::Pad => 'P',
            sam::alignment::record::cigar::op::Kind::SequenceMatch => '=',
            sam::alignment::record::cigar::op::Kind::SequenceMismatch => 'X',
        });
    }

    Ok(Some(s))
}

fn soft_clip_len(ops: &[sam::alignment::record::cigar::Op]) -> usize {
    use sam::alignment::record::cigar::op::Kind;

    let mut len = 0;
    let mut left = 0;

    if matches!(ops.first().map(|op| op.kind()), Some(Kind::SoftClip)) {
        len += ops[0].len();
    } else if ops.len() > 1
        && matches!(ops[0].kind(), Kind::HardClip)
        && matches!(ops[1].kind(), Kind::SoftClip)
    {
        left = 1;
        len += ops[1].len();
    }

    if ops.len() > left + 1 && matches!(ops.last().map(|op| op.kind()), Some(Kind::SoftClip)) {
        len += ops[ops.len() - 1].len();
    } else if ops.len() > left + 2
        && matches!(ops[ops.len() - 1].kind(), Kind::HardClip)
        && matches!(ops[ops.len() - 2].kind(), Kind::SoftClip)
    {
        len += ops[ops.len() - 2].len();
    }

    len
}

fn hard_clip_len(ops: &[sam::alignment::record::cigar::Op]) -> usize {
    use sam::alignment::record::cigar::op::Kind;

    let left = ops
        .first()
        .filter(|op| matches!(op.kind(), Kind::HardClip))
        .map_or(0, |op| op.len());
    let right = ops
        .last()
        .filter(|op| matches!(op.kind(), Kind::HardClip))
        .map_or(0, |op| op.len());

    left + right
}

fn sequence_string<R>(record: &R) -> Option<String>
where
    R: sam::alignment::Record + ?Sized,
{
    let sequence = record.sequence();

    Some(String::from_utf8_lossy(&sequence.iter().collect::<Vec<_>>()).into_owned())
}

fn quality_scores_string<R>(record: &R) -> io::Result<Option<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    let quality_scores = record.quality_scores();

    let scores = quality_scores.iter().collect::<io::Result<Vec<_>>>()?;

    String::from_utf8(scores)
        .map(Some)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn fastq_quality_scores_string<R>(record: &R) -> io::Result<String>
where
    R: sam::alignment::Record + ?Sized,
{
    let scores = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<_>>>()?;
    let bytes = scores
        .into_iter()
        .map(|score| {
            score.checked_add(b'!').ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FASTQ quality score overflow")
            })
        })
        .collect::<io::Result<Vec<_>>>()?;

    String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_group_library<R>(header: &Header, record: &R) -> io::Result<Option<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::data::field::{Tag, Value};
    use sam::header::record::value::map::read_group::tag;

    let data = record.data();
    let Some(result) = data.get(&Tag::READ_GROUP) else {
        return Ok(None);
    };

    let Value::String(read_group_id) = result? else {
        return Ok(None);
    };

    Ok(header
        .read_groups()
        .iter()
        .find(|(id, _)| AsRef::<[u8]>::as_ref(id) == AsRef::<[u8]>::as_ref(read_group_id))
        .and_then(|(_, read_group)| read_group.other_fields().get(&tag::LIBRARY))
        .map(|library| String::from_utf8_lossy(library).into_owned()))
}

fn data_values<R>(record: &R) -> io::Result<Vec<([u8; 2], SamFilterFieldValue)>>
where
    R: sam::alignment::Record + ?Sized,
{
    record
        .data()
        .iter()
        .filter_map(|result| match result {
            Ok((tag, value)) => sam_filter_field_value(value).map(|value| Ok((tag.into(), value))),
            Err(e) => Some(Err(e)),
        })
        .collect()
}

fn sam_filter_field_value(
    value: sam::alignment::record::data::field::Value<'_>,
) -> Option<SamFilterFieldValue> {
    use sam::alignment::record::data::field::Value;

    match value {
        Value::Character(b) => Some(SamFilterFieldValue::String(char::from(b).to_string())),
        Value::Int8(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::UInt8(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::Int16(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::UInt16(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::Int32(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::UInt32(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::Float(n) => Some(SamFilterFieldValue::Number(f64::from(n))),
        Value::String(s) | Value::Hex(s) => Some(SamFilterFieldValue::String(
            String::from_utf8_lossy(s).into_owned(),
        )),
        Value::Array(_) => None,
    }
}

fn parse_bracketed_tag(symbol: &str) -> Option<([u8; 2], usize)> {
    let bytes = symbol.as_bytes();

    (bytes.len() >= 4 && bytes[0] == b'[' && bytes[3] == b']').then(|| ([bytes[1], bytes[2]], 4))
}

/// Reports base modifications in the style of HTSlib's `test_mod` helper.
pub fn base_modification_report_from_sam_path<P>(
    src: P,
    report_unchecked: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut report = String::new();

    for result in reader.records() {
        let record = result?;
        push_base_modification_report(&mut report, &header, &record, report_unchecked, false)?;
    }

    Ok(report)
}

/// Reports base modifications with extended type metadata in the style of HTSlib's `test_mod -x`.
pub fn extended_base_modification_report_from_sam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut report = String::new();

    for result in reader.records() {
        let record = result?;
        push_base_modification_report(&mut report, &header, &record, false, true)?;
    }

    Ok(report)
}

/// Reports a SAM pileup with base modifications in the style of HTSlib's `pileup_mod`.
pub fn base_modification_pileup_report_from_sam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut pileup_records = Vec::new();

    for result in reader.records() {
        let record = result?;
        if let Some(record) = PileupRecord::try_from_record(&header, &record)? {
            pileup_records.push(record);
        }
    }

    push_base_modification_pileup_report(&pileup_records)
}

/// Reports pileup columns in the style of HTSlib's `test/pileup.c` helper.
pub fn pileup_report_from_sam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result?;
        if let Some(record) = TestPileupRecord::try_from_record(&header, &record)? {
            records.push(record);
        }
    }

    push_test_pileup_report(&records)
}

/// Reports BAM pileup columns in the style of HTSlib's `test/pileup.c` helper.
pub fn pileup_report_from_bam_path<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result?;
        if let Some(record) = TestPileupRecord::try_from_record(&header, &record)? {
            records.push(record);
        }
    }

    push_test_pileup_report(&records)
}

/// Builds synchronized pileup columns across multiple SAM/BAM inputs.
pub fn synchronized_pileup_from_alignment_paths<P>(
    paths: &[P],
) -> io::Result<Vec<SynchronizedPileupColumn>>
where
    P: AsRef<Path>,
{
    let inputs = paths
        .iter()
        .map(read_test_pileup_records_from_alignment_path)
        .collect::<io::Result<Vec<_>>>()?;
    let mut sites: BTreeMap<SynchronizedPileupSite, Vec<SynchronizedPileupEntry>> = BTreeMap::new();

    for (input_index, records) in inputs.iter().enumerate() {
        for (record_index, record) in records.iter().enumerate() {
            for (column_index, column) in record.columns.iter().enumerate() {
                sites
                    .entry((record.reference_name.clone(), column.reference_position))
                    .or_default()
                    .push((input_index, record_index, column_index));
            }
        }
    }

    sites
        .into_iter()
        .map(|((reference_name, zero_based_position), entries)| {
            let mut depths_by_input = vec![0; inputs.len()];
            let mut bases_by_input = vec![String::new(); inputs.len()];
            let mut qualities_by_input = vec![String::new(); inputs.len()];

            for input_index in 0..inputs.len() {
                let input_entries = entries
                    .iter()
                    .filter(|(entry_input_index, _, _)| *entry_input_index == input_index)
                    .map(|(_, record_index, column_index)| {
                        let record = &inputs[input_index][*record_index];
                        let column = &record.columns[*column_index];
                        (*record_index, record, column)
                    })
                    .collect::<Vec<_>>();

                depths_by_input[input_index] = input_entries.len();
                let adjusted_qualities = adjusted_test_pileup_qualities(&input_entries)?;

                for ((_, record, column), quality) in input_entries.iter().zip(adjusted_qualities) {
                    push_test_pileup_bases(&mut bases_by_input[input_index], record, column)?;
                    push_test_pileup_quality_score(&mut qualities_by_input[input_index], quality);
                }
            }

            Ok(SynchronizedPileupColumn {
                reference_name,
                position: zero_based_position + 1,
                total_depth: depths_by_input.iter().sum(),
                depths_by_input,
                bases_by_input,
                qualities_by_input,
            })
        })
        .collect()
}

fn read_test_pileup_records_from_alignment_path<P>(src: P) -> io::Result<Vec<TestPileupRecord>>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();
    if src
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
    {
        let mut reader = File::open(src).map(bam::io::Reader::new)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result?;
            if let Some(record) = TestPileupRecord::try_from_record(&header, &record)? {
                records.push(record);
            }
        }

        Ok(records)
    } else {
        let mut reader = File::open(src)
            .map(BufReader::new)
            .map(sam::io::Reader::new)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result?;
            if let Some(record) = TestPileupRecord::try_from_record(&header, &record)? {
                records.push(record);
            }
        }

        Ok(records)
    }
}

/// Queries BGZF-compressed SAM records from a local file using its associated BAI or CSI index.
pub fn query_sam_records_from_path<P>(src: P, region: &Region) -> io::Result<Vec<sam::Record>>
where
    P: AsRef<Path>,
{
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = sam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let query = reader.query(&header, region)?;

    query.records().collect()
}

/// Writes a SAM text view of indexed BGZF-compressed SAM records for regions in request order.
pub fn view_sam_regions_as_text_from_path<P>(src: P, regions: &[Region]) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let index = read_associated_bam_index(&src)?;
    let raw_header = read_raw_bgzf_sam_header_text(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = sam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(raw_header.into_bytes());
    let mut seen_records = HashSet::new();

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;
            let mut line_writer = sam::io::Writer::new(Vec::new());
            line_writer.write_record(&header, &record)?;
            let line = line_writer.into_inner();

            if seen_records.insert(line.clone()) {
                writer.get_mut().extend(line);
            }
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_raw_bgzf_sam_header_text<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    let mut reader = File::open(data_path)
        .map(bgzf::io::Reader::new)
        .map(BufReader::new)?;
    let mut raw_header = String::new();
    let mut line = String::new();

    loop {
        line.clear();

        if reader.read_line(&mut line)? == 0 || !line.starts_with('@') {
            break;
        }

        raw_header.push_str(&line);
    }

    Ok(raw_header)
}

/// Reads a BAM header from a reader.
pub fn read_bam_header<R>(reader: R) -> io::Result<Header>
where
    R: Read,
{
    let mut reader = bam::io::Reader::new(reader);

    reader.read_header()
}

/// Reads a BAM header from a local file.
pub fn read_bam_header_from_path<P>(src: P) -> io::Result<Header>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(read_bam_header)
}

/// Counts BAM records from a reader.
pub fn count_bam_records<R>(reader: R) -> io::Result<usize>
where
    R: Read,
{
    let mut reader = bam::io::Reader::new(reader);
    reader.read_header()?;

    reader
        .records()
        .try_fold(0, |n, result| result.map(|_| n + 1))
}

/// Counts BAM records from a local file.
pub fn count_bam_records_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(count_bam_records)
}

/// Reads BAM input without producing output and returns the number of records seen.
pub fn benchmark_bam_view_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    count_bam_records_from_path(src)
}

/// Reads BAM records from a local file into format-neutral summaries.
pub fn summarize_bam_records_from_path<P>(src: P) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    reader
        .records()
        .map(|result| result.and_then(|record| summarize_alignment_record(&header, &record)))
        .collect()
}

/// Writes a FASTQ view of the first `limit` BAM records.
pub fn view_bam_as_fastq_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fastq_record(&mut writer, &record)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of the first `limit` BAM records, optionally appending read numbers.
pub fn view_bam_as_fastq_text_from_path_with_limit_and_suffix<P>(
    src: P,
    limit: Option<usize>,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fastq_record_with_suffix(&mut writer, &record, append_read_number)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of BAM records passing flag filters.
pub fn view_bam_as_fastq_text_from_path_with_flag_filter<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fastq_record(&mut writer, &record)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTQ view of filtered BAM records, optionally appending read numbers.
pub fn view_bam_as_fastq_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fastq_record_with_suffix(&mut writer, &record, append_read_number)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of the first `limit` BAM records.
pub fn view_bam_as_fasta_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fasta_record(&mut writer, &record)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of the first `limit` BAM records, optionally appending read numbers.
pub fn view_bam_as_fasta_text_from_path_with_limit_and_suffix<P>(
    src: P,
    limit: Option<usize>,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_fasta_record_with_suffix(&mut writer, &record, append_read_number)?;
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of BAM records passing flag filters.
pub fn view_bam_as_fasta_text_from_path_with_flag_filter<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fasta_record(&mut writer, &record)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a FASTA view of filtered BAM records, optionally appending read numbers.
pub fn view_bam_as_fasta_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_fasta_record_with_suffix(&mut writer, &record, append_read_number)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes split FASTQ views of BAM records passing flag filters.
pub fn view_bam_as_fastq_split_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut split = FastxSplitBuffers::default();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_split_fastq_record(&mut split, &record, append_read_number)?;
        }
    }

    split.into_text()
}

/// Writes split FASTA views of BAM records passing flag filters.
pub fn view_bam_as_fasta_split_text_from_path_with_flag_filter_and_suffix<P>(
    src: P,
    require_flags: u16,
    exclude_flags: u16,
    exclude_all_flags: u16,
    append_read_number: bool,
) -> io::Result<FastxSplitText>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;
    let mut split = FastxSplitBuffers::default();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, require_flags, exclude_flags, exclude_all_flags)? {
            write_split_fasta_record(&mut split, &record, append_read_number)?;
        }
    }

    split.into_text()
}

/// Writes BAM input as BAM, including the header and all records.
pub fn write_bam_from_path<P, W>(src: P, dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    File::open(src).and_then(|reader| write_bam(reader, dst))
}

/// Writes BAM input from a reader as BAM, including the header and all records.
pub fn write_bam<R, W>(reader: R, dst: W) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes BAM input records matching an HTSlib-style filter expression to BAM output.
pub fn write_bam_matching_filter_from_path<P, W>(src: P, filter: &str, dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| write_bam_matching_filter(reader, filter, dst))
}

/// Writes BAM records from a reader matching an HTSlib-style filter expression to BAM output.
pub fn write_bam_matching_filter<R, W>(reader: R, filter: &str, dst: W) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes BAM input as CRAM using a FASTA reference.
pub fn write_cram_from_bam_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    File::open(src).and_then(|reader| {
        write_cram_from_bam_reader_with_reference_repository(
            reader,
            reference_sequence_repository,
            writer,
        )
    })
}

/// Writes BAM input from a reader as CRAM using a FASTA reference.
pub fn write_cram_from_bam_reader_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_from_bam_reader_with_reference_repository(
        reader,
        reference_sequence_repository,
        writer,
    )
}

fn write_cram_from_bam_reader_with_reference_repository<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes BAM input records matching an HTSlib-style filter expression to CRAM output.
pub fn write_cram_matching_filter_from_bam_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    File::open(src).and_then(|reader| {
        write_cram_matching_filter_from_bam_reader_with_reference_repository(
            reader,
            reference_sequence_repository,
            filter,
            writer,
        )
    })
}

/// Writes BAM records from a reader matching a filter expression to CRAM output.
pub fn write_cram_matching_filter_from_bam_reader_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_matching_filter_from_bam_reader_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
        writer,
    )
}

fn write_cram_matching_filter_from_bam_reader_with_reference_repository<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes indexed BAM records overlapping the given regions to CRAM output using a FASTA reference.
pub fn write_bam_regions_as_cram_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes indexed BAM records overlapping the given regions and matching a filter to CRAM output.
pub fn write_bam_regions_matching_filter_as_cram_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                writer.write_alignment_record(&header, &record)?;
            }
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes indexed BAM records overlapping the given regions to BAM output.
///
/// Records are emitted in request order and duplicates are preserved across
/// overlapping or repeated regions, matching the behavior needed by
/// `samtools view -P -b`.
pub fn write_bam_regions_from_path<P, W>(src: P, regions: &[Region], dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;
            writer.write_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes indexed BAM records overlapping the given regions and matching a filter to BAM output.
pub fn write_bam_regions_matching_filter_from_path<P, W>(
    src: P,
    regions: &[Region],
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                writer.write_record(&header, &record)?;
            }
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes BAM input records with all required flag bits set to BAM output.
pub fn write_bam_records_with_required_flags_from_path<P, W>(
    src: P,
    required_flags: u16,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        let flags = u16::from(record.flags());
        if flags & required_flags == required_flags {
            writer.write_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes SAM input as BAM, including the header and all records.
pub fn write_bam_from_sam_path<P, W>(src: P, dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    write_bam_from_sam_path_with_compression_level(
        src,
        dst,
        bgzf::io::writer::CompressionLevel::default(),
    )
}

/// Writes SAM input records matching an HTSlib-style filter expression to BAM output.
pub fn write_bam_matching_filter_from_sam_path<P, W>(src: P, filter: &str, dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    write_bam_matching_filter_from_sam_reader(&mut reader, filter, dst)
}

pub fn write_bam_matching_filter_from_sam_reader<R, W>(
    reader: &mut sam::io::Reader<R>,
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let bgzf_writer = bgzf::io::writer::Builder::default().build_from_writer(dst);
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes SAM input as BAM using an explicit BGZF compression level.
pub fn write_bam_from_sam_path_with_compression_level<P, W>(
    src: P,
    dst: W,
    compression_level: bgzf::io::writer::CompressionLevel,
) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    write_bam_from_sam_reader_with_compression_level(&mut reader, dst, compression_level)
}

/// Writes SAM input as BAM from any buffered reader, including the header and all records.
pub fn write_bam_from_sam_reader<R, W>(reader: R, dst: W) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    let mut reader = sam::io::Reader::new(reader);
    write_bam_from_sam_reader_with_compression_level(
        &mut reader,
        dst,
        bgzf::io::writer::CompressionLevel::default(),
    )
}

/// Writes SAM input as BAM from any buffered reader using an explicit BGZF compression level.
pub fn write_bam_from_sam_reader_with_compression_level<R, W>(
    reader: &mut sam::io::Reader<R>,
    dst: W,
    compression_level: bgzf::io::writer::CompressionLevel,
) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let header = reader.read_header()?;
    let bgzf_writer = bgzf::io::writer::Builder::default()
        .set_compression_level(compression_level)
        .build_from_writer(dst);
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes SAM input as BAM and writes a BAI index for the generated BAM file.
pub fn write_bam_from_sam_path_with_bai<P, Q, R>(src: P, bam_dst: Q, bai_dst: R) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let encoded = write_bam_from_sam_path(src, Vec::new())?;

    std::fs::write(&bam_dst, encoded)?;

    let index = build_bai(&bam_dst)?;
    write_bai(bai_dst, &index)
}

/// Adds NUL padding to the raw SAM header inside a BGZF-compressed BAM stream.
pub fn add_bam_header_nul_padding(data: &[u8], extra_nuls: usize) -> io::Result<Vec<u8>> {
    let mut uncompressed = crate::bgzf_compat::read_all(data)?;

    add_uncompressed_bam_header_nul_padding(&mut uncompressed, extra_nuls)?;

    crate::bgzf_compat::write_all(Vec::new(), &uncompressed)
}

fn add_uncompressed_bam_header_nul_padding(
    data: &mut Vec<u8>,
    extra_nuls: usize,
) -> io::Result<()> {
    const MAGIC_NUMBER_LEN: usize = 4;
    const HEADER_LEN_OFFSET: usize = MAGIC_NUMBER_LEN;
    const HEADER_TEXT_OFFSET: usize = HEADER_LEN_OFFSET + std::mem::size_of::<i32>();
    const MAGIC_NUMBER: &[u8; 4] = b"BAM\x01";

    if extra_nuls == 0 {
        return Ok(());
    }

    if data.len() < HEADER_TEXT_OFFSET || &data[..MAGIC_NUMBER_LEN] != MAGIC_NUMBER {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ));
    }

    let l_text = i32::from_le_bytes(
        data[HEADER_LEN_OFFSET..HEADER_TEXT_OFFSET]
            .try_into()
            .expect("slice length is fixed"),
    );

    if l_text < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "negative BAM header length",
        ));
    }

    let header_len =
        usize::try_from(l_text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let insert_at = HEADER_TEXT_OFFSET
        .checked_add(header_len)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BAM header length overflow"))?;

    if insert_at > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "BAM header is truncated",
        ));
    }

    let new_header_len = header_len
        .checked_add(extra_nuls)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "BAM header length overflow"))?;
    let new_l_text = i32::try_from(new_header_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    data[HEADER_LEN_OFFSET..HEADER_TEXT_OFFSET].copy_from_slice(&new_l_text.to_le_bytes());
    data.splice(insert_at..insert_at, std::iter::repeat_n(0, extra_nuls));

    Ok(())
}

/// Writes a SAM text view of BAM records, including the header.
pub fn view_bam_as_sam_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let data_path = associated_data_path(&src);
    let mut reader = File::open(data_path).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let raw_header = read_raw_bam_header_text_for_view(&src, &header)?;
    let mut writer = sam::io::Writer::new(raw_header.into_bytes());

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a SAM text view of BAM records from a reader with an optional record limit.
pub fn view_bam_as_sam_text<R>(reader: R, limit: Option<usize>) -> io::Result<String>
where
    R: Read,
{
    use sam::alignment::io::Write as _;

    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Queries BAM records from a local file using its associated BAI or CSI index.
pub fn query_bam_records_from_path<P>(src: P, region: &Region) -> io::Result<Vec<bam::Record>>
where
    P: AsRef<Path>,
{
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let query = reader.query(&header, region)?;

    query.records().collect()
}

/// Returns an owning iterator over BAM records from an indexed local file.
pub fn iter_bam_records_from_path<P>(
    src: P,
    region: &Region,
) -> io::Result<std::vec::IntoIter<bam::Record>>
where
    P: AsRef<Path>,
{
    query_bam_records_from_path(src, region).map(Vec::into_iter)
}

/// Queries BAM records for multiple regions, preserving region order and duplicates.
pub fn query_bam_regions_from_path<P>(src: P, regions: &[Region]) -> io::Result<Vec<bam::Record>>
where
    P: AsRef<Path>,
{
    let mut records = Vec::new();

    for region in regions {
        records.extend(query_bam_records_from_path(&src, region)?);
    }

    Ok(records)
}

/// Writes a SAM text view of indexed BAM records for regions in request order.
pub fn view_bam_regions_as_sam_text_from_path<P>(src: P, regions: &[Region]) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_bam_regions_as_sam_text_from_path_with_dedup(src, regions, false)
}

/// Writes a SAM text view of indexed BAM records for regions in request order.
pub fn view_bam_regions_as_sam_text_from_path_with_dedup<P>(
    src: P,
    regions: &[Region],
    deduplicate: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_bam_regions_as_sam_text_from_path_with_options(src, regions, deduplicate, None, None)
}

/// Writes a SAM text view of indexed BAM records with an optional record limit.
pub fn view_bam_regions_as_sam_text_from_path_with_limit<P>(
    src: P,
    regions: &[Region],
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_bam_regions_as_sam_text_from_path_with_options(src, regions, false, limit, None)
}

/// Writes a SAM text view of indexed BAM records matching an HTSlib-style filter expression.
pub fn view_bam_regions_as_sam_text_matching_filter_from_path<P>(
    src: P,
    regions: &[Region],
    filter: &str,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let filter = Filter::new(filter);

    view_bam_regions_as_sam_text_from_path_with_options(src, regions, false, None, Some(&filter))
}

fn view_bam_regions_as_sam_text_from_path_with_options<P>(
    src: P,
    regions: &[Region],
    deduplicate: bool,
    limit: Option<usize>,
    filter: Option<&Filter>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut raw_header = read_raw_bam_header_text_for_view(&src, &header)?;

    if deduplicate {
        raw_header = trim_header_line_trailing_spaces(&raw_header);
    }

    let mut writer = sam::io::Writer::new(raw_header.into_bytes());
    let mut seen_records = HashSet::new();
    let mut remaining = limit.unwrap_or(usize::MAX);

    for region in regions {
        if remaining == 0 {
            break;
        }

        let query = reader.query(&header, region)?;

        for result in query.records() {
            if remaining == 0 {
                break;
            }

            let record = result?;

            if let Some(filter) = filter
                && !record_matches_filter(&header, &record, filter)?
            {
                continue;
            }

            let mut line_writer = sam::io::Writer::new(Vec::new());
            line_writer.write_alignment_record(&header, &record)?;
            let line = line_writer.into_inner();

            if !deduplicate || seen_records.insert(line.clone()) {
                writer.get_mut().extend(line);
                remaining -= 1;
            }
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_raw_bam_header_text<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    let mut reader = File::open(data_path).map(bam::io::Reader::new)?;
    let mut header_reader = reader.header_reader();

    header_reader.read_magic_number()?;

    let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;
    let mut raw_header = String::new();
    raw_sam_header_reader.read_to_string(&mut raw_header)?;
    raw_sam_header_reader.discard_to_end()?;

    Ok(raw_header)
}

fn read_raw_bam_header_text_for_view<P>(src: P, header: &Header) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut raw_header = read_raw_bam_header_text(src)?;

    if !raw_header.lines().any(|line| line.starts_with("@SQ\t")) {
        for (name, reference_sequence) in header.reference_sequences() {
            raw_header.push_str(&format!(
                "@SQ\tSN:{name}\tLN:{}\n",
                usize::from(reference_sequence.length())
            ));
        }
    }

    Ok(raw_header)
}

/// Counts BAM records from a local indexed file that intersect a region.
pub fn count_bam_records_in_region_from_path<P>(src: P, region: &Region) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    query_bam_records_from_path(src, region).map(|records| records.len())
}

/// Counts indexed BAM records overlapping the given regions and matching an HTSlib-style filter expression.
pub fn count_bam_records_in_regions_matching_filter_from_path<P>(
    src: P,
    regions: &[Region],
    filter: &str,
) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    let filter = Filter::new(filter);
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut count = 0;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                count += 1;
            }
        }
    }

    Ok(count)
}

/// Reads a CRAM header from a reader.
pub fn read_cram_header<R>(reader: R) -> io::Result<Header>
where
    R: Read,
{
    let mut reader = cram::io::Reader::new(reader);

    reader.read_header()
}

/// Reads a CRAM header from a local file.
pub fn read_cram_header_from_path<P>(src: P) -> io::Result<Header>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(read_cram_header)
}

/// Reads CRAM records from a local file into format-neutral summaries.
pub fn summarize_cram_records_from_path<P>(src: P) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(cram::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut records: Vec<_> = reader
        .records(&header)
        .map(|result| result.and_then(|record| summarize_alignment_record(&header, &record)))
        .collect::<io::Result<_>>()?;

    apply_htslib_cram_template_lengths(&mut records);

    Ok(records)
}

/// Reads CRAM records from a local file with a FASTA reference into format-neutral summaries.
pub fn summarize_cram_records_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    summarize_cram_records_from_path_with_reference_repository(src, reference_sequence_repository)
}

/// Reads CRAM input without producing output and returns the number of records seen.
pub fn benchmark_cram_view_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
) -> io::Result<usize>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    summarize_cram_records_from_path_with_reference(src, reference_src).map(|records| records.len())
}

fn summarize_cram_records_from_path_with_reference_repository<P>(
    src: P,
    reference_sequence_repository: fasta::Repository,
) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;

    let mut records: Vec<_> = reader
        .records(&header)
        .map(|result| result.and_then(|record| summarize_alignment_record(&header, &record)))
        .collect::<io::Result<_>>()?;

    apply_htslib_cram_template_lengths(&mut records);

    Ok(records)
}

/// Writes CRAM records decoded from a local CRAM file using a FASTA reference.
pub fn write_cram_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        write_cram_from_reader_with_reference_repository(
            reader,
            reference_sequence_repository,
            writer,
        )
    })
}

/// Writes CRAM records decoded from a reader using a FASTA reference.
pub fn write_cram_from_reader_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_from_reader_with_reference_repository(reader, reference_sequence_repository, writer)
}

fn write_cram_from_reader_with_reference_repository<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository.clone())
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes CRAM input records matching an HTSlib-style filter expression to CRAM output.
pub fn write_cram_matching_filter_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        write_cram_matching_filter_from_reader_with_reference_repository(
            reader,
            reference_sequence_repository,
            filter,
            writer,
        )
    })
}

/// Writes CRAM records from a reader matching a filter expression to CRAM output.
pub fn write_cram_matching_filter_from_reader_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_matching_filter_from_reader_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
        writer,
    )
}

fn write_cram_matching_filter_from_reader_with_reference_repository<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository.clone())
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes indexed CRAM records overlapping the given regions to BAM output.
pub fn write_cram_regions_as_bam_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes indexed CRAM records overlapping the given regions and matching a filter to BAM output.
pub fn write_cram_regions_matching_filter_as_bam_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                writer.write_alignment_record(&header, &record)?;
            }
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes indexed CRAM records overlapping the given regions to CRAM output using a FASTA reference.
pub fn write_cram_regions_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_writer(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes indexed CRAM records overlapping the given regions and matching a filter to CRAM output.
pub fn write_cram_regions_matching_filter_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_writer(dst);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                writer.write_alignment_record(&header, &record)?;
            }
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes CRAM records with all required flag bits set to BAM output.
pub fn write_cram_records_with_required_flags_as_bam_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    required_flags: u16,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        write_cram_records_with_required_flags_as_bam_with_reference_repository(
            reader,
            repository,
            required_flags,
            dst,
        )
    })
}

/// Writes CRAM records from a reader with all required flag bits set to BAM output.
pub fn write_cram_records_with_required_flags_as_bam_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    required_flags: u16,
    dst: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_records_with_required_flags_as_bam_with_reference_repository(
        reader,
        repository,
        required_flags,
        dst,
    )
}

fn write_cram_records_with_required_flags_as_bam_with_reference_repository<R, W>(
    reader: R,
    repository: fasta::Repository,
    required_flags: u16,
    dst: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags().bits();
        if flags & required_flags == required_flags {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes CRAM records matching an HTSlib-style filter expression to BAM output.
pub fn write_cram_records_matching_filter_as_bam_from_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        write_cram_records_matching_filter_as_bam_with_reference_repository(
            reader, repository, filter, dst,
        )
    })
}

/// Writes CRAM records from a reader matching a filter expression to BAM output.
pub fn write_cram_records_matching_filter_as_bam_with_reference<R, Q, W>(
    reader: R,
    reference_src: Q,
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    R: Read,
    Q: AsRef<Path>,
    W: Write,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_records_matching_filter_as_bam_with_reference_repository(
        reader, repository, filter, dst,
    )
}

fn write_cram_records_matching_filter_as_bam_with_reference_repository<R, W>(
    reader: R,
    repository: fasta::Repository,
    filter: &str,
    dst: W,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = bam::io::Writer::new(dst);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes CRAM records decoded from a local SAM file using a FASTA reference.
pub fn write_cram_from_sam_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    write_cram_from_sam_reader_with_reference_repository(
        &mut reader,
        reference_sequence_repository,
        writer,
    )
}

/// Writes SAM input from a reader as CRAM using a FASTA reference.
pub fn write_cram_from_sam_reader_with_reference<R, Q, W>(
    reader: &mut sam::io::Reader<R>,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    R: BufRead,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_from_sam_reader_with_reference_repository(
        reader,
        reference_sequence_repository,
        writer,
    )
}

fn write_cram_from_sam_reader_with_reference_repository<R, W>(
    reader: &mut sam::io::Reader<R>,
    reference_sequence_repository: fasta::Repository,
    writer: W,
) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes SAM input records matching an HTSlib-style filter expression to CRAM output.
pub fn write_cram_matching_filter_from_sam_path_with_reference<P, Q, W>(
    src: P,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    write_cram_matching_filter_from_sam_reader_with_reference_repository(
        &mut reader,
        reference_sequence_repository,
        filter,
        writer,
    )
}

/// Writes SAM input records from a reader matching a filter expression to CRAM output.
pub fn write_cram_matching_filter_from_sam_reader_with_reference<R, Q, W>(
    reader: &mut sam::io::Reader<R>,
    reference_src: Q,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: BufRead,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_matching_filter_from_sam_reader_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
        writer,
    )
}

fn write_cram_matching_filter_from_sam_reader_with_reference_repository<R, W>(
    reader: &mut sam::io::Reader<R>,
    reference_sequence_repository: fasta::Repository,
    filter: &str,
    writer: W,
) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .build_from_writer(writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if record_matches_filter(&header, &record, &filter)? {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

/// Writes SAM input as CRAM and writes a CRAI index for the generated CRAM file.
pub fn write_cram_from_sam_path_with_reference_and_crai<P, Q, R, S>(
    src: P,
    reference_src: Q,
    cram_dst: R,
    crai_dst: S,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
    S: AsRef<Path>,
{
    let encoded = write_cram_from_sam_path_with_reference(src, reference_src, Vec::new())?;

    std::fs::write(&cram_dst, encoded)?;

    let index = build_cram_crai(&cram_dst)?;
    write_cram_crai(crai_dst, &index)
}

/// Queries CRAM records from a local file using its associated CRAI index.
pub fn query_cram_records_from_path<P>(
    src: P,
    region: &Region,
) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
{
    query_cram_records_from_path_with_reference_repository(
        src,
        region,
        fasta::Repository::default(),
    )
}

/// Builds a noodles FASTA reference repository for CRAM decoding from a local FASTA file.
pub fn cram_reference_repository_from_fasta_path<P>(
    reference_src: P,
) -> io::Result<fasta::Repository>
where
    P: AsRef<Path>,
{
    let reference_src = reference_src.as_ref();
    let primary_index_src = reference_src.with_extension("fa.fai");
    let fallback_index_src = append_fai_extension(reference_src);
    let reference_index = fasta::fai::fs::read(&primary_index_src).or_else(|primary_err| {
        fasta::fai::fs::read(&fallback_index_src).map_err(|fallback_err| {
            io_error_with_path(
                io::ErrorKind::NotFound,
                "read FASTA index",
                &fallback_index_src,
                io::Error::new(
                    fallback_err.kind(),
                    format!(
                        "{}; also failed to read {}: {}",
                        fallback_err,
                        primary_index_src.display(),
                        primary_err
                    ),
                ),
            )
        })
    })?;
    let reference_reader = File::open(reference_src)
        .map(BufReader::new)
        .with_io_path_context("open FASTA reference", reference_src)?;
    let reference_reader = fasta::io::IndexedReader::new(reference_reader, reference_index);

    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reference_reader),
    ))
}

/// Queries CRAM records from a local file using its associated CRAI index and a FASTA reference.
pub fn query_cram_records_from_path_with_reference<P, Q>(
    src: P,
    region: &Region,
    reference_src: Q,
) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;

    query_cram_records_from_path_with_reference_repository(src, region, repository)
}

/// Returns an owning iterator over CRAM records from an indexed local file and FASTA reference.
pub fn iter_cram_records_from_path_with_reference<P, Q>(
    src: P,
    region: &Region,
    reference_src: Q,
) -> io::Result<std::vec::IntoIter<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    query_cram_records_from_path_with_reference(src, region, reference_src).map(Vec::into_iter)
}

/// Counts indexed CRAM records overlapping the given regions and matching an HTSlib-style filter expression.
pub fn count_cram_records_in_regions_matching_filter_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    filter: &str,
) -> io::Result<usize>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let filter = Filter::new(filter);
    let mut count = 0;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;

            if record_matches_filter(&header, &record, &filter)? {
                count += 1;
            }
        }
    }

    Ok(count)
}

/// Writes a SAM text view of indexed CRAM records for regions in request order.
pub fn view_cram_regions_as_sam_text_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    deduplicate: bool,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequences = read_fasta_sequences(&reference_src)?;
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;

    view_cram_regions_as_sam_text_from_path_with_reference_repository(
        src,
        repository,
        &reference_sequences,
        regions,
        deduplicate,
        None,
        None,
    )
}

/// Writes a SAM text view of indexed CRAM records with an optional record limit.
pub fn view_cram_regions_as_sam_text_from_path_with_reference_and_limit<P, Q>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequences = read_fasta_sequences(&reference_src)?;
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;

    view_cram_regions_as_sam_text_from_path_with_reference_repository(
        src,
        repository,
        &reference_sequences,
        regions,
        false,
        limit,
        None,
    )
}

/// Writes a SAM text view of CRAM records with an optional record limit.
pub fn view_cram_as_sam_text_from_path_with_reference_and_limit<P, Q>(
    src: P,
    reference_src: Q,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let reference_sequences = read_fasta_sequences(&reference_src)?;
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let raw_header = read_raw_cram_header_text(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(raw_header.into_bytes());

    for result in reader.records(&header).take(limit.unwrap_or(usize::MAX)) {
        let mut record = result?;

        add_md_and_nm_to_record(&reference_sequences, &header, &mut record)?;
        writer.write_alignment_record(&header, &record)?;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a SAM text view of CRAM records from a reader with an optional record limit.
pub fn view_cram_as_sam_text_with_reference<R, Q>(
    reader: R,
    reference_src: Q,
    limit: Option<usize>,
) -> io::Result<String>
where
    R: Read,
    Q: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let reference_sequences = read_fasta_sequences(&reference_src)?;
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(Vec::new());

    writer.write_header(&header)?;

    for result in reader.records(&header).take(limit.unwrap_or(usize::MAX)) {
        let mut record = result?;

        add_md_and_nm_to_record(&reference_sequences, &header, &mut record)?;
        writer.write_alignment_record(&header, &record)?;
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes a SAM text view of indexed CRAM records matching an HTSlib-style filter expression.
pub fn view_cram_regions_as_sam_text_matching_filter_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    filter: &str,
) -> io::Result<String>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reference_sequences = read_fasta_sequences(&reference_src)?;
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let filter = Filter::new(filter);

    view_cram_regions_as_sam_text_from_path_with_reference_repository(
        src,
        repository,
        &reference_sequences,
        regions,
        false,
        None,
        Some(&filter),
    )
}

fn query_cram_records_from_path_with_reference_repository<P>(
    src: P,
    region: &Region,
    reference_sequence_repository: fasta::Repository,
) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
{
    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let query = reader.query(&header, region)?;

    query.collect()
}

fn view_cram_regions_as_sam_text_from_path_with_reference_repository<P>(
    src: P,
    reference_sequence_repository: fasta::Repository,
    reference_sequences: &HashMap<Vec<u8>, Vec<u8>>,
    regions: &[Region],
    deduplicate: bool,
    limit: Option<usize>,
    filter: Option<&Filter>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let mut raw_header = read_raw_cram_header_text(&src)?;

    if deduplicate {
        raw_header = trim_header_line_trailing_spaces(&raw_header);
    }

    let index = read_associated_cram_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository)
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(raw_header.into_bytes());
    let mut seen_records = HashSet::new();
    let mut remaining = limit.unwrap_or(usize::MAX);

    for region in regions {
        if remaining == 0 {
            break;
        }

        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid reference sequence name: {}",
                        String::from_utf8_lossy(region.name())
                    ),
                )
            })?;
        let query = reader.query(&header, region)?;

        for result in query {
            if remaining == 0 {
                break;
            }

            let mut record = result?;

            if !record_intersects_region(&record, reference_sequence_id, region) {
                continue;
            }

            if let Some(filter) = filter
                && !record_matches_filter(&header, &record, filter)?
            {
                continue;
            }

            add_md_and_nm_to_record(reference_sequences, &header, &mut record)?;

            let mut line_writer = sam::io::Writer::new(Vec::new());
            line_writer.write_alignment_record(&header, &record)?;
            let line = line_writer.into_inner();

            if !deduplicate || seen_records.insert(line.clone()) {
                writer.get_mut().extend(line);
                remaining -= 1;
            }
        }
    }

    String::from_utf8(writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Reads the embedded raw SAM header text from a local CRAM file.
pub fn read_raw_cram_header_text<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    let mut reader = File::open(data_path).map(cram::io::Reader::new)?;
    let mut header_reader = reader.header_reader();

    header_reader.read_magic_number()?;
    header_reader.read_format_version()?;
    header_reader.read_file_id()?;

    let mut container_reader = header_reader.container_reader()?;
    let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    let mut raw_header = String::new();
    raw_sam_header_reader.read_to_string(&mut raw_header)?;
    raw_sam_header_reader.discard_to_end()?;

    Ok(raw_header)
}

fn trim_header_line_trailing_spaces(header: &str) -> String {
    header
        .lines()
        .map(|line| format!("{}\n", line.trim_end_matches(' ')))
        .collect()
}

fn read_fasta_sequences<P>(src: P) -> io::Result<HashMap<Vec<u8>, Vec<u8>>>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;
    let mut sequences = HashMap::new();

    for result in reader.records() {
        let record = result?;
        sequences.insert(record.name().to_vec(), record.sequence().as_ref().to_vec());
    }

    Ok(sequences)
}

fn record_intersects_region(
    record: &sam::alignment::RecordBuf,
    reference_sequence_id: usize,
    region: &Region,
) -> bool {
    if record.reference_sequence_id() != Some(reference_sequence_id) {
        return false;
    }

    match (record.alignment_start(), record.alignment_end()) {
        (Some(start), Some(end)) => region.interval().intersects((start..=end).into()),
        _ => false,
    }
}

fn add_md_and_nm_to_record(
    reference_sequences: &HashMap<Vec<u8>, Vec<u8>>,
    header: &Header,
    record: &mut sam::alignment::RecordBuf,
) -> io::Result<()> {
    use sam::alignment::record::data::field::Tag;
    use sam::alignment::record_buf::data::field::Value;

    let reference_sequence_id = record
        .reference_sequence_id()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing reference sequence"))?;
    let (reference_sequence_name, _) = header
        .reference_sequences()
        .get_index(reference_sequence_id)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
        })?;
    let reference_sequence = reference_sequences
        .get(reference_sequence_name.as_ref() as &[u8])
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "missing reference sequence: {}",
                    String::from_utf8_lossy(reference_sequence_name)
                ),
            )
        })?;
    let alignment_start = record
        .alignment_start()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing alignment start"))?;
    let (md, nm) = calculate_md_and_nm(
        reference_sequence,
        usize::from(alignment_start) - 1,
        record.cigar().as_ref(),
        record.sequence().as_ref(),
    )?;

    let data = record.data_mut();
    let read_group = data.remove(&Tag::READ_GROUP);
    data.insert(Tag::MISMATCHED_POSITIONS, Value::from(md));
    data.insert(Tag::EDIT_DISTANCE, Value::from(nm));

    if let Some((tag, value)) = read_group {
        data.insert(tag, value);
    }

    Ok(())
}

fn calculate_md_and_nm(
    reference_sequence: &[u8],
    mut reference_index: usize,
    cigar: &[sam::alignment::record::cigar::Op],
    read_sequence: &[u8],
) -> io::Result<(String, i32)> {
    use sam::alignment::record::cigar::op::Kind;

    let mut read_index = 0;
    let mut match_count = 0;
    let mut edit_distance = 0;
    let mut md = String::new();

    for op in cigar {
        match op.kind() {
            Kind::Match => {
                for _ in 0..op.len() {
                    let reference_base =
                        reference_sequence.get(reference_index).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "reference position out of range",
                            )
                        })?;
                    let read_base = read_sequence.get(read_index).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "read position out of range")
                    })?;

                    if reference_base.eq_ignore_ascii_case(read_base) {
                        match_count += 1;
                    } else {
                        md.push_str(&match_count.to_string());
                        md.push(reference_base.to_ascii_uppercase() as char);
                        match_count = 0;
                        edit_distance += 1;
                    }

                    reference_index += 1;
                    read_index += 1;
                }
            }
            Kind::SequenceMatch => {
                match_count += op.len();
                reference_index += op.len();
                read_index += op.len();
            }
            Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    let reference_base =
                        reference_sequence.get(reference_index).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "reference position out of range",
                            )
                        })?;

                    md.push_str(&match_count.to_string());
                    md.push(reference_base.to_ascii_uppercase() as char);
                    match_count = 0;
                    edit_distance += 1;
                    reference_index += 1;
                    read_index += 1;
                }
            }
            Kind::Insertion => {
                read_index += op.len();
                edit_distance += i32::try_from(op.len())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            }
            Kind::Deletion => {
                md.push_str(&match_count.to_string());
                md.push('^');

                for _ in 0..op.len() {
                    let reference_base =
                        reference_sequence.get(reference_index).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "reference position out of range",
                            )
                        })?;
                    md.push(reference_base.to_ascii_uppercase() as char);
                    reference_index += 1;
                }

                match_count = 0;
                edit_distance += i32::try_from(op.len())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            }
            Kind::Skip => {
                reference_index += op.len();
            }
            Kind::SoftClip => {
                read_index += op.len();
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    md.push_str(&match_count.to_string());

    Ok((md, edit_distance))
}

fn append_fai_extension(src: &Path) -> PathBuf {
    let mut path = src.as_os_str().to_os_string();
    path.push(".fai");
    PathBuf::from(path)
}

/// Returns the number of reference sequences declared in an alignment header.
pub fn reference_sequence_count(header: &Header) -> usize {
    header.reference_sequences().len()
}

fn summarize_alignment_record<R>(header: &Header, record: &R) -> io::Result<AlignmentRecordSummary>
where
    R: sam::alignment::Record + ?Sized,
{
    let name = record.name().map(|name| name.to_vec());
    let flags = record.flags()?;
    let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    let alignment_start = record.alignment_start().transpose()?;
    let mapping_quality = record
        .mapping_quality()
        .transpose()?
        .map(|mapping_quality| mapping_quality.get());
    let cigar = record.cigar().iter().collect::<io::Result<_>>()?;
    let mate_reference_sequence_id = record.mate_reference_sequence_id(header).transpose()?;
    let mate_alignment_start = record.mate_alignment_start().transpose()?;
    let template_length = record.template_length()?;
    let sequence = record.sequence().iter().collect();
    let quality_scores = record.quality_scores().iter().collect::<io::Result<_>>()?;

    Ok(AlignmentRecordSummary {
        name,
        flags,
        reference_sequence_id,
        alignment_start,
        mapping_quality,
        cigar,
        mate_reference_sequence_id,
        mate_alignment_start,
        template_length,
        sequence,
        quality_scores,
    })
}

fn push_base_modification_report<R>(
    dst: &mut String,
    _header: &Header,
    record: &R,
    report_unchecked: bool,
    extended: bool,
) -> io::Result<()>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::{
        alignment::record::Flags, record::data::field::value::base_modifications::group::Status,
    };

    let sequence = record.sequence().iter().collect::<Vec<_>>();
    let flags = record.flags()?;
    let is_reverse_complemented = flags.contains(Flags::REVERSE_COMPLEMENTED);
    let data = record.data();
    let base_modifications =
        parse_record_base_modifications(&data, is_reverse_complemented, &sequence)?;
    let calls = record_base_modification_calls(record, &base_modifications)?;

    for (i, base) in sequence.iter().enumerate() {
        dst.push_str(&format!("{i}\t{}", char::from(*base)));

        let mut sep = '\t';
        for call in &calls {
            if call.position == i {
                dst.push(sep);
                call.push_to(dst, extended);
                sep = ' ';
            }
        }

        if report_unchecked {
            push_unchecked_base_modification_calls(
                dst,
                &base_modifications,
                is_reverse_complemented,
                i,
                *base,
                sep,
            );
        }

        dst.push('\n');
    }

    dst.push_str("---\nPresent:");

    for group in base_modifications.as_ref() {
        let status = match group.status() {
            Some(Status::Explicit) => '?',
            Some(Status::Implicit) | None => '.',
        };

        for &modification in group.modifications() {
            dst.push(' ');
            dst.push_str(&format_present_modification(modification));
            dst.push(status);
        }
    }

    dst.push('\n');

    for (i, base) in sequence.iter().enumerate() {
        let mut line = String::new();
        let mut sep = '\t';

        for call in &calls {
            if call.position == i {
                line.push(sep);
                call.push_to(&mut line, false);
                sep = ' ';
            }
        }

        if report_unchecked {
            push_unchecked_base_modification_calls(
                &mut line,
                &base_modifications,
                is_reverse_complemented,
                i,
                *base,
                sep,
            );
        }

        if !line.is_empty() {
            dst.push_str(&format!("{i}\t{}", char::from(*base)));
            dst.push_str(&line);
            dst.push('\n');
        }
    }

    dst.push_str("\n===\n\n");

    Ok(())
}

fn parse_record_base_modifications(
    data: &dyn sam::alignment::record::data::Data,
    is_reverse_complemented: bool,
    sequence: &[u8],
) -> io::Result<sam::record::data::field::value::BaseModifications> {
    use sam::{
        alignment::{
            record::data::field::{Tag, Value},
            record_buf::Sequence,
        },
        record::data::field::value::BaseModifications,
    };

    validate_base_modification_sequence_length(data, sequence.len())?;

    let mm_value = match data.get(&Tag::BASE_MODIFICATIONS).transpose()? {
        Some(value) => Some(value),
        None => data.get(&Tag::new(b'M', b'm')).transpose()?,
    };

    match mm_value {
        Some(Value::String(s)) => {
            let sequence_buf = Sequence::from(sequence.to_vec());

            match BaseModifications::parse(s.as_ref(), is_reverse_complemented, &sequence_buf) {
                Ok(base_modifications) => Ok(base_modifications),
                Err(_) => {
                    parse_htslib_base_modifications(s.as_ref(), is_reverse_complemented, sequence)
                }
            }
        }
        Some(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "MM tag is not a string",
        )),
        None => Ok(BaseModifications::from(Vec::new())),
    }
}

fn validate_base_modification_sequence_length(
    data: &dyn sam::alignment::record::data::Data,
    sequence_len: usize,
) -> io::Result<()> {
    use sam::alignment::record::data::field::Tag;

    let Some(value) = data
        .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
        .transpose()?
    else {
        return Ok(());
    };

    let Some(n) = value.as_int() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "MN tag is not an integer",
        ));
    };
    let Ok(n) = usize::try_from(n) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "MN tag is negative",
        ));
    };

    if n == sequence_len {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "MN tag does not match sequence length",
        ))
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct PileupRecord {
    reference_name: String,
    start: usize,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    qpos_by_reference_position: Vec<Option<usize>>,
    modification_calls: Vec<BaseModificationCall>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct TestPileupRecord {
    name: Option<Vec<u8>>,
    reference_name: String,
    start: usize,
    mapping_quality: u8,
    is_reverse: bool,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    columns: Vec<TestPileupColumn>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct TestPileupColumn {
    reference_position: usize,
    qpos: usize,
    base: Option<u8>,
    is_deletion: bool,
    is_refskip: bool,
    is_head: bool,
    is_tail: bool,
    insertion_after: Vec<u8>,
    deletion_after: usize,
}

impl TestPileupRecord {
    fn try_from_record<R>(header: &Header, record: &R) -> io::Result<Option<Self>>
    where
        R: sam::alignment::Record + ?Sized,
    {
        use sam::alignment::record::Flags;

        let flags = record.flags()?;

        if flags.intersects(Flags::UNMAPPED | Flags::SECONDARY | Flags::QC_FAIL | Flags::DUPLICATE)
        {
            return Ok(None);
        }

        let Some(reference_sequence_id) = record.reference_sequence_id(header).transpose()? else {
            return Ok(None);
        };
        let Some(alignment_start) = record.alignment_start().transpose()? else {
            return Ok(None);
        };

        let name = record.name().map(|name| name.to_vec());
        let reference_name = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "missing reference sequence")
            })?;
        let sequence = record.sequence().iter().collect::<Vec<_>>();
        let mut quality_scores = record
            .quality_scores()
            .iter()
            .collect::<io::Result<Vec<_>>>()?;

        if quality_scores.is_empty() {
            quality_scores.resize(sequence.len(), u8::MAX);
        }

        let mapping_quality = record
            .mapping_quality()
            .transpose()?
            .map_or(255, |q| q.get());
        let start = usize::from(alignment_start) - 1;
        let mut columns = test_pileup_columns(record, start, &sequence)?;

        if let Some(first) = columns.first_mut() {
            first.is_head = true;
        }

        if let Some(last) = columns.last_mut() {
            last.is_tail = true;
        }

        Ok(Some(Self {
            name,
            reference_name,
            start,
            mapping_quality,
            is_reverse: flags.is_reverse_complemented(),
            sequence,
            quality_scores,
            columns,
        }))
    }
}

fn test_pileup_columns<R>(
    record: &R,
    start: usize,
    sequence: &[u8],
) -> io::Result<Vec<TestPileupColumn>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::cigar::op::Kind;

    let mut columns = Vec::new();
    let mut reference_position = start;
    let mut qpos = 0;
    let mut last_query_column = None;

    for result in record.cigar().iter() {
        let op = result?;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    let base = sequence.get(qpos).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR qpos is out of range")
                    })?;
                    columns.push(TestPileupColumn {
                        reference_position,
                        qpos,
                        base: Some(base),
                        is_deletion: false,
                        is_refskip: false,
                        is_head: false,
                        is_tail: false,
                        insertion_after: Vec::new(),
                        deletion_after: 0,
                    });
                    last_query_column = Some(columns.len() - 1);
                    reference_position += 1;
                    qpos += 1;
                }
            }
            Kind::Insertion => {
                if let Some(i) = last_query_column {
                    let end = qpos.checked_add(op.len()).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR qpos overflow")
                    })?;
                    let bases = sequence.get(qpos..end).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR qpos is out of range")
                    })?;
                    columns[i].insertion_after.extend_from_slice(bases);
                }

                qpos += op.len();
            }
            Kind::SoftClip => qpos += op.len(),
            Kind::Deletion | Kind::Skip => {
                if op.kind() == Kind::Deletion
                    && let Some(i) = last_query_column
                {
                    columns[i].deletion_after = op.len();
                }

                for _ in 0..op.len() {
                    columns.push(TestPileupColumn {
                        reference_position,
                        qpos,
                        base: None,
                        is_deletion: op.kind() == Kind::Deletion,
                        is_refskip: op.kind() == Kind::Skip,
                        is_head: false,
                        is_tail: false,
                        insertion_after: Vec::new(),
                        deletion_after: 0,
                    });
                    last_query_column = Some(columns.len() - 1);
                    reference_position += 1;
                }
            }
            Kind::Pad => {
                if let Some(i) = last_query_column {
                    columns[i]
                        .insertion_after
                        .extend(std::iter::repeat_n(b'*', op.len()));
                }
            }
            Kind::HardClip => {}
        }
    }

    Ok(columns)
}

impl PileupRecord {
    fn try_from_record<R>(header: &Header, record: &R) -> io::Result<Option<Self>>
    where
        R: sam::alignment::Record + ?Sized,
    {
        use sam::alignment::record::Flags;

        let flags = record.flags()?;

        if flags.intersects(Flags::UNMAPPED | Flags::SECONDARY | Flags::QC_FAIL | Flags::DUPLICATE)
        {
            return Ok(None);
        }

        let Some(reference_sequence_id) = record.reference_sequence_id(header).transpose()? else {
            return Ok(None);
        };
        let Some(alignment_start) = record.alignment_start().transpose()? else {
            return Ok(None);
        };

        let reference_name = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "missing reference sequence")
            })?;
        let sequence = record.sequence().iter().collect::<Vec<_>>();
        let mut quality_scores = record
            .quality_scores()
            .iter()
            .collect::<io::Result<Vec<_>>>()?;

        if quality_scores.is_empty() {
            quality_scores.resize(sequence.len(), u8::MAX);
        }

        let data = record.data();
        let base_modifications =
            parse_record_base_modifications(&data, flags.is_reverse_complemented(), &sequence)?;
        let modification_calls = record_base_modification_calls(record, &base_modifications)?;
        let qpos_by_reference_position = qpos_by_reference_position(record)?;

        Ok(Some(Self {
            reference_name,
            start: usize::from(alignment_start) - 1,
            sequence,
            quality_scores,
            qpos_by_reference_position,
            modification_calls,
        }))
    }
}

fn qpos_by_reference_position<R>(record: &R) -> io::Result<Vec<Option<usize>>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::cigar::op::Kind;

    let mut qpos_by_reference_position = Vec::new();
    let mut qpos = 0;

    for result in record.cigar().iter() {
        let op = result?;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    qpos_by_reference_position.push(Some(qpos));
                    qpos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => qpos += op.len(),
            Kind::Deletion | Kind::Skip => {
                qpos_by_reference_position.extend(std::iter::repeat_n(None, op.len()));
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(qpos_by_reference_position)
}

fn push_base_modification_pileup_report(records: &[PileupRecord]) -> io::Result<String> {
    let Some(first_record) = records.first() else {
        return Ok(String::new());
    };

    let reference_name = &first_record.reference_name;
    let start = records.iter().map(|record| record.start).min().unwrap_or(0);
    let end = records
        .iter()
        .map(|record| record.start + record.qpos_by_reference_position.len())
        .max()
        .unwrap_or(start);
    let mut report = String::new();

    for reference_position in start..end {
        let mut bases = String::new();
        let mut qualities = String::new();

        for record in records {
            if record.reference_name != *reference_name || reference_position < record.start {
                continue;
            }

            let offset = reference_position - record.start;
            let Some(qpos) = record
                .qpos_by_reference_position
                .get(offset)
                .copied()
                .flatten()
            else {
                continue;
            };

            let base = record.sequence.get(qpos).copied().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CIGAR qpos is out of range")
            })?;
            bases.push(char::from(base));
            push_pileup_modification_calls(&mut bases, &record.modification_calls, qpos);

            let quality_score = record.quality_scores.get(qpos).copied().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "quality score qpos is out of range",
                )
            })?;
            qualities.push(char::from(b'!'.saturating_add(quality_score).min(b'~')));
        }

        if !bases.is_empty() {
            report.push_str(&format!(
                "{reference_name}\t{reference_position}\t{bases}\t{qualities}\n"
            ));
        }
    }

    Ok(report)
}

fn push_test_pileup_report(records: &[TestPileupRecord]) -> io::Result<String> {
    let Some(first_record) = records.first() else {
        return Ok(String::new());
    };

    let reference_name = &first_record.reference_name;
    let start = records
        .iter()
        .flat_map(|record| {
            record
                .columns
                .iter()
                .map(|column| column.reference_position)
        })
        .min()
        .unwrap_or(0);
    let end = records
        .iter()
        .flat_map(|record| {
            record
                .columns
                .iter()
                .map(|column| column.reference_position + 1)
        })
        .max()
        .unwrap_or(start);
    let mut report = String::new();

    for reference_position in start..end {
        let mut bases = String::new();
        let mut depth = 0;
        let mut entries = Vec::new();

        for (record_index, record) in records.iter().enumerate() {
            if record.reference_name != *reference_name || reference_position < record.start {
                continue;
            }

            let Some(column) = record
                .columns
                .iter()
                .find(|column| column.reference_position == reference_position)
            else {
                continue;
            };

            depth += 1;
            entries.push((record_index, record, column));
        }

        let adjusted_qualities = adjusted_test_pileup_qualities(&entries)?;
        let mut qualities = String::with_capacity(entries.len());

        for ((_, record, column), quality) in entries.iter().zip(adjusted_qualities) {
            push_test_pileup_bases(&mut bases, record, column)?;
            push_test_pileup_quality_score(&mut qualities, quality);
        }

        if depth > 0 {
            report.push_str(&format!(
                "{reference_name}\t{}\t{depth}\t{bases}\t{qualities}\n",
                reference_position + 1
            ));
        }
    }

    Ok(report)
}

fn adjusted_test_pileup_qualities(
    entries: &[(usize, &TestPileupRecord, &TestPileupColumn)],
) -> io::Result<Vec<u8>> {
    let mut qualities = entries
        .iter()
        .map(|(_, record, column)| {
            Ok(record
                .quality_scores
                .get(column.qpos)
                .copied()
                .unwrap_or(u8::MAX))
        })
        .collect::<io::Result<Vec<_>>>()?;

    for i in 0..entries.len() {
        let (_, left_record, left_column) = entries[i];

        for j in (i + 1)..entries.len() {
            let (_, right_record, right_column) = entries[j];

            if !test_pileup_records_are_mates(left_record, right_record) {
                continue;
            }

            let same_base = test_pileup_columns_match(left_column, right_column);

            if left_column.is_deletion != right_column.is_deletion {
                let (left_multiplier, right_multiplier) =
                    overlap_quality_multipliers(left_record.name.as_deref().unwrap_or_default());
                let quality = qualities[i].saturating_add(qualities[j]).min(200);

                if left_column.is_deletion {
                    qualities[i] = quality * left_multiplier;
                    qualities[j] = ((f32::from(qualities[j]) * 0.8) as u8) * right_multiplier;
                } else {
                    qualities[j] = quality * right_multiplier;
                    qualities[i] = ((f32::from(qualities[i]) * 0.8) as u8) * left_multiplier;
                }
            } else if same_base {
                let (left_multiplier, right_multiplier) =
                    overlap_quality_multipliers(left_record.name.as_deref().unwrap_or_default());
                let quality = qualities[i].saturating_add(qualities[j]).min(200);
                qualities[i] = quality * left_multiplier;
                qualities[j] = quality * right_multiplier;
            } else if qualities[i] > qualities[j] {
                qualities[i] = (f32::from(qualities[i]) * 0.8) as u8;
                qualities[j] = 0;
            } else if qualities[i] < qualities[j] {
                qualities[j] = (f32::from(qualities[j]) * 0.8) as u8;
                qualities[i] = 0;
            } else {
                let (left_multiplier, right_multiplier) =
                    overlap_quality_multipliers(left_record.name.as_deref().unwrap_or_default());
                let quality = (f32::from(qualities[i]) * 0.8) as u8;
                qualities[i] = quality * left_multiplier;
                qualities[j] = quality * right_multiplier;
            };
        }
    }

    Ok(qualities)
}

fn overlap_quality_multipliers(name: &[u8]) -> (u8, u8) {
    if wang_hash(x31_hash_string(name)) & 1 == 1 {
        (1, 0)
    } else {
        (0, 1)
    }
}

fn x31_hash_string(name: &[u8]) -> u32 {
    let Some((&first, rest)) = name.split_first() else {
        return 0;
    };
    let mut hash = u32::from(first);

    for &byte in rest {
        hash = hash
            .wrapping_shl(5)
            .wrapping_sub(hash)
            .wrapping_add(u32::from(byte));
    }

    hash
}

fn wang_hash(mut key: u32) -> u32 {
    key = key.wrapping_add(!(key.wrapping_shl(15)));
    key ^= key >> 10;
    key = key.wrapping_add(key.wrapping_shl(3));
    key ^= key >> 6;
    key = key.wrapping_add(!(key.wrapping_shl(11)));
    key ^= key >> 16;
    key
}

fn test_pileup_records_are_mates(left: &TestPileupRecord, right: &TestPileupRecord) -> bool {
    left.name.is_some() && left.name == right.name && left.is_reverse != right.is_reverse
}

fn test_pileup_columns_match(
    left_column: &TestPileupColumn,
    right_column: &TestPileupColumn,
) -> bool {
    match (left_column.base, right_column.base) {
        (Some(left), Some(right)) => left.eq_ignore_ascii_case(&right),
        (None, None) => {
            left_column.is_deletion == right_column.is_deletion
                && left_column.is_refskip == right_column.is_refskip
        }
        _ => false,
    }
}

fn push_test_pileup_bases(
    dst: &mut String,
    record: &TestPileupRecord,
    column: &TestPileupColumn,
) -> io::Result<()> {
    if column.is_head {
        dst.push('^');
        dst.push(char::from(
            b'!'.saturating_add(record.mapping_quality.min(93)),
        ));
    }

    if column.is_deletion {
        dst.push('*');
    } else if column.is_refskip {
        dst.push(if record.is_reverse { '<' } else { '>' });
    } else {
        let base = column
            .base
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing pileup base"))?;
        let base = if record.is_reverse {
            base.to_ascii_lowercase()
        } else {
            base.to_ascii_uppercase()
        };
        dst.push(char::from(base));
    }

    if !column.insertion_after.is_empty() {
        dst.push_str(&format!("+{}(", column.insertion_after.len()));

        for base in &column.insertion_after {
            let base = if record.is_reverse {
                base.to_ascii_lowercase()
            } else {
                base.to_ascii_uppercase()
            };
            dst.push(char::from(base));
        }

        dst.push(')');
    }

    if column.deletion_after > 0 {
        dst.push_str(&format!("-{}()", column.deletion_after));
    }

    if column.is_tail {
        dst.push('$');
    }

    Ok(())
}

fn push_test_pileup_quality_score(dst: &mut String, quality_score: u8) {
    let q = quality_score
        .checked_add(33)
        .filter(|q| *q < b'~')
        .unwrap_or(b'~');

    dst.push(char::from(q));
}

fn push_pileup_modification_calls(dst: &mut String, calls: &[BaseModificationCall], qpos: usize) {
    let mut calls = calls.iter().filter(|call| call.position == qpos).peekable();

    if calls.peek().is_none() {
        return;
    }

    dst.push('[');

    for call in calls {
        dst.push(match call.strand {
            sam::record::data::field::value::base_modifications::group::Strand::Forward => '+',
            sam::record::data::field::value::base_modifications::group::Strand::Reverse => '-',
        });
        dst.push_str(&format_call_modification(call.modification));
        dst.push_str(&call.probability.unwrap_or(u8::MAX).to_string());
    }

    dst.push(']');
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct BaseModificationCall {
    position: usize,
    canonical_base: u8,
    strand: sam::record::data::field::value::base_modifications::group::Strand,
    modification: sam::record::data::field::value::base_modifications::group::Modification,
    status: Option<sam::record::data::field::value::base_modifications::group::Status>,
    probability: Option<u8>,
}

impl BaseModificationCall {
    fn push_to(&self, dst: &mut String, extended: bool) {
        dst.push(char::from(self.canonical_base));
        dst.push(match self.strand {
            sam::record::data::field::value::base_modifications::group::Strand::Forward => '+',
            sam::record::data::field::value::base_modifications::group::Strand::Reverse => '-',
        });
        dst.push_str(&format_call_modification(self.modification));

        if extended {
            dst.push(match self.status {
                Some(
                    sam::record::data::field::value::base_modifications::group::Status::Explicit,
                ) => '?',
                Some(
                    sam::record::data::field::value::base_modifications::group::Status::Implicit,
                )
                | None => '.',
            });
        }

        match self.probability {
            Some(probability) => dst.push_str(&probability.to_string()),
            None => dst.push('.'),
        }
    }
}

fn record_base_modification_probabilities<R>(record: &R) -> io::Result<Option<Vec<u8>>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::data::field::{Tag, Value, value::array::Array};

    let data = record.data();
    let ml_value = match data
        .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
        .transpose()?
    {
        Some(value) => Some(value),
        None => data.get(&Tag::new(b'M', b'l')).transpose()?,
    };

    match ml_value {
        Some(Value::Array(Array::UInt8(values))) => {
            values.iter().collect::<io::Result<_>>().map(Some)
        }
        Some(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "ML tag is not a B:C array",
        )),
        None => Ok(None),
    }
}

fn base_modification_calls(
    base_modifications: &sam::record::data::field::value::BaseModifications,
    probabilities: Option<&[u8]>,
) -> io::Result<Vec<BaseModificationCall>> {
    let mut calls = Vec::new();
    let mut probability_iter = probabilities.into_iter().flatten();

    for group in base_modifications.as_ref() {
        for &position in group.positions() {
            for &modification in group.modifications() {
                let probability = probability_iter.next().copied();
                if probabilities.is_some() && probability.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "ML tag has fewer probabilities than MM calls",
                    ));
                }

                calls.push(BaseModificationCall {
                    position,
                    canonical_base: u8::from(group.unmodified_base()),
                    strand: group.strand(),
                    modification,
                    status: group.status(),
                    probability,
                });
            }
        }
    }

    if probability_iter.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "ML tag has more probabilities than MM calls",
        ));
    }

    Ok(calls)
}

fn record_base_modification_calls<R>(
    record: &R,
    base_modifications: &sam::record::data::field::value::BaseModifications,
) -> io::Result<Vec<BaseModificationCall>>
where
    R: sam::alignment::Record + ?Sized,
{
    let probabilities = record_base_modification_probabilities(record)?;

    base_modification_calls(base_modifications, probabilities.as_deref())
}

fn parse_htslib_base_modifications(
    src: &[u8],
    is_reverse_complemented: bool,
    sequence: &[u8],
) -> io::Result<sam::record::data::field::value::BaseModifications> {
    use sam::record::data::field::value::{
        BaseModifications,
        base_modifications::{
            Group,
            group::{Status, Strand, UnmodifiedBase},
        },
    };

    let mut groups = Vec::new();

    for raw_group in src.split(|b| *b == b';').filter(|group| !group.is_empty()) {
        let (&raw_base, rest) = raw_group
            .split_first()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing MM base"))?;
        let unmodified_base = UnmodifiedBase::try_from(raw_base)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let (&raw_strand, rest) = rest
            .split_first()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing MM strand"))?;
        let strand = Strand::try_from(raw_strand)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let (raw_modifications, rest) = split_modifications(rest);
        let modifications = parse_htslib_modification_codes(raw_modifications)?;
        let (status, rest) = match rest.split_first() {
            Some((b'.', rest)) => (Some(Status::Implicit), rest),
            Some((b'?', rest)) => (Some(Status::Explicit), rest),
            _ => (None, rest),
        };
        let skip_counts = parse_htslib_skip_counts(rest)?;
        let positions = decode_htslib_modification_positions(
            &skip_counts,
            is_reverse_complemented,
            sequence,
            unmodified_base,
        )?;

        groups.push(Group::new(
            unmodified_base,
            strand,
            modifications,
            status,
            positions,
        ));
    }

    Ok(BaseModifications::from(groups))
}

fn split_modifications(src: &[u8]) -> (&[u8], &[u8]) {
    let end = src
        .iter()
        .position(|b| matches!(*b, b'.' | b'?' | b','))
        .unwrap_or(src.len());

    src.split_at(end)
}

fn parse_htslib_modification_codes(
    src: &[u8],
) -> io::Result<Vec<sam::record::data::field::value::base_modifications::group::Modification>> {
    use sam::record::data::field::value::base_modifications::group::Modification;

    if src.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing MM modification",
        ));
    }

    if src.iter().all(u8::is_ascii_digit) {
        let id = std::str::from_utf8(src)
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid ChEBI ID"))?;
        Ok(vec![Modification::ChebiId(id)])
    } else {
        src.iter()
            .copied()
            .map(|code| {
                Modification::try_from(code)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .collect()
    }
}

fn parse_htslib_skip_counts(src: &[u8]) -> io::Result<Vec<usize>> {
    if src.is_empty() {
        return Ok(Vec::new());
    }

    if !src.starts_with(b",") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid MM skip count list",
        ));
    }

    src[1..]
        .split(|b| *b == b',')
        .map(|raw_count| {
            std::str::from_utf8(raw_count)
                .ok()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid MM skip count"))
        })
        .collect()
}

fn decode_htslib_modification_positions(
    skip_counts: &[usize],
    is_reverse_complemented: bool,
    sequence: &[u8],
    unmodified_base: sam::record::data::field::value::base_modifications::group::UnmodifiedBase,
) -> io::Result<Vec<usize>> {
    let mut positions = Vec::with_capacity(skip_counts.len());
    let canonical_base = if is_reverse_complemented {
        unmodified_base.complement()
    } else {
        unmodified_base
    };
    let canonical_base = u8::from(canonical_base);
    let candidate_positions = || {
        let iter = sequence.iter().enumerate().filter_map(move |(i, base)| {
            (canonical_base == b'N' || *base == canonical_base).then_some(i)
        });

        if is_reverse_complemented {
            Box::new(iter.rev()) as Box<dyn Iterator<Item = usize>>
        } else {
            Box::new(iter) as Box<dyn Iterator<Item = usize>>
        }
    };
    let mut candidates = candidate_positions();

    for &count in skip_counts {
        let position = candidates.nth(count).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "MM skip count is out of range")
        })?;
        positions.push(position);
    }

    Ok(positions)
}

fn push_unchecked_base_modification_calls(
    dst: &mut String,
    base_modifications: &sam::record::data::field::value::BaseModifications,
    is_reverse_complemented: bool,
    position: usize,
    base: u8,
    mut sep: char,
) {
    for group in base_modifications.as_ref() {
        if group.status()
            != Some(sam::record::data::field::value::base_modifications::group::Status::Explicit)
            || group.positions().contains(&position)
            || !is_group_base_at(group.unmodified_base(), is_reverse_complemented, base)
        {
            continue;
        }

        for &modification in group.modifications() {
            dst.push(sep);
            dst.push(char::from(u8::from(group.unmodified_base())));
            dst.push(match group.strand() {
                sam::record::data::field::value::base_modifications::group::Strand::Forward => '+',
                sam::record::data::field::value::base_modifications::group::Strand::Reverse => '-',
            });
            dst.push_str(&format_call_modification(modification));
            dst.push('#');
            sep = ' ';
        }
    }
}

fn is_group_base_at(
    base: sam::record::data::field::value::base_modifications::group::UnmodifiedBase,
    is_reverse_complemented: bool,
    actual: u8,
) -> bool {
    let expected = if is_reverse_complemented {
        base.complement()
    } else {
        base
    };

    u8::from(expected) == actual
}

fn format_call_modification(
    modification: sam::record::data::field::value::base_modifications::group::Modification,
) -> String {
    use sam::record::data::field::value::base_modifications::group::Modification;

    match modification {
        Modification::Code(code) => char::from(code).to_string(),
        Modification::ChebiId(id) => format!("({id})"),
    }
}

fn format_present_modification(
    modification: sam::record::data::field::value::base_modifications::group::Modification,
) -> String {
    use sam::record::data::field::value::base_modifications::group::Modification;

    match modification {
        Modification::Code(code) => char::from(code).to_string(),
        Modification::ChebiId(id) => format!("#-{id}"),
    }
}

fn apply_htslib_cram_template_lengths(records: &mut [AlignmentRecordSummary]) {
    let mut start = 0;

    while start < records.len() {
        let name = records[start].name.clone();
        let mut end = start + 1;

        while end < records.len() && records[end].name == name {
            end += 1;
        }

        apply_htslib_cram_template_lengths_to_group(&mut records[start..end]);
        start = end;
    }
}

fn apply_htslib_cram_template_lengths_to_group(records: &mut [AlignmentRecordSummary]) {
    if records.len() < 2 {
        return;
    }

    let Some(reference_sequence_id) = records[0].reference_sequence_id else {
        return;
    };

    if !records
        .iter()
        .all(|record| record.reference_sequence_id == Some(reference_sequence_id))
    {
        for record in records {
            record.template_length = 0;
        }

        return;
    }

    let Some(left) = records
        .iter()
        .filter_map(|record| record.alignment_start.map(usize::from))
        .min()
    else {
        return;
    };
    let Some(right) = records.iter().filter_map(alignment_summary_end).max() else {
        return;
    };

    let left_count = records
        .iter()
        .filter(|record| record.alignment_start.map(usize::from) == Some(left))
        .count();
    let right_count = records
        .iter()
        .filter(|record| alignment_summary_end(record) == Some(right))
        .count();
    let Ok(len) = i32::try_from(right - left + 1) else {
        return;
    };

    let mut next_len = len;

    for (i, record) in records.iter_mut().enumerate() {
        let start = record.alignment_start.map(usize::from);
        let end = alignment_summary_end(record);

        record.template_length = if i == 0
            && start == Some(left)
            && (end.is_some_and(|end| end < right) || left_count <= 1)
        {
            next_len = -len;
            len
        } else if i == 0
            && start == Some(left)
            && end == Some(right)
            && left_count > 1
            && right_count > 1
        {
            if record.flags.is_first_segment() {
                next_len = -len;
                len
            } else {
                -len
            }
        } else if i == 0 {
            -len
        } else {
            next_len
        };
    }
}

fn alignment_summary_end(record: &AlignmentRecordSummary) -> Option<usize> {
    let start = record.alignment_start.map(usize::from)?;
    let span = record
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum::<usize>();

    Some(start + span.saturating_sub(1))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::PathBuf;

    use super::{
        count_bam_records_from_path, count_bam_records_in_region_from_path,
        count_sam_records_from_path, mpileup_baq_from_alignment, mpileup_indel_alignment_score,
        query_bam_regions_from_path, read_bam_header_from_path, read_cram_header_from_path,
        read_sam_header_from_path, reference_sequence_count,
        synchronized_pileup_from_alignment_paths,
        view_sam_as_fastq_split_text_from_reader_with_flag_filter_and_suffix,
        write_bam_from_sam_reader, write_bam_regions_from_path,
    };

    fn fixture(path: &str) -> PathBuf {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        if path.starts_with("bcftools/") {
            base.join("..").join(path)
        } else {
            base.join(path)
        }
    }

    #[test]
    fn test_read_sam_header_and_records() {
        let path = fixture("htslib/test/xx#minimal.sam");
        let header = read_sam_header_from_path(&path).unwrap();

        assert_eq!(reference_sequence_count(&header), 2);
        assert_eq!(count_sam_records_from_path(path).unwrap(), 8);
    }

    #[test]
    fn test_write_bam_from_sam_reader() {
        let sam = b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\n";
        let bam_data = write_bam_from_sam_reader(Cursor::new(sam), Vec::new()).unwrap();
        let mut reader = crate::bam::io::Reader::new(Cursor::new(bam_data));
        let header = reader.read_header().unwrap();
        assert_eq!(reference_sequence_count(&header), 0);
        let mut count = 0;
        for result in reader.records() {
            result.unwrap();
            count += 1;
        }
        assert_eq!(count, 1);
    }

    #[test]
    fn test_view_sam_as_fastq_split_text_from_reader() {
        let sam = b"@HD\tVN:1.6\nr1\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\nr1\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\n";
        let mut reader = crate::sam::io::Reader::new(Cursor::new(sam));
        let split = view_sam_as_fastq_split_text_from_reader_with_flag_filter_and_suffix(
            &mut reader,
            0,
            0,
            0,
            false,
        )
        .unwrap();

        assert_eq!(split.read1, "@r1\nAC\n+\n!!\n");
        assert_eq!(split.read2, "@r1\nTG\n+\n##\n");
        assert_eq!(split.singleton, "");
    }

    #[test]
    fn test_read_bam_header_and_records() {
        let path = fixture("htslib/test/range.bam");
        let header = read_bam_header_from_path(&path).unwrap();

        assert!(reference_sequence_count(&header) > 0);
        assert!(count_bam_records_from_path(path).unwrap() > 0);
    }

    #[test]
    fn test_query_bam_records() {
        let path = fixture("htslib/test/range.bam");
        let region = "CHROMOSOME_II:2980-2980".parse().unwrap();

        assert_eq!(
            count_bam_records_in_region_from_path(&path, &region).unwrap(),
            1
        );

        let regions = [
            "CHROMOSOME_II:2980-2980".parse().unwrap(),
            "CHROMOSOME_IV:1500-1500".parse().unwrap(),
            "CHROMOSOME_II:2980-2980".parse().unwrap(),
            "CHROMOSOME_I:1000-1100".parse().unwrap(),
        ];

        assert_eq!(
            query_bam_regions_from_path(path, &regions).unwrap().len(),
            7
        );
    }

    #[test]
    fn test_write_bam_regions_from_path() {
        let path = fixture("htslib/test/range.bam");
        let regions = [
            "CHROMOSOME_II:2980-2980".parse().unwrap(),
            "CHROMOSOME_IV:1500-1500".parse().unwrap(),
            "CHROMOSOME_II:2980-2980".parse().unwrap(),
        ];

        let bam_data = write_bam_regions_from_path(&path, &regions, Vec::new()).unwrap();
        let mut reader = crate::bam::io::Reader::new(Cursor::new(bam_data));
        let header = reader.read_header().unwrap();
        assert!(reference_sequence_count(&header) > 0);

        let mut count = 0;
        for result in reader.records() {
            result.unwrap();
            count += 1;
        }

        assert_eq!(count, 4);
    }

    #[test]
    fn test_read_cram_header_and_records() {
        let path = fixture("htslib/test/range.cram");
        let header = read_cram_header_from_path(&path).unwrap();

        assert!(reference_sequence_count(&header) > 0);
    }

    #[test]
    fn test_mpileup_baq_helper_uses_probaln_path() {
        let reference = b"ACGTACGTACGT";
        let sequence = b"ACGTTCGT";
        let qualities = vec![30; sequence.len()];

        let baq = mpileup_baq_from_alignment(sequence, &qualities, reference, 2, "8M", false)
            .unwrap()
            .unwrap();
        let extended = mpileup_baq_from_alignment(sequence, &qualities, reference, 2, "8M", true)
            .unwrap()
            .unwrap();

        assert_eq!(baq.len(), sequence.len());
        assert_eq!(extended.len(), sequence.len());
        assert!(baq.bytes().all(|b| b >= 64));
        assert!(extended.bytes().all(|b| b >= 64));
    }

    #[test]
    fn test_mpileup_baq_helper_skips_reference_skips() {
        let reference = b"ACGTACGTACGT";
        let sequence = b"ACGT";
        let qualities = vec![30; sequence.len()];

        let baq = mpileup_baq_from_alignment(sequence, &qualities, reference, 2, "2M2N2M", false)
            .unwrap();

        assert_eq!(baq, None);
    }

    #[test]
    fn test_mpileup_indel_alignment_score_uses_bam2bcf_parameters() {
        let exact =
            mpileup_indel_alignment_score(b"ACGT", b"ACGT", Some(&[30, 30, 30, 30])).unwrap();
        let inserted =
            mpileup_indel_alignment_score(b"ACGT", b"ACGTT", Some(&[30, 30, 30, 30, 30])).unwrap();

        assert_eq!(exact, 5);
        assert!(inserted > exact);
    }

    #[test]
    fn test_synchronized_pileup_reports_per_input_columns() {
        let left = std::env::temp_dir().join(format!(
            "htslib-rs-sync-pileup-{}-left.sam",
            std::process::id()
        ));
        let right = std::env::temp_dir().join(format!(
            "htslib-rs-sync-pileup-{}-right.sam",
            std::process::id()
        ));

        std::fs::write(
            &left,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:20\n\
             left\t0\tsq0\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();
        std::fs::write(
            &right,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:20\n\
             right\t0\tsq0\t2\t60\t4M\t*\t0\t0\tCGTA\tJJJJ\n",
        )
        .unwrap();

        let columns = synchronized_pileup_from_alignment_paths(&[left.clone(), right.clone()])
            .inspect(|_| {
                let _ = std::fs::remove_file(&left);
                let _ = std::fs::remove_file(&right);
            })
            .unwrap();

        assert!(!columns.is_empty());
        assert!(
            columns
                .iter()
                .all(|column| column.depths_by_input.len() == 2)
        );
        assert!(
            columns
                .iter()
                .all(|column| column.bases_by_input.len() == 2)
        );
        assert!(
            columns
                .iter()
                .all(|column| column.qualities_by_input.len() == 2)
        );
        assert!(columns.iter().any(|column| {
            column.depths_by_input[0] > 0
                && column.depths_by_input[1] > 0
                && column.total_depth == column.depths_by_input.iter().sum::<usize>()
        }));
    }
}
