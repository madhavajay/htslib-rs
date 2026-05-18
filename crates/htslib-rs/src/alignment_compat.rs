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
    aux_values: Vec<([u8; 2], Vec<u8>)>,
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

    /// Returns the 1-based alignment start position, if any.
    pub fn alignment_start(&self) -> Option<usize> {
        self.alignment_start.map(usize::from)
    }

    /// Returns the 1-based mate alignment start position, if any.
    pub fn mate_alignment_start(&self) -> Option<usize> {
        self.mate_alignment_start.map(usize::from)
    }

    /// Returns the read name bytes (without the trailing NUL), if any.
    pub fn name_bytes(&self) -> Option<&[u8]> {
        self.name.as_deref()
    }

    /// Returns the ASCII sequence bytes.
    pub fn sequence_bytes(&self) -> &[u8] {
        &self.sequence
    }

    /// Query length derived from the CIGAR: the sum of query-consuming
    /// operations (`M`/`I`/`S`/`=`/`X`), matching `samtools view -m` and the
    /// upstream test harness's `querylen`. An empty CIGAR (e.g. an unmapped
    /// record with `*`) yields 0, so such records are excluded by `-m INT`.
    pub fn cigar_query_len(&self) -> usize {
        use sam::alignment::record::cigar::op::Kind;

        self.cigar
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::Insertion
                        | Kind::SoftClip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                )
            })
            .map(|op| op.len())
            .sum()
    }

    /// Returns the raw phred quality-score bytes.
    pub fn quality_score_bytes(&self) -> &[u8] {
        &self.quality_scores
    }

    /// Returns the SAM-text payload for an auxiliary tag, if present.
    ///
    /// The payload is the value after the `TAG:T:` prefix, matching
    /// `samtools view -d TAG[:VAL]` value comparison.
    pub fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
        self.aux_values
            .iter()
            .find_map(|(stored, value)| (*stored == tag).then_some(value.as_slice()))
    }

    /// Returns the `RG:Z:` read-group id, if present.
    pub fn read_group_id(&self) -> Option<&[u8]> {
        self.aux_value(*b"RG")
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

/// Writes CRAM records matching an HTSlib-style filter expression as SAM text
/// using a synthetic all-`N` reference built from the CRAM header.
///
/// This is intended for reference-independent filters on CRAMs that can be
/// decoded without reconstructing reference bases exactly. The decoded
/// sequence bytes are not meaningful for reference-compressed CRAM records.
pub fn view_cram_as_sam_text_matching_filter_from_path_synthesizing_reference<P>(
    src: P,
    filter: &str,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(&src);
    let header = read_cram_header_from_path(&data_path)?;
    let reference_sequence_repository = synthetic_cram_reference_repository_from_header(&header);
    let reader = File::open(data_path)?;
    view_cram_as_sam_text_matching_filter_with_reference_repository(
        reader,
        reference_sequence_repository,
        filter,
    )
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
    refid: Option<usize>,
    pos: Option<usize>,
    flag: u16,
    mapq: Option<u8>,
    mrname: Option<String>,
    mrefid: Option<usize>,
    mpos: Option<usize>,
    tlen: i32,
    qlen: usize,
    cigar: Option<String>,
    ncigar: usize,
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
        let refid = record.reference_sequence_id(header).transpose()?;
        let pos = record.alignment_start().transpose()?.map(usize::from);
        let flag = record.flags()?.bits();
        let mapq = record.mapping_quality().transpose()?.map(|mapq| mapq.get());
        let mrname = record
            .mate_reference_sequence(header)
            .transpose()?
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned());
        let mrefid = record.mate_reference_sequence_id(header).transpose()?;
        let mpos = record.mate_alignment_start().transpose()?.map(usize::from);
        let tlen = record.template_length()?;
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
        let ncigar = cigar_ops.len();
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
            refid,
            pos,
            flag,
            mapq,
            mrname,
            mrefid,
            mpos,
            tlen,
            qlen,
            cigar,
            ncigar,
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
        } else if symbol.starts_with("rnext") {
            (ExprValue::string(self.mrname.as_deref().unwrap_or("*")), 5)
        } else if symbol.starts_with("mrname") {
            (ExprValue::string(self.mrname.as_deref().unwrap_or("*")), 6)
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
        } else if symbol.starts_with("endpos") {
            (number_or_undefined(self.endpos().map(|n| n as f64)), 6)
        } else if symbol.starts_with("pnext") {
            (ExprValue::number(self.mpos.unwrap_or(0) as f64), 5)
        } else if symbol.starts_with("mpos") {
            (ExprValue::number(self.mpos.unwrap_or(0) as f64), 4)
        } else if let Some((value, len)) = flag_expr_value(symbol, self.flag) {
            (value, len)
        } else if symbol.starts_with("flag") {
            (ExprValue::number(f64::from(self.flag)), 4)
        } else if symbol.starts_with("mapq") {
            (number_or_undefined(self.mapq.map(f64::from)), 4)
        } else if symbol.starts_with("mrefid") {
            (ExprValue::number(ref_id_value(self.mrefid)), 6)
        } else if symbol.starts_with("refid") {
            (ExprValue::number(ref_id_value(self.refid)), 5)
        } else if symbol.starts_with("qlen") {
            (ExprValue::number(self.qlen as f64), 4)
        } else if symbol.starts_with("rlen") {
            (ExprValue::number(self.rlen as f64), 4)
        } else if symbol.starts_with("sclen") {
            (ExprValue::number(self.sclen as f64), 5)
        } else if symbol.starts_with("hclen") {
            (ExprValue::number(self.hclen as f64), 5)
        } else if symbol.starts_with("ncigar") {
            (ExprValue::number(self.ncigar as f64), 6)
        } else if symbol.starts_with("tlen") {
            (ExprValue::number(f64::from(self.tlen)), 4)
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

    fn endpos(&self) -> Option<usize> {
        self.pos.map(|pos| pos + self.rlen.saturating_sub(1))
    }
}

fn ref_id_value(id: Option<usize>) -> f64 {
    id.map_or(-1.0, |id| id as f64)
}

fn flag_expr_value(symbol: &str, flag: u16) -> Option<(ExprValue, usize)> {
    let suffix = symbol.strip_prefix("flag.")?;
    let (mask, len): (u16, usize) = if suffix.starts_with("paired") {
        (0x1, "paired".len())
    } else if suffix.starts_with("proper_pair") {
        (0x2, "proper_pair".len())
    } else if suffix.starts_with("unmap") {
        (0x4, "unmap".len())
    } else if suffix.starts_with("munmap") {
        (0x8, "munmap".len())
    } else if suffix.starts_with("reverse") {
        (0x10, "reverse".len())
    } else if suffix.starts_with("mreverse") {
        (0x20, "mreverse".len())
    } else if suffix.starts_with("read1") {
        (0x40, "read1".len())
    } else if suffix.starts_with("read2") {
        (0x80, "read2".len())
    } else if suffix.starts_with("secondary") {
        (0x100, "secondary".len())
    } else if suffix.starts_with("qcfail") {
        (0x200, "qcfail".len())
    } else if suffix.starts_with("dup") {
        (0x400, "dup".len())
    } else if suffix.starts_with("supplementary") {
        (0x800, "supplementary".len())
    } else {
        return None;
    };

    Some((
        ExprValue::number(f64::from(flag & mask)),
        "flag.".len() + len,
    ))
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
    let options = PileupOptions::default();
    let inputs = paths
        .iter()
        .map(|path| read_test_pileup_records_from_alignment_path(path, &options))
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

/// Record-selection options for the pileup iterators.
///
/// `Default` mirrors HTSlib's pileup default: exclude unmapped (`0x4`),
/// secondary (`0x100`), QC-fail (`0x200`) and duplicate (`0x400`) records.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PileupOptions {
    /// Skip a record if it has any of these flag bits set.
    pub exclude_flags: u16,
    /// Skip a record unless it has all of these flag bits set.
    pub require_flags: u16,
    /// Skip a record whose mapping quality is below this value.
    pub min_mapping_quality: u8,
    /// Apply HTSlib's smart overlap removal (zero one mate's overlapping
    /// base qualities) — HTSlib's `MPLP_SMART_OVERLAPS` default.
    pub detect_overlaps: bool,
    /// Discard "orphan"/anomalous reads: paired but not in a proper pair —
    /// HTSlib mpileup's `MPLP_NO_ORPHAN` default (cleared by `-A`).
    pub discard_orphans: bool,
    /// Apply default BAQ realignment to base qualities when a reference is
    /// provided.
    pub apply_baq: bool,
}

impl Default for PileupOptions {
    fn default() -> Self {
        Self {
            exclude_flags: 0x4 | 0x100 | 0x200 | 0x400,
            require_flags: 0,
            min_mapping_quality: 0,
            detect_overlaps: true,
            discard_orphans: false,
            apply_baq: false,
        }
    }
}

/// Ports HTSlib's `tweak_overlap_quality` / `overlap_push`: for proper-pair
/// mates that overlap on the reference, one mate's overlapping base qualities
/// are zeroed (and the surviving mate's boosted) so the duplicate coverage is
/// not double-counted. Operates per input, mutating `quality_scores` in place.
fn apply_overlap_correction(records: &mut [TestPileupRecord]) {
    use std::collections::HashMap;

    let mut groups: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        let Some(name) = rec.name.as_ref() else {
            continue;
        };
        // mate mapped (0x8 clear) and a proper pair (0x2 set)
        if rec.flags & 0x8 != 0 || rec.flags & 0x2 == 0 {
            continue;
        }
        if let (Some(mtid), Some(_)) = (rec.mate_reference_sequence_id, rec.mate_start) {
            if mtid != rec.reference_sequence_id {
                continue;
            }
        } else {
            continue;
        }
        groups.entry(name.clone()).or_default().push(i);
    }

    for indices in groups.values() {
        if indices.len() != 2 {
            continue;
        }
        let (i, j) = (indices[0], indices[1]);

        // `a` is the mate stored first (its mate lies to the right:
        // mpos >= pos); fall back to the leftmost start.
        let a_first = |r: &TestPileupRecord| r.mate_start.is_some_and(|m| m >= r.start);
        let (ai, bi) = if a_first(&records[i]) && !a_first(&records[j]) {
            (i, j)
        } else if a_first(&records[j]) && !a_first(&records[i]) {
            (j, i)
        } else if records[i].start <= records[j].start {
            (i, j)
        } else {
            (j, i)
        };

        // Wild-cigar guard mirrors overlap_push.
        {
            let b = &records[bi];
            let b_end = b.start + b.columns.len();
            if (b.template_length.unsigned_abs() as usize) >= 2 * b.sequence.len().max(1)
                && b.mate_start.is_some_and(|m| m >= b_end)
            {
                continue;
            }
        }

        let name = records[ai].name.clone().unwrap_or_default();
        let modify_b = wang_hash(x31_hash_string(&name)) & 1 == 1;
        let (amul, bmul): (u32, u32) = if modify_b { (1, 0) } else { (0, 1) };

        // Shared reference positions with real bases in both mates.
        let map_of = |r: &TestPileupRecord| {
            r.columns
                .iter()
                .filter_map(|c| c.base.map(|b| (c.reference_position, (c.qpos, b))))
                .collect::<HashMap<usize, (usize, u8)>>()
        };
        let a_map = map_of(&records[ai]);
        let b_map = map_of(&records[bi]);

        let mut updates: Vec<(bool, usize, u8)> = Vec::new();
        for (refpos, &(a_qpos, a_base)) in &a_map {
            let Some(&(b_qpos, b_base)) = b_map.get(refpos) else {
                continue;
            };
            let qa = u32::from(records[ai].quality_scores.get(a_qpos).copied().unwrap_or(0));
            let qb = u32::from(records[bi].quality_scores.get(b_qpos).copied().unwrap_or(0));
            let (new_a, new_b): (u8, u8) = if a_base.eq_ignore_ascii_case(&b_base) {
                let q = (qa + qb).min(200);
                ((amul * q) as u8, (bmul * q) as u8)
            } else if qa > qb {
                ((qa as f64 * 0.8) as u8, 0)
            } else if qa < qb {
                (0, (qb as f64 * 0.8) as u8)
            } else {
                (
                    (amul as f64 * 0.8 * qa as f64) as u8,
                    (bmul as f64 * 0.8 * qb as f64) as u8,
                )
            };
            updates.push((true, a_qpos, new_a));
            updates.push((false, b_qpos, new_b));
        }

        for (is_a, qpos, value) in updates {
            let rec = if is_a {
                &mut records[ai]
            } else {
                &mut records[bi]
            };
            if let Some(slot) = rec.quality_scores.get_mut(qpos) {
                *slot = value;
            }
        }
    }
}

fn apply_baq_to_pileup_records(
    records: &mut [TestPileupRecord],
    reference_sequences: &HashMap<Vec<u8>, Vec<u8>>,
) -> io::Result<()> {
    for record in records {
        let Some(reference) = reference_sequences.get(record.reference_name.as_bytes()) else {
            continue;
        };
        let cigar = pileup_cigar_string(&record.cigar);
        let Some(baq) = mpileup_baq_from_alignment(
            &record.sequence,
            &record.quality_scores,
            reference,
            record.start,
            &cigar,
            false,
        )?
        else {
            continue;
        };

        for (quality, baq) in record.quality_scores.iter_mut().zip(baq.bytes()) {
            *quality = quality.saturating_sub(baq.saturating_sub(64));
        }
    }

    Ok(())
}

fn pileup_cigar_string(cigar: &[(u8, usize)]) -> String {
    let mut out = String::new();
    for &(op, len) in cigar {
        out.push_str(&len.to_string());
        out.push(match op {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => 'M',
        });
    }
    out
}

/// A single read's contribution to a pileup column (HTSlib `bam_pileup1_t`-shaped).
#[derive(Clone, Debug, PartialEq)]
pub struct PileupRead {
    /// Read name, when present.
    pub name: Option<Vec<u8>>,
    /// Mapping quality (255 when unavailable).
    pub mapping_quality: u8,
    /// Whether the read is reverse-complemented.
    pub is_reverse: bool,
    /// Query base at this column; `None` for a deletion or reference skip.
    pub base: Option<u8>,
    /// Base quality (Phred) at this column; `None` for a deletion or reference skip.
    pub quality: Option<u8>,
    /// Raw read quality at `qpos` regardless of deletion/refskip (`0` when
    /// `qpos` is past the sequence). Mirrors HTSlib's
    /// `qpos < l_qseq ? qual[qpos] : 0` used by mpileup's base-quality gate.
    pub qpos_quality: u8,
    /// 0-based offset into the read sequence aligned at this column.
    pub qpos: usize,
    /// This column is a deletion in the read (CIGAR `D`).
    pub is_deletion: bool,
    /// This column is a reference skip in the read (CIGAR `N`).
    pub is_refskip: bool,
    /// This is the first aligned column of the read.
    pub is_head: bool,
    /// This is the last aligned column of the read.
    pub is_tail: bool,
    /// Indel immediately following this column: `>0` insertion length, `<0`
    /// deletion length, `0` none (HTSlib `bam_pileup1_t::indel`).
    pub indel: i32,
    /// Inserted bases immediately following this column (when `indel > 0`).
    pub insertion: Vec<u8>,
    /// Consensus Bayesian `poly_len` at this column (homopolymer run
    /// length), from the per-read `nm_init` precompute with default
    /// `nm_halo`/`sc_cost`.
    pub bayes_poly: i32,
    /// Consensus Bayesian `nm_local` at this column (local-NM score),
    /// same precompute.
    pub bayes_nm_local: f64,
}

/// A pileup column synchronized across one or more inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct PileupColumn {
    /// Reference sequence name.
    pub reference_name: String,
    /// 1-based reference position.
    pub position: usize,
    /// Reads overlapping this column, grouped by input in argument order.
    pub reads_by_input: Vec<Vec<PileupRead>>,
}

impl PileupColumn {
    /// Total depth across all inputs at this column.
    pub fn total_depth(&self) -> usize {
        self.reads_by_input.iter().map(Vec::len).sum()
    }

    /// Depth contributed by a single input index.
    pub fn depth_of_input(&self, input_index: usize) -> usize {
        self.reads_by_input.get(input_index).map_or(0, Vec::len)
    }
}

/// Faithful port of `samtools/bam_consensus.c` `nm_init`'s per-read
/// precompute (default options: `adj_qual` on, `homopoly_fix` off,
/// non-`BAYES_116`). Produces the packed `local_nm` array: high 8 bits
/// = homopolymer run length, low 24 bits = local-NM score. `seq` is
/// ASCII bases, `qual` Phred, `cigar` as `(bam_op,len)`, `md` the MD
/// aux bytes. `poly_len`/`nm_local` index this array.
///
/// `homopoly_fix` (opt-in, no test fixtures) is not modelled.
#[allow(clippy::needless_range_loop)]
pub fn compute_local_nm(
    seq: &[u8],
    qual: &[u8],
    cigar: &[(u8, usize)],
    md: Option<&[u8]>,
    nm_halo: i64,
    sc_cost: i32,
    mode_bayes116: bool,
) -> Vec<i32> {
    let qlen = seq.len();
    let mut local_nm = vec![0i32; qlen];
    if qlen == 0 {
        return local_nm;
    }
    let q = |i: usize| qual.get(i).copied().unwrap_or(0) as i32;
    let poly_adj = 1.0f64; // homopoly_fix off

    // --- adj_qual: accumulate the local quality-deficit into local_nm
    let qhalo = 8usize;
    let qhalop = 2usize;
    let mut qmin = q(0);
    let mut qminp = q(0);
    let base0 = seq[0];
    for i in 1..qlen {
        if seq[i] != base0 {
            break;
        }
        if i < qhalop && qminp > q(i) {
            qminp = q(i);
        }
    }
    let mut i = 0usize;
    while i < qlen && i < qhalo {
        if qmin > q(i) {
            qmin = q(i);
        }
        i += 1;
    }
    while i + qhalo < qlen {
        // homopoly_fix off => polyl == polyr == 0, pl == 0
        let t = if mode_bayes116 {
            (q(i) + 5 * qmin) / 4
        } else {
            q(i) / 3 + (qminp as f64 * poly_adj) as i32
        };
        if t < q(i) {
            local_nm[i] += q(i) - t;
        }
        qminp = q(i);
        // inner k-loop over [max(0,i-qhalop)..=min(0,i+qhalop)] is empty
        // for i>qhalop (polyl=polyr=0), matching upstream.
        if qmin > q(i + qhalo) {
            qmin = q(i + qhalo);
        } else if qmin <= q(i - qhalo) {
            qmin = 99;
            for j in (i - qhalo + 1)..=(i + qhalo) {
                if qmin > q(j) {
                    qmin = q(j);
                }
            }
        }
        i += 1;
    }
    while i < qlen {
        let t = if mode_bayes116 {
            (q(i) + 5 * qmin) / 4
        } else {
            q(i) / 3 + (qminp as f64 * poly_adj) as i32
        };
        if t < q(i) {
            local_nm[i] += q(i) - t;
        }
        i += 1;
    }

    // --- homopolymer run length into the high 8 bits
    let mut i = 0usize;
    while i < qlen {
        let b = seq[i];
        let mut j = i + 1;
        while j < qlen && seq[j] == b {
            j += 1;
        }
        let mut poly = (j - i - 1) as i32;
        if poly > 100 {
            poly = 100;
        }
        // HALO == 0 => k in [i, j)
        for k in i..j {
            let cur = local_nm[k];
            local_nm[k] = (poly.max(cur >> 24) << 24) | (cur & ((1 << 24) - 1));
        }
        i = j;
    }

    // --- soft-clip cost at the read ends
    let is_sc = |idx: usize| -> bool {
        cigar
            .get(idx)
            .map(|&(op, _)| {
                op == 4 // S
                    || (op == 5 // H then S
                        && cigar.len() > 1
                        && cigar.get(idx + 1).is_some_and(|&(o, _)| o == 4))
            })
            .unwrap_or(false)
    };
    let halo = nm_halo;
    if !cigar.is_empty() && (cigar[0].0 == 4 || is_sc(0)) {
        let mut k = 0i64;
        while k < halo && (k as usize) < qlen {
            local_nm[k as usize] += sc_cost;
            k += 1;
        }
        while k < halo * 2 && (k as usize) < qlen {
            local_nm[k as usize] += sc_cost >> 1;
            k += 1;
        }
    }
    let last = cigar.len().wrapping_sub(1);
    let tail_sc = !cigar.is_empty()
        && (cigar[last].0 == 4
            || (cigar[last].0 == 5 && cigar.len() > 1 && cigar[last - 1].0 == 4));
    if tail_sc {
        let mut k = qlen as i64 - 1;
        while k >= qlen as i64 - halo && k >= 0 {
            local_nm[k as usize] += sc_cost;
            k -= 1;
        }
        while k >= qlen as i64 - halo * 2 && k >= 0 {
            local_nm[k as usize] += sc_cost >> 1;
            k -= 1;
        }
    }

    // --- MD walk: haloed mismatch cost (upstream does NOT advance `pos`
    // past the mismatch base; deletions `^...` are skipped).
    if let Some(md) = md {
        let mut pos: i64 = 0;
        let mut p = 0usize;
        while p < md.len() {
            let c = md[p];
            if c.is_ascii_digit() {
                let mut n: i64 = 0;
                while p < md.len() && md[p].is_ascii_digit() {
                    n = n * 10 + (md[p] - b'0') as i64;
                    p += 1;
                }
                pos += n;
                continue;
            }
            if c == b'^' {
                p += 1;
                while p < md.len() && !md[p].is_ascii_digit() {
                    p += 1;
                }
                continue;
            }
            let h = halo;
            let mut k = if pos - h * 2 >= 0 { pos - h * 2 } else { 0 };
            while k < pos - h && (k as usize) < qlen {
                local_nm[k as usize] += 5;
                k += 1;
            }
            while k < pos + h && (k as usize) < qlen {
                local_nm[k as usize] += 10;
                k += 1;
            }
            while k < pos + h * 2 && (k as usize) < qlen {
                local_nm[k as usize] += 5;
                k += 1;
            }
            p += 1;
        }
    }

    local_nm
}

/// `poly_len(pos)` / `nm_local(pos)` over a precomputed `local_nm`,
/// where `idx` is the upstream `pos - b->core.pos` (= `seq_offset+1`).
pub fn local_nm_poly(local_nm: &[i32], idx: i64) -> i32 {
    if idx >= 0 && (idx as usize) < local_nm.len() {
        local_nm[idx as usize] >> 24
    } else {
        0
    }
}

pub fn local_nm_score(local_nm: &[i32], idx: i64) -> f64 {
    if local_nm.is_empty() {
        return 0.0;
    }
    let mask = (1 << 24) - 1;
    if idx < 0 {
        (local_nm[0] & mask) as f64
    } else if idx as usize >= local_nm.len() {
        (local_nm[local_nm.len() - 1] & mask) as f64
    } else {
        (local_nm[idx as usize] & mask) as f64 / 10.0
    }
}

fn pileup_read_from_column(
    record: &TestPileupRecord,
    column: &TestPileupColumn,
    local_nm: &[i32],
) -> PileupRead {
    let qpos_quality = record.quality_scores.get(column.qpos).copied().unwrap_or(0);
    let quality = if column.base.is_some() {
        record.quality_scores.get(column.qpos).copied()
    } else {
        None
    };
    let (indel, insertion) = if !column.insertion_after.is_empty() {
        (
            i32::try_from(column.insertion_after.len()).unwrap_or(i32::MAX),
            column.insertion_after.clone(),
        )
    } else if column.deletion_after > 0 {
        (
            -i32::try_from(column.deletion_after).unwrap_or(i32::MAX),
            Vec::new(),
        )
    } else {
        (0, Vec::new())
    };

    PileupRead {
        name: record.name.clone(),
        mapping_quality: record.mapping_quality,
        is_reverse: record.is_reverse,
        base: column.base,
        quality,
        qpos_quality,
        qpos: column.qpos,
        is_deletion: column.is_deletion,
        is_refskip: column.is_refskip,
        is_head: column.is_head,
        is_tail: column.is_tail,
        indel,
        insertion,
        // Upstream indexes nm[] at `pos - b->core.pos` = seq_offset+1.
        bayes_poly: local_nm_poly(local_nm, column.qpos as i64 + 1),
        bayes_nm_local: local_nm_score(local_nm, column.qpos as i64 + 1),
    }
}

fn pileup_columns_from_inputs(inputs: &[Vec<TestPileupRecord>]) -> Vec<PileupColumn> {
    let mut sites: BTreeMap<(String, usize), Vec<Vec<PileupRead>>> = BTreeMap::new();

    for (input_index, records) in inputs.iter().enumerate() {
        for record in records {
            // Per-read consensus Bayesian precompute (default nm
            // params); reused across this record's columns.
            let local_nm = compute_local_nm(
                &record.sequence,
                &record.quality_scores,
                &record.cigar,
                record.md.as_deref(),
                50,
                60,
                false,
            );
            for column in &record.columns {
                let entry = sites
                    .entry((record.reference_name.clone(), column.reference_position))
                    .or_insert_with(|| vec![Vec::new(); inputs.len()]);
                entry[input_index].push(pileup_read_from_column(record, column, &local_nm));
            }
        }
    }

    sites
        .into_iter()
        .map(
            |((reference_name, zero_based_position), reads_by_input)| PileupColumn {
                reference_name,
                position: zero_based_position + 1,
                reads_by_input,
            },
        )
        .collect()
}

/// Builds synchronized pileup columns across multiple SAM/BAM inputs.
///
/// Each column reports, per input, the [`PileupRead`] entries overlapping a
/// reference position (HTSlib `bam_plp`/`sam_pileup`-shaped). Unmapped,
/// secondary, QC-fail and duplicate records are excluded, matching HTSlib's
/// default pileup behavior.
pub fn pileup_from_alignment_paths<P>(paths: &[P]) -> io::Result<Vec<PileupColumn>>
where
    P: AsRef<Path>,
{
    pileup_from_alignment_paths_with_options(paths, &PileupOptions::default())
}

/// Like [`pileup_from_alignment_paths`] with explicit record-selection options.
pub fn pileup_from_alignment_paths_with_options<P>(
    paths: &[P],
    options: &PileupOptions,
) -> io::Result<Vec<PileupColumn>>
where
    P: AsRef<Path>,
{
    let mut inputs = paths
        .iter()
        .map(|path| read_test_pileup_records_from_alignment_path(path, options))
        .collect::<io::Result<Vec<_>>>()?;

    if options.detect_overlaps {
        for records in &mut inputs {
            apply_overlap_correction(records);
        }
    }

    Ok(pileup_columns_from_inputs(&inputs))
}

/// Like [`pileup_from_alignment_paths`] but also accepts CRAM inputs, decoding
/// them against the supplied FASTA reference.
pub fn pileup_from_alignment_paths_with_reference<P, Q>(
    paths: &[P],
    reference_src: Q,
) -> io::Result<Vec<PileupColumn>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    pileup_from_alignment_paths_with_reference_and_options(
        paths,
        reference_src,
        &PileupOptions::default(),
    )
}

/// Like [`pileup_from_alignment_paths_with_reference`] with explicit
/// record-selection options.
pub fn pileup_from_alignment_paths_with_reference_and_options<P, Q>(
    paths: &[P],
    reference_src: Q,
    options: &PileupOptions,
) -> io::Result<Vec<PileupColumn>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let repository = cram_reference_repository_from_fasta_path(&reference_src)?;
    let reference_sequences = if options.apply_baq {
        Some(read_fasta_sequences_maybe_compressed(reference_src)?)
    } else {
        None
    };
    let mut inputs = paths
        .iter()
        .map(|path| {
            read_test_pileup_records_from_alignment_path_with_reference(
                path,
                repository.clone(),
                options,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;

    if let Some(reference_sequences) = reference_sequences.as_ref() {
        for records in &mut inputs {
            apply_baq_to_pileup_records(records, reference_sequences)?;
        }
    }

    if options.detect_overlaps {
        for records in &mut inputs {
            apply_overlap_correction(records);
        }
    }

    Ok(pileup_columns_from_inputs(&inputs))
}

/// Iterator form of [`pileup_from_alignment_paths`].
pub fn iter_pileup_from_alignment_paths<P>(
    paths: &[P],
) -> io::Result<std::vec::IntoIter<PileupColumn>>
where
    P: AsRef<Path>,
{
    pileup_from_alignment_paths(paths).map(|columns| columns.into_iter())
}

fn read_test_pileup_records_from_alignment_path_with_reference<P>(
    src: P,
    reference_sequence_repository: fasta::Repository,
    options: &PileupOptions,
) -> io::Result<Vec<TestPileupRecord>>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();
    if src
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("cram"))
    {
        let data_path = associated_data_path(src);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(reference_sequence_repository)
            .build_from_path(data_path)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        for result in reader.records(&header) {
            let record = result?;
            if let Some(record) =
                TestPileupRecord::try_from_record_with_options(&header, &record, options)?
            {
                records.push(record);
            }
        }

        Ok(records)
    } else {
        read_test_pileup_records_from_alignment_path(src, options)
    }
}

fn read_test_pileup_records_from_alignment_path<P>(
    src: P,
    options: &PileupOptions,
) -> io::Result<Vec<TestPileupRecord>>
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
            if let Some(record) =
                TestPileupRecord::try_from_record_with_options(&header, &record, options)?
            {
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
            if let Some(record) =
                TestPileupRecord::try_from_record_with_options(&header, &record, options)?
            {
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

/// Copies a BAM file to BAM output, but first rewrites the **header
/// text** through `transform` (header serialized to SAM text, the
/// callback returns the replacement text, which is parsed back). Used
/// for binary `@PG` insertion on BAM-input → BAM-output paths where
/// the records stay binary. Record bodies are streamed unchanged.
pub fn write_bam_from_path_transforming_header<P, W, F>(
    src: P,
    dst: W,
    transform: F,
) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
    F: FnOnce(&str) -> io::Result<String>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut header_writer = sam::io::Writer::new(Vec::new());
    header_writer.write_header(&header)?;
    let header_text = String::from_utf8(header_writer.into_inner())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let new_text = transform(&header_text)?;
    let new_header = sam::io::Reader::new(io::Cursor::new(new_text.into_bytes())).read_header()?;

    let mut writer = bam::io::Writer::new(dst);
    writer.write_header(&new_header)?;
    for result in reader.records() {
        let record = result?;
        writer.write_record(&new_header, &record)?;
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

/// CRAM encoder knobs forwarded from `samtools view -O cram,...`.
///
/// `records_per_slice` / `slices_per_container` correspond to the
/// `seqs_per_slice` / `slices_per_slice` output options; `None` keeps
/// the noodles default. `embed_reference` matches `embed_ref=1`.
#[derive(Clone, Copy, Debug, Default)]
pub struct CramWriteOptions {
    pub embed_reference: bool,
    pub records_per_slice: Option<usize>,
    pub slices_per_container: Option<usize>,
}

impl CramWriteOptions {
    /// Returns options that only embed the reference (the previous
    /// `embed_ref`-only behavior).
    pub fn embedded() -> Self {
        Self {
            embed_reference: true,
            ..Self::default()
        }
    }

    fn configure(&self, builder: cram::io::writer::Builder) -> cram::io::writer::Builder {
        let mut builder = builder.set_embed_reference(self.embed_reference);
        if let Some(n) = self.records_per_slice {
            builder = builder.set_records_per_slice(n);
        }
        if let Some(n) = self.slices_per_container {
            builder = builder.set_slices_per_container(n);
        }
        builder
    }
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

/// As [`write_cram_from_bam_path_with_reference`], with explicit CRAM
/// encoder options (`seqs_per_slice` / `slices_per_slice` /
/// `embed_ref`).
pub fn write_cram_from_bam_path_with_reference_and_options<P, Q, W>(
    src: P,
    reference_src: Q,
    options: CramWriteOptions,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;
    File::open(src).and_then(|reader| {
        write_cram_from_bam_reader_with_reference_repository_opts(
            reader,
            reference_sequence_repository,
            writer,
            options,
        )
    })
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
    write_cram_from_bam_reader_with_reference_repository_opts(
        reader,
        reference_sequence_repository,
        writer,
        CramWriteOptions::default(),
    )
}

fn write_cram_from_bam_reader_with_reference_repository_opts<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    writer: W,
    options: CramWriteOptions,
) -> io::Result<W>
where
    R: Read,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let mut reader = bam::io::Reader::new(reader);
    let header = reader.read_header()?;
    let mut writer = options
        .configure(
            cram::io::writer::Builder::default()
                .set_reference_sequence_repository(reference_sequence_repository),
        )
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

/// Writes indexed BAM records overlapping the given regions to BAM output using
/// noodles' multithreaded BGZF writer.
pub fn write_bam_regions_from_path_with_worker_count<P>(
    src: P,
    regions: &[Region],
    worker_count: NonZero<usize>,
) -> io::Result<Vec<u8>>
where
    P: AsRef<Path>,
{
    let index = read_associated_bam_index(&src)?;
    let data_path = associated_data_path(&src);
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(worker_count, Vec::new());
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query.records() {
            let record = result?;
            writer.write_record(&header, &record)?;
        }
    }

    let mut bgzf_writer = writer.into_inner();
    bgzf_writer.finish()
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

/// Writes BAM input records with all required flag bits set to BAM output using
/// noodles' multithreaded BGZF writer.
pub fn write_bam_records_with_required_flags_from_path_with_worker_count<P>(
    src: P,
    required_flags: u16,
    worker_count: NonZero<usize>,
) -> io::Result<Vec<u8>>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(worker_count, Vec::new());
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        let flags = u16::from(record.flags());
        if flags & required_flags == required_flags {
            writer.write_record(&header, &record)?;
        }
    }

    let mut bgzf_writer = writer.into_inner();
    bgzf_writer.finish()
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

/// Reads CRAM records into summaries **without an external reference**,
/// by synthesizing an all-`N` repository sized from the CRAM header's
/// `@SQ` lines.
///
/// CRAM stores each record's BAM flags, reference id, position, mapping
/// quality and read name in core/external data series independently of
/// the reference; only the *sequence* (and NM/MD-derived values) is
/// reconstructed against it. So for consumers that only need those
/// reference-independent fields — `idxstats`, `flagstat` — a synthetic
/// reference yields byte-identical counts while letting the noodles
/// CRAM decoder run without erroring on a missing reference. (The
/// decoded `sequence` bytes are *not* meaningful with this path.)
pub fn summarize_cram_records_from_path_synthesizing_reference<P>(
    src: P,
) -> io::Result<Vec<AlignmentRecordSummary>>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(&src);
    let header = read_cram_header_from_path(&data_path)?;
    let repository = synthetic_cram_reference_repository_from_header(&header);
    summarize_cram_records_from_path_with_reference_repository(src, repository)
}

fn synthetic_cram_reference_repository_from_header(header: &Header) -> fasta::Repository {
    let records: Vec<fasta::Record> = header
        .reference_sequences()
        .iter()
        .map(|(name, reference_sequence)| {
            let len = usize::from(reference_sequence.length());
            let definition = fasta::record::Definition::new(name.clone(), None);
            let sequence = fasta::record::Sequence::from(vec![b'N'; len]);
            fasta::Record::new(definition, sequence)
        })
        .collect();

    fasta::Repository::new(records)
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

/// As [`write_cram_from_path_with_reference`], with explicit CRAM
/// encoder options (`seqs_per_slice` / `slices_per_slice` /
/// `embed_ref`).
pub fn write_cram_from_path_with_reference_and_options<P, Q, W>(
    src: P,
    reference_src: Q,
    options: CramWriteOptions,
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
        write_cram_from_reader_with_reference_repository_opts(
            reader,
            reference_sequence_repository,
            writer,
            options,
        )
    })
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
    write_cram_from_reader_with_reference_repository_opts(
        reader,
        reference_sequence_repository,
        writer,
        CramWriteOptions::default(),
    )
}

fn write_cram_from_reader_with_reference_repository_opts<R, W>(
    reader: R,
    reference_sequence_repository: fasta::Repository,
    writer: W,
    options: CramWriteOptions,
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
    let mut writer = options
        .configure(
            cram::io::writer::Builder::default()
                .set_reference_sequence_repository(reference_sequence_repository),
        )
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

/// Writes indexed CRAM records overlapping the given regions to BAM output using
/// noodles' multithreaded BGZF writer.
pub fn write_cram_regions_as_bam_from_path_with_reference_with_worker_count<P, Q>(
    src: P,
    reference_src: Q,
    regions: &[Region],
    worker_count: NonZero<usize>,
) -> io::Result<Vec<u8>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
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
    let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(worker_count, Vec::new());
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for region in regions {
        let query = reader.query(&header, region)?;

        for result in query {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
    }

    let mut bgzf_writer = writer.into_inner();
    bgzf_writer.finish()
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

/// Writes CRAM records with all required flag bits set to BAM output using
/// noodles' multithreaded BGZF writer.
pub fn write_cram_records_with_required_flags_as_bam_from_path_with_reference_with_worker_count<
    P,
    Q,
>(
    src: P,
    reference_src: Q,
    required_flags: u16,
    worker_count: NonZero<usize>,
) -> io::Result<Vec<u8>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    File::open(data_path).and_then(|reader| {
        write_cram_records_with_required_flags_as_bam_with_reference_repository_worker_count(
            reader,
            repository,
            required_flags,
            worker_count,
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

fn write_cram_records_with_required_flags_as_bam_with_reference_repository_worker_count<R>(
    reader: R,
    repository: fasta::Repository,
    required_flags: u16,
    worker_count: NonZero<usize>,
) -> io::Result<Vec<u8>>
where
    R: Read,
{
    use sam::alignment::io::Write as _;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(reader);
    let header = reader.read_header()?;
    let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(worker_count, Vec::new());
    let mut writer = bam::io::Writer::from(bgzf_writer);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags().bits();
        if flags & required_flags == required_flags {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    let mut bgzf_writer = writer.into_inner();
    bgzf_writer.finish()
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

/// As [`write_cram_from_sam_reader_with_reference`], but **embeds**
/// each mapped slice's reference span in-container
/// (`samtools view -O cram,embed_ref=1`) so the CRAM decodes with no
/// external reference.
pub fn write_cram_from_sam_reader_with_reference_embedded<R, Q, W>(
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

    write_cram_from_sam_reader_with_reference_repository_opts(
        reader,
        reference_sequence_repository,
        writer,
        CramWriteOptions::embedded(),
    )
}

/// As [`write_cram_from_sam_reader_with_reference`], with explicit
/// CRAM encoder options (`seqs_per_slice` / `slices_per_slice` /
/// `embed_ref`).
pub fn write_cram_from_sam_reader_with_reference_and_options<R, Q, W>(
    reader: &mut sam::io::Reader<R>,
    reference_src: Q,
    options: CramWriteOptions,
    writer: W,
) -> io::Result<W>
where
    R: BufRead,
    Q: AsRef<Path>,
    W: Write,
{
    let reference_sequence_repository = cram_reference_repository_from_fasta_path(reference_src)?;

    write_cram_from_sam_reader_with_reference_repository_opts(
        reader,
        reference_sequence_repository,
        writer,
        options,
    )
}

/// As [`write_cram_from_sam_path_with_reference`], with explicit CRAM
/// encoder options.
pub fn write_cram_from_sam_path_with_reference_and_options<P, Q, W>(
    src: P,
    reference_src: Q,
    options: CramWriteOptions,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let file = File::open(src)?;
    let mut reader = sam::io::Reader::new(io::BufReader::new(file));
    write_cram_from_sam_reader_with_reference_and_options(
        &mut reader,
        reference_src,
        options,
        writer,
    )
}

/// As [`write_cram_from_sam_path_with_reference`], but embeds the
/// reference in-container (`embed_ref=1`).
pub fn write_cram_from_sam_path_with_reference_embedded<P, Q, W>(
    src: P,
    reference_src: Q,
    writer: W,
) -> io::Result<W>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    W: Write,
{
    let file = File::open(src)?;
    let mut reader = sam::io::Reader::new(io::BufReader::new(file));
    write_cram_from_sam_reader_with_reference_embedded(&mut reader, reference_src, writer)
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
    write_cram_from_sam_reader_with_reference_repository_opts(
        reader,
        reference_sequence_repository,
        writer,
        CramWriteOptions::default(),
    )
}

fn write_cram_from_sam_reader_with_reference_repository_opts<R, W>(
    reader: &mut sam::io::Reader<R>,
    reference_sequence_repository: fasta::Repository,
    writer: W,
    options: CramWriteOptions,
) -> io::Result<W>
where
    R: BufRead,
    W: Write,
{
    use sam::alignment::io::Write as _;

    let header = reader.read_header()?;
    let mut writer = options
        .configure(
            cram::io::writer::Builder::default()
                .set_reference_sequence_repository(reference_sequence_repository),
        )
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

/// Queries CRAM records from a local file using its associated CRAI index and
/// a synthetic all-`N` reference repository sized from the CRAM header.
///
/// This is suitable only for consumers that do not inspect decoded read bases
/// or reference-derived tags. Flags, coordinates, mapping qualities, CIGAR,
/// read names, and quality scores are reference-independent CRAM fields.
pub fn query_cram_records_from_path_synthesizing_reference<P>(
    src: P,
    region: &Region,
) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(&src);
    let header = read_cram_header_from_path(&data_path)?;
    let repository = synthetic_cram_reference_repository_from_header(&header);

    query_cram_records_from_path_with_reference_repository(src, region, repository)
}

/// Builds a noodles FASTA reference repository for CRAM decoding from a local FASTA file.
pub fn cram_reference_repository_from_fasta_path<P>(
    reference_src: P,
) -> io::Result<fasta::Repository>
where
    P: AsRef<Path>,
{
    let reference_src = reference_src.as_ref();
    if path_is_bgzf_like(reference_src) {
        let records = read_fasta_sequences_maybe_compressed(reference_src)?
            .into_iter()
            .map(|(name, sequence)| {
                fasta::Record::new(
                    fasta::record::Definition::new(name, None),
                    fasta::record::Sequence::from(sequence),
                )
            })
            .collect::<Vec<_>>();
        return Ok(fasta::Repository::new(records));
    }

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

/// Reads **every** record of a CRAM file (no region/index required), decoding
/// against a FASTA reference and returning owned [`sam::alignment::RecordBuf`]
/// values with full sequence/quality/aux/flags preserved.
///
/// This is the non-region complement of
/// [`query_cram_records_from_path_with_reference`]; the `summarize_*` path
/// only yields coordinate summaries and discards per-record sequence/quality.
pub fn query_cram_records_all_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let repository = cram_reference_repository_from_fasta_path(reference_src)?;
    let data_path = associated_data_path(src);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;

    reader
        .records(&header)
        .map(|result| {
            result.and_then(|record| {
                sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)
            })
        })
        .collect()
}

/// Owning-iterator form of [`query_cram_records_all_from_path_with_reference`].
pub fn iter_cram_records_all_from_path_with_reference<P, Q>(
    src: P,
    reference_src: Q,
) -> io::Result<std::vec::IntoIter<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    query_cram_records_all_from_path_with_reference(src, reference_src).map(Vec::into_iter)
}

/// Reads **every** record of a CRAM file with **no external reference**,
/// using an empty FASTA repository. This decodes correctly for CRAMs
/// built with an embedded reference (`embed_ref`), where the reference
/// bases travel inside the container; reference-compressed CRAMs that
/// need an external reference will error (use the
/// `*_with_reference` variant for those).
///
/// Returns full [`sam::alignment::RecordBuf`] values (sequence,
/// quality, CIGAR, aux, flags) — the non-region, no-reference
/// complement needed by `samtools reference`'s MD path on CRAM input.
pub fn query_cram_records_all_from_path<P>(src: P) -> io::Result<Vec<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
{
    let data_path = associated_data_path(src);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(fasta::Repository::default())
        .build_from_path(data_path)?;
    let header = reader.read_header()?;

    reader
        .records(&header)
        .map(|result| {
            result.and_then(|record| {
                sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)
            })
        })
        .collect()
}

/// Owning-iterator form of [`query_cram_records_all_from_path`].
pub fn iter_cram_records_all_from_path<P>(
    src: P,
) -> io::Result<std::vec::IntoIter<sam::alignment::RecordBuf>>
where
    P: AsRef<Path>,
{
    query_cram_records_all_from_path(src).map(Vec::into_iter)
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

/// Writes a SAM text view of CRAM records using a synthetic all-`N`
/// reference built from the CRAM header.
///
/// This mirrors [`summarize_cram_records_from_path_synthesizing_reference`]:
/// it lets callers render CRAMs without an external reference when their
/// downstream behavior does not depend on true reference bases. It deliberately
/// does not synthesize MD/NM tags from the all-`N` reference.
pub fn view_cram_as_sam_text_from_path_synthesizing_reference_and_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    use sam::alignment::io::Write as _;

    let raw_header = read_raw_cram_header_text(&src)?;
    let data_path = associated_data_path(&src);
    let header = read_cram_header_from_path(&data_path)?;
    let repository = synthetic_cram_reference_repository_from_header(&header);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(raw_header.into_bytes());

    for result in reader.records(&header).take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
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

    // noodles' CRAM `query` yields every record in the slices that overlap the
    // region (slice-granular), so callers must still filter each record to the
    // requested interval — mirroring the SAM-output path. Without this the
    // count/metric callers over-count records near slice boundaries.
    let mut records = Vec::new();
    for result in query {
        let record = result?;
        if record_intersects_region(&record, reference_sequence_id, region) {
            records.push(record);
        }
    }
    Ok(records)
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
        sequences.insert(
            fasta_record_primary_name(record.name()).to_vec(),
            record.sequence().as_ref().to_vec(),
        );
    }

    Ok(sequences)
}

fn read_fasta_sequences_maybe_compressed<P>(src: P) -> io::Result<HashMap<Vec<u8>, Vec<u8>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    if !path_is_bgzf_like(path) {
        return read_fasta_sequences(path);
    }

    let mut bytes = Vec::new();
    bgzf::io::Reader::new(File::open(path)?).read_to_end(&mut bytes)?;
    let mut sequences = HashMap::new();
    let mut name = None;
    let mut sequence = Vec::new();

    for line in bytes.split(|&b| b == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some(name) = name.take() {
                sequences.insert(name, std::mem::take(&mut sequence));
            }
            name = Some(fasta_record_primary_name(&line[1..]).to_vec());
        } else if line.first() != Some(&b';') {
            sequence.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }

    if let Some(name) = name {
        sequences.insert(name, sequence);
    }

    Ok(sequences)
}

fn path_is_bgzf_like(path: &Path) -> bool {
    path.extension().is_some_and(|extension| {
        extension.eq_ignore_ascii_case("gz") || extension.eq_ignore_ascii_case("bgz")
    })
}

fn fasta_record_primary_name(name: &[u8]) -> &[u8] {
    let end = name
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(name.len());
    &name[..end]
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

    let Some(reference_sequence_id) = record.reference_sequence_id() else {
        return Ok(());
    };
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
    let aux_values = record
        .data()
        .iter()
        .map(|result| {
            let (tag, value) = result?;
            Ok((tag.into(), summary_aux_payload(value)?.into_bytes()))
        })
        .collect::<io::Result<_>>()?;

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
        aux_values,
    })
}

fn summary_aux_payload(
    value: sam::alignment::record::data::field::Value<'_>,
) -> io::Result<String> {
    use sam::alignment::record::data::field::{Value, value::Array};

    let payload = match value {
        Value::Character(n) => char::from(n).to_string(),
        Value::Int8(n) => n.to_string(),
        Value::UInt8(n) => n.to_string(),
        Value::Int16(n) => n.to_string(),
        Value::UInt16(n) => n.to_string(),
        Value::Int32(n) => n.to_string(),
        Value::UInt32(n) => n.to_string(),
        Value::Float(n) => n.to_string(),
        Value::String(s) | Value::Hex(s) => String::from_utf8_lossy(s).into_owned(),
        Value::Array(Array::Int8(values)) => format!("c,{}", join_summary_array(values.iter())?),
        Value::Array(Array::UInt8(values)) => format!("C,{}", join_summary_array(values.iter())?),
        Value::Array(Array::Int16(values)) => format!("s,{}", join_summary_array(values.iter())?),
        Value::Array(Array::UInt16(values)) => format!("S,{}", join_summary_array(values.iter())?),
        Value::Array(Array::Int32(values)) => format!("i,{}", join_summary_array(values.iter())?),
        Value::Array(Array::UInt32(values)) => format!("I,{}", join_summary_array(values.iter())?),
        Value::Array(Array::Float(values)) => {
            format!("f,{}", join_summary_array(values.iter())?)
        }
    };

    Ok(payload)
}

fn join_summary_array<N>(iter: Box<dyn Iterator<Item = io::Result<N>> + '_>) -> io::Result<String>
where
    N: std::fmt::Display,
{
    let mut values = Vec::new();
    for result in iter {
        values.push(result?.to_string());
    }
    Ok(values.join(","))
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
    reference_sequence_id: usize,
    start: usize,
    mapping_quality: u8,
    is_reverse: bool,
    flags: u16,
    mate_reference_sequence_id: Option<usize>,
    mate_start: Option<usize>,
    template_length: i32,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    /// CIGAR ops as `(bam_op_code, len)` — M0 I1 D2 N3 S4 H5 P6 =7 X8.
    /// Needed by the consensus Bayesian `nm_init` precompute
    /// (soft-clip cost + the MD reference walk).
    cigar: Vec<(u8, usize)>,
    /// `MD` aux tag bytes, when present (drives the per-base local-NM).
    md: Option<Vec<u8>>,
    /// `NM` aux tag value, when present.
    nm: Option<i64>,
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
        Self::try_from_record_with_options(header, record, &PileupOptions::default())
    }

    fn try_from_record_with_options<R>(
        header: &Header,
        record: &R,
        options: &PileupOptions,
    ) -> io::Result<Option<Self>>
    where
        R: sam::alignment::Record + ?Sized,
    {
        let flags = record.flags()?;
        let flag_bits = u16::from(flags);

        if flag_bits & options.exclude_flags != 0
            || flag_bits & options.require_flags != options.require_flags
        {
            return Ok(None);
        }

        // HTSlib mpileup MPLP_NO_ORPHAN: drop paired-but-not-proper reads.
        if options.discard_orphans && flag_bits & 0x1 != 0 && flag_bits & 0x2 == 0 {
            return Ok(None);
        }

        let Some(reference_sequence_id) = record.reference_sequence_id(header).transpose()? else {
            return Ok(None);
        };
        let Some(alignment_start) = record.alignment_start().transpose()? else {
            return Ok(None);
        };

        let mapping_quality = record
            .mapping_quality()
            .transpose()?
            .map_or(255, |q| q.get());

        if mapping_quality < options.min_mapping_quality {
            return Ok(None);
        }

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

        // Mate info only feeds the heuristic overlap correction; never fail
        // the whole pileup if a record's mate fields are unreadable.
        let mate_reference_sequence_id = record
            .mate_reference_sequence_id(header)
            .and_then(Result::ok);
        let mate_start = record
            .mate_alignment_start()
            .and_then(Result::ok)
            .map(|p| usize::from(p) - 1);
        let template_length = record.template_length().unwrap_or(0);
        let start = usize::from(alignment_start) - 1;
        let mut columns = test_pileup_columns(record, start, &sequence)?;

        if let Some(first) = columns.first_mut() {
            first.is_head = true;
        }

        if let Some(last) = columns.last_mut() {
            last.is_tail = true;
        }

        // Capture CIGAR + MD/NM for the consensus Bayesian `nm_init`
        // precompute (additive; unused by the existing pileup paths).
        use sam::alignment::record::cigar::op::Kind as CigKind;
        let mut cigar = Vec::new();
        for result in record.cigar().iter() {
            let op = result?;
            let code = match op.kind() {
                CigKind::Match => 0u8,
                CigKind::Insertion => 1,
                CigKind::Deletion => 2,
                CigKind::Skip => 3,
                CigKind::SoftClip => 4,
                CigKind::HardClip => 5,
                CigKind::Pad => 6,
                CigKind::SequenceMatch => 7,
                CigKind::SequenceMismatch => 8,
            };
            cigar.push((code, op.len()));
        }
        use sam::alignment::record::data::field::{Tag, Value};
        let data = record.data();
        let md = match data.get(&Tag::MISMATCHED_POSITIONS).transpose()? {
            Some(Value::String(s)) => Some(s.to_string().into_bytes()),
            _ => None,
        };
        let nm = match data.get(&Tag::EDIT_DISTANCE).transpose()? {
            Some(v) => v.as_int(),
            None => None,
        };

        Ok(Some(Self {
            name,
            reference_name,
            reference_sequence_id,
            start,
            mapping_quality,
            is_reverse: flags.is_reverse_complemented(),
            flags: flag_bits,
            mate_reference_sequence_id,
            mate_start,
            template_length,
            sequence,
            quality_scores,
            cigar,
            md,
            nm,
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
                    let base = match sequence.get(qpos).copied() {
                        Some(base) => Some(base),
                        None if sequence.is_empty() => None,
                        None => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CIGAR qpos is out of range",
                            ));
                        }
                    };
                    columns.push(TestPileupColumn {
                        reference_position,
                        qpos,
                        base,
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
    use std::io::{Cursor, Write};
    use std::num::NonZero;
    use std::path::PathBuf;

    use crate::bgzf;

    use super::{
        PileupColumn, compute_local_nm, count_bam_records, count_bam_records_from_path,
        count_bam_records_in_region_from_path, count_bam_records_matching_filter,
        count_sam_records_from_path, count_sam_records_matching_filter,
        cram_reference_repository_from_fasta_path, local_nm_poly, local_nm_score,
        mpileup_baq_from_alignment, mpileup_indel_alignment_score, pileup_from_alignment_paths,
        pileup_from_alignment_paths_with_reference, query_bam_regions_from_path,
        read_bam_header_from_path, read_cram_header_from_path, read_sam_header_from_path,
        reference_sequence_count, synchronized_pileup_from_alignment_paths, view_bam_as_sam_text,
        view_sam_as_fastq_split_text_from_reader_with_flag_filter_and_suffix,
        write_bam_from_sam_reader,
        write_bam_records_with_required_flags_from_path_with_worker_count,
        write_bam_regions_from_path, write_bam_regions_from_path_with_worker_count,
    };

    fn column_at(columns: &[PileupColumn], position: usize) -> &PileupColumn {
        columns
            .iter()
            .find(|column| column.position == position)
            .unwrap_or_else(|| panic!("no pileup column at position {position}"))
    }

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
    fn test_write_bam_from_sam_reader_resolves_reference_alias() {
        let sam =
            b"@HD\tVN:1.6\n@SQ\tSN:r3\tLN:50\tAN:ref3\nr1\t0\tref3\t1\t30\t1M\t*\t0\t0\tA\t!\n";
        let bam_data = write_bam_from_sam_reader(Cursor::new(sam), Vec::new()).unwrap();

        let text = view_bam_as_sam_text(Cursor::new(bam_data), None).unwrap();

        assert!(text.contains("\n@SQ\tSN:r3\tLN:50\tAN:ref3\n"));
        assert!(text.contains("\nr1\t0\tr3\t1\t30\t1M\t*\t0\t0\tA\t!"));
    }

    #[test]
    fn test_filter_expression_flag_names_match_htslib_symbols() {
        let sam = concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "proper\t99\tref\t1\t30\t4M\t=\t20\t19\tACGT\t!!!!\n",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t!!\n",
            "plain\t0\tref\t2\t30\t4M\t*\t0\t0\tTGCA\t####\n",
        )
        .as_bytes();

        assert_eq!(
            count_sam_records_matching_filter(Cursor::new(sam), "flag.proper_pair").unwrap(),
            1
        );
        assert_eq!(
            count_sam_records_matching_filter(Cursor::new(sam), "flag.unmap").unwrap(),
            1
        );

        let bam_data = write_bam_from_sam_reader(Cursor::new(sam), Vec::new()).unwrap();
        assert_eq!(
            count_bam_records_matching_filter(Cursor::new(bam_data), "flag.proper_pair").unwrap(),
            1
        );
    }

    #[test]
    fn test_filter_expression_cigar_derived_symbols_match_htslib_symbols() {
        let sam = concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "plain\t0\tref\t2\t30\t4M\t*\t0\t0\tTGCA\t####\n",
            "soft\t0\tref\t5\t30\t2S4M\t*\t0\t0\tAATGCA\t!!!!!!\n",
            "hard\t0\tref\t10\t30\t3H4M\t*\t0\t0\tACGT\t!!!!\n",
        )
        .as_bytes();

        let bam_data = write_bam_from_sam_reader(Cursor::new(sam), Vec::new()).unwrap();
        assert_eq!(
            count_sam_records_matching_filter(Cursor::new(sam), "rlen>=4").unwrap(),
            3
        );
        assert_eq!(
            count_bam_records_matching_filter(Cursor::new(&bam_data), "rlen>=4").unwrap(),
            3
        );
        assert_eq!(
            count_sam_records_matching_filter(Cursor::new(sam), "endpos>=13").unwrap(),
            1
        );
        assert_eq!(
            count_bam_records_matching_filter(Cursor::new(&bam_data), "endpos>=13").unwrap(),
            1
        );
        assert_eq!(
            count_bam_records_matching_filter(Cursor::new(&bam_data), "sclen>0").unwrap(),
            1
        );
        assert_eq!(
            count_bam_records_matching_filter(Cursor::new(bam_data), "hclen>0").unwrap(),
            1
        );
    }

    #[test]
    fn test_filter_expression_mate_and_reference_symbols_match_htslib_symbols() {
        let sam = concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "@SQ\tSN:alt\tLN:100\n",
            "pair\t99\tref\t2\t30\t4M\t=\t20\t22\tTGCA\t####\n",
            "pair\t147\tref\t20\t30\t4M\t=\t2\t-22\tACGT\t!!!!\n",
            "other\t0\talt\t5\t30\t4M\t*\t0\t0\tNNNN\t!!!!\n",
        )
        .as_bytes();

        let bam_data = write_bam_from_sam_reader(Cursor::new(sam), Vec::new()).unwrap();

        for expr in [
            "mpos>0",
            "pnext>0",
            "tlen!=0",
            "rnext==\"ref\"",
            "mrname==\"ref\"",
            "refid==0",
            "mrefid==0",
            "ncigar==1",
        ] {
            let expected = if expr == "ncigar==1" { 3 } else { 2 };
            assert_eq!(
                count_sam_records_matching_filter(Cursor::new(sam), expr).unwrap(),
                expected,
                "SAM expression {expr}"
            );
            assert_eq!(
                count_bam_records_matching_filter(Cursor::new(&bam_data), expr).unwrap(),
                expected,
                "BAM expression {expr}"
            );
        }
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
    fn test_write_bam_regions_from_path_with_worker_count() {
        let path = fixture("htslib/test/range.bam");
        let regions = ["CHROMOSOME_II:2980-2980".parse().unwrap()];

        let bam_data = write_bam_regions_from_path_with_worker_count(
            &path,
            &regions,
            NonZero::new(2).unwrap(),
        )
        .unwrap();
        let mut reader = crate::bam::io::Reader::new(Cursor::new(bam_data));
        let header = reader.read_header().unwrap();
        assert!(reference_sequence_count(&header) > 0);

        let mut count = 0;
        for result in reader.records() {
            result.unwrap();
            count += 1;
        }

        assert_eq!(count, 1);
    }

    #[test]
    fn test_write_bam_records_with_required_flags_with_worker_count() {
        let path = fixture("htslib/test/range.bam");

        let bam_data = write_bam_records_with_required_flags_from_path_with_worker_count(
            &path,
            0,
            NonZero::new(2).unwrap(),
        )
        .unwrap();

        assert!(count_bam_records(Cursor::new(bam_data)).unwrap() > 0);
    }

    #[test]
    fn test_read_cram_header_and_records() {
        let path = fixture("htslib/test/range.cram");
        let header = read_cram_header_from_path(&path).unwrap();

        assert!(reference_sequence_count(&header) > 0);
    }

    #[test]
    fn test_cram_reference_repository_reads_bgzf_fasta_primary_name() {
        let path =
            std::env::temp_dir().join(format!("htslib-rs-bgzf-ref-{}.fa.gz", std::process::id()));
        {
            let file = std::fs::File::create(&path).unwrap();
            let mut writer = bgzf::io::Writer::new(file);
            writer.write_all(b">sq0 sq0:1-12\nACGTACGTACGT\n").unwrap();
            writer.finish().unwrap();
        }

        let repository = cram_reference_repository_from_fasta_path(&path).unwrap();
        let sequence = repository.get(b"sq0").transpose().unwrap().unwrap();
        assert_eq!(sequence.as_ref().as_ref(), b"ACGTACGTACGT");
        assert!(repository.get(b"sq0 sq0:1-12").is_none());

        let _ = std::fs::remove_file(path);
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

    #[test]
    fn test_pileup_iterator_reports_base_quality_indel_and_refskip() {
        let path =
            std::env::temp_dir().join(format!("htslib-rs-plp-cigar-{}-in.sam", std::process::id()));
        std::fs::write(
            &path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:30\n\
             del\t0\tsq0\t1\t60\t4M2D4M\t*\t0\t0\tACGTACGT\tIIIIJJJJ\n\
             ins\t0\tsq0\t1\t60\t4M2I4M\t*\t0\t0\tAAAACCGGGG\tIIIIKKLLLL\n\
             skip\t0\tsq0\t1\t60\t4M3N4M\t*\t0\t0\tTTTTGGGG\tIIIIJJJJ\n",
        )
        .unwrap();

        let columns = pileup_from_alignment_paths(std::slice::from_ref(&path))
            .inspect(|_| {
                let _ = std::fs::remove_file(&path);
            })
            .unwrap();

        // Column 4 (zero-based ref index 3): last matched base before each indel.
        let c4 = column_at(&columns, 4);
        assert_eq!(c4.total_depth(), 3);
        let reads = &c4.reads_by_input[0];

        let del = reads.iter().find(|r| r.indel < 0).unwrap();
        assert_eq!(del.base, Some(b'T'));
        assert_eq!(del.quality, Some(b'I' - 33));
        assert_eq!(del.qpos, 3);
        assert_eq!(del.indel, -2);
        assert!(del.insertion.is_empty());

        let ins = reads.iter().find(|r| r.indel > 0).unwrap();
        assert_eq!(ins.base, Some(b'A'));
        assert_eq!(ins.indel, 2);
        assert_eq!(ins.insertion, b"CC");

        let plain = reads.iter().find(|r| r.indel == 0).unwrap();
        assert_eq!(plain.base, Some(b'T'));
        assert_eq!(plain.indel, 0);

        // Column 5 (zero-based 4): deletion in `del`, refskip in `skip`.
        let c5 = column_at(&columns, 5);
        let c5_reads = &c5.reads_by_input[0];
        let deleted = c5_reads.iter().find(|r| r.is_deletion).unwrap();
        assert_eq!(deleted.base, None);
        assert_eq!(deleted.quality, None);
        let skipped = c5_reads.iter().find(|r| r.is_refskip).unwrap();
        assert_eq!(skipped.base, None);
        assert_eq!(skipped.quality, None);

        // Heads at the first column, tails at the read's final aligned column.
        let c1 = column_at(&columns, 1);
        assert!(c1.reads_by_input[0].iter().all(|r| r.is_head));
        let last = column_at(&columns, 10);
        assert!(last.reads_by_input[0].iter().any(|r| r.is_tail));
    }

    #[test]
    fn test_pileup_iterator_merges_multiple_inputs() {
        let path =
            std::env::temp_dir().join(format!("htslib-rs-plp-merge-{}-in.sam", std::process::id()));
        std::fs::write(
            &path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:20\n\
             a\t0\tsq0\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let columns = pileup_from_alignment_paths(&[path.clone(), path.clone()])
            .inspect(|_| {
                let _ = std::fs::remove_file(&path);
            })
            .unwrap();

        assert_eq!(columns.len(), 4);
        for column in &columns {
            assert_eq!(column.reads_by_input.len(), 2);
            assert_eq!(column.depth_of_input(0), 1);
            assert_eq!(column.depth_of_input(1), 1);
            assert_eq!(column.total_depth(), 2);
            assert_eq!(column.reads_by_input[0], column.reads_by_input[1]);
        }
    }

    #[test]
    fn test_pileup_options_filter_by_flags_and_mapq() {
        use super::{PileupOptions, pileup_from_alignment_paths_with_options};

        let path =
            std::env::temp_dir().join(format!("htslib-rs-plp-opts-{}-in.sam", std::process::id()));
        std::fs::write(
            &path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:20\n\
             keep\t0\tsq0\t1\t40\t4M\t*\t0\t0\tACGT\tIIII\n\
             dup\t1024\tsq0\t1\t40\t4M\t*\t0\t0\tACGT\tIIII\n\
             lowmq\t0\tsq0\t1\t5\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        // Default: excludes the duplicate, keeps both mapq-40 and mapq-5.
        let default_cols = pileup_from_alignment_paths_with_options(
            std::slice::from_ref(&path),
            &PileupOptions::default(),
        )
        .unwrap();
        assert_eq!(column_at(&default_cols, 1).total_depth(), 2);

        // Raise the mapq floor: only the mapq-40 read survives.
        let mq_cols = pileup_from_alignment_paths_with_options(
            std::slice::from_ref(&path),
            &PileupOptions {
                min_mapping_quality: 10,
                ..PileupOptions::default()
            },
        )
        .unwrap();
        assert_eq!(column_at(&mq_cols, 1).total_depth(), 1);

        // Empty exclude mask: the duplicate is now included.
        let cols = pileup_from_alignment_paths_with_options(
            std::slice::from_ref(&path),
            &PileupOptions {
                exclude_flags: 0,
                ..PileupOptions::default()
            },
        )
        .inspect(|_| {
            let _ = std::fs::remove_file(&path);
        })
        .unwrap();
        assert_eq!(column_at(&cols, 1).total_depth(), 3);
    }

    #[test]
    fn test_pileup_overlap_removal_and_orphan_filter() {
        use super::{PileupOptions, pileup_from_alignment_paths_with_options};

        let path =
            std::env::temp_dir().join(format!("htslib-rs-plp-olap-{}-in.sam", std::process::id()));
        // `p`: proper FR pair whose mates overlap over sq0:5-8.
        // `o`: paired but not a proper pair (an orphan).
        std::fs::write(
            &path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:sq0\tLN:50\n\
             p\t99\tsq0\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII\n\
             o\t65\tsq0\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
             p\t147\tsq0\t5\t60\t8M\t=\t1\t-12\tACGTACGT\tIIIIIIII\n",
        )
        .unwrap();

        // Overlap removal on (default): at sq0:5 the two `p` mates collapse —
        // one quality zeroed, the survivor boosted to the capped sum (80).
        let cols = pileup_from_alignment_paths_with_options(
            std::slice::from_ref(&path),
            &PileupOptions {
                discard_orphans: true,
                ..PileupOptions::default()
            },
        )
        .unwrap();
        let c5 = column_at(&cols, 5);
        let p_reads: Vec<_> = c5.reads_by_input[0]
            .iter()
            .filter(|r| r.name.as_deref() == Some(b"p"))
            .collect();
        assert_eq!(p_reads.len(), 2);
        let mut quals: Vec<u8> = p_reads.iter().map(|r| r.qpos_quality).collect();
        quals.sort_unstable();
        assert_eq!(quals, vec![0, 80]);
        // The orphan `o` is dropped everywhere when discard_orphans is set.
        assert!(cols.iter().all(|c| {
            c.reads_by_input[0]
                .iter()
                .all(|r| r.name.as_deref() != Some(b"o"))
        }));

        // Orphan retained when discard_orphans is cleared.
        let cols = pileup_from_alignment_paths_with_options(
            std::slice::from_ref(&path),
            &PileupOptions {
                discard_orphans: false,
                ..PileupOptions::default()
            },
        )
        .inspect(|_| {
            let _ = std::fs::remove_file(&path);
        })
        .unwrap();
        assert!(cols.iter().any(|c| {
            c.reads_by_input[0]
                .iter()
                .any(|r| r.name.as_deref() == Some(b"o"))
        }));
    }

    #[test]
    fn test_cram_all_records_match_bam_equivalent() {
        use super::query_cram_records_all_from_path_with_reference;
        use crate::sam::alignment::RecordBuf;

        let cram = fixture("htslib/test/range.cram");
        let reference = fixture("htslib/test/ce.fa");
        let bam = fixture("htslib/test/range.bam");

        let cram_records =
            query_cram_records_all_from_path_with_reference(&cram, &reference).unwrap();

        let mut reader = crate::bam::io::Reader::new(std::fs::File::open(&bam).unwrap());
        let bam_header = reader.read_header().unwrap();
        let bam_records: Vec<RecordBuf> = reader
            .records()
            .map(|r| RecordBuf::try_from_alignment_record(&bam_header, &r.unwrap()).unwrap())
            .collect();

        assert!(!cram_records.is_empty());
        assert_eq!(cram_records.len(), bam_records.len());

        for (c, b) in cram_records.iter().zip(&bam_records) {
            assert_eq!(c.name(), b.name());
            assert_eq!(c.flags(), b.flags());
            assert_eq!(c.alignment_start(), b.alignment_start());
            assert_eq!(c.sequence().as_ref(), b.sequence().as_ref());
            assert_eq!(c.quality_scores().as_ref(), b.quality_scores().as_ref());
            // NM is reference-derived and not stored in CRAM; noodles does
            // not synthesize it on decode, so it is recomputed by callers
            // (e.g. stats/reference) rather than asserted here.
        }
    }

    #[test]
    fn write_bam_from_path_transforming_header_rewrites_header_keeps_records() {
        use super::{summarize_bam_records_from_path, write_bam_from_path_transforming_header};
        let bam = fixture("htslib/test/range.bam");
        let dir = std::env::temp_dir().join(format!("htslib-rs-bamhdr-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("out.bam");

        let dst = std::fs::File::create(&out).unwrap();
        write_bam_from_path_transforming_header(&bam, dst, |text| {
            // Append a @PG-shaped line before the first non-@ line.
            Ok(format!("{text}@PG\tID:probe\tPN:probe\n"))
        })
        .unwrap();

        let before = summarize_bam_records_from_path(&bam).unwrap();
        let after = summarize_bam_records_from_path(&out).unwrap();
        assert_eq!(before.len(), after.len());
        for (a, b) in before.iter().zip(&after) {
            assert_eq!(a.flags_u16(), b.flags_u16());
            assert_eq!(a.reference_sequence_id(), b.reference_sequence_id());
            assert_eq!(a.alignment_start(), b.alignment_start());
        }
        let hdr = read_bam_header_from_path(&out).unwrap();
        assert!(
            hdr.programs().as_ref().keys().any(|k| {
                let id: &[u8] = k.as_ref();
                id == b"probe"
            }),
            "transformed header must carry the injected @PG"
        );
    }

    #[test]
    fn cram_summaries_without_reference_match_bam_flags_and_tids() {
        // The synthetic-reference path must yield the same record
        // count and per-record (flags, reference id, position) as the
        // BAM equivalent — the only fields idxstats/flagstat need —
        // without an external reference.
        use super::summarize_cram_records_from_path_synthesizing_reference;
        use super::{summarize_bam_records_from_path, summarize_cram_records_from_path};

        let cram = fixture("htslib/test/range.cram");
        let bam = fixture("htslib/test/range.bam");

        // The plain no-reference path errors on reference-compressed
        // CRAM (noodles eagerly resolves), which is exactly why the
        // synthesizing variant exists.
        assert!(summarize_cram_records_from_path(&cram).is_err());

        let cram_records = summarize_cram_records_from_path_synthesizing_reference(&cram).unwrap();
        let bam_records = summarize_bam_records_from_path(&bam).unwrap();

        assert!(!cram_records.is_empty());
        assert_eq!(cram_records.len(), bam_records.len());
        for (c, b) in cram_records.iter().zip(&bam_records) {
            assert_eq!(c.flags_u16(), b.flags_u16());
            assert_eq!(c.reference_sequence_id(), b.reference_sequence_id());
            assert_eq!(c.alignment_start(), b.alignment_start());
        }
    }

    #[test]
    fn query_cram_region_synthesizing_reference_matches_reference_backed_core_fields() {
        use super::{
            query_cram_records_from_path_synthesizing_reference,
            query_cram_records_from_path_with_reference,
        };
        use crate::core::Region;
        use crate::sam::alignment::record::Cigar as _;

        let cram = fixture("htslib/test/range.cram");
        let reference = fixture("htslib/test/ce.fa");
        let region: Region = "CHROMOSOME_II:2980-2980".parse().unwrap();

        let synthetic =
            query_cram_records_from_path_synthesizing_reference(&cram, &region).unwrap();
        let reference_backed =
            query_cram_records_from_path_with_reference(&cram, &region, &reference).unwrap();

        assert!(!synthetic.is_empty());
        assert_eq!(synthetic.len(), reference_backed.len());
        for (synthetic, reference_backed) in synthetic.iter().zip(&reference_backed) {
            assert_eq!(synthetic.name(), reference_backed.name());
            assert_eq!(synthetic.flags(), reference_backed.flags());
            assert_eq!(
                synthetic.reference_sequence_id(),
                reference_backed.reference_sequence_id()
            );
            assert_eq!(
                synthetic.alignment_start(),
                reference_backed.alignment_start()
            );
            assert_eq!(
                synthetic.mapping_quality(),
                reference_backed.mapping_quality()
            );
            assert_eq!(
                synthetic.quality_scores().as_ref(),
                reference_backed.quality_scores().as_ref()
            );
            let synthetic_cigar = synthetic
                .cigar()
                .iter()
                .map(|op| {
                    let op = op.unwrap();
                    (op.kind(), op.len())
                })
                .collect::<Vec<_>>();
            let reference_cigar = reference_backed
                .cigar()
                .iter()
                .map(|op| {
                    let op = op.unwrap();
                    (op.kind(), op.len())
                })
                .collect::<Vec<_>>();
            assert_eq!(synthetic_cigar, reference_cigar);
        }
    }

    #[test]
    fn query_cram_records_all_from_path_errors_on_reference_compressed_cram() {
        // The no-external-reference all-record reader is for
        // embed_ref CRAM; a reference-compressed CRAM (range.cram)
        // must error cleanly rather than silently mis-decoding.
        // (The positive embed_ref path is proven by the samtools-rs
        // `reference` CRAM integration test, whose fixture is an
        // embed_ref CRAM not shipped in htslib-rs/htslib/test.)
        use super::query_cram_records_all_from_path;
        let cram = fixture("htslib/test/range.cram");
        assert!(query_cram_records_all_from_path(&cram).is_err());
    }

    #[test]
    fn test_pileup_iterator_cram_matches_bam() {
        let bam = fixture("htslib/test/range.bam");
        let cram = fixture("htslib/test/range.cram");
        let reference = fixture("htslib/test/ce.fa");

        let mut from_bam = pileup_from_alignment_paths(std::slice::from_ref(&bam)).unwrap();
        let mut from_cram =
            pileup_from_alignment_paths_with_reference(std::slice::from_ref(&cram), &reference)
                .unwrap();

        // The consensus-Bayesian `bayes_*` fields derive from the `MD`
        // tag, which CRAM regenerates and BAM carries verbatim, so they
        // can legitimately differ by container. They are not part of
        // "does the CRAM pileup match the BAM pileup"; normalise them.
        let strip = |cols: &mut Vec<PileupColumn>| {
            for c in cols {
                for input in &mut c.reads_by_input {
                    for r in input {
                        r.bayes_poly = 0;
                        r.bayes_nm_local = 0.0;
                    }
                }
            }
        };
        strip(&mut from_bam);
        strip(&mut from_cram);

        assert!(!from_bam.is_empty());
        assert_eq!(from_bam, from_cram);
    }

    #[test]
    fn compute_local_nm_packs_poly_and_md_mismatch_cost() {
        // 10bp read, qual 40, no clips, MD "4A5" -> one mismatch at
        // reference offset 4. Defaults: nm_halo=50, sc_cost=60.
        let seq = b"ACGTACGTAC";
        let qual = [40u8; 10];
        let cigar = [(0u8, 10usize)]; // 10M
        let lnm = compute_local_nm(seq, &qual, &cigar, Some(b"4A5"), 50, 60, false);
        assert_eq!(lnm.len(), 10);

        // Homopolymer high-8-bits: this seq has no runs (poly 0 every
        // base), so >>24 == 0 everywhere.
        for &v in &lnm {
            assert_eq!(v >> 24, 0);
        }
        // The single MD mismatch at pos 4 adds the +10 inner halo to
        // every base (halo=50 spans the whole 10bp read), plus the
        // adj_qual deficit. So all low-24 scores are > 0.
        for &v in &lnm {
            assert!(v & ((1 << 24) - 1) > 0, "mismatch halo + adj_qual");
        }
        // local_nm_score divides the masked value by 10; idx clamping.
        assert_eq!(local_nm_score(&lnm, -1), (lnm[0] & 0xff_ffff) as f64);
        assert_eq!(local_nm_score(&lnm, 99), (lnm[9] & 0xff_ffff) as f64);
        assert_eq!(local_nm_score(&lnm, 4), (lnm[4] & 0xff_ffff) as f64 / 10.0);

        // A homopolymer read: AAAAA -> first base sees a run of 4
        // following, poly=4 in the high bits for the whole run.
        let hp = compute_local_nm(b"AAAAA", &[30u8; 5], &[(0u8, 5)], None, 50, 60, false);
        assert_eq!(hp[0] >> 24, 4);
        assert_eq!(local_nm_poly(&hp, 0), 4);
        assert_eq!(local_nm_poly(&hp, 10), 0); // out of range -> 0

        // Empty read is handled.
        assert!(compute_local_nm(b"", &[], &[], None, 50, 60, false).is_empty());
    }

    #[test]
    fn cram_write_options_records_per_slice_partitions_containers() {
        use super::{
            CramWriteOptions, write_cram_from_sam_path_with_reference,
            write_cram_from_sam_path_with_reference_and_options,
        };
        use crate::cram;

        fn count_containers(buf: &[u8]) -> usize {
            let mut reader = cram::io::reader::Builder::default().build_from_reader(buf);
            reader.read_header().unwrap();

            let mut containers = 0;
            let mut container = cram::io::reader::Container::default();
            while reader.read_container(&mut container).unwrap() != 0 {
                containers += 1;
            }
            containers
        }

        let sam = fixture("htslib/test/ce#1000.sam");
        let reference = fixture("htslib/test/ce.fa");

        // Default: all ~1000 records collapse into one container.
        let default_buf =
            write_cram_from_sam_path_with_reference(&sam, &reference, Vec::new()).unwrap();
        assert_eq!(count_containers(&default_buf), 1);

        // seqs_per_slice=100 must cut a new slice/container every 100
        // records, so a 1000-record file yields multiple containers.
        let options = CramWriteOptions {
            records_per_slice: Some(100),
            ..CramWriteOptions::default()
        };
        let chunked_buf = write_cram_from_sam_path_with_reference_and_options(
            &sam,
            &reference,
            options,
            Vec::new(),
        )
        .unwrap();
        assert!(
            count_containers(&chunked_buf) > 1,
            "seqs_per_slice=100 should produce more than one container"
        );
    }
}
