//! HTSlib-compatible FASTQ index helpers backed by noodles FASTQ.

use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};

use noodles::fastq;

/// A FASTQ index.
pub type Index = fastq::fai::Index;

/// Builds a FASTQ index from a reader.
pub fn build_index<R>(reader: R) -> io::Result<Index>
where
    R: BufRead,
{
    let mut indexer = fastq::io::Indexer::new(reader);
    let mut index = Vec::new();

    while let Some(record) = indexer.index_record()? {
        index.push(record);
    }

    Ok(index)
}

/// Reads a FASTQ index from a reader.
pub fn read_index<R>(mut reader: R) -> io::Result<Index>
where
    R: BufRead,
{
    let mut buf = String::new();
    let mut index = Vec::new();

    loop {
        buf.clear();

        match fastq::fai::io::Reader::new(&mut reader).read_record(&mut buf)? {
            0 => break,
            _ => {
                let record = buf
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                index.push(record);
            }
        }
    }

    Ok(index)
}

/// Returns whether the index has a sequence with `name`.
pub fn has_sequence(index: &Index, name: &str) -> bool {
    find_record(index, name).is_some()
}

/// Returns a sequence name by index.
pub fn sequence_name(index: &Index, i: usize) -> Option<&str> {
    index.get(i).map(|record| record.name())
}

/// Returns the sequence length for `name`.
pub fn sequence_len(index: &Index, name: &str) -> Option<u64> {
    find_record(index, name).map(|record| record.length())
}

/// Returns the number of bases per sequence line for `name`.
pub fn line_length(index: &Index, name: &str) -> Option<u64> {
    find_record(index, name).map(|record| record.line_bases())
}

/// Fetches the indexed sequence and quality string for a single-line FASTQ record.
///
/// This is intentionally scoped to the behavior currently exposed by
/// `noodles-fastq`'s FAI indexer. Multiline FASTQ parity remains tracked by the
/// scripted FASTQ TODOs.
pub fn fetch_sequence_and_quality<R>(
    mut reader: R,
    index: &Index,
    name: &str,
) -> io::Result<Option<(Vec<u8>, Vec<u8>)>>
where
    R: Read + Seek,
{
    let Some(record) = find_record(index, name) else {
        return Ok(None);
    };

    let len = usize::try_from(record.length())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let sequence = read_at(&mut reader, record.sequence_offset(), len)?;
    let quality = read_at(&mut reader, record.quality_scores_offset(), len)?;

    Ok(Some((sequence, quality)))
}

/// Options for converting FASTA/FASTQ records to SAM.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FastxToSamOptions {
    /// Preserve tab-delimited auxiliary fields from the FASTA/FASTQ definition line.
    pub include_aux: bool,
    /// Optional auxiliary tag allow-list used when `include_aux` is enabled.
    pub aux_tags: Option<Vec<String>>,
    /// Parse CASAVA 1.8 read number, filter flag, and barcode fields.
    pub casava: bool,
    /// Auxiliary tag used for CASAVA barcodes.
    pub barcode_tag: Option<String>,
    /// Auxiliary tag used for barcode qualities.
    pub barcode_quality_tag: Option<String>,
    /// Use the second whitespace-delimited definition token as the read name when present.
    pub name2: bool,
    /// Extract an inline UMI barcode into this auxiliary tag.
    pub umi_tag: Option<String>,
    /// Attach this read group ID as an `RG:Z` auxiliary tag.
    pub read_group_id: Option<String>,
}

/// Options for converting SAM records to FASTA/FASTQ.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SamToFastxOptions {
    /// Preserve SAM auxiliary fields on the FASTA/FASTQ definition line.
    pub include_aux: bool,
    /// Render CASAVA 1.8 read number, filter flag, and barcode fields.
    pub casava: bool,
    /// Auxiliary tag used for CASAVA barcodes.
    pub barcode_tag: Option<String>,
    /// Append UMI barcodes from this auxiliary tag to the read name.
    pub umi_tag: Option<String>,
    /// Append `/1` or `/2` to paired read names.
    pub append_read_number: bool,
}

/// Converts FASTQ records to SAM text.
pub fn write_sam_from_fastq<R, W>(
    reader: R,
    writer: &mut W,
    options: &FastxToSamOptions,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
    for record in read_fastq_records(reader)? {
        write_fastx_sam_record(writer, &record, options)?;
    }

    Ok(())
}

/// Converts paired FASTQ records to SAM text.
pub fn write_sam_from_paired_fastq<R1, R2, W>(
    read1: R1,
    read2: R2,
    writer: &mut W,
    options: &FastxToSamOptions,
) -> io::Result<()>
where
    R1: BufRead,
    R2: BufRead,
    W: Write,
{
    let read1_records = read_fastq_records(read1)?;
    let read2_records = read_fastq_records(read2)?;

    if read1_records.len() != read2_records.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "paired FASTQ inputs have different record counts",
        ));
    }

    for (mut r1, mut r2) in read1_records.into_iter().zip(read2_records) {
        ensure_read_number(&mut r1.definition, 1);
        ensure_read_number(&mut r2.definition, 2);
        write_fastx_sam_record(writer, &r1, options)?;
        write_fastx_sam_record(writer, &r2, options)?;
    }

    Ok(())
}

/// Converts paired FASTQ records with optional index FASTQ records to SAM text.
pub fn write_sam_from_paired_fastq_with_indexes<R1, R2, I1, I2, W>(
    read1: R1,
    read2: R2,
    index1: Option<I1>,
    index2: Option<I2>,
    writer: &mut W,
    options: &FastxToSamOptions,
    index_on_both_reads: bool,
) -> io::Result<()>
where
    R1: BufRead,
    R2: BufRead,
    I1: BufRead,
    I2: BufRead,
    W: Write,
{
    let read1_records = read_fastq_records(read1)?;
    let read2_records = read_fastq_records(read2)?;
    let index1_records = index1.map(read_fastq_records).transpose()?;
    let index2_records = index2.map(read_fastq_records).transpose()?;

    let read_count = read1_records.len();
    if read_count != read2_records.len()
        || index1_records
            .as_ref()
            .is_some_and(|records| records.len() != read_count)
        || index2_records
            .as_ref()
            .is_some_and(|records| records.len() != read_count)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "input FASTQ files have different record counts",
        ));
    }

    for i in 0..read_count {
        let mut r1 = read1_records[i].clone();
        let mut r2 = read2_records[i].clone();
        ensure_read_number(&mut r1.definition, 1);
        ensure_read_number(&mut r2.definition, 2);

        let index_aux = index_aux_fields(
            index1_records.as_ref().and_then(|records| records.get(i)),
            index2_records.as_ref().and_then(|records| records.get(i)),
            options,
        );

        write_fastx_sam_record_with_extra_aux(writer, &r1, options, &index_aux)?;
        if index_on_both_reads {
            write_fastx_sam_record_with_extra_aux(writer, &r2, options, &index_aux)?;
        } else {
            write_fastx_sam_record(writer, &r2, options)?;
        }
    }

    Ok(())
}

/// Converts FASTQ records with optional index FASTQ records to SAM text.
pub fn write_sam_from_fastq_with_indexes<R, I1, I2, W>(
    reader: R,
    index1: Option<I1>,
    index2: Option<I2>,
    writer: &mut W,
    options: &FastxToSamOptions,
) -> io::Result<()>
where
    R: BufRead,
    I1: BufRead,
    I2: BufRead,
    W: Write,
{
    let records = read_fastq_records(reader)?;
    let index1_records = index1.map(read_fastq_records).transpose()?;
    let index2_records = index2.map(read_fastq_records).transpose()?;

    let read_count = records.len();
    if index1_records
        .as_ref()
        .is_some_and(|records| records.len() != read_count)
        || index2_records
            .as_ref()
            .is_some_and(|records| records.len() != read_count)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "input FASTQ files have different record counts",
        ));
    }

    for (i, record) in records.iter().enumerate() {
        let index_aux = index_aux_fields(
            index1_records.as_ref().and_then(|records| records.get(i)),
            index2_records.as_ref().and_then(|records| records.get(i)),
            options,
        );
        write_fastx_sam_record_with_extra_aux(writer, record, options, &index_aux)?;
    }

    Ok(())
}

/// Converts FASTA records to SAM text.
pub fn write_sam_from_fasta<R, W>(
    reader: R,
    writer: &mut W,
    options: &FastxToSamOptions,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
    for record in read_fasta_records(reader)? {
        write_fastx_sam_record(writer, &record, options)?;
    }

    Ok(())
}

/// Converts SAM records to FASTQ text.
pub fn write_fastq_from_sam<R, W>(
    reader: R,
    writer: &mut W,
    options: &SamToFastxOptions,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
    for record in read_sam_records(reader)? {
        let name = sam_fastx_name(&record, options);
        let aux = sam_fastx_aux(&record, options);

        writeln!(writer, "@{name}{aux}")?;
        writeln!(writer, "{}", record.sequence)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", record.quality)?;
    }

    Ok(())
}

/// Converts SAM records to FASTA text.
pub fn write_fasta_from_sam<R, W>(
    reader: R,
    writer: &mut W,
    options: &SamToFastxOptions,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
    for record in read_sam_records(reader)? {
        let name = sam_fastx_name(&record, options);
        let aux = sam_fastx_aux(&record, options);

        writeln!(writer, ">{name}{aux}")?;
        writeln!(writer, "{}", record.sequence)?;
    }

    Ok(())
}

fn find_record<'a>(index: &'a Index, name: &str) -> Option<&'a fastq::fai::Record> {
    index.iter().find(|record| record.name() == name)
}

fn read_at<R>(reader: &mut R, offset: u64, len: usize) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    let mut buf = vec![0; len];
    reader.seek(SeekFrom::Start(offset))?;
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct FastxRecord {
    definition: String,
    sequence: String,
    quality: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SamRecord {
    name: String,
    flags: u16,
    sequence: String,
    quality: String,
    aux: Vec<String>,
}

fn read_fastq_records<R>(mut reader: R) -> io::Result<Vec<FastxRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();
    let mut line = String::new();

    loop {
        line.clear();

        if reader.read_line(&mut line)? == 0 {
            break;
        }

        trim_newline(&mut line);

        if !line.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid FASTQ definition",
            ));
        }

        let definition = line[1..].to_string();
        let mut sequence = String::new();

        loop {
            line.clear();

            if reader.read_line(&mut line)? == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "missing FASTQ plus line",
                ));
            }

            trim_newline(&mut line);

            if line.starts_with('+') {
                break;
            }

            sequence.push_str(&line);
        }

        let mut quality = String::new();

        while quality.len() < sequence.len() {
            line.clear();

            if reader.read_line(&mut line)? == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "short FASTQ quality string",
                ));
            }

            trim_newline(&mut line);
            quality.push_str(&line);
        }

        if quality.len() != sequence.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ quality length differs from sequence length",
            ));
        }

        records.push(FastxRecord {
            definition,
            sequence,
            quality: Some(quality),
        });
    }

    Ok(records)
}

fn read_fasta_records<R>(mut reader: R) -> io::Result<Vec<FastxRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();
    let mut line = String::new();
    let mut current: Option<FastxRecord> = None;

    while reader.read_line(&mut line)? != 0 {
        trim_newline(&mut line);

        if let Some(definition) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }

            current = Some(FastxRecord {
                definition: definition.to_string(),
                sequence: String::new(),
                quality: None,
            });
        } else if let Some(record) = current.as_mut() {
            record.sequence.push_str(&line);
        } else if !line.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid FASTA sequence before definition",
            ));
        }

        line.clear();
    }

    if let Some(record) = current {
        records.push(record);
    }

    Ok(records)
}

fn write_fastx_sam_record<W>(
    writer: &mut W,
    record: &FastxRecord,
    options: &FastxToSamOptions,
) -> io::Result<()>
where
    W: Write,
{
    write_fastx_sam_record_with_extra_aux(writer, record, options, &[])
}

fn write_fastx_sam_record_with_extra_aux<W>(
    writer: &mut W,
    record: &FastxRecord,
    options: &FastxToSamOptions,
    extra_aux: &[String],
) -> io::Result<()>
where
    W: Write,
{
    let ParsedDefinition {
        name,
        read_number,
        filtered,
        aux,
    } = parse_definition(&record.definition, options);
    let mut flags = match read_number {
        Some(1) => 77,
        Some(2) => 141,
        _ => 4,
    };

    if filtered {
        flags |= 0x200;
    }

    let quality = record.quality.as_deref().unwrap_or("*");

    write!(
        writer,
        "{name}\t{flags}\t*\t0\t0\t*\t*\t0\t0\t{}\t{quality}",
        record.sequence
    )?;

    for field in aux {
        write!(writer, "\t{field}")?;
    }

    for field in extra_aux {
        write!(writer, "\t{field}")?;
    }

    if let Some(read_group_id) = options.read_group_id.as_deref() {
        write!(writer, "\tRG:Z:{read_group_id}")?;
    }

    writeln!(writer)
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ParsedDefinition {
    name: String,
    read_number: Option<u8>,
    filtered: bool,
    aux: Vec<String>,
}

fn parse_definition(definition: &str, options: &FastxToSamOptions) -> ParsedDefinition {
    let tokens: Vec<&str> = definition.split_whitespace().collect();
    let raw_name = if options.name2 {
        tokens
            .get(1)
            .or_else(|| tokens.first())
            .copied()
            .unwrap_or("")
    } else {
        tokens.first().copied().unwrap_or("")
    };

    let mut name = raw_name.to_string();
    let read_number = strip_read_number(&mut name);
    let mut casava_read_number = None;
    let mut filtered = false;
    let mut aux = if options.include_aux {
        fastx_aux_fields(&tokens, options)
    } else {
        Vec::new()
    };

    if options.casava
        && let Some(casava) = tokens.get(1).and_then(|field| parse_casava_field(field))
    {
        casava_read_number = casava.read_number;
        filtered = casava.filtered;

        if let Some(barcode) = casava.barcode.filter(|barcode| barcode != "0")
            && let Some(tag) = barcode_tag(options.barcode_tag.as_deref())
        {
            aux.push(format!("{tag}:Z:{barcode}"));
        }
    }

    if let Some(tag) = options.umi_tag.as_deref().filter(|tag| tag.len() == 2)
        && let Some((umi_name, umi)) = extract_umi(&name)
    {
        name = umi_name;
        aux.push(format!("{tag}:Z:{umi}"));
    }

    ParsedDefinition {
        name,
        read_number: read_number.or(casava_read_number),
        filtered,
        aux,
    }
}

fn fastx_aux_fields(tokens: &[&str], options: &FastxToSamOptions) -> Vec<String> {
    tokens
        .iter()
        .skip(1)
        .filter(|field| field.matches(':').count() >= 2)
        .filter(|field| {
            let tag = field.get(..2);
            let has_tag_separator = field.as_bytes().get(2) == Some(&b':');
            match (tag, has_tag_separator, options.aux_tags.as_ref()) {
                (Some(tag), true, Some(tags)) => tags.iter().any(|wanted| wanted == tag),
                (_, _, Some(_)) => false,
                _ => true,
            }
        })
        .map(|field| normalize_aux_field(field))
        .collect()
}

fn normalize_aux_field(field: &str) -> String {
    let mut parts = field.splitn(3, ':');
    let Some(tag) = parts.next() else {
        return field.to_string();
    };
    let Some(kind) = parts.next() else {
        return field.to_string();
    };
    let Some(value) = parts.next() else {
        return field.to_string();
    };

    if kind == "f" {
        return format!("{tag}:{kind}:{}", normalize_float_exponent(value));
    }

    field.to_string()
}

fn normalize_float_exponent(value: &str) -> String {
    let Some(e_pos) = value.find(['e', 'E']) else {
        return value.to_string();
    };
    let mantissa = &value[..=e_pos];
    let exponent = &value[e_pos + 1..];
    if exponent.starts_with(['+', '-']) {
        value.to_string()
    } else {
        format!("{mantissa}+{exponent}")
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct CasavaField {
    read_number: Option<u8>,
    filtered: bool,
    barcode: Option<String>,
}

fn parse_casava_field(field: &str) -> Option<CasavaField> {
    let mut parts = field.split(':');
    let read_number = match parts.next()? {
        "1" => Some(1),
        "2" => Some(2),
        _ => None,
    };
    let filtered = matches!(parts.next(), Some("Y"));
    let _control_number = parts.next();
    let barcode = parts.next().map(String::from);

    Some(CasavaField {
        read_number,
        filtered,
        barcode,
    })
}

fn barcode_tag(tag: Option<&str>) -> Option<&str> {
    match tag {
        Some(tag) if tag.len() == 2 => Some(tag),
        Some(_) => None,
        None => Some("BC"),
    }
}

fn barcode_quality_tag(tag: Option<&str>) -> Option<&str> {
    match tag {
        Some(tag) if tag.len() == 2 => Some(tag),
        Some(_) => None,
        None => Some("QT"),
    }
}

fn index_aux_fields(
    index1: Option<&FastxRecord>,
    index2: Option<&FastxRecord>,
    options: &FastxToSamOptions,
) -> Vec<String> {
    let indexes = [index1, index2].into_iter().flatten().collect::<Vec<_>>();
    if indexes.is_empty() {
        return Vec::new();
    }

    let sequence = indexes
        .iter()
        .map(|record| record.sequence.as_str())
        .collect::<Vec<_>>()
        .join("-");
    let quality = indexes
        .iter()
        .map(|record| record.quality.as_deref().unwrap_or("*"))
        .collect::<Vec<_>>()
        .join(" ");

    let mut aux = Vec::new();
    if let Some(tag) = barcode_tag(options.barcode_tag.as_deref()) {
        aux.push(format!("{tag}:Z:{sequence}"));
    }
    if let Some(tag) = barcode_quality_tag(options.barcode_quality_tag.as_deref()) {
        aux.push(format!("{tag}:Z:{quality}"));
    }
    aux
}

fn strip_read_number(name: &mut String) -> Option<u8> {
    match name.as_bytes() {
        bytes if bytes.ends_with(b"/1") => {
            name.truncate(name.len() - 2);
            Some(1)
        }
        bytes if bytes.ends_with(b"/2") => {
            name.truncate(name.len() - 2);
            Some(2)
        }
        _ => None,
    }
}

fn ensure_read_number(definition: &mut String, read_number: u8) {
    let mut tokens = definition.splitn(2, char::is_whitespace);
    let mut name = tokens.next().unwrap_or("").to_string();
    let rest = tokens.next();

    strip_read_number(&mut name);
    name.push('/');
    name.push(char::from(b'0' + read_number));

    *definition = if let Some(rest) = rest {
        format!("{name} {rest}")
    } else {
        name
    };
}

fn extract_umi(name: &str) -> Option<(String, String)> {
    let hash = name.rfind('#');
    let colon = name[..hash.unwrap_or(name.len())].rfind(':')?;
    let umi_end = hash.unwrap_or(name.len());

    if colon + 1 >= umi_end {
        return None;
    }

    let umi = &name[colon + 1..umi_end];

    if umi.as_bytes().iter().all(u8::is_ascii_digit) {
        return None;
    }

    let mut output_name = String::from(&name[..colon]);

    if let Some(hash) = hash {
        output_name.push_str(&name[hash..]);
    }

    let umi = umi.replace('+', "-");

    Some((output_name, umi))
}

fn read_sam_records<R>(mut reader: R) -> io::Result<Vec<SamRecord>>
where
    R: BufRead,
{
    let mut records = Vec::new();
    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        trim_newline(&mut line);

        if line.is_empty() || line.starts_with('@') {
            line.clear();
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 11 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SAM record has fewer than 11 fields",
            ));
        }

        let flags = fields[1]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        records.push(SamRecord {
            name: fields[0].to_string(),
            flags,
            sequence: fields[9].to_string(),
            quality: fields[10].to_string(),
            aux: fields[11..]
                .iter()
                .map(|field| (*field).to_string())
                .collect(),
        });

        line.clear();
    }

    Ok(records)
}

fn sam_fastx_name(record: &SamRecord, options: &SamToFastxOptions) -> String {
    let mut name = record.name.clone();

    if let Some(umi) = sam_umi(record, options.umi_tag.as_deref()) {
        insert_umi(&mut name, umi);
    }

    if options.append_read_number {
        if record.flags & 0x40 != 0 {
            name.push_str("/1");
        } else if record.flags & 0x80 != 0 {
            name.push_str("/2");
        }
    }

    name
}

fn sam_fastx_aux(record: &SamRecord, options: &SamToFastxOptions) -> String {
    if options.casava {
        let read_number = if record.flags & 0x80 != 0 { 2 } else { 1 };
        let filter = if record.flags & 0x200 != 0 { "Y" } else { "N" };
        let barcode = sam_barcode(record, options.barcode_tag.as_deref()).unwrap_or("0");

        return format!(" {read_number}:{filter}:0:{barcode}");
    }

    if options.include_aux && !record.aux.is_empty() {
        format!("\t{}", record.aux.join("\t"))
    } else {
        String::new()
    }
}

fn sam_barcode<'a>(record: &'a SamRecord, tag: Option<&str>) -> Option<&'a str> {
    let tag = barcode_tag(tag)?;
    let prefix = format!("{tag}:Z:");

    record
        .aux
        .iter()
        .find_map(|field| field.strip_prefix(&prefix))
}

fn sam_umi<'a>(record: &'a SamRecord, tag: Option<&str>) -> Option<&'a str> {
    let tag = tag.filter(|tag| tag.len() == 2).unwrap_or("RX");
    let prefix = format!("{tag}:Z:");

    record
        .aux
        .iter()
        .find_map(|field| field.strip_prefix(&prefix))
}

fn insert_umi(name: &mut String, umi: &str) {
    let umi = umi.replace('-', "+");

    if let Some(hash) = name.rfind('#') {
        name.insert_str(hash, &format!(":{umi}"));
    } else {
        name.push(':');
        name.push_str(&umi);
    }
}

fn trim_newline(line: &mut String) {
    if line.ends_with('\n') {
        line.pop();

        if line.ends_with('\r') {
            line.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{
        FastxToSamOptions, SamToFastxOptions, build_index, fetch_sequence_and_quality,
        has_sequence, line_length, read_index, sequence_len, sequence_name, write_fasta_from_sam,
        write_fastq_from_sam, write_sam_from_fasta, write_sam_from_fastq,
        write_sam_from_fastq_with_indexes, write_sam_from_paired_fastq,
        write_sam_from_paired_fastq_with_indexes,
    };

    const MINIMAL_FASTQ: &[u8] = include_bytes!("../../../htslib/test/fastq/minimal.fq");
    const SINGLE_FASTQ: &[u8] = include_bytes!("../../../htslib/test/fastq/single.fq");

    #[test]
    fn test_build_and_read_index() {
        let index = build_index(Cursor::new(MINIMAL_FASTQ)).unwrap();
        assert_eq!(index.len(), 1);

        let text = b"x\t1\t3\t1\t2\t7\n";
        assert_eq!(read_index(Cursor::new(text)).unwrap(), index);
    }

    #[test]
    fn test_index_lookup_helpers() {
        let index = build_index(Cursor::new(SINGLE_FASTQ)).unwrap();

        assert!(has_sequence(&index, "HS25_09827:2:1201:1505:59795#49"));
        assert!(!has_sequence(&index, "missing"));
        assert_eq!(
            sequence_name(&index, 0),
            Some("HS25_09827:2:1201:1505:59795#49")
        );
        assert_eq!(
            sequence_len(&index, "HS25_09827:2:1201:1505:59795#49"),
            Some(100)
        );
        assert_eq!(
            line_length(&index, "HS25_09827:2:1201:1505:59795#49"),
            Some(100)
        );
    }

    #[test]
    fn test_fetch_sequence_and_quality() {
        let index = build_index(Cursor::new(MINIMAL_FASTQ)).unwrap();
        let (sequence, quality) =
            fetch_sequence_and_quality(Cursor::new(MINIMAL_FASTQ), &index, "x")
                .unwrap()
                .unwrap();

        assert_eq!(sequence, b"A");
        assert_eq!(quality, b"+");

        assert!(
            fetch_sequence_and_quality(Cursor::new(MINIMAL_FASTQ), &index, "missing")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn test_write_sam_from_fastq() {
        let mut out = Vec::new();
        write_sam_from_fastq(
            Cursor::new(MINIMAL_FASTQ),
            &mut out,
            &FastxToSamOptions::default(),
        )
        .unwrap();

        assert_eq!(out, b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t+\n");
    }

    #[test]
    fn test_write_sam_from_fasta() {
        let mut out = Vec::new();
        write_sam_from_fasta(
            Cursor::new(b">x\nA\n"),
            &mut out,
            &FastxToSamOptions::default(),
        )
        .unwrap();

        assert_eq!(out, b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n");
    }

    #[test]
    fn test_write_sam_from_fastq_with_aux_tag_filter() {
        let mut out = Vec::new();
        let options = FastxToSamOptions {
            include_aux: true,
            aux_tags: Some(vec![String::from("XZ"), String::from("AA")]),
            ..Default::default()
        };

        write_sam_from_fastq(
            Cursor::new(b"@x\tXX:i:10\tXZ:i:20\tAA:Z:ok\nA\n+\n+\n"),
            &mut out,
            &options,
        )
        .unwrap();

        assert_eq!(out, b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t+\tXZ:i:20\tAA:Z:ok\n");
    }

    #[test]
    fn test_write_sam_from_fastq_normalizes_float_aux_exponent() {
        let mut out = Vec::new();
        let options = FastxToSamOptions {
            include_aux: true,
            ..Default::default()
        };

        write_sam_from_fastq(
            Cursor::new(b"@x\tFF:f:-1e20\nA\n+\n+\n"),
            &mut out,
            &options,
        )
        .unwrap();

        assert_eq!(out, b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t+\tFF:f:-1e+20\n");
    }

    #[test]
    fn test_write_sam_from_fastq_with_read_group() {
        let mut out = Vec::new();
        let options = FastxToSamOptions {
            read_group_id: Some(String::from("rg1")),
            ..Default::default()
        };

        write_sam_from_fastq(Cursor::new(MINIMAL_FASTQ), &mut out, &options).unwrap();

        assert_eq!(out, b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t+\tRG:Z:rg1\n");
    }

    #[test]
    fn test_write_sam_from_paired_fastq() {
        let mut out = Vec::new();
        write_sam_from_paired_fastq(
            Cursor::new(b"@x\nA\n+\n+\n"),
            Cursor::new(b"@x\nT\n+\n-\n"),
            &mut out,
            &FastxToSamOptions::default(),
        )
        .unwrap();

        assert_eq!(
            out,
            b"x\t77\t*\t0\t0\t*\t*\t0\t0\tA\t+\nx\t141\t*\t0\t0\t*\t*\t0\t0\tT\t-\n"
        );
    }

    #[test]
    fn test_write_sam_from_paired_fastq_with_indexes() {
        let mut out = Vec::new();
        write_sam_from_paired_fastq_with_indexes(
            Cursor::new(b"@x\nA\n+\n+\n"),
            Cursor::new(b"@x\nT\n+\n-\n"),
            Some(Cursor::new(b"@x\nCG\n+\n12\n")),
            Some(Cursor::new(b"@x\nTA\n+\n34\n")),
            &mut out,
            &FastxToSamOptions::default(),
            false,
        )
        .unwrap();

        assert_eq!(
            out,
            b"x\t77\t*\t0\t0\t*\t*\t0\t0\tA\t+\tBC:Z:CG-TA\tQT:Z:12 34\nx\t141\t*\t0\t0\t*\t*\t0\t0\tT\t-\n"
        );
    }

    #[test]
    fn test_write_sam_from_paired_fastq_with_custom_index_tags_on_both_reads() {
        let mut out = Vec::new();
        let options = FastxToSamOptions {
            barcode_tag: Some(String::from("OX")),
            barcode_quality_tag: Some(String::from("BZ")),
            ..Default::default()
        };

        write_sam_from_paired_fastq_with_indexes(
            Cursor::new(b"@x\nA\n+\n+\n"),
            Cursor::new(b"@x\nT\n+\n-\n"),
            Some(Cursor::new(b"@x\nCG\n+\n12\n")),
            None::<Cursor<&[u8]>>,
            &mut out,
            &options,
            true,
        )
        .unwrap();

        assert_eq!(
            out,
            b"x\t77\t*\t0\t0\t*\t*\t0\t0\tA\t+\tOX:Z:CG\tBZ:Z:12\nx\t141\t*\t0\t0\t*\t*\t0\t0\tT\t-\tOX:Z:CG\tBZ:Z:12\n"
        );
    }

    #[test]
    fn test_write_sam_from_fastq_with_indexes() {
        let mut out = Vec::new();
        write_sam_from_fastq_with_indexes(
            Cursor::new(b"@x/1\nA\n+\n+\n@y/2\nT\n+\n-\n"),
            Some(Cursor::new(b"@x\nCG\n+\n12\n@y\nTA\n+\n34\n")),
            None::<Cursor<&[u8]>>,
            &mut out,
            &FastxToSamOptions::default(),
        )
        .unwrap();

        assert_eq!(
            out,
            b"x\t77\t*\t0\t0\t*\t*\t0\t0\tA\t+\tBC:Z:CG\tQT:Z:12\ny\t141\t*\t0\t0\t*\t*\t0\t0\tT\t-\tBC:Z:TA\tQT:Z:34\n"
        );
    }

    #[test]
    fn test_write_fastq_and_fasta_from_sam() {
        let sam = b"x\t4\t*\t0\t0\t*\t*\t0\t0\tA\t+\n";
        let mut out = Vec::new();
        write_fastq_from_sam(Cursor::new(sam), &mut out, &SamToFastxOptions::default()).unwrap();
        assert_eq!(out, MINIMAL_FASTQ);

        let mut out = Vec::new();
        write_fasta_from_sam(Cursor::new(sam), &mut out, &SamToFastxOptions::default()).unwrap();
        assert_eq!(out, b">x\nA\n");
    }
}
