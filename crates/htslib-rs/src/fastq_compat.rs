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
    /// Parse CASAVA 1.8 read number, filter flag, and barcode fields.
    pub casava: bool,
    /// Auxiliary tag used for CASAVA barcodes.
    pub barcode_tag: Option<String>,
    /// Use the second whitespace-delimited definition token as the read name when present.
    pub name2: bool,
    /// Extract an inline UMI barcode into this auxiliary tag.
    pub umi_tag: Option<String>,
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
        tokens
            .iter()
            .skip(1)
            .filter(|field| field.matches(':').count() >= 2)
            .map(|field| (*field).to_string())
            .collect()
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
