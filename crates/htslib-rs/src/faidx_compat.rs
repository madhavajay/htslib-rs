//! HTSlib-compatible FASTA/FASTQ index helpers backed by noodles.

use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};

use noodles::{
    core::{Position, Region},
    fasta, fastq,
};

/// A FASTA index.
pub type Index = fasta::fai::Index;

/// A FASTQ index.
pub type FastqIndex = fastq::fai::Index;

/// Builds a FASTA index from a reader.
pub fn build_index<R>(reader: R) -> io::Result<Index>
where
    R: BufRead,
{
    let mut reader = reader;
    let mut src = Vec::new();
    reader.read_to_end(&mut src)?;

    let mut lines = LinesWithOffsets::new(&src);
    let mut records = Vec::new();

    while let Some(record) = index_fasta_record(&mut lines)? {
        records.push(record);
    }

    Ok(Index::from(records))
}

/// Builds a FASTQ index from a reader.
pub fn build_fastq_index<R>(reader: R) -> io::Result<FastqIndex>
where
    R: BufRead,
{
    let mut reader = reader;
    let mut src = Vec::new();
    reader.read_to_end(&mut src)?;

    let mut lines = LinesWithOffsets::new(&src);
    let mut records = Vec::new();

    while let Some(record) = index_fastq_record(&mut lines)? {
        records.push(record);
    }

    Ok(records)
}

/// Reads a FASTA index from a reader.
pub fn read_index<R>(reader: R) -> io::Result<Index>
where
    R: BufRead,
{
    fasta::fai::io::Reader::new(reader).read_index()
}

/// Reads a FASTQ index from a reader.
pub fn read_fastq_index<R>(reader: R) -> io::Result<FastqIndex>
where
    R: BufRead,
{
    let mut reader = fastq::fai::io::Reader::new(reader);
    let mut records = Vec::new();
    let mut buf = String::new();

    while reader.read_record(&mut buf)? != 0 {
        let record = buf
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        records.push(record);
        buf.clear();
    }

    Ok(records)
}

/// Writes a FASTA index to a writer.
pub fn write_index<W>(writer: W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    let mut writer = fasta::fai::io::Writer::new(writer);
    writer.write_index(index)
}

/// Writes a FASTQ index to a writer.
pub fn write_fastq_index<W>(writer: W, index: &FastqIndex) -> io::Result<()>
where
    W: Write,
{
    let mut writer = fastq::fai::io::Writer::new(writer);

    for record in index {
        writer.write_record(record)?;
    }

    Ok(())
}

/// Returns whether the index has a sequence with `name`.
pub fn has_sequence(index: &Index, name: &str) -> bool {
    find_record(index, name).is_some()
}

/// Returns a sequence name by index.
pub fn sequence_name(index: &Index, i: usize) -> Option<&[u8]> {
    index.as_ref().get(i).map(|record| record.name().as_ref())
}

/// Returns the sequence length for `name`.
pub fn sequence_len(index: &Index, name: &str) -> Option<u64> {
    find_record(index, name).map(|record| record.length())
}

/// Returns the number of bases per line for `name`.
pub fn line_length(index: &Index, name: &str) -> Option<u64> {
    find_record(index, name).map(|record| record.line_bases())
}

/// Returns whether the FASTQ index has a sequence with `name`.
pub fn has_fastq_sequence(index: &FastqIndex, name: &str) -> bool {
    find_fastq_record(index, name).is_some()
}

/// Returns a FASTQ sequence name by index.
pub fn fastq_sequence_name(index: &FastqIndex, i: usize) -> Option<&str> {
    index.get(i).map(|record| record.name())
}

/// Returns the FASTQ sequence length for `name`.
pub fn fastq_sequence_len(index: &FastqIndex, name: &str) -> Option<u64> {
    find_fastq_record(index, name).map(|record| record.length())
}

/// Returns the number of bases per FASTQ sequence line for `name`.
pub fn fastq_line_length(index: &FastqIndex, name: &str) -> Option<u64> {
    find_fastq_record(index, name).map(|record| record.line_bases())
}

/// Fetches a FASTA sequence for a noodles region.
pub fn fetch_sequence<R>(reader: R, index: &Index, region: &Region) -> io::Result<Vec<u8>>
where
    R: BufRead + Seek,
{
    let mut reader = fasta::io::Reader::new(reader);
    let record = reader.query(index, region)?;
    Ok(record.sequence().as_ref().to_vec())
}

/// Fetches a FASTA sequence for an HTSlib-style region string.
pub fn fetch_region_sequence<R>(reader: &mut R, index: &Index, region: &str) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    let region = parse_region(region)?;
    let name = region_name(&region)?;
    let record = find_record(index, name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "reference sequence not found"))?;
    let (start, end) = resolve_interval(&region, record.length())?;

    read_indexed_range(
        reader,
        record.offset(),
        record.line_bases(),
        record.line_width(),
        start,
        end,
    )
}

/// Fetches a FASTQ sequence for an HTSlib-style region string.
pub fn fetch_fastq_region_sequence<R>(
    reader: &mut R,
    index: &FastqIndex,
    region: &str,
) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    let region = parse_region(region)?;
    let name = region_name(&region)?;
    let record = find_fastq_record(index, name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "reference sequence not found"))?;
    let (start, end) = resolve_interval(&region, record.length())?;

    read_indexed_range(
        reader,
        record.sequence_offset(),
        record.line_bases(),
        record.line_width(),
        start,
        end,
    )
}

/// Fetches FASTQ qualities for an HTSlib-style region string.
pub fn fetch_fastq_region_quality<R>(
    reader: &mut R,
    index: &FastqIndex,
    region: &str,
) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    let region = parse_region(region)?;
    let name = region_name(&region)?;
    let record = find_fastq_record(index, name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "reference sequence not found"))?;
    let (start, end) = resolve_interval(&region, record.length())?;

    read_indexed_range(
        reader,
        record.quality_scores_offset(),
        record.line_bases(),
        record.line_width(),
        start,
        end,
    )
}

/// Renders FASTA retrieval output using HTSlib test_faidx formatting.
pub fn write_fasta_retrieval<R, W>(
    reader: &mut R,
    index: &Index,
    regions: &[&str],
    writer: &mut W,
) -> io::Result<()>
where
    R: Read + Seek,
    W: Write,
{
    for region in regions {
        let sequence = fetch_region_sequence(reader, index, region)?;
        write_retrieval_record(writer, b'>', region, &sequence, None)?;
    }

    Ok(())
}

/// Renders FASTQ retrieval output using HTSlib test_faidx formatting.
pub fn write_fastq_retrieval<R, W>(
    reader: &mut R,
    index: &FastqIndex,
    regions: &[&str],
    writer: &mut W,
) -> io::Result<()>
where
    R: Read + Seek,
    W: Write,
{
    for region in regions {
        let sequence = fetch_fastq_region_sequence(reader, index, region)?;
        let quality = fetch_fastq_region_quality(reader, index, region)?;
        write_retrieval_record(writer, b'@', region, &sequence, Some(&quality))?;
    }

    Ok(())
}

/// Renders FASTQ sequence retrieval as FASTA output.
pub fn write_fastq_as_fasta_retrieval<R, W>(
    reader: &mut R,
    index: &FastqIndex,
    regions: &[&str],
    writer: &mut W,
) -> io::Result<()>
where
    R: Read + Seek,
    W: Write,
{
    for region in regions {
        let sequence = fetch_fastq_region_sequence(reader, index, region)?;
        write_retrieval_record(writer, b'>', region, &sequence, None)?;
    }

    Ok(())
}

fn find_record<'a>(index: &'a Index, name: &str) -> Option<&'a fasta::fai::Record> {
    index
        .as_ref()
        .iter()
        .find(|record| record.name() == name.as_bytes())
}

fn index_fasta_record(lines: &mut LinesWithOffsets<'_>) -> io::Result<Option<fasta::fai::Record>> {
    let Some((_, definition)) = lines.next() else {
        return Ok(None);
    };

    if !definition.starts_with(b">") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTA definition",
        ));
    }

    let name = parse_fasta_name(definition);

    let Some((sequence_offset, first_line)) = lines.peek() else {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "missing FASTA sequence",
        ));
    };

    if first_line.starts_with(b">") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty FASTA sequence",
        ));
    }

    let line_bases = trimmed_len(first_line) as u64;
    let line_width = first_line.len() as u64;
    let mut length = 0;

    while let Some((_, line)) = lines.peek() {
        if line.starts_with(b">") {
            break;
        }

        length += trimmed_len(line) as u64;
        lines.next();
    }

    Ok(Some(fasta::fai::Record::new(
        name,
        length,
        sequence_offset,
        line_bases,
        line_width,
    )))
}

fn parse_fasta_name(definition: &[u8]) -> Vec<u8> {
    definition[1..]
        .trim_ascii_start()
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default()
        .to_vec()
}

fn find_fastq_record<'a>(index: &'a FastqIndex, name: &str) -> Option<&'a fastq::fai::Record> {
    index.iter().find(|record| record.name() == name)
}

fn index_fastq_record(lines: &mut LinesWithOffsets<'_>) -> io::Result<Option<fastq::fai::Record>> {
    let Some((_, definition)) = lines.next() else {
        return Ok(None);
    };

    if !definition.starts_with(b"@") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTQ definition",
        ));
    }

    let name = parse_fastq_name(definition)?;
    let (sequence_offset, line_bases, line_width, length) = index_fastq_sequence(lines)?;
    let quality_scores_offset = index_fastq_quality(lines, length)?;

    Ok(Some(fastq::fai::Record::new(
        name,
        length,
        sequence_offset,
        line_bases,
        line_width,
        quality_scores_offset,
    )))
}

fn parse_fastq_name(definition: &[u8]) -> io::Result<String> {
    let name = definition[1..]
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default();

    std::str::from_utf8(name)
        .map(String::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn index_fastq_sequence(lines: &mut LinesWithOffsets<'_>) -> io::Result<(u64, u64, u64, u64)> {
    let Some((sequence_offset, first_line)) = lines.next() else {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "missing FASTQ sequence",
        ));
    };

    if first_line.starts_with(b"+") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty FASTQ sequence",
        ));
    }

    let line_bases = trimmed_len(first_line) as u64;
    let line_width = first_line.len() as u64;
    let mut length = line_bases;

    loop {
        let Some((_, line)) = lines.peek() else {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "missing FASTQ plus line",
            ));
        };

        if line.starts_with(b"+") {
            lines.next();
            break;
        }

        length += trimmed_len(line) as u64;
        lines.next();
    }

    Ok((sequence_offset, line_bases, line_width, length))
}

fn index_fastq_quality(lines: &mut LinesWithOffsets<'_>, expected_len: u64) -> io::Result<u64> {
    let Some((quality_scores_offset, first_line)) = lines.next() else {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "missing FASTQ quality scores",
        ));
    };

    let mut length = trimmed_len(first_line) as u64;

    while length < expected_len {
        let Some((_, line)) = lines.next() else {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "short FASTQ quality scores",
            ));
        };

        length += trimmed_len(line) as u64;
    }

    if length == expected_len {
        Ok(quality_scores_offset)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ quality scores longer than sequence",
        ))
    }
}

#[derive(Debug)]
struct LinesWithOffsets<'a> {
    src: &'a [u8],
    pos: usize,
    peeked: Option<(u64, &'a [u8])>,
}

impl<'a> LinesWithOffsets<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            pos: 0,
            peeked: None,
        }
    }

    fn peek(&mut self) -> Option<(u64, &'a [u8])> {
        if self.peeked.is_none() {
            self.peeked = self.read_next();
        }

        self.peeked
    }

    fn next(&mut self) -> Option<(u64, &'a [u8])> {
        self.peeked.take().or_else(|| self.read_next())
    }

    fn read_next(&mut self) -> Option<(u64, &'a [u8])> {
        if self.pos >= self.src.len() {
            return None;
        }

        let start = self.pos;
        let rel_end = self.src[start..]
            .iter()
            .position(|b| *b == b'\n')
            .map(|i| start + i + 1)
            .unwrap_or(self.src.len());
        self.pos = rel_end;

        Some((start as u64, &self.src[start..rel_end]))
    }
}

fn trimmed_len(line: &[u8]) -> usize {
    match line.iter().rposition(|b| !b.is_ascii_whitespace()) {
        Some(i) => i + 1,
        None => 0,
    }
}

fn parse_region(region: &str) -> io::Result<Region> {
    region
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn region_name(region: &Region) -> io::Result<&str> {
    std::str::from_utf8(region.name()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn resolve_interval(region: &Region, len: u64) -> io::Result<(u64, u64)> {
    let interval = region.interval();
    let start = interval
        .start()
        .map(position_to_u64)
        .transpose()?
        .unwrap_or(1);
    let end = interval
        .end()
        .map(position_to_u64)
        .transpose()?
        .unwrap_or(len)
        .min(len);

    if start == 0 || start > end || start > len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid faidx query interval",
        ));
    }

    Ok((start, end))
}

fn position_to_u64(position: Position) -> io::Result<u64> {
    u64::try_from(usize::from(position)).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn read_indexed_range<R>(
    reader: &mut R,
    offset: u64,
    line_bases: u64,
    line_width: u64,
    start: u64,
    end: u64,
) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    let start = start - 1;
    let pos = offset + start / line_bases * line_width + start % line_bases;
    let len =
        usize::try_from(end - start).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    reader.seek(SeekFrom::Start(pos))?;

    let mut sequence = Vec::with_capacity(len);
    let mut buf = [0; 8192];

    while sequence.len() < len {
        let n = reader.read(&mut buf)?;

        if n == 0 {
            break;
        }

        for b in &buf[..n] {
            if !b.is_ascii_whitespace() {
                sequence.push(*b);

                if sequence.len() == len {
                    break;
                }
            }
        }
    }

    if sequence.len() == len {
        Ok(sequence)
    } else {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "short faidx sequence read",
        ))
    }
}

fn write_retrieval_record<W>(
    writer: &mut W,
    marker: u8,
    region: &str,
    sequence: &[u8],
    quality: Option<&[u8]>,
) -> io::Result<()>
where
    W: Write,
{
    writeln!(
        writer,
        "{}{} length: {}",
        char::from(marker),
        region,
        sequence.len()
    )?;
    write_wrapped(writer, sequence)?;

    if let Some(quality) = quality {
        if quality.len() != sequence.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sequence and quality lengths differ",
            ));
        }

        writeln!(writer, "+")?;
        write_wrapped(writer, quality)?;
    }

    Ok(())
}

fn write_wrapped<W>(writer: &mut W, data: &[u8]) -> io::Result<()>
where
    W: Write,
{
    for chunk in data.chunks(50) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use noodles::core::Region;

    use super::{
        build_fastq_index, build_index, fastq_line_length, fastq_sequence_len, fastq_sequence_name,
        fetch_fastq_region_quality, fetch_fastq_region_sequence, fetch_region_sequence,
        fetch_sequence, has_fastq_sequence, has_sequence, line_length, read_fastq_index,
        read_index, sequence_len, sequence_name, write_fasta_retrieval,
        write_fastq_as_fasta_retrieval, write_fastq_index, write_fastq_retrieval, write_index,
    };

    const FASTA: &[u8] = include_bytes!("../../../repos/htslib/test/xx.fa");
    const FAI: &[u8] = include_bytes!("../../../repos/htslib/test/xx.fa.fai");
    const FASTQ: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/fastqs.fq");
    const FASTQ_FAI: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/fastqs.fq.expected.fai");

    #[test]
    fn test_read_and_build_index() {
        let expected = read_index(Cursor::new(FAI)).unwrap();
        let actual = build_index(Cursor::new(FASTA)).unwrap();

        assert_eq!(actual, expected);

        let mut buf = Vec::new();
        write_index(&mut buf, &actual).unwrap();
        assert_eq!(buf, FAI);
    }

    #[test]
    fn test_read_and_build_fastq_index() {
        let expected = read_fastq_index(Cursor::new(FASTQ_FAI)).unwrap();
        let actual = build_fastq_index(Cursor::new(FASTQ)).unwrap();

        assert_eq!(actual, expected);

        let mut buf = Vec::new();
        write_fastq_index(&mut buf, &actual).unwrap();
        assert_eq!(buf, FASTQ_FAI);
    }

    #[test]
    fn test_index_lookup_helpers() {
        let index = read_index(Cursor::new(FAI)).unwrap();

        assert!(has_sequence(&index, "xx"));
        assert!(!has_sequence(&index, "missing"));
        assert_eq!(sequence_name(&index, 0), Some(&b"xx"[..]));
        assert_eq!(sequence_name(&index, 99), None);
        assert_eq!(sequence_len(&index, "zz"), Some(30));
        assert_eq!(line_length(&index, "xx"), Some(20));
    }

    #[test]
    fn test_fastq_index_lookup_helpers() {
        let index = read_fastq_index(Cursor::new(FASTQ_FAI)).unwrap();

        assert!(has_fastq_sequence(&index, "FAKE0005_1"));
        assert!(!has_fastq_sequence(&index, "missing"));
        assert_eq!(fastq_sequence_name(&index, 0), Some("FAKE0005_1"));
        assert_eq!(fastq_sequence_name(&index, 999), None);
        assert_eq!(fastq_sequence_len(&index, "FSRRS4401CM938_1"), Some(453));
        assert_eq!(fastq_line_length(&index, "FAKE0005_3"), Some(63));
    }

    #[test]
    fn test_fetch_sequence() {
        let index = read_index(Cursor::new(FAI)).unwrap();

        let region: Region = "xx:1-10".parse().unwrap();
        assert_eq!(
            fetch_sequence(Cursor::new(FASTA), &index, &region).unwrap(),
            b"AAAAAAAAAA"
        );

        let region: Region = "zz:21-30".parse().unwrap();
        assert_eq!(
            fetch_sequence(Cursor::new(FASTA), &index, &region).unwrap(),
            b"CCCCCCCCCC"
        );
    }

    #[test]
    fn test_fetch_region_sequence() {
        let index = read_index(Cursor::new(FAI)).unwrap();
        let mut reader = Cursor::new(FASTA);

        assert_eq!(
            fetch_region_sequence(&mut reader, &index, "xx:1-10").unwrap(),
            b"AAAAAAAAAA"
        );
        assert_eq!(
            fetch_region_sequence(&mut reader, &index, "zz:21-30").unwrap(),
            b"CCCCCCCCCC"
        );
    }

    #[test]
    fn test_fetch_fastq_region_sequence_and_quality() {
        let index = read_fastq_index(Cursor::new(FASTQ_FAI)).unwrap();
        let mut reader = Cursor::new(FASTQ);

        assert_eq!(
            fetch_fastq_region_sequence(&mut reader, &index, "FAKE0006_1:4-12").unwrap(),
            b"TGCATGCAT"
        );
        assert_eq!(
            fetch_fastq_region_quality(&mut reader, &index, "FAKE0006_1:4-12").unwrap(),
            b"{zyxwvuts"
        );
    }

    #[test]
    fn test_write_retrieval_output() {
        let fasta_index = read_index(Cursor::new(FAI)).unwrap();
        let mut fasta_reader = Cursor::new(FASTA);
        let mut out = Vec::new();
        write_fasta_retrieval(&mut fasta_reader, &fasta_index, &["xx:1-10"], &mut out).unwrap();
        assert_eq!(out, b">xx:1-10 length: 10\nAAAAAAAAAA\n");

        let fastq_index = read_fastq_index(Cursor::new(FASTQ_FAI)).unwrap();
        let mut fastq_reader = Cursor::new(FASTQ);
        let mut out = Vec::new();
        write_fastq_retrieval(
            &mut fastq_reader,
            &fastq_index,
            &["FAKE0006_1:4-12"],
            &mut out,
        )
        .unwrap();
        assert_eq!(
            out,
            b"@FAKE0006_1:4-12 length: 9\nTGCATGCAT\n+\n{zyxwvuts\n"
        );

        let mut fastq_reader = Cursor::new(FASTQ);
        let mut out = Vec::new();
        write_fastq_as_fasta_retrieval(
            &mut fastq_reader,
            &fastq_index,
            &["FAKE0006_1:4-12"],
            &mut out,
        )
        .unwrap();
        assert_eq!(out, b">FAKE0006_1:4-12 length: 9\nTGCATGCAT\n");
    }
}
