//! HTSlib-compatible tabix helpers backed by noodles tabix and CSI readers.

use std::{
    fs::File,
    io::{self, BufRead, Write},
    path::Path,
};

use crate::{
    bgzf,
    core::{Position, Region},
    csi,
    index_compat::{associated_data_path, read_associated_vcf_index},
    tabix,
};

/// A tabix-indexed text format preset.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextFormat {
    /// BED, using 0-based start and 1-based end coordinates.
    Bed,
    /// GFF/GTF, using 1-based closed coordinates.
    Gff,
    /// VCF, using 1-based variant positions.
    Vcf,
}

/// Builds a BGZF-compressed BED stream and a matching TBI index.
pub fn build_bed_bgzf_and_index<R, W>(reader: R, writer: W) -> io::Result<(W, tabix::Index)>
where
    R: BufRead,
    W: Write,
{
    build_bgzf_and_index(reader, writer, TextFormat::Bed)
}

/// Builds a BGZF-compressed text stream and a matching TBI index.
pub fn build_bgzf_and_index<R, W>(
    mut reader: R,
    writer: W,
    format: TextFormat,
) -> io::Result<(W, tabix::Index)>
where
    R: BufRead,
    W: Write,
{
    let mut writer = bgzf::io::Writer::new(writer);
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(format.header_builder().build());

    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        let start_position = writer.virtual_position();
        writer.write_all(line.as_bytes())?;
        let end_position = writer.virtual_position();

        if let Some((reference_sequence_name, start, end)) = parse_index_fields(format, &line)? {
            let chunk = csi::binning_index::index::reference_sequence::bin::Chunk::new(
                start_position,
                end_position,
            );
            indexer.add_record(reference_sequence_name, start, end, chunk)?;
        }

        line.clear();
    }

    let writer = writer.finish()?;
    let index = indexer.build();

    Ok((writer, index))
}

/// Builds a BGZF-compressed text stream and a matching CSI index.
pub fn build_bgzf_and_csi<R, W>(
    reader: R,
    writer: W,
    format: TextFormat,
) -> io::Result<(W, csi::Index)>
where
    R: BufRead,
    W: Write,
{
    build_bgzf_and_csi_with_min_shift(reader, writer, format, 14)
}

/// Builds a BGZF-compressed text stream and a matching CSI index with a custom min_shift.
pub fn build_bgzf_and_csi_with_min_shift<R, W>(
    mut reader: R,
    writer: W,
    format: TextFormat,
    min_shift: u8,
) -> io::Result<(W, csi::Index)>
where
    R: BufRead,
    W: Write,
{
    use csi::binning_index::{
        Indexer,
        index::{
            header::ReferenceSequenceNames, reference_sequence::bin::Chunk,
            reference_sequence::index::BinnedIndex,
        },
    };

    let mut writer = bgzf::io::Writer::new(writer);
    let mut reference_sequence_names = ReferenceSequenceNames::new();
    let mut records = Vec::new();
    let mut max_position = 0;
    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        let start_position = writer.virtual_position();
        writer.write_all(line.as_bytes())?;
        let end_position = writer.virtual_position();

        if let Some((reference_sequence_name, start, end)) = parse_index_fields(format, &line)? {
            let (reference_sequence_id, _) =
                reference_sequence_names.insert_full(reference_sequence_name.into());
            let chunk = Chunk::new(start_position, end_position);
            max_position = max_position.max(usize::from(end));
            records.push((reference_sequence_id, start, end, chunk));
        }

        line.clear();
    }

    let writer = writer.finish()?;
    let mut indexer =
        Indexer::<BinnedIndex>::new(min_shift, csi_depth_for_position(min_shift, max_position));

    for (reference_sequence_id, start, end, chunk) in records {
        indexer.add_record(Some((reference_sequence_id, start, end, true)), chunk)?;
    }

    let header = format
        .header_builder()
        .set_reference_sequence_names(reference_sequence_names.clone())
        .build();
    let index = indexer
        .set_header(header)
        .build(reference_sequence_names.len());

    Ok((writer, index))
}

fn csi_depth_for_position(min_shift: u8, max_position: usize) -> u8 {
    const DEFAULT_DEPTH: u8 = 6;

    let mut depth = DEFAULT_DEPTH;

    while (max_position as u128) > csi_max_position(min_shift, depth) {
        depth += 1;
    }

    depth
}

fn csi_max_position(min_shift: u8, depth: u8) -> u128 {
    let bit_count = u32::from(min_shift) + 3 * u32::from(depth);
    1_u128 << bit_count
}

/// Writes a BGZF-compressed BED stream and a matching TBI index to local files.
pub fn write_bed_bgzf_and_index<R, P, Q>(reader: R, bgzf_dst: P, tbi_dst: Q) -> io::Result<()>
where
    R: BufRead,
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let bgzf_file = File::create(bgzf_dst)?;
    let (_, index) = build_bed_bgzf_and_index(reader, bgzf_file)?;
    tabix::fs::write(tbi_dst, &index)
}

/// Writes a BGZF-compressed text stream and a matching TBI index to local files.
pub fn write_bgzf_and_index<R, P, Q>(
    reader: R,
    bgzf_dst: P,
    tbi_dst: Q,
    format: TextFormat,
) -> io::Result<()>
where
    R: BufRead,
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let bgzf_file = File::create(bgzf_dst)?;
    let (_, index) = build_bgzf_and_index(reader, bgzf_file, format)?;
    tabix::fs::write(tbi_dst, &index)
}

/// Writes a BGZF-compressed text stream and a matching CSI index to local files.
pub fn write_bgzf_and_csi<R, P, Q>(
    reader: R,
    bgzf_dst: P,
    csi_dst: Q,
    format: TextFormat,
) -> io::Result<()>
where
    R: BufRead,
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let bgzf_file = File::create(bgzf_dst)?;
    let (_, index) = build_bgzf_and_csi(reader, bgzf_file, format)?;
    csi::fs::write(csi_dst, &index)
}

/// Queries tabix-indexed text records from a local BGZF file.
pub fn query_records_from_path<P>(
    src: P,
    index: tabix::Index,
    region: &Region,
) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    let file = File::open(src)?;
    let mut reader = csi::io::IndexedReader::new(file, index);
    let query = reader.query(region)?;

    query
        .map(|result| result.map(|record| record.as_ref().to_string()))
        .collect()
}

/// Queries multiple tabix regions and formats output like `tabix --separate-regions`.
pub fn query_records_from_path_separate_regions<P>(
    src: P,
    index: tabix::Index,
    regions: &[(&str, Region)],
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut out = String::new();

    for (label, region) in regions {
        out.push('#');
        out.push_str(label);
        out.push('\n');

        for record in query_records_from_path(&src, index.clone(), region)? {
            out.push_str(&record);
            out.push('\n');
        }
    }

    Ok(out)
}

/// Queries CSI-indexed text records from a local BGZF file.
pub fn query_csi_records_from_path<P>(
    src: P,
    index: csi::Index,
    region: &Region,
) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    let file = File::open(src)?;
    let mut reader = csi::io::IndexedReader::new(file, index);
    let query = reader.query(region)?;

    query
        .map(|result| result.map(|record| record.as_ref().to_string()))
        .collect()
}

/// Queries CSI-indexed VCF text records, including HTSlib-style `INFO/END` spans.
pub fn query_vcf_csi_records_from_path<P>(
    src: P,
    index: csi::Index,
    region: &Region,
) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    drop(index);

    query_vcf_text_records_from_path(src, region)
}

fn query_vcf_text_records_from_path<P>(src: P, region: &Region) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    let file = File::open(src)?;
    let mut reader = io::BufReader::new(bgzf::io::Reader::new(file));
    let mut records = Vec::new();
    let mut line = String::new();

    while reader.read_line(&mut line)? != 0 {
        if vcf_line_intersects_region(&line, region)? {
            records.push(line.trim_end_matches(['\r', '\n']).to_string());
        }

        line.clear();
    }

    Ok(records)
}

/// Queries VCF text records from a local BGZF file using its associated TBI or CSI index.
pub fn query_vcf_records_from_associated_path<P>(src: P, region: &Region) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    let index = read_associated_vcf_index(&src)?;
    drop(index);

    query_vcf_text_records_from_path(associated_data_path(src), region)
}

impl TextFormat {
    fn header_builder(self) -> csi::binning_index::index::header::Builder {
        match self {
            Self::Bed => csi::binning_index::index::header::Builder::bed(),
            Self::Gff => csi::binning_index::index::header::Builder::gff(),
            Self::Vcf => csi::binning_index::index::header::Builder::vcf(),
        }
    }
}

fn parse_index_fields(
    format: TextFormat,
    line: &str,
) -> io::Result<Option<(&str, Position, Position)>> {
    const COMMENT_PREFIX: char = '#';

    let line = line.trim_end_matches(['\r', '\n']);

    if line.is_empty() || line.starts_with(COMMENT_PREFIX) {
        return Ok(None);
    }

    match format {
        TextFormat::Bed => parse_bed_index_fields(line),
        TextFormat::Gff => parse_gff_index_fields(line),
        TextFormat::Vcf => parse_vcf_index_fields(line),
    }
}

fn parse_bed_index_fields(line: &str) -> io::Result<Option<(&str, Position, Position)>> {
    let reference_sequence_name = get_field(line, 0, "reference name")?;

    let start = parse_position(get_field(line, 1, "BED start")?, 1)?;
    let end = parse_position(get_field(line, 2, "BED end")?, 0)?;

    Ok(Some((reference_sequence_name, start, end)))
}

fn parse_gff_index_fields(line: &str) -> io::Result<Option<(&str, Position, Position)>> {
    let reference_sequence_name = get_field(line, 0, "reference name")?;
    let start = parse_position(get_field(line, 3, "GFF start")?, 0)?;
    let end = parse_position(get_field(line, 4, "GFF end")?, 0)?;

    Ok(Some((reference_sequence_name, start, end)))
}

fn parse_vcf_index_fields(line: &str) -> io::Result<Option<(&str, Position, Position)>> {
    let reference_sequence_name = get_field(line, 0, "reference name")?;
    let start = parse_position(get_field(line, 1, "VCF position")?, 0)?;
    let end = parse_vcf_info_end(get_field(line, 7, "VCF INFO")?)?.unwrap_or(start);
    let end = end.max(start);

    Ok(Some((reference_sequence_name, start, end)))
}

fn parse_vcf_info_end(info: &str) -> io::Result<Option<Position>> {
    for field in info.split(';') {
        let Some(value) = field.strip_prefix("END=") else {
            continue;
        };

        return parse_position(value, 0).map(Some);
    }

    Ok(None)
}

fn vcf_line_intersects_region(line: &str, region: &Region) -> io::Result<bool> {
    let Some((reference_sequence_name, start, end)) = parse_index_fields(TextFormat::Vcf, line)?
    else {
        return Ok(false);
    };

    if reference_sequence_name.as_bytes() != region.name() {
        return Ok(false);
    }

    let record_start = usize::from(start);
    let record_end = usize::from(end);
    let interval = region.interval();
    let region_start = interval.start().map(usize::from).unwrap_or(1);
    let region_end = interval.end().map(usize::from).unwrap_or(usize::MAX);

    Ok(record_start <= region_end && region_start <= record_end)
}

fn get_field<'a>(line: &'a str, i: usize, name: &str) -> io::Result<&'a str> {
    line.split('\t').nth(i).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("missing {name} field at index {i}"),
        )
    })
}

fn parse_position(s: &str, addend: usize) -> io::Result<Position> {
    s.parse::<usize>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|n| {
            Position::try_from(n + addend)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{
        TextFormat, build_bed_bgzf_and_index, build_bgzf_and_csi, build_bgzf_and_index,
        query_csi_records_from_path, query_records_from_path,
    };

    #[test]
    fn test_build_and_query_bed_index() {
        let data = b"sq0\t7\t13\na\nsq0\t20\t34\n";
        assert!(build_bed_bgzf_and_index(Cursor::new(data), Vec::new()).is_err());

        let data = b"sq0\t7\t13\nsq0\t20\t34\n";
        let (encoded, index) = build_bed_bgzf_and_index(Cursor::new(data), Vec::new()).unwrap();

        let path = std::env::temp_dir().join(format!(
            "htslib-rs-tabix-{}-{}.bed.gz",
            std::process::id(),
            "unit"
        ));

        std::fs::write(&path, encoded).unwrap();

        let region = "sq0:8-13".parse().unwrap();
        let records = query_records_from_path(&path, index, &region).unwrap();

        std::fs::remove_file(path).unwrap();

        assert_eq!(records, ["sq0\t7\t13"]);
    }

    #[test]
    fn test_build_and_query_vcf_index() {
        let data = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nsq0\t8\t.\tA\tC\t.\t.\t.\nsq0\t21\t.\tA\tG\t.\t.\t.\n";
        let (encoded, index) =
            build_bgzf_and_index(Cursor::new(data), Vec::new(), TextFormat::Vcf).unwrap();

        let path = std::env::temp_dir().join(format!(
            "htslib-rs-tabix-{}-{}.vcf.gz",
            std::process::id(),
            "vcf-unit"
        ));

        std::fs::write(&path, encoded).unwrap();

        let region = "sq0:8-8".parse().unwrap();
        let records = query_records_from_path(&path, index, &region).unwrap();

        std::fs::remove_file(path).unwrap();

        assert_eq!(records, ["sq0\t8\t.\tA\tC\t.\t.\t."]);
    }

    #[test]
    fn test_build_and_query_csi_index() {
        let data = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nsq0\t8\t.\tA\tC\t.\t.\t.\nsq0\t21\t.\tA\tG\t.\t.\t.\n";
        let (encoded, index) =
            build_bgzf_and_csi(Cursor::new(data), Vec::new(), TextFormat::Vcf).unwrap();

        let path = std::env::temp_dir().join(format!(
            "htslib-rs-tabix-{}-{}.vcf.gz",
            std::process::id(),
            "csi-unit"
        ));

        std::fs::write(&path, encoded).unwrap();

        let region = "sq0:21-21".parse().unwrap();
        let records = query_csi_records_from_path(&path, index, &region).unwrap();

        std::fs::remove_file(path).unwrap();

        assert_eq!(records, ["sq0\t21\t.\tA\tG\t.\t.\t."]);
    }
}
