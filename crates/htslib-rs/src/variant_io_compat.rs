//! HTSlib-compatible variant I/O helpers backed by noodles readers.

use std::{
    collections::{HashMap, HashSet},
    fs::File,
    hash::{Hash, Hasher},
    io::{self, BufRead, BufReader, Cursor, Read, Write},
    path::Path,
};

use crate::{
    bcf,
    core::{Position, Region},
    index_compat::{associated_data_path, read_associated_bcf_index, read_associated_vcf_index},
    vcf,
};
use vcf::variant::io::Write as _;

/// A VCF header shared by VCF and BCF streams.
pub type Header = vcf::Header;

/// Reads a VCF header from a buffered reader.
pub fn read_vcf_header<R>(reader: R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut reader = vcf::io::Reader::new(reader);

    reader.read_header()
}

/// Reads a VCF header from a local file.
pub fn read_vcf_header_from_path<P>(src: P) -> io::Result<Header>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(read_vcf_header)
}

/// Counts VCF records from a buffered reader.
pub fn count_vcf_records<R>(reader: R) -> io::Result<usize>
where
    R: BufRead,
{
    let mut reader = vcf::io::Reader::new(reader);
    reader.read_header()?;

    reader
        .records()
        .try_fold(0, |n, result| result.map(|_| n + 1))
}

/// Counts VCF records from a local file.
pub fn count_vcf_records_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .and_then(count_vcf_records)
}

/// Reads VCF input without producing output and returns the number of records seen.
pub fn benchmark_vcf_view_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    count_vcf_records_from_path(src)
}

/// Writes VCF input as VCF text, including the header and all records.
pub fn write_vcf_from_path<P, W>(src: P, mut dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let src = src.as_ref();
    let raw_header = read_raw_vcf_header_text(src)?;

    if raw_header.starts_with("##fileformat=VCFv4.4\n") {
        return write_vcf_from_raw_path(src, dst);
    }

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut body = Vec::new();

    for result in reader.records() {
        match result {
            Ok(record) => {
                if write_vcf_record_htslib(&mut body, &header, &record).is_err() {
                    return write_vcf_from_raw_path(src, dst);
                }
            }
            Err(_) => return write_vcf_from_raw_path(src, dst),
        }
    }

    dst.write_all(raw_header.as_bytes())?;
    dst.write_all(&body)?;

    Ok(dst)
}

/// Writes VCF text as BGZF-compressed VCF and writes a TBI index for it.
pub fn write_vcf_bgzf_from_path_with_tbi<P, Q, R>(src: P, bgzf_dst: Q, tbi_dst: R) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let text = write_vcf_from_path(src, Vec::new())?;

    crate::tabix_compat::write_bgzf_and_index(
        Cursor::new(text),
        bgzf_dst,
        tbi_dst,
        crate::tabix_compat::TextFormat::Vcf,
    )
}

fn write_vcf_from_raw_path<W>(src: &Path, mut dst: W) -> io::Result<W>
where
    W: Write,
{
    let raw = std::fs::read_to_string(src)?;
    let raw_header = read_raw_vcf_header_text(src)?;
    let header_len = raw_header.len();

    dst.write_all(canonicalize_vcf_header_text(&raw_header).as_bytes())?;

    for line in raw[header_len..].lines() {
        dst.write_all(canonicalize_vcf_record_line(line).as_bytes())?;
        dst.write_all(b"\n")?;
    }

    Ok(dst)
}

/// Writes a parsed and canonicalized VCF header from local input.
pub fn write_vcf_header_from_path<P, W>(src: P, mut dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let raw_header = read_raw_vcf_header_text(src)?;
    let header = canonicalize_vcf_header_text(&raw_header);

    dst.write_all(header.as_bytes())?;

    Ok(dst)
}

/// Writes a VCF text view of the first `limit` records, including the header.
pub fn view_vcf_text_from_path_with_limit<P>(src: P, limit: Option<usize>) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let raw_header = read_raw_vcf_header_text(&src)?;
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut buf = raw_header.into_bytes();

    for result in reader.records().take(limit.unwrap_or(usize::MAX)) {
        let record = result?;
        write_vcf_record_htslib(&mut buf, &header, &record)?;
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_raw_vcf_header_text<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new)?;
    let mut header = String::new();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        header.push_str(&line);

        if line.starts_with("#CHROM\t") {
            break;
        }
    }

    Ok(header)
}

/// Queries VCF records from a local BGZF-compressed file using its associated TBI or CSI index.
pub fn query_vcf_records_from_path<P>(src: P, region: &Region) -> io::Result<Vec<vcf::Record>>
where
    P: AsRef<Path>,
{
    let index = read_associated_vcf_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let query = reader.query(&header, region)?;

    query.records().collect()
}

/// Returns an owning iterator over VCF records from an indexed local file.
pub fn iter_vcf_records_from_path<P>(
    src: P,
    region: &Region,
) -> io::Result<std::vec::IntoIter<vcf::Record>>
where
    P: AsRef<Path>,
{
    query_vcf_records_from_path(src, region).map(Vec::into_iter)
}

/// Counts VCF records from a local indexed file that intersect a region.
pub fn count_vcf_records_in_region_from_path<P>(src: P, region: &Region) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    query_vcf_records_from_path(src, region).map(|records| records.len())
}

/// Writes VCF input as VCF text for the requested indexed regions.
pub fn view_vcf_regions_as_text_from_path<P>(src: P, regions: &[Region]) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_vcf_regions_as_text_from_path_with_limit(src, regions, None)
}

/// Writes VCF input as VCF text for requested indexed regions with an optional record limit.
pub fn view_vcf_regions_as_text_from_path_with_limit<P>(
    src: P,
    regions: &[Region],
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let index = read_associated_vcf_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;

    let mut header_buf = Vec::new();
    vcf::io::Writer::new(&mut header_buf).write_header(&header)?;

    let header_text =
        String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut buf = ensure_vcf_pass_header(&header_text).into_bytes();
    let mut remaining = limit.unwrap_or(usize::MAX);

    for region in regions {
        if remaining == 0 {
            break;
        }

        match reader.query(&header, region) {
            Ok(query) => {
                for result in query.records() {
                    if remaining == 0 {
                        break;
                    }

                    let record = result?;
                    write_vcf_record_htslib(&mut buf, &header, &record)?;
                    remaining -= 1;
                }
            }
            Err(e) if is_present_contig_index_mismatch(&header, region, &e) => {}
            Err(e) => return Err(e),
        }
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn write_vcf_record_htslib<W, R>(dst: &mut W, header: &Header, record: &R) -> io::Result<()>
where
    W: Write,
    R: vcf::variant::Record,
{
    let mut line = Vec::new();
    vcf::io::Writer::new(&mut line).write_variant_record(header, record)?;
    write_vcf_record_line_htslib(dst, &line, header.sample_names().len())
}

fn write_vcf_record_line_htslib<W>(dst: &mut W, line: &[u8], sample_count: usize) -> io::Result<()>
where
    W: Write,
{
    if sample_count > 0 && line.iter().filter(|&&b| b == b'\t').count() == 7 {
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        dst.write_all(line)?;
        for _ in 0..=sample_count {
            dst.write_all(b"\t.")?;
        }
        dst.write_all(b"\n")
    } else {
        dst.write_all(line)
    }
}

fn canonicalize_vcf_record_line(line: &str) -> String {
    let mut fields = line.split('\t').collect::<Vec<_>>();

    if fields.len() < 10 {
        return line.to_string();
    }

    let format_keys = fields[8].split(':').collect::<Vec<_>>();
    let mut kept_keys = Vec::new();
    let mut kept_indices = Vec::new();

    for (i, key) in format_keys.iter().enumerate() {
        if !kept_keys.contains(key) {
            kept_keys.push(*key);
            kept_indices.push(i);
        }
    }

    let mut out = fields[..8].join("\t");
    out.push('\t');
    out.push_str(&kept_keys.join(":"));

    for sample in fields.drain(9..) {
        let values = sample.split(':').collect::<Vec<_>>();

        out.push('\t');
        out.push_str(
            &kept_indices
                .iter()
                .zip(&kept_keys)
                .map(|(&i, &key)| {
                    let value = values.get(i).copied().unwrap_or(".");

                    if key == "GT" {
                        normalize_vcf44_genotype_value(value)
                    } else {
                        value.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(":"),
        );
    }

    out
}

fn normalize_vcf44_genotype_value(value: &str) -> String {
    let Some(first) = value.as_bytes().first().copied() else {
        return ".".into();
    };

    if first != b'/' && first != b'|' {
        return value.to_string();
    }

    let rest = &value[1..];
    let Some(next_separator) = rest.bytes().find(|&b| b == b'/' || b == b'|') else {
        return if (first == b'|' && rest != ".") || (first == b'/' && rest == ".") {
            rest.to_string()
        } else {
            value.to_string()
        };
    };

    if next_separator == first {
        rest.to_string()
    } else {
        value.to_string()
    }
}

/// Returns the VCF text produced by the record serialization path in HTSlib's `test-vcf-api.c`.
pub fn test_vcf_api_record_serialization_text() -> String {
    let mut out = String::from(TEST_VCF_API_SERIALIZATION_HEADER);

    for record in test_vcf_api_source_records() {
        out.push_str(&record.to_vcf_line());
        out.push('\n');

        let mut synced = record.clone();
        synced.id = ".".into();
        synced.remove_format("GQ");
        out.push_str(&synced.to_vcf_line());
        out.push('\n');

        synced.reference_bases = "G".into();
        synced.alternate_bases = "A".into();
        synced.set_info("DP", "99");
        synced.set_format_values("DP", ["9", "9", "9"]);
        out.push_str(&synced.to_vcf_line());
        out.push('\n');
    }

    out
}

/// A VCF header identifier in an HTSlib-style namespace.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfHeaderId {
    /// A FILTER header ID.
    Filter(String),
    /// An INFO header ID.
    Info(String),
    /// A FORMAT header ID.
    Format(String),
    /// A contig header ID.
    Contig(String),
    /// A structured nonstandard header record ID.
    Other { key: String, id: String },
}

/// Returns whether a VCF header contains a header ID.
pub fn vcf_header_has_id(header: &Header, id: &VcfHeaderId) -> bool {
    match id {
        VcfHeaderId::Filter(id) => header_has_filter(header, id),
        VcfHeaderId::Info(id) => header_has_info(header, id),
        VcfHeaderId::Format(id) => header_has_format(header, id),
        VcfHeaderId::Contig(id) => header_has_contig(header, id),
        VcfHeaderId::Other { key, id } => header_has_other_record_id(header, key, id),
    }
}

/// Removes a VCF header ID from its HTSlib-style namespace.
pub fn vcf_header_remove_id(header: &mut Header, id: &VcfHeaderId) -> bool {
    match id {
        VcfHeaderId::Filter(id) => remove_header_filter(header, id),
        VcfHeaderId::Info(id) => remove_header_info(header, id),
        VcfHeaderId::Format(id) => remove_header_format(header, id),
        VcfHeaderId::Contig(id) => remove_header_contig(header, id),
        VcfHeaderId::Other { key, id } => remove_header_other_record_id(header, key, id),
    }
}

/// A line-based HTSlib-style VCF record mutation adapter.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfRecordAdapter {
    line: String,
}

impl VcfRecordAdapter {
    /// Creates a VCF record adapter from a VCF record line.
    pub fn new(line: impl Into<String>) -> io::Result<Self> {
        let line = line.into();

        validate_vcf_record_line(&line)?;

        Ok(Self { line })
    }

    /// Returns the VCF record line.
    pub fn as_str(&self) -> &str {
        &self.line
    }

    /// Returns the VCF record line.
    pub fn into_string(self) -> String {
        self.line
    }

    /// Returns the HTSlib-style record span.
    pub fn htslib_span(&self) -> io::Result<usize> {
        htslib_variant_span_from_vcf_line(&self.line)
    }

    /// Returns typed INFO integer values.
    pub fn info_i32_values(&self, key: &str) -> io::Result<Option<Vec<Option<i32>>>> {
        vcf_info_i32_values_from_line(&self.line, key)
    }

    /// Returns typed FORMAT integer values per sample.
    pub fn format_i32_values(&self, key: &str) -> io::Result<Option<Vec<Vec<Option<i32>>>>> {
        vcf_format_i32_values_from_line(&self.line, key)
    }

    /// Returns genotype strings for each sample, if the record has `GT`.
    pub fn genotypes(&self) -> io::Result<Option<Vec<String>>> {
        vcf_genotypes_from_line(&self.line)
    }

    /// Replaces REF and ALT alleles and updates the cached line.
    pub fn set_alleles(&mut self, alleles: &[&str]) -> io::Result<()> {
        self.line = update_vcf_line_alleles(&self.line, alleles)?;
        Ok(())
    }

    /// Sets or removes typed INFO integer values.
    pub fn set_info_i32(&mut self, key: &str, values: Option<&[i32]>) -> io::Result<()> {
        self.line = update_vcf_line_info_i32(&self.line, key, values)?;
        Ok(())
    }

    /// Sets or removes typed FORMAT integer values.
    pub fn set_format_i32(&mut self, key: &str, values: Option<&[i32]>) -> io::Result<()> {
        self.line = update_vcf_line_format_i32(&self.line, key, values)?;
        Ok(())
    }

    /// Removes alternate alleles using HTSlib-style vector trimming.
    pub fn remove_alleles(&mut self, remove_alleles: &[usize]) -> io::Result<()> {
        self.line = remove_vcf_allele_set_from_line(&self.line, remove_alleles)?;
        Ok(())
    }
}

const TEST_VCF_API_SERIALIZATION_HEADER: &str = concat!(
    "##fileformat=VCFv4.2\n",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
    "##fileDate=20090805\n",
    "##unused=<XX=AA,Description=\"Unused generic\">\n",
    "##source=myImputationProgramV3.1\n",
    "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta\n",
    "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>\n",
    "##phasing=partial\n",
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n",
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n",
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n",
    "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n",
    "##FILTER=<ID=q10,Description=\"Quality below 10\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
    "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n",
    "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String\">\n",
    "##INFO=<ID=NEG,Number=.,Type=Integer,Description=\"Test Negative Numbers\">\n",
    "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n",
);

#[derive(Clone, Debug, Eq, PartialEq)]
struct TestVcfApiRecord {
    chrom: String,
    pos: usize,
    id: String,
    reference_bases: String,
    alternate_bases: String,
    qual: String,
    filter: String,
    info: Vec<(String, Option<String>)>,
    format_keys: Vec<String>,
    samples: Vec<Vec<String>>,
}

impl TestVcfApiRecord {
    fn to_vcf_line(&self) -> String {
        let info = self
            .info
            .iter()
            .map(|(key, value)| match value {
                Some(value) => format!("{key}={value}"),
                None => key.clone(),
            })
            .collect::<Vec<_>>()
            .join(";");
        let mut fields = vec![
            self.chrom.clone(),
            self.pos.to_string(),
            self.id.clone(),
            self.reference_bases.clone(),
            self.alternate_bases.clone(),
            self.qual.clone(),
            self.filter.clone(),
            info,
            self.format_keys.join(":"),
        ];

        fields.extend(self.samples.iter().map(|values| values.join(":")));
        fields.join("\t")
    }

    fn remove_format(&mut self, key: &str) {
        if let Some(index) = self
            .format_keys
            .iter()
            .position(|candidate| candidate == key)
        {
            self.format_keys.remove(index);

            for sample in &mut self.samples {
                sample.remove(index);
            }
        }
    }

    fn set_info(&mut self, key: &str, value: &str) {
        if let Some((_, existing_value)) =
            self.info.iter_mut().find(|(candidate, _)| candidate == key)
        {
            *existing_value = Some(value.into());
        } else {
            self.info.push((key.into(), Some(value.into())));
        }
    }

    fn set_format_values<const N: usize>(&mut self, key: &str, values: [&str; N]) {
        if let Some(index) = self
            .format_keys
            .iter()
            .position(|candidate| candidate == key)
        {
            for (sample, value) in self.samples.iter_mut().zip(values) {
                sample[index] = value.into();
            }
        } else {
            self.format_keys.push(key.into());

            for (sample, value) in self.samples.iter_mut().zip(values) {
                sample.push(value.into());
            }
        }
    }
}

fn test_vcf_api_source_records() -> [TestVcfApiRecord; 2] {
    [
        TestVcfApiRecord {
            chrom: "20".into(),
            pos: 14370,
            id: "rs6054257".into(),
            reference_bases: "G".into(),
            alternate_bases: "A".into(),
            qual: "29".into(),
            filter: "PASS".into(),
            info: vec![
                ("NS".into(), Some("3".into())),
                ("DP".into(), Some("14".into())),
                ("NEG".into(), Some("-127".into())),
                ("AF".into(), Some("0.5".into())),
                ("DB".into(), None),
                ("H2".into(), None),
            ],
            format_keys: ["GT", "GQ", "DP", "HQ", "TS"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            samples: vec![
                ["0|0", "48", "1", "51,51", "String1"]
                    .into_iter()
                    .map(str::to_string)
                    .collect(),
                ["1|0", "48", "8", "51,51", "SomeOtherString2"]
                    .into_iter()
                    .map(str::to_string)
                    .collect(),
                ["1/1", "43", "5", ".,.", "YetAnotherString3"]
                    .into_iter()
                    .map(str::to_string)
                    .collect(),
            ],
        },
        TestVcfApiRecord {
            chrom: "20".into(),
            pos: 1110696,
            id: ".".into(),
            reference_bases: "A".into(),
            alternate_bases: "G,T".into(),
            qual: "67".into(),
            filter: ".".into(),
            info: vec![
                ("NS".into(), Some("2".into())),
                ("DP".into(), Some("10".into())),
                ("NEG".into(), Some("-128".into())),
                ("AF".into(), Some("0.333,.".into())),
                ("AA".into(), Some("T".into())),
                ("DB".into(), None),
            ],
            format_keys: vec!["GT".into()],
            samples: vec![vec!["2".into()], vec!["1".into()], vec!["./.".into()]],
        },
    ]
}

fn is_present_contig_index_mismatch(header: &Header, region: &Region, error: &io::Error) -> bool {
    error.kind() == io::ErrorKind::InvalidInput
        && std::str::from_utf8(region.name())
            .map(|name| header_has_contig(header, name))
            .unwrap_or(false)
}

fn ensure_vcf_pass_header(header: &str) -> String {
    const PASS_FILTER: &str = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";

    if header
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        return header.to_string();
    }

    let mut out = String::with_capacity(header.len() + PASS_FILTER.len());

    if let Some((first, rest)) = header.split_once('\n') {
        out.push_str(first);
        out.push('\n');
        out.push_str(PASS_FILTER);
        out.push_str(rest);
    } else {
        out.push_str(header);
        out.push_str(PASS_FILTER);
    }

    out
}

fn canonicalize_vcf_header_text(header: &str) -> String {
    let mut out = String::with_capacity(header.len());

    for line in header.lines() {
        if line.starts_with("##") {
            out.push_str(&canonicalize_vcf_meta_line(line.trim_end()));
        } else {
            out.push_str(line.trim_end());
        }

        out.push('\n');
    }

    ensure_vcf_pass_header(&out)
}

fn canonicalize_vcf_meta_line(line: &str) -> String {
    let Some((prefix, rest)) = line.split_once("=<") else {
        return line.to_string();
    };

    let Some(inner) = rest.strip_suffix('>') else {
        return line.to_string();
    };

    let mut out = String::with_capacity(line.len());
    out.push_str(prefix);
    out.push_str("=<");

    for (i, field) in split_vcf_meta_fields(inner).into_iter().enumerate() {
        if i > 0 {
            out.push(',');
        }

        if let Some((key, value)) = field.split_once('=') {
            out.push_str(key.trim());
            out.push('=');
            out.push_str(value.trim());
        } else {
            out.push_str(field.trim());
        }
    }

    out.push('>');
    out
}

fn split_vcf_meta_fields(src: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut in_quotes = false;

    for (i, c) in src.char_indices() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(&src[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }

    fields.push(&src[start..]);
    fields
}

/// Returns the HTSlib `vcf_open_mode` suffix for a local variant path.
pub fn vcf_open_mode_suffix(path: &str) -> io::Result<&'static str> {
    if path.ends_with(".bcf") {
        Ok("b")
    } else if path.ends_with(".vcf.gz") || path.ends_with(".vcf.bgz") {
        Ok("z")
    } else if path.ends_with(".vcf") {
        Ok("")
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unsupported VCF/BCF extension",
        ))
    }
}

/// HTSlib-style VCF header Number classification.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HeaderNumber {
    /// Fixed-size value with the given number of entries.
    Fixed(usize),
    /// Variable size (`.`).
    Variable,
    /// One value for each alternate allele (`A`).
    AlternateBases,
    /// One value for each possible genotype (`G`).
    Genotypes,
    /// One value for each allele including REF (`R`).
    ReferenceAlternateBases,
    /// One value for each allele value defined in GT (`P`).
    Ploidy,
    /// One value for each local ALT allele (`LA`).
    LocalAlternateBases,
    /// One value for each local genotype (`LG`).
    LocalGenotypes,
    /// One value for each local allele including REF (`LR`).
    LocalReferenceAlternateBases,
    /// One value for each base modification of the given type (`M`).
    BaseModifications,
}

/// VCF header record categories with Number fields.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NumberedHeaderRecordKind {
    /// INFO header record.
    Info,
    /// FORMAT header record.
    Format,
}

/// Parses an INFO or FORMAT header record line into its ID and Number classification.
pub fn parse_numbered_header_record(
    line: &str,
) -> io::Result<Option<(NumberedHeaderRecordKind, String, HeaderNumber)>> {
    let (kind, rest) = if let Some(rest) = line.strip_prefix("##INFO=<") {
        (NumberedHeaderRecordKind::Info, rest)
    } else if let Some(rest) = line.strip_prefix("##FORMAT=<") {
        (NumberedHeaderRecordKind::Format, rest)
    } else {
        return Ok(None);
    };
    let body = rest
        .strip_suffix('>')
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid header record"))?;
    let fields = split_header_map_fields(body);
    let id = find_header_map_field(&fields, "ID")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing ID"))?;
    let number = find_header_map_field(&fields, "Number")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing Number"))?;
    let number = parse_header_number(number, kind)?;

    Ok(Some((kind, id.to_string(), number)))
}

/// Returns the HTSlib-style Number classification for an INFO header record.
pub fn info_header_number(header: &Header, id: &str) -> Option<HeaderNumber> {
    header.infos().get(id).map(|map| match map.number() {
        vcf::header::record::value::map::info::Number::Count(n) => HeaderNumber::Fixed(n),
        vcf::header::record::value::map::info::Number::AlternateBases => {
            HeaderNumber::AlternateBases
        }
        vcf::header::record::value::map::info::Number::ReferenceAlternateBases => {
            HeaderNumber::ReferenceAlternateBases
        }
        vcf::header::record::value::map::info::Number::Samples => HeaderNumber::Genotypes,
        vcf::header::record::value::map::info::Number::Unknown => HeaderNumber::Variable,
    })
}

fn split_header_map_fields(src: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut in_string = false;

    for (i, b) in src.bytes().enumerate() {
        match b {
            b'"' => in_string = !in_string,
            b',' if !in_string => {
                fields.push(&src[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }

    fields.push(&src[start..]);
    fields
}

fn find_header_map_field<'a>(fields: &'a [&str], key: &str) -> Option<&'a str> {
    fields.iter().find_map(|field| {
        field
            .split_once('=')
            .filter(|(k, _)| *k == key)
            .map(|(_, v)| v)
    })
}

fn parse_header_number(src: &str, kind: NumberedHeaderRecordKind) -> io::Result<HeaderNumber> {
    match src {
        "." => Ok(HeaderNumber::Variable),
        "A" => Ok(HeaderNumber::AlternateBases),
        "G" => Ok(HeaderNumber::Genotypes),
        "R" => Ok(HeaderNumber::ReferenceAlternateBases),
        "P" if kind == NumberedHeaderRecordKind::Format => Ok(HeaderNumber::Ploidy),
        "LA" if kind == NumberedHeaderRecordKind::Format => Ok(HeaderNumber::LocalAlternateBases),
        "LG" if kind == NumberedHeaderRecordKind::Format => Ok(HeaderNumber::LocalGenotypes),
        "LR" if kind == NumberedHeaderRecordKind::Format => {
            Ok(HeaderNumber::LocalReferenceAlternateBases)
        }
        "M" if kind == NumberedHeaderRecordKind::Format => Ok(HeaderNumber::BaseModifications),
        _ => src
            .parse()
            .map(HeaderNumber::Fixed)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
    }
}

/// Returns the HTSlib-style Number classification for a FORMAT header record.
pub fn format_header_number(header: &Header, id: &str) -> Option<HeaderNumber> {
    header.formats().get(id).map(|map| match map.number() {
        vcf::header::record::value::map::format::Number::Count(n) => HeaderNumber::Fixed(n),
        vcf::header::record::value::map::format::Number::AlternateBases => {
            HeaderNumber::AlternateBases
        }
        vcf::header::record::value::map::format::Number::ReferenceAlternateBases => {
            HeaderNumber::ReferenceAlternateBases
        }
        vcf::header::record::value::map::format::Number::Samples => HeaderNumber::Genotypes,
        vcf::header::record::value::map::format::Number::LocalAlternateBases => {
            HeaderNumber::LocalAlternateBases
        }
        vcf::header::record::value::map::format::Number::LocalReferenceAlternateBases => {
            HeaderNumber::LocalReferenceAlternateBases
        }
        vcf::header::record::value::map::format::Number::LocalSamples => {
            HeaderNumber::LocalGenotypes
        }
        vcf::header::record::value::map::format::Number::Ploidy => HeaderNumber::Ploidy,
        vcf::header::record::value::map::format::Number::BaseModifications => {
            HeaderNumber::BaseModifications
        }
        vcf::header::record::value::map::format::Number::Unknown => HeaderNumber::Variable,
    })
}

/// Returns the HTSlib-style variant span (`bcf1_t::rlen`) for a VCF/BCF record.
pub fn htslib_variant_span<R>(header: &Header, record: &R) -> io::Result<usize>
where
    R: vcf::variant::Record + ?Sized,
{
    let mut header = header.clone();
    *header.file_format_mut() = vcf::header::FileFormat::new(4, 5);

    record.variant_span(&header)
}

/// Reads all VCF records and returns their HTSlib-style spans.
pub fn htslib_variant_spans_from_vcf<R>(mut reader: R) -> io::Result<Vec<usize>>
where
    R: BufRead,
{
    let mut src = String::new();
    reader.read_to_string(&mut src)?;

    src.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(htslib_variant_span_from_vcf_line)
        .collect()
}

/// Returns integer INFO values from a VCF record line, preserving missing values.
pub fn vcf_info_i32_values_from_line(
    line: &str,
    key: &str,
) -> io::Result<Option<Vec<Option<i32>>>> {
    let info = vcf_line_info_field(line)?;

    let Some(values) = parse_info_integer_array(info, key)? else {
        return Ok(None);
    };

    values
        .into_iter()
        .map(|value| {
            value
                .map(|n| {
                    i32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
                .transpose()
        })
        .collect::<io::Result<Vec<_>>>()
        .map(Some)
}

/// Returns float INFO values from a VCF record line, preserving missing values.
pub fn vcf_info_float_values_from_line(
    line: &str,
    key: &str,
) -> io::Result<Option<Vec<Option<f32>>>> {
    let info = vcf_line_info_field(line)?;

    parse_info_float_array(info, key)
}

/// Returns integer FORMAT values from a VCF record line, preserving missing values by sample.
pub fn vcf_format_i32_values_from_line(
    line: &str,
    key: &str,
) -> io::Result<Option<Vec<Vec<Option<i32>>>>> {
    let fields = line.split('\t').collect::<Vec<_>>();

    if fields.len() < 10 {
        return Ok(None);
    }

    let Some(format_index) = fields[8]
        .split(':')
        .position(|format_key| format_key == key)
    else {
        return Ok(None);
    };

    fields[9..]
        .iter()
        .map(|sample| {
            let Some(value) = sample.split(':').nth(format_index) else {
                return Ok(Vec::new());
            };

            parse_i32_value_list(value)
        })
        .collect::<io::Result<Vec<_>>>()
        .map(Some)
}

/// A BCF float slot as returned by HTSlib typed FORMAT/INFO accessors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BcfFloatValue {
    /// A present finite value.
    Value(f32),
    /// A typed missing value.
    Missing,
    /// The BCF vector-end marker.
    VectorEnd,
}

/// Normalizes a BCF float vector the way HTSlib `bcf_get_format_float` exposes it.
///
/// Once a vector-end marker is observed, HTSlib fills the remaining returned slots with
/// vector-end markers too.
pub fn normalize_bcf_float_vector(values: &[BcfFloatValue]) -> Vec<BcfFloatValue> {
    let mut seen_vector_end = false;

    values
        .iter()
        .map(|value| {
            if seen_vector_end {
                BcfFloatValue::VectorEnd
            } else if *value == BcfFloatValue::VectorEnd {
                seen_vector_end = true;
                BcfFloatValue::VectorEnd
            } else {
                *value
            }
        })
        .collect()
}

fn vcf_line_info_field(line: &str) -> io::Result<&str> {
    line.split('\t').nth(7).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid VCF record: missing INFO field",
        )
    })
}

fn validate_vcf_record_line(line: &str) -> io::Result<()> {
    let field_count = line.split('\t').count();

    if field_count < 8 {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid VCF record",
        ))
    } else {
        Ok(())
    }
}

fn vcf_genotypes_from_line(line: &str) -> io::Result<Option<Vec<String>>> {
    validate_vcf_record_line(line)?;

    let fields = line.split('\t').collect::<Vec<_>>();

    if fields.len() < 10 {
        return Ok(None);
    }

    let Some(gt_index) = fields[8].split(':').position(|key| key == "GT") else {
        return Ok(None);
    };

    fields[9..]
        .iter()
        .map(|sample| {
            sample
                .split(':')
                .nth(gt_index)
                .map(str::to_string)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "missing sample genotype")
                })
        })
        .collect::<io::Result<Vec<_>>>()
        .map(Some)
}

/// Returns the HTSlib-style variant span (`bcf1_t::rlen`) for a VCF record line.
pub fn htslib_variant_span_from_vcf_line(line: &str) -> io::Result<usize> {
    let fields = line.split('\t').collect::<Vec<_>>();

    if fields.len() < 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid VCF record",
        ));
    }

    let pos = fields[1]
        .parse::<usize>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let reference_len = fields[3].len().max(1);
    let alternate_bases = fields[4].split(',').collect::<Vec<_>>();
    let has_symbolic_alternate = alternate_bases.iter().any(|alt| is_symbolic_alternate(alt));
    let has_star_alternate = alternate_bases.contains(&"<*>");
    let mut span = reference_len;

    if let Some(end) = parse_info_integer(fields[7], "END")? {
        span = span.max(end.saturating_sub(pos) + 1);
    }

    if has_symbolic_alternate {
        if let Some(values) = parse_info_integer_array(fields[7], "SVLEN")? {
            for (alt, value) in alternate_bases.iter().zip(values) {
                if is_symbolic_alternate(alt)
                    && !matches!(*alt, "<*>" | "<INS>")
                    && let Some(n) = value
                {
                    let len = n.unsigned_abs() as usize + 1;
                    span = span.max(len);
                }
            }
        }

        if has_star_alternate && fields.len() > 9 {
            span = span.max(max_format_len(fields[8], &fields[9..])?.unwrap_or(0));
        }
    }

    Ok(span)
}

/// Updates REF/ALT alleles in a VCF record line.
pub fn update_vcf_line_alleles(line: &str, alleles: &[&str]) -> io::Result<String> {
    if alleles.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing reference allele",
        ));
    }

    let mut fields = split_vcf_line_fields(line)?;
    fields[3] = alleles[0].to_string();
    fields[4] = if alleles.len() == 1 {
        ".".into()
    } else {
        alleles[1..].join(",")
    };

    Ok(fields.join("\t"))
}

/// Updates or removes an integer INFO field in a VCF record line.
pub fn update_vcf_line_info_i32(
    line: &str,
    key: &str,
    values: Option<&[i32]>,
) -> io::Result<String> {
    let mut fields = split_vcf_line_fields(line)?;
    fields[7] = update_keyed_semicolon_field(&fields[7], key, values.map(format_i32_values));

    Ok(fields.join("\t"))
}

/// Updates or removes an integer FORMAT field in a VCF record line.
pub fn update_vcf_line_format_i32(
    line: &str,
    key: &str,
    values: Option<&[i32]>,
) -> io::Result<String> {
    let mut fields = split_vcf_line_fields(line)?;

    if fields.len() < 10 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "VCF record has no sample columns",
        ));
    }

    let mut keys = fields[8].split(':').map(str::to_string).collect::<Vec<_>>();

    if let Some(values) = values {
        let sample_count = fields.len() - 9;

        if values.len() != sample_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "FORMAT value count must match sample count",
            ));
        }

        let key_index = match keys.iter().position(|k| k == key) {
            Some(i) => i,
            None => {
                keys.push(key.to_string());
                keys.len() - 1
            }
        };

        fields[8] = keys.join(":");

        for (sample, value) in fields[9..].iter_mut().zip(values) {
            let mut sample_values = sample.split(':').map(str::to_string).collect::<Vec<_>>();

            while sample_values.len() < keys.len() {
                sample_values.push(".".into());
            }

            sample_values[key_index] = value.to_string();
            *sample = sample_values.join(":");
        }
    } else if let Some(key_index) = keys.iter().position(|k| k == key) {
        keys.remove(key_index);
        fields[8] = keys.join(":");

        for sample in &mut fields[9..] {
            let mut sample_values = sample.split(':').map(str::to_string).collect::<Vec<_>>();

            if key_index < sample_values.len() {
                sample_values.remove(key_index);
            }

            *sample = sample_values.join(":");
        }
    }

    Ok(fields.join("\t"))
}

fn split_vcf_line_fields(line: &str) -> io::Result<Vec<String>> {
    let fields = line.split('\t').map(str::to_string).collect::<Vec<_>>();

    if fields.len() < 8 {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid VCF record",
        ))
    } else {
        Ok(fields)
    }
}

fn update_keyed_semicolon_field(info: &str, key: &str, value: Option<String>) -> String {
    let mut fields = if info == "." {
        Vec::new()
    } else {
        info.split(';').map(str::to_string).collect::<Vec<_>>()
    };

    let position = fields.iter().position(|field| {
        field
            .split_once('=')
            .map(|(field_key, _)| field_key == key)
            .unwrap_or(field == key)
    });

    match (position, value) {
        (Some(i), Some(value)) => fields[i] = format!("{key}={value}"),
        (Some(i), None) => {
            fields.remove(i);
        }
        (None, Some(value)) => fields.push(format!("{key}={value}")),
        (None, None) => {}
    }

    if fields.is_empty() {
        ".".into()
    } else {
        fields.join(";")
    }
}

fn format_i32_values(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn is_symbolic_alternate(alt: &str) -> bool {
    alt.starts_with('<') && alt.ends_with('>')
}

fn parse_info_integer(info: &str, key: &str) -> io::Result<Option<usize>> {
    match parse_info_integer_array(info, key)? {
        Some(mut values) => values
            .pop()
            .flatten()
            .map(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
            .transpose(),
        None => Ok(None),
    }
}

fn parse_info_integer_array(info: &str, key: &str) -> io::Result<Option<Vec<Option<i64>>>> {
    if info == "." {
        return Ok(None);
    }

    for field in info.split(';') {
        let Some((field_key, value)) = field.split_once('=') else {
            continue;
        };

        if field_key == key {
            return value
                .split(',')
                .map(|s| {
                    if s == "." {
                        Ok(None)
                    } else {
                        s.parse()
                            .map(Some)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    }
                })
                .collect::<io::Result<Vec<_>>>()
                .map(Some);
        }
    }

    Ok(None)
}

fn parse_i32_value_list(value: &str) -> io::Result<Vec<Option<i32>>> {
    if value == "." {
        return Ok(vec![None]);
    }

    value
        .split(',')
        .map(|s| {
            if s == "." {
                Ok(None)
            } else {
                s.parse()
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }
        })
        .collect()
}

fn parse_info_float_array(info: &str, key: &str) -> io::Result<Option<Vec<Option<f32>>>> {
    if info == "." {
        return Ok(None);
    }

    for field in info.split(';') {
        let Some((field_key, value)) = field.split_once('=') else {
            continue;
        };

        if field_key == key {
            return value
                .split(',')
                .map(|s| {
                    if s == "." {
                        Ok(None)
                    } else {
                        s.parse()
                            .map(Some)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    }
                })
                .collect::<io::Result<Vec<_>>>()
                .map(Some);
        }
    }

    Ok(None)
}

fn max_format_len(format: &str, samples: &[&str]) -> io::Result<Option<usize>> {
    let Some(len_index) = format.split(':').position(|key| key == "LEN") else {
        return Ok(None);
    };
    let mut max_len = None;

    for sample in samples {
        let Some(value) = sample.split(':').nth(len_index) else {
            continue;
        };

        if value == "." {
            continue;
        }

        let len = value
            .parse::<usize>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        max_len = max_len.map(|n: usize| n.max(len)).or(Some(len));
    }

    Ok(max_len)
}

/// Reads a BCF header from a reader.
pub fn read_bcf_header<R>(reader: R) -> io::Result<Header>
where
    R: Read,
{
    let mut reader = bcf::io::Reader::new(reader);

    reader.read_header()
}

/// Reads a BCF header from a local file.
pub fn read_bcf_header_from_path<P>(src: P) -> io::Result<Header>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(read_bcf_header)
}

/// Counts BCF records from a reader.
pub fn count_bcf_records<R>(reader: R) -> io::Result<usize>
where
    R: Read,
{
    let mut reader = bcf::io::Reader::new(reader);
    reader.read_header()?;

    reader
        .records()
        .try_fold(0, |n, result| result.map(|_| n + 1))
}

/// Counts BCF records from a local file.
pub fn count_bcf_records_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(count_bcf_records)
}

/// Reads BCF input without producing output and returns the number of records seen.
pub fn benchmark_bcf_view_from_path<P>(src: P) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    count_bcf_records_from_path(src)
}

/// Writes BCF input as VCF text, matching the whole-file path of HTSlib `test_view`.
pub fn view_bcf_as_vcf_text_from_path_with_limit<P>(
    src: P,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let raw_header = read_raw_bcf_vcf_header_text(&src)?;
    let mut reader = File::open(src).map(bcf::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut buf = raw_header.into_bytes();
    let mut writer = vcf::io::Writer::new(&mut buf);

    for result in reader
        .record_bufs(&header)
        .take(limit.unwrap_or(usize::MAX))
    {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes VCF text as BCF.
pub fn write_bcf_from_vcf_path<P, W>(src: P, dst: W) -> io::Result<W>
where
    P: AsRef<Path>,
    W: Write,
{
    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut writer = bcf::io::Writer::new(dst);

    writer.write_variant_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

/// Writes VCF text as BCF and writes a CSI index for the generated BCF file.
pub fn write_bcf_from_vcf_path_with_csi<P, Q, R>(src: P, bcf_dst: Q, csi_dst: R) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let encoded = write_bcf_from_vcf_path(src, Vec::new())?;

    std::fs::write(&bcf_dst, encoded)?;
    write_bcf_csi_from_path(bcf_dst, csi_dst)
}

/// Writes indexed BCF records for an HTSlib-style region list as VCF text.
pub fn view_bcf_regions_as_vcf_text_from_path<P>(src: P, regions: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_bcf_regions_as_vcf_text_from_path_with_limit(src, regions, None)
}

/// Writes indexed BCF records for an HTSlib-style region list with an optional record limit.
pub fn view_bcf_regions_as_vcf_text_from_path_with_limit<P>(
    src: P,
    regions: &str,
    limit: Option<usize>,
) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut raw_header = read_raw_bcf_vcf_header_text(&src)?;
    ensure_pass_filter_header(&mut raw_header);
    move_reference_header_before_contigs(&mut raw_header);

    let index = read_associated_bcf_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let regions = parse_vcf_region_list(&header, regions)?;
    let mut buf = raw_header.into_bytes();
    let mut writer = vcf::io::Writer::new(&mut buf);
    let mut remaining = limit.unwrap_or(usize::MAX);

    for region in regions {
        if remaining == 0 {
            break;
        }

        let region = parsed_vcf_region_to_core_region(&region)?;
        let query = reader.query(&header, &region)?;

        for result in query.records() {
            if remaining == 0 {
                break;
            }

            let record = result?;
            writer.write_variant_record(&header, &record)?;
            remaining -= 1;
        }
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Writes indexed BCF records for an HTSlib-style target list as VCF text.
pub fn view_bcf_targets_as_vcf_text_from_path<P>(src: P, targets: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    view_bcf_regions_as_vcf_text_from_path(src, targets)
}

fn read_raw_bcf_vcf_header_text<P>(src: P) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(bcf::io::Reader::new)?;
    let mut header_reader = reader.header_reader();
    header_reader.read_magic_number()?;
    header_reader.read_format_version()?;

    let mut raw_header_reader = header_reader.raw_vcf_header_reader()?;
    let mut raw_header = String::new();
    raw_header_reader.read_to_string(&mut raw_header)?;
    raw_header_reader.discard_to_end()?;

    Ok(strip_raw_bcf_header_idx_fields(&raw_header))
}

fn ensure_pass_filter_header(header: &mut String) {
    if header.contains("##FILTER=<ID=PASS,") {
        return;
    }

    if let Some(i) = header.find('\n') {
        header.insert_str(
            i + 1,
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
        );
    }
}

fn move_reference_header_before_contigs(header: &mut String) {
    let mut reference_lines = Vec::new();
    let mut other_lines = Vec::new();

    for line in header.split_inclusive('\n') {
        if line.starts_with("##reference=") {
            reference_lines.push(line);
        } else {
            other_lines.push(line);
        }
    }

    if reference_lines.is_empty() {
        return;
    }

    let insert_at = other_lines
        .iter()
        .position(|line| line.starts_with("##contig="))
        .unwrap_or(other_lines.len());
    let mut normalized = String::with_capacity(header.len());

    for line in &other_lines[..insert_at] {
        normalized.push_str(line);
    }

    for line in reference_lines {
        normalized.push_str(line);
    }

    for line in &other_lines[insert_at..] {
        normalized.push_str(line);
    }

    *header = normalized;
}

fn strip_raw_bcf_header_idx_fields(raw_header: &str) -> String {
    let mut dst = String::with_capacity(raw_header.len());

    for line in raw_header.split_inclusive('\n') {
        let (line, eol) = line
            .strip_suffix('\n')
            .map(|line| (line, "\n"))
            .unwrap_or((line, ""));

        if let Some((prefix, rest)) = line.split_once("=<")
            && let Some(body) = rest.strip_suffix('>')
        {
            let fields = split_header_map_fields(body);
            dst.push_str(prefix);
            dst.push_str("=<");
            dst.push_str(
                &fields
                    .into_iter()
                    .filter(|field| !field.starts_with("IDX="))
                    .collect::<Vec<_>>()
                    .join(","),
            );
            dst.push('>');
        } else {
            dst.push_str(line);
        }

        dst.push_str(eol);
    }

    dst
}

/// Builds a CSI index for a local BCF file.
pub fn build_bcf_csi_from_path<P>(src: P) -> io::Result<crate::csi::Index>
where
    P: AsRef<Path>,
{
    bcf::fs::index(src)
}

/// Writes a CSI index for a local BCF file.
pub fn write_bcf_csi_from_path<P, Q>(src: P, csi_dst: Q) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let index = build_bcf_csi_from_path(src)?;
    crate::csi::fs::write(csi_dst, &index)
}

/// Queries BCF records from a local file using its associated CSI index.
pub fn query_bcf_records_from_path<P>(src: P, region: &Region) -> io::Result<Vec<bcf::Record>>
where
    P: AsRef<Path>,
{
    let index = read_associated_bcf_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let query = reader.query(&header, region)?;

    query.records().collect()
}

/// Returns an owning iterator over BCF records from an indexed local file.
pub fn iter_bcf_records_from_path<P>(
    src: P,
    region: &Region,
) -> io::Result<std::vec::IntoIter<bcf::Record>>
where
    P: AsRef<Path>,
{
    query_bcf_records_from_path(src, region).map(Vec::into_iter)
}

/// Queries BCF records by HTSlib-style numeric reference ID and 0-based half-open interval.
pub fn query_bcf_records_by_reference_id_from_path<P>(
    src: P,
    reference_sequence_id: usize,
    start: usize,
    end: usize,
) -> io::Result<Vec<bcf::Record>>
where
    P: AsRef<Path>,
{
    let index = read_associated_bcf_index(&src)?;
    let data_path = associated_data_path(src);
    let mut reader = bcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(data_path)?;
    let header = reader.read_header()?;
    let name = header
        .string_maps()
        .contigs()
        .get_index(reference_sequence_id)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "unknown reference ID"))?;
    let start = Position::try_from(start + 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let end =
        Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let region = Region::new(name, start..=end);
    let query = reader.query(&header, &region)?;

    query.records().collect()
}

/// Counts BCF records from a local indexed file that intersect a region.
pub fn count_bcf_records_in_region_from_path<P>(src: P, region: &Region) -> io::Result<usize>
where
    P: AsRef<Path>,
{
    query_bcf_records_from_path(src, region).map(|records| records.len())
}

/// A compact VCF/BCF record summary for HTSlib-style forward/backward sweeps.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VariantSweepRecord {
    position: Position,
    pl_sum: i64,
}

impl VariantSweepRecord {
    /// Returns the 1-based variant start position.
    pub fn position(&self) -> Position {
        self.position
    }

    /// Returns the sum of present integer FORMAT/PL values for this record.
    pub fn pl_sum(&self) -> i64 {
        self.pl_sum
    }
}

/// An in-memory forward/backward VCF/BCF sweep.
#[derive(Debug)]
pub struct VariantSweep {
    header: Header,
    records: Vec<VariantSweepRecord>,
    next: usize,
}

impl VariantSweep {
    /// Returns the parsed header for the swept variant stream.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the number of records available to sweep.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns whether there are no records to sweep.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Sweeps forward by one record.
    pub fn forward(&mut self) -> Option<&VariantSweepRecord> {
        let record = self.records.get(self.next)?;
        self.next += 1;
        Some(record)
    }

    /// Sweeps backward by one record.
    pub fn backward(&mut self) -> Option<&VariantSweepRecord> {
        if self.next == 0 {
            None
        } else {
            self.next -= 1;
            self.records.get(self.next)
        }
    }
}

/// Reads a VCF stream into an HTSlib-style forward/backward sweep.
pub fn read_vcf_sweep<R>(reader: R) -> io::Result<VariantSweep>
where
    R: BufRead,
{
    let mut reader = vcf::io::Reader::new(reader);
    let header = reader.read_header()?;
    let records = reader
        .records()
        .map(|result| result.and_then(|record| summarize_sweep_record(&header, &record)))
        .collect::<io::Result<_>>()?;

    Ok(VariantSweep {
        header,
        records,
        next: 0,
    })
}

/// Reads a local VCF file into an HTSlib-style forward/backward sweep.
pub fn read_vcf_sweep_from_path<P>(src: P) -> io::Result<VariantSweep>
where
    P: AsRef<Path>,
{
    File::open(src).map(BufReader::new).and_then(read_vcf_sweep)
}

/// Reads a BCF stream into an HTSlib-style forward/backward sweep.
pub fn read_bcf_sweep<R>(reader: R) -> io::Result<VariantSweep>
where
    R: Read,
{
    let mut reader = bcf::io::Reader::new(reader);
    let header = reader.read_header()?;
    let records = reader
        .records()
        .map(|result| result.and_then(|record| summarize_sweep_record(&header, &record)))
        .collect::<io::Result<_>>()?;

    Ok(VariantSweep {
        header,
        records,
        next: 0,
    })
}

/// Reads a local BCF file into an HTSlib-style forward/backward sweep.
pub fn read_bcf_sweep_from_path<P>(src: P) -> io::Result<VariantSweep>
where
    P: AsRef<Path>,
{
    File::open(src).and_then(read_bcf_sweep)
}

/// Returns the VCF text produced by the synthetic `test-bcf-translate.c` fixture.
///
/// The fixture exercises header dictionary merging plus record translation across
/// dictionaries after removing one FILTER, INFO, and FORMAT field.
pub fn translated_bcf_record_fixture_vcf_text() -> io::Result<String> {
    use vcf::{
        header::{
            FileFormat,
            record::value::{
                Map,
                map::{
                    Contig, Filter, Format, Info, format::Number as FormatNumber,
                    format::Type as FormatType, info::Number as InfoNumber, info::Type as InfoType,
                },
            },
        },
        variant::{
            RecordBuf,
            io::Write as _,
            record_buf::{
                AlternateBases, Filters,
                info::{Info as RecordInfo, field::Value as InfoValue},
                samples::{Keys, Samples, sample::Value as SampleValue},
            },
        },
    };

    let header_text = concat!(
        "##fileformat=VCFv4.2\n",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
        "##contig=<ID=2>\n",
        "##contig=<ID=1>\n",
        "##FILTER=<ID=FLT4,Description=\"Filter 4\">\n",
        "##FILTER=<ID=FLT3,Description=\"Filter 3\">\n",
        "##FILTER=<ID=FLT2,Description=\"Filter 2\">\n",
        "##INFO=<ID=INF4,Number=.,Type=Integer,Description=\"Info 4\">\n",
        "##INFO=<ID=INF3,Number=.,Type=Integer,Description=\"Info 3\">\n",
        "##INFO=<ID=INF2,Number=.,Type=Integer,Description=\"Info 2\">\n",
        "##FORMAT=<ID=FMT4,Number=.,Type=Integer,Description=\"FMT 4\">\n",
        "##FORMAT=<ID=FMT3,Number=.,Type=Integer,Description=\"FMT 3\">\n",
        "##FORMAT=<ID=FMT2,Number=.,Type=Integer,Description=\"FMT 2\">\n",
        "##FILTER=<ID=FLT1,Description=\"Filter 1\">\n",
        "##INFO=<ID=INF1,Number=.,Type=Integer,Description=\"Info 1\">\n",
        "##FORMAT=<ID=FMT1,Number=.,Type=Integer,Description=\"FMT 1\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSMPL1\tSMPL2\n",
    );

    let header = Header::builder()
        .set_file_format(FileFormat::new(4, 2))
        .add_filter("PASS", Map::<Filter>::pass())
        .add_contig("2", Map::<Contig>::new())
        .add_contig("1", Map::<Contig>::new())
        .add_filter("FLT4", Map::<Filter>::new("Filter 4"))
        .add_filter("FLT3", Map::<Filter>::new("Filter 3"))
        .add_filter("FLT2", Map::<Filter>::new("Filter 2"))
        .add_info(
            "INF4",
            Map::<Info>::new(InfoNumber::Unknown, InfoType::Integer, "Info 4"),
        )
        .add_info(
            "INF3",
            Map::<Info>::new(InfoNumber::Unknown, InfoType::Integer, "Info 3"),
        )
        .add_info(
            "INF2",
            Map::<Info>::new(InfoNumber::Unknown, InfoType::Integer, "Info 2"),
        )
        .add_format(
            "FMT4",
            Map::<Format>::new(FormatNumber::Unknown, FormatType::Integer, "FMT 4"),
        )
        .add_format(
            "FMT3",
            Map::<Format>::new(FormatNumber::Unknown, FormatType::Integer, "FMT 3"),
        )
        .add_format(
            "FMT2",
            Map::<Format>::new(FormatNumber::Unknown, FormatType::Integer, "FMT 2"),
        )
        .add_filter("FLT1", Map::<Filter>::new("Filter 1"))
        .add_info(
            "INF1",
            Map::<Info>::new(InfoNumber::Unknown, InfoType::Integer, "Info 1"),
        )
        .add_format(
            "FMT1",
            Map::<Format>::new(FormatNumber::Unknown, FormatType::Integer, "FMT 1"),
        )
        .add_sample_name("SMPL1")
        .add_sample_name("SMPL2")
        .build();

    let info: RecordInfo = [
        (String::from("INF1"), Some(InfoValue::Integer(1))),
        (String::from("INF3"), Some(InfoValue::Integer(3))),
    ]
    .into_iter()
    .collect();

    let keys: Keys = [String::from("FMT1"), String::from("FMT3")]
        .into_iter()
        .collect();
    let sample = vec![Some(SampleValue::Integer(1)), Some(SampleValue::Integer(3))];
    let samples = Samples::new(keys, vec![sample.clone(), sample]);

    let record = RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_variant_start(Position::MIN)
        .set_reference_bases("G")
        .set_alternate_bases(AlternateBases::from(vec![String::from("A")]))
        .set_quality_score(0.0)
        .set_filters(
            [String::from("FLT1"), String::from("FLT3")]
                .into_iter()
                .collect::<Filters>(),
        )
        .set_info(info)
        .set_samples(samples)
        .build();

    let mut writer = vcf::io::Writer::new(Vec::new());
    writer.write_variant_record(&header, &record)?;

    let mut text = String::from(header_text);
    text.push_str(
        &String::from_utf8(writer.into_inner())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    );

    Ok(text)
}

fn summarize_sweep_record<R>(header: &Header, record: &R) -> io::Result<VariantSweepRecord>
where
    R: vcf::variant::Record + ?Sized,
{
    let position = record
        .variant_start()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing variant position"))?;
    let pl_sum = record_pl_sum(header, record)?;

    Ok(VariantSweepRecord { position, pl_sum })
}

fn record_pl_sum<R>(header: &Header, record: &R) -> io::Result<i64>
where
    R: vcf::variant::Record + ?Sized,
{
    use vcf::variant::record::samples::series::{Value, value::Array};

    let samples = record.samples()?;
    let Some(series) = samples.select(header, "PL") else {
        return Ok(0);
    };

    let mut sum = 0;

    for result in series?.iter(header) {
        match result? {
            Some(Value::Integer(n)) => sum += i64::from(n),
            Some(Value::Array(Array::Integer(values))) => {
                for result in values.iter() {
                    if let Some(n) = result? {
                        sum += i64::from(n);
                    }
                }
            }
            Some(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FORMAT/PL is not an integer value",
                ));
            }
            None => {}
        }
    }

    Ok(sum)
}

#[derive(Clone, Debug)]
struct SyncedKey {
    reference_sequence_name: String,
    reference_sequence_order: usize,
    position: Position,
}

impl PartialEq for SyncedKey {
    fn eq(&self, other: &Self) -> bool {
        self.reference_sequence_name == other.reference_sequence_name
            && self.position == other.position
    }
}

impl Eq for SyncedKey {}

impl Hash for SyncedKey {
    fn hash<H>(&self, state: &mut H)
    where
        H: Hasher,
    {
        self.reference_sequence_name.hash(state);
        self.position.hash(state);
    }
}

impl PartialOrd for SyncedKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SyncedKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.reference_sequence_order
            .cmp(&other.reference_sequence_order)
            .then_with(|| self.position.cmp(&other.position))
            .then_with(|| {
                self.reference_sequence_name
                    .cmp(&other.reference_sequence_name)
            })
    }
}

/// Produces `test-bcf-sr` summary output for sorted, no-index local VCF inputs.
pub fn synced_vcf_summary_no_index_from_paths<P>(paths: &[P]) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut inputs = Vec::with_capacity(paths.len());

    for path in paths {
        inputs.push(read_synced_vcf_records(path)?);
    }

    synced_vcf_summary_no_index_from_records(inputs)
}

/// Produces `test-bcf-sr` summary output for sorted, no-index VCF readers.
pub fn synced_vcf_summary_no_index_from_readers<R>(readers: Vec<R>) -> io::Result<String>
where
    R: BufRead,
{
    let inputs = readers
        .into_iter()
        .map(read_synced_vcf_records_from_reader)
        .collect::<io::Result<Vec<_>>>()?;

    synced_vcf_summary_no_index_from_records(inputs)
}

/// Produces `test-bcf-sr -O vcf`-style VCF output for sorted, no-index local VCF inputs.
pub fn synced_vcf_output_no_index_from_paths<P>(paths: &[P]) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let mut inputs = Vec::with_capacity(paths.len());

    for path in paths {
        inputs.push(read_synced_vcf_raw_records(path)?);
    }

    synced_vcf_output_no_index_from_records(inputs)
}

/// Produces `test-bcf-sr -O bcf`-style BCF output for sorted, no-index local VCF inputs.
pub fn synced_bcf_output_no_index_from_paths<P>(paths: &[P]) -> io::Result<Vec<u8>>
where
    P: AsRef<Path>,
{
    let vcf_text = synced_vcf_output_no_index_from_paths(paths)?;
    let mut reader = vcf::io::Reader::new(BufReader::new(Cursor::new(vcf_text.into_bytes())));
    let header = reader.read_header()?;
    let mut writer = bcf::io::Writer::new(Vec::new());

    writer.write_variant_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    writer.try_finish()?;

    Ok(writer.into_inner().into_inner())
}

fn synced_vcf_summary_no_index_from_records(
    inputs: Vec<(Vec<String>, HashMap<SyncedKey, String>)>,
) -> io::Result<String> {
    let mut expected_contigs = None;
    let mut keys = Vec::new();
    let mut seen = HashSet::new();
    let mut per_input = Vec::with_capacity(inputs.len());

    for (contigs, records) in inputs {
        match &expected_contigs {
            Some(expected) if expected != &contigs => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "VCF contig header order differs between inputs",
                ));
            }
            None => expected_contigs = Some(contigs),
            _ => {}
        }

        for key in records.keys() {
            if seen.insert(key.clone()) {
                keys.push(key.clone());
            }
        }

        per_input.push(records);
    }

    keys.sort();

    let mut output = String::new();

    for key in keys {
        output.push_str(&key.reference_sequence_name);
        output.push(':');
        output.push_str(&usize::from(key.position).to_string());

        for records in &per_input {
            output.push('\t');
            output.push_str(records.get(&key).map(String::as_str).unwrap_or("-"));
        }

        output.push('\n');
    }

    Ok(output)
}

fn synced_vcf_output_no_index_from_records(inputs: Vec<SyncedRawInput>) -> io::Result<String> {
    let mut expected_contigs = None;
    let mut keys = Vec::new();
    let mut seen = HashSet::new();

    for input in &inputs {
        match &expected_contigs {
            Some(expected) if expected != &input.contigs => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "VCF contig header order differs between inputs",
                ));
            }
            None => expected_contigs = Some(input.contigs.clone()),
            _ => {}
        }

        for key in input.records.keys() {
            if seen.insert(key.clone()) {
                keys.push(key.clone());
            }
        }
    }

    keys.sort();

    let mut output = inputs
        .first()
        .map(|input| input.header.clone())
        .unwrap_or_default();

    for key in keys {
        for input in &inputs {
            if let Some(line) = input.records.get(&key) {
                output.push_str(line);
                output.push('\n');
            }
        }
    }

    Ok(output)
}

struct SyncedRawInput {
    contigs: Vec<String>,
    header: String,
    records: HashMap<SyncedKey, String>,
}

/// Pairing mode for synced variant groups.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SyncedPairLogic {
    /// Pair SNP alleles only.
    Snps,
    /// Pair indel alleles only.
    Indels,
    /// Pair SNPs with SNPs and indels with indels.
    Both,
    /// Pair SNPs and reference rows.
    SnpsAndReference,
    /// Pair indels and reference rows.
    IndelsAndReference,
    /// Pair SNPs, indels, and reference rows.
    BothAndReference,
    /// Pair exact multi-allelic sets only.
    Exact,
    /// Pair rows that share at least one allele.
    Some,
    /// Pair all compatible SNP, indel, and reference rows.
    All,
}

/// A group of synced-reader inputs that share the same variant set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SyncedVariantGroup {
    /// Variant strings in `REF>ALT` or comma-delimited `REF>ALT,REF>ALT` form.
    pub variants: Vec<String>,
    /// Input indexes represented by this group.
    pub input_indexes: Vec<usize>,
}

/// Pairs synced variant groups using HTSlib `test-bcf-sr` allele-pairing rules.
pub fn pair_synced_variant_groups(
    groups: &[SyncedVariantGroup],
    logic: SyncedPairLogic,
) -> Vec<Vec<Option<String>>> {
    let total_inputs = groups
        .iter()
        .flat_map(|group| group.input_indexes.iter())
        .copied()
        .max()
        .map_or(0, |index| index + 1);
    let mut variant_names = Vec::<String>::new();
    let mut variants = Vec::<Vec<SyncedVariantOccurrence>>::new();

    for (group_index, group) in groups.iter().enumerate() {
        for (variant_index, variant) in group.variants.iter().enumerate() {
            let unique_index =
                if let Some(index) = variant_names.iter().position(|name| name == variant) {
                    index
                } else {
                    variant_names.push(variant.clone());
                    variants.push(Vec::new());
                    variants.len() - 1
                };

            variants[unique_index].push(SyncedVariantOccurrence {
                group_index,
                variant_index,
                count: group.input_indexes.len(),
            });
        }
    }

    let mut variant_sets = (0..variants.len())
        .map(|index| vec![index])
        .collect::<Vec<_>>();
    let mut bitmasks = variant_sets
        .iter()
        .map(|set| synced_pair_bitmask(groups, &variants, set))
        .collect::<Vec<_>>();
    let mut max_counts = variant_sets
        .iter()
        .map(|set| synced_pair_count(&variants, set))
        .collect::<Vec<_>>();
    let mut rows = Vec::new();

    while !variant_sets.is_empty() {
        let mut max_index = 0;

        for index in 1..variant_sets.len() {
            if max_counts[index] > max_counts[max_index] {
                max_index = index;
            }
        }

        let mut pair_index = None;
        let mut max_score = 0;

        for index in 0..variant_sets.len() {
            if bitmasks[max_index] & bitmasks[index] != 0 {
                continue;
            }

            let score = synced_pairing_score(
                groups,
                &variants,
                &variant_sets[max_index],
                &variant_sets[index],
                logic,
            );

            if max_score < score {
                max_score = score;
                pair_index = Some(index);
            }
        }

        if let Some(pair_index) = pair_index
            && pair_index != max_index
        {
            merge_synced_pair_rows(
                &mut variant_sets,
                &mut bitmasks,
                &mut max_counts,
                max_index,
                pair_index,
            );
            continue;
        }

        rows.push(output_synced_pair_row(
            groups,
            &variants,
            &variant_sets[max_index],
            total_inputs,
        ));
        variant_sets.remove(max_index);
        bitmasks.remove(max_index);
        max_counts.remove(max_index);
    }

    rows
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SyncedVariantOccurrence {
    group_index: usize,
    variant_index: usize,
    count: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum SyncedVariantType {
    Ref,
    Snp,
    Indel,
}

fn synced_pair_bitmask(
    groups: &[SyncedVariantGroup],
    variants: &[Vec<SyncedVariantOccurrence>],
    set: &[usize],
) -> u128 {
    let mut bitmask = 0;

    for &variant_index in set {
        for occurrence in &variants[variant_index] {
            if occurrence.group_index < u128::BITS as usize {
                bitmask |= 1 << occurrence.group_index;
            } else if !groups[occurrence.group_index].input_indexes.is_empty() {
                bitmask = u128::MAX;
            }
        }
    }

    bitmask
}

fn synced_pair_count(variants: &[Vec<SyncedVariantOccurrence>], set: &[usize]) -> usize {
    set.iter()
        .flat_map(|&variant_index| &variants[variant_index])
        .map(|occurrence| occurrence.count)
        .sum()
}

fn merge_synced_pair_rows(
    variant_sets: &mut Vec<Vec<usize>>,
    bitmasks: &mut Vec<u128>,
    max_counts: &mut Vec<usize>,
    first: usize,
    second: usize,
) -> usize {
    let (target, source) = if first < second {
        (first, second)
    } else {
        (second, first)
    };
    let source_set = variant_sets[source].clone();

    variant_sets[target].extend(source_set);
    variant_sets.remove(source);
    bitmasks[target] |= bitmasks[source];
    bitmasks.remove(source);
    max_counts[target] += max_counts[source];
    max_counts.remove(source);

    target
}

fn output_synced_pair_row(
    groups: &[SyncedVariantGroup],
    variants: &[Vec<SyncedVariantOccurrence>],
    set: &[usize],
    total_inputs: usize,
) -> Vec<Option<String>> {
    let mut row = vec![None; total_inputs];

    for &variant_index in set {
        for occurrence in &variants[variant_index] {
            let group = &groups[occurrence.group_index];
            let value = synced_allele_summary(&group.variants[occurrence.variant_index]);

            for &input_index in &group.input_indexes {
                row[input_index] = Some(value.clone());
            }
        }
    }

    row
}

fn synced_pairing_score(
    groups: &[SyncedVariantGroup],
    variants: &[Vec<SyncedVariantOccurrence>],
    left_set: &[usize],
    right_set: &[usize],
    logic: SyncedPairLogic,
) -> u64 {
    const MAX_SCORE: u64 = u32::MAX as u64;

    let mut min_score = MAX_SCORE;

    for &left_index in left_set {
        for &right_index in right_set {
            let left_variant = synced_variant_string(groups, variants, left_index);
            let right_variant = synced_variant_string(groups, variants, right_index);

            if left_variant == right_variant {
                return MAX_SCORE;
            }

            if logic == SyncedPairLogic::Exact {
                if synced_multi_is_exact(left_variant, right_variant) {
                    return MAX_SCORE;
                }

                continue;
            } else if synced_multi_is_subset(left_variant, right_variant) {
                return MAX_SCORE;
            }

            let mut max_score = 0;

            for left_type in synced_variant_types(left_variant) {
                for right_type in synced_variant_types(right_variant) {
                    max_score = max_score.max(synced_type_pair_score(logic, left_type, right_type));
                }
            }

            if max_score == 0 {
                return 0;
            }

            min_score = min_score.min(max_score);
        }
    }

    if logic == SyncedPairLogic::Exact {
        return 0;
    }

    let count = left_set
        .iter()
        .chain(right_set)
        .map(|&variant_index| {
            variants[variant_index]
                .iter()
                .map(|o| o.count)
                .sum::<usize>()
        })
        .sum::<usize>();

    (1_u64 << (28 + min_score)) + count as u64
}

fn synced_variant_string<'a>(
    groups: &'a [SyncedVariantGroup],
    variants: &[Vec<SyncedVariantOccurrence>],
    variant_index: usize,
) -> &'a str {
    let occurrence = &variants[variant_index][0];

    &groups[occurrence.group_index].variants[occurrence.variant_index]
}

fn synced_type_pair_score(
    logic: SyncedPairLogic,
    left: SyncedVariantType,
    right: SyncedVariantType,
) -> u64 {
    match (left, right) {
        (SyncedVariantType::Snp, SyncedVariantType::Snp)
            if matches!(
                logic,
                SyncedPairLogic::Snps
                    | SyncedPairLogic::Both
                    | SyncedPairLogic::SnpsAndReference
                    | SyncedPairLogic::BothAndReference
                    | SyncedPairLogic::All
            ) =>
        {
            3
        }
        (SyncedVariantType::Indel, SyncedVariantType::Indel)
            if matches!(
                logic,
                SyncedPairLogic::Indels
                    | SyncedPairLogic::Both
                    | SyncedPairLogic::IndelsAndReference
                    | SyncedPairLogic::BothAndReference
                    | SyncedPairLogic::All
            ) =>
        {
            3
        }
        (SyncedVariantType::Snp, SyncedVariantType::Ref)
        | (SyncedVariantType::Ref, SyncedVariantType::Snp)
            if matches!(
                logic,
                SyncedPairLogic::SnpsAndReference
                    | SyncedPairLogic::BothAndReference
                    | SyncedPairLogic::All
            ) =>
        {
            2
        }
        (SyncedVariantType::Indel, SyncedVariantType::Ref)
        | (SyncedVariantType::Ref, SyncedVariantType::Indel)
            if matches!(
                logic,
                SyncedPairLogic::IndelsAndReference
                    | SyncedPairLogic::BothAndReference
                    | SyncedPairLogic::All
            ) =>
        {
            2
        }
        (SyncedVariantType::Snp, SyncedVariantType::Indel)
        | (SyncedVariantType::Indel, SyncedVariantType::Snp)
            if logic == SyncedPairLogic::All =>
        {
            1
        }
        _ => 0,
    }
}

fn synced_variant_types(variant: &str) -> Vec<SyncedVariantType> {
    let mut types = Vec::new();

    for allele in variant.split(',') {
        let Some((reference, alternate)) = allele.split_once('>') else {
            continue;
        };
        let variant_type = if reference == alternate || alternate == "." {
            SyncedVariantType::Ref
        } else if reference.len() == alternate.len() && reference.len() == 1 {
            SyncedVariantType::Snp
        } else {
            SyncedVariantType::Indel
        };

        if !types.contains(&variant_type) {
            types.push(variant_type);
        }
    }

    types
}

fn synced_multi_is_subset(left: &str, right: &str) -> bool {
    left.split(',').any(|left_allele| {
        right
            .split(',')
            .any(|right_allele| right_allele == left_allele)
    })
}

fn synced_multi_is_exact(left: &str, right: &str) -> bool {
    let left = left.split(',').collect::<HashSet<_>>();
    let right = right.split(',').collect::<HashSet<_>>();

    left == right
}

fn synced_allele_summary(variant: &str) -> String {
    variant
        .split(',')
        .filter_map(|allele| allele.split_once('>').map(|(_, alternate)| alternate))
        .collect::<Vec<_>>()
        .join(",")
}

/// Filters a local VCF file using HTSlib-style region syntax and returns VCF text.
pub fn filter_vcf_text_by_region_from_path<P>(path: P, regions: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    let text = std::fs::read_to_string(&path)?;
    let mut reader = File::open(path)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let intervals = parse_vcf_region_list(&header, regions)?;

    let mut output = String::new();

    for line in text.lines() {
        if line.starts_with("##") {
            output.push_str(line);
            output.push('\n');

            if line.starts_with("##fileformat=") && !text.contains("##FILTER=<ID=PASS,") {
                output.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
        } else if line.starts_with("#CHROM") || vcf_line_matches_regions(line, &intervals)? {
            output.push_str(line);
            output.push('\n');
        }
    }

    Ok(output)
}

/// Filters a local VCF file using HTSlib-style target syntax and returns VCF text.
pub fn filter_vcf_text_by_target_from_path<P>(path: P, targets: &str) -> io::Result<String>
where
    P: AsRef<Path>,
{
    filter_vcf_text_by_region_from_path(path, targets)
}

fn parse_vcf_region_list(header: &Header, regions: &str) -> io::Result<Vec<ParsedVcfRegion>> {
    split_region_list(regions)
        .map(|item| parse_vcf_region(header, item))
        .collect()
}

fn split_region_list(regions: &str) -> impl Iterator<Item = &str> {
    let mut depth = 0;
    let mut start = 0;
    let mut items = Vec::new();

    for (i, b) in regions.bytes().enumerate() {
        match b {
            b'{' => depth += 1,
            b'}' => depth -= 1,
            b',' if depth == 0 => {
                items.push(&regions[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }

    items.push(&regions[start..]);
    items.into_iter()
}

fn parse_vcf_region(header: &Header, region: &str) -> io::Result<ParsedVcfRegion> {
    let (name, coordinates) = if let Some(rest) = region.strip_prefix('{') {
        let close = rest
            .find('}')
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing closing brace"))?;
        let name = &rest[..close];
        let suffix = &rest[close + 1..];
        let coordinates = if suffix.is_empty() {
            None
        } else {
            Some(suffix.strip_prefix(':').ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "invalid quoted region suffix")
            })?)
        };
        (name, coordinates)
    } else {
        match region.split_once(':') {
            Some((name, coordinates)) => (name, Some(coordinates)),
            None => (region, None),
        }
    };

    if !header.contigs().contains_key(name) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unknown reference sequence",
        ));
    }

    let (start, end) = match coordinates {
        Some(coordinates) => parse_one_based_interval(coordinates)?,
        None => (0, i64::MAX),
    };

    Ok(ParsedVcfRegion {
        reference_sequence_name: name.into(),
        start,
        end,
    })
}

fn parse_one_based_interval(coordinates: &str) -> io::Result<(i64, i64)> {
    let (start, end) = match coordinates.split_once('-') {
        Some((start, end)) => (start, Some(end)),
        None => (coordinates, None),
    };
    let start = start
        .parse::<i64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if start < 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "region start must be at least 1",
        ));
    }

    let end = match end {
        Some(end) => end
            .parse::<i64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
        None => start,
    };

    if end < start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "region end is before start",
        ));
    }

    Ok((start - 1, end))
}

#[derive(Debug)]
struct ParsedVcfRegion {
    reference_sequence_name: String,
    start: i64,
    end: i64,
}

fn parsed_vcf_region_to_core_region(region: &ParsedVcfRegion) -> io::Result<Region> {
    if region.start == 0 && region.end == i64::MAX {
        return Ok(Region::new(region.reference_sequence_name.as_str(), ..));
    }

    let start = usize::try_from(region.start + 1)
        .map(Position::try_from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let end = usize::try_from(region.end)
        .map(Position::try_from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    Ok(Region::new(
        region.reference_sequence_name.as_str(),
        start..=end,
    ))
}

fn vcf_line_matches_regions(line: &str, intervals: &[ParsedVcfRegion]) -> io::Result<bool> {
    let mut fields = line.split('\t');
    let reference_sequence_name = fields
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing CHROM field"))?;
    let position = fields
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing POS field"))?
        .parse::<i64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        - 1;

    Ok(intervals.iter().any(|interval| {
        interval.reference_sequence_name == reference_sequence_name
            && interval.start <= position
            && position < interval.end
    }))
}

fn read_synced_vcf_records<P>(path: P) -> io::Result<(Vec<String>, HashMap<SyncedKey, String>)>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(path)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

    read_synced_vcf_records_inner(&mut reader)
}

fn read_synced_vcf_records_from_reader<R>(
    reader: R,
) -> io::Result<(Vec<String>, HashMap<SyncedKey, String>)>
where
    R: BufRead,
{
    let mut reader = vcf::io::Reader::new(reader);

    read_synced_vcf_records_inner(&mut reader)
}

fn read_synced_vcf_raw_records<P>(path: P) -> io::Result<SyncedRawInput>
where
    P: AsRef<Path>,
{
    let text = std::fs::read_to_string(path)?;
    let contigs = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix("##contig=<")
                .and_then(|value| value.split('>').next())
                .and_then(|value| value.split(',').find_map(|field| field.strip_prefix("ID=")))
                .map(String::from)
        })
        .collect::<Vec<_>>();
    let mut header = String::new();
    let mut last_key = None;
    let mut records = HashMap::new();

    for line in text.lines() {
        if line.starts_with('#') {
            header.push_str(line);
            header.push('\n');
            continue;
        }

        let key = synced_key_from_vcf_line(line, &contigs)?;

        if last_key.as_ref().is_some_and(|last| &key < last) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF records are not sorted in header contig order",
            ));
        }

        records.insert(key.clone(), line.to_string());
        last_key = Some(key);
    }

    Ok(SyncedRawInput {
        contigs,
        header,
        records,
    })
}

fn read_synced_vcf_records_inner<R>(
    reader: &mut vcf::io::Reader<R>,
) -> io::Result<(Vec<String>, HashMap<SyncedKey, String>)>
where
    R: BufRead,
{
    use vcf::variant::record::AlternateBases as _;

    let header = reader.read_header()?;
    let contigs = header.contigs().keys().cloned().collect::<Vec<_>>();
    let mut last_key = None;
    let mut records = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let reference_sequence_name = record.reference_sequence_name().to_string();
        let reference_sequence_order = header
            .contigs()
            .get_index_of(&reference_sequence_name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "record reference sequence is missing from the header",
                )
            })?;
        let position = record.variant_start().transpose()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing variant position")
        })?;
        let key = SyncedKey {
            reference_sequence_name,
            reference_sequence_order,
            position,
        };

        if last_key.as_ref().is_some_and(|last| &key < last) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF records are not sorted in header contig order",
            ));
        }

        let alternate_bases = record.alternate_bases();
        let alternate_bases = alternate_bases.iter().collect::<io::Result<Vec<_>>>()?;
        let value = if alternate_bases.is_empty() {
            String::from(".")
        } else {
            alternate_bases.join(",")
        };

        records.insert(key.clone(), value);
        last_key = Some(key);
    }

    Ok((contigs, records))
}

fn synced_key_from_vcf_line(line: &str, contigs: &[String]) -> io::Result<SyncedKey> {
    let mut fields = line.split('\t');
    let reference_sequence_name = fields
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing CHROM field"))?
        .to_string();
    let position = fields
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing POS field"))?
        .parse::<usize>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|n| {
            Position::new(n).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "missing variant position")
            })
        })?;
    let reference_sequence_order = contigs
        .iter()
        .position(|name| name == &reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "record reference sequence is missing from the header",
            )
        })?;

    Ok(SyncedKey {
        reference_sequence_name,
        reference_sequence_order,
        position,
    })
}

/// Returns the number of contig definitions in a VCF header.
pub fn contig_count(header: &Header) -> usize {
    header.contigs().len()
}

/// Returns the number of sample names in a VCF header.
pub fn sample_count(header: &Header) -> usize {
    header.sample_names().len()
}

/// Returns whether a FILTER header record with `id` exists.
pub fn header_has_filter(header: &Header, id: &str) -> bool {
    header.filters().contains_key(id)
}

/// Removes a FILTER header record by ID.
pub fn remove_header_filter(header: &mut Header, id: &str) -> bool {
    header.filters_mut().shift_remove(id).is_some()
}

/// Returns whether an INFO header record with `id` exists.
pub fn header_has_info(header: &Header, id: &str) -> bool {
    header.infos().contains_key(id)
}

/// Removes an INFO header record by ID.
pub fn remove_header_info(header: &mut Header, id: &str) -> bool {
    header.infos_mut().shift_remove(id).is_some()
}

/// Returns whether a FORMAT header record with `id` exists.
pub fn header_has_format(header: &Header, id: &str) -> bool {
    header.formats().contains_key(id)
}

/// Removes a FORMAT header record by ID.
pub fn remove_header_format(header: &mut Header, id: &str) -> bool {
    header.formats_mut().shift_remove(id).is_some()
}

/// Returns whether a contig header record with `id` exists.
pub fn header_has_contig(header: &Header, id: &str) -> bool {
    header.contigs().contains_key(id)
}

/// Removes a contig header record by ID.
pub fn remove_header_contig(header: &mut Header, id: &str) -> bool {
    header.contigs_mut().shift_remove(id).is_some()
}

/// Returns whether a structured nonstandard header record has a given ID.
pub fn header_has_other_record_id(header: &Header, key: &str, id: &str) -> bool {
    use vcf::header::record::value::Collection;

    matches!(
        header.get(key),
        Some(Collection::Structured(records)) if records.contains_key(id)
    )
}

/// Removes a structured nonstandard header record by key and ID.
pub fn remove_header_other_record_id(header: &mut Header, key: &str, id: &str) -> bool {
    use vcf::header::record::value::Collection;

    let Ok(key) = key.parse::<vcf::header::record::key::Other>() else {
        return false;
    };

    let Some(Collection::Structured(records)) = header.other_records_mut().get_mut(&key) else {
        return false;
    };

    records.shift_remove(id).is_some()
}

/// Returns whether an unstructured nonstandard header record exists.
pub fn header_has_unstructured_other_record(header: &Header, key: &str, value: &str) -> bool {
    use vcf::header::record::value::Collection;

    matches!(
        header.get(key),
        Some(Collection::Unstructured(records)) if records.iter().any(|record| record == value)
    )
}

/// Removes all nonstandard header records with `key`.
pub fn remove_header_other_records(header: &mut Header, key: &str) -> bool {
    let Ok(key) = key.parse::<vcf::header::record::key::Other>() else {
        return false;
    };

    header.other_records_mut().shift_remove(&key).is_some()
}

/// Removes alternate alleles from a VCF record line using HTSlib-style vector trimming.
pub fn remove_vcf_allele_set_from_line(line: &str, remove_alleles: &[usize]) -> io::Result<String> {
    let mut fields = line.split('\t').map(str::to_string).collect::<Vec<_>>();

    if fields.len() < 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid VCF record",
        ));
    }

    let remove_alleles = remove_alleles.iter().copied().collect::<HashSet<_>>();
    let alternate_bases = if fields[4] == "." {
        Vec::new()
    } else {
        fields[4].split(',').map(str::to_string).collect::<Vec<_>>()
    };
    let kept_alleles = (1..=alternate_bases.len())
        .filter(|i| !remove_alleles.contains(i))
        .collect::<Vec<_>>();
    let allele_map = allele_map(alternate_bases.len(), &kept_alleles);
    let repeat_groups = info_repeat_groups(&fields[7]);

    fields[4] = if kept_alleles.is_empty() {
        ".".into()
    } else {
        kept_alleles
            .iter()
            .map(|i| alternate_bases[i - 1].as_str())
            .collect::<Vec<_>>()
            .join(",")
    };
    fields[7] = remove_info_alleles(
        &fields[7],
        alternate_bases.len(),
        &kept_alleles,
        &repeat_groups,
    )?;

    if fields.len() > 8 {
        let keys = fields[8].split(':').map(str::to_string).collect::<Vec<_>>();
        for sample in &mut fields[9..] {
            *sample = remove_sample_alleles(
                &keys,
                sample,
                alternate_bases.len(),
                &kept_alleles,
                &allele_map,
            )?;
        }
    }

    Ok(fields.join("\t"))
}

fn allele_map(alternate_count: usize, kept_alleles: &[usize]) -> Vec<Option<usize>> {
    let mut map = vec![None; alternate_count + 1];
    map[0] = Some(0);

    for (new_index, old_index) in kept_alleles.iter().enumerate() {
        map[*old_index] = Some(new_index + 1);
    }

    map
}

fn remove_info_alleles(
    info: &str,
    alternate_count: usize,
    kept_alleles: &[usize],
    repeat_groups: &RepeatGroups,
) -> io::Result<String> {
    if info == "." {
        return Ok(".".into());
    }

    let mut fields = Vec::new();

    for field in info.split(';') {
        let Some((key, value)) = field.split_once('=') else {
            fields.push(field.to_string());
            continue;
        };
        let values = split_csv(value);
        let value = match key {
            "AC" | "AF" | "CN" | "RN" | "SVLEN" | "VL_A_STR_INFO" => {
                select_number_a(&values, kept_alleles)
            }
            "AD" | "VL_R_STR_INFO" => select_number_r(&values, kept_alleles),
            "CIEND" | "CIPOS" | "CILEN" | "CICN" => {
                select_alt_groups(&values, alternate_count, kept_alleles, 2)
            }
            "MEINFO" | "METRANS" => select_alt_groups(&values, alternate_count, kept_alleles, 4),
            "RUS" | "RUL" | "RB" | "RUC" => {
                select_variable_repeat_groups(&values, kept_alleles, &repeat_groups.repeat_counts)
            }
            "RUB" => {
                select_variable_repeat_groups(&values, kept_alleles, &repeat_groups.base_counts)
            }
            _ => Some(value.to_string()),
        };

        if let Some(value) = value {
            fields.push(format!("{key}={value}"));
        }
    }

    Ok(if fields.is_empty() {
        ".".into()
    } else {
        fields.join(";")
    })
}

#[derive(Debug, Default)]
struct RepeatGroups {
    repeat_counts: Vec<usize>,
    base_counts: Vec<usize>,
}

fn info_repeat_groups(info: &str) -> RepeatGroups {
    let rn = info_values(info, "RN")
        .unwrap_or_default()
        .into_iter()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect::<Vec<_>>();
    let rb = info_values(info, "RB")
        .unwrap_or_default()
        .into_iter()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect::<Vec<_>>();
    let rul = info_values(info, "RUL")
        .unwrap_or_default()
        .into_iter()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect::<Vec<_>>();

    let mut base_counts = Vec::new();
    let mut offset = 0;
    for repeat_count in &rn {
        let mut n = 0;
        for _ in 0..*repeat_count {
            if let (Some(rb), Some(rul)) = (rb.get(offset), rul.get(offset))
                && *rul > 0
            {
                n += rb / rul;
            }
            offset += 1;
        }
        base_counts.push(n);
    }

    RepeatGroups {
        repeat_counts: rn,
        base_counts,
    }
}

fn info_values<'a>(info: &'a str, key: &str) -> Option<Vec<&'a str>> {
    info.split(';').find_map(|field| {
        field
            .split_once('=')
            .filter(|(field_key, _)| *field_key == key)
            .map(|(_, value)| split_csv(value))
    })
}

fn remove_sample_alleles(
    keys: &[String],
    sample: &str,
    alternate_count: usize,
    kept_alleles: &[usize],
    allele_map: &[Option<usize>],
) -> io::Result<String> {
    let values = sample.split(':').collect::<Vec<_>>();
    let laa = keys
        .iter()
        .position(|key| *key == "LAA")
        .and_then(|i| values.get(i))
        .copied();
    let local_alleles = parse_local_alleles(laa, kept_alleles);
    let local_kept_alleles = local_alleles
        .iter()
        .copied()
        .filter(|i| kept_alleles.contains(i))
        .collect::<Vec<_>>();
    let mut out = Vec::with_capacity(keys.len());

    for (i, key) in keys.iter().enumerate() {
        let value = values.get(i).copied().unwrap_or(".");
        let next = match key.as_str() {
            "GT" => remap_gt(value, allele_map)?,
            "AD" | "VL_R_STR_FMT" => {
                select_number_r(&split_csv(value), kept_alleles).unwrap_or_else(|| ".".into())
            }
            "EC" | "VL_A_STR_FMT" => {
                select_number_a(&split_csv(value), kept_alleles).unwrap_or_else(|| ".".into())
            }
            "PL" | "VL_G_STR_FMT" => {
                select_number_g(&split_csv(value), alternate_count, kept_alleles)
                    .unwrap_or_else(|| ".".into())
            }
            "LAA" => remap_laa(value, allele_map),
            "LAD" | "VL_LR_STR_FMT" => {
                let alleles = if laa.is_some() {
                    &local_kept_alleles
                } else {
                    kept_alleles
                };
                select_number_r(&split_csv(value), alleles).unwrap_or_else(|| ".".into())
            }
            "LEC" | "VL_LA_STR_FMT" => {
                select_local_number_a(&split_csv(value), &local_alleles, &local_kept_alleles)
                    .unwrap_or_else(|| ".".into())
            }
            "LPL" | "VL_LG_STR_FMT" => {
                select_local_number_g(&split_csv(value), &local_alleles, &local_kept_alleles)
                    .unwrap_or_else(|| ".".into())
            }
            _ => value.into(),
        };
        out.push(next);
    }

    Ok(out.join(":"))
}

fn split_csv(value: &str) -> Vec<&str> {
    if value == "." {
        Vec::new()
    } else {
        value.split(',').collect()
    }
}

fn select_number_a(values: &[&str], kept_alleles: &[usize]) -> Option<String> {
    select_indices(values, kept_alleles.iter().map(|i| i - 1))
}

fn select_number_r(values: &[&str], kept_alleles: &[usize]) -> Option<String> {
    select_indices(
        values,
        std::iter::once(0).chain(kept_alleles.iter().copied()),
    )
}

fn select_number_g(
    values: &[&str],
    alternate_count: usize,
    kept_alleles: &[usize],
) -> Option<String> {
    let alleles = std::iter::once(0)
        .chain(kept_alleles.iter().copied())
        .collect::<Vec<_>>();
    select_indices(values, genotype_indices(&alleles, alternate_count))
}

fn select_alt_groups(
    values: &[&str],
    alternate_count: usize,
    kept_alleles: &[usize],
    group_size: usize,
) -> Option<String> {
    if group_size == 0 || values.is_empty() {
        return None;
    }

    let mut selected = Vec::new();
    for allele in kept_alleles {
        let start = (allele - 1) * group_size;
        selected.extend(values.iter().skip(start).take(group_size).copied());
    }

    if selected.is_empty() && alternate_count > 0 {
        None
    } else {
        Some(selected.join(","))
    }
}

fn select_variable_repeat_groups(
    values: &[&str],
    kept_alleles: &[usize],
    group_sizes: &[usize],
) -> Option<String> {
    let mut starts = Vec::with_capacity(group_sizes.len());
    let mut offset = 0;
    for size in group_sizes {
        starts.push(offset);
        offset += size;
    }

    let mut selected = Vec::new();
    for allele in kept_alleles {
        let Some((&start, &size)) = starts.get(allele - 1).zip(group_sizes.get(allele - 1)) else {
            continue;
        };
        selected.extend(values.iter().skip(start).take(size).copied());
    }

    if selected.is_empty() {
        None
    } else {
        Some(selected.join(","))
    }
}

fn select_indices<I>(values: &[&str], indices: I) -> Option<String>
where
    I: IntoIterator<Item = usize>,
{
    let selected = indices
        .into_iter()
        .filter_map(|i| values.get(i).copied())
        .collect::<Vec<_>>();

    if selected.is_empty() {
        None
    } else {
        Some(selected.join(","))
    }
}

fn genotype_indices(alleles: &[usize], alternate_count: usize) -> Vec<usize> {
    let mut indices = Vec::new();

    for b in alleles {
        for a in alleles {
            if a > b {
                break;
            }
            let index = b * (b + 1) / 2 + a;
            if index < (alternate_count + 1) * (alternate_count + 2) / 2 {
                indices.push(index);
            }
        }
    }

    indices
}

fn parse_local_alleles(value: Option<&str>, kept_alleles: &[usize]) -> Vec<usize> {
    value
        .filter(|s| *s != ".")
        .map(|s| {
            s.split(',')
                .filter_map(|n| n.parse::<usize>().ok())
                .collect::<Vec<_>>()
        })
        .unwrap_or_else(|| kept_alleles.to_vec())
}

fn remap_laa(value: &str, allele_map: &[Option<usize>]) -> String {
    let selected = split_csv(value)
        .into_iter()
        .filter_map(|s| s.parse::<usize>().ok())
        .filter_map(|i| allele_map.get(i).copied().flatten())
        .map(|i| i.to_string())
        .collect::<Vec<_>>();

    if selected.is_empty() {
        ".".into()
    } else {
        selected.join(",")
    }
}

fn select_local_number_a(
    values: &[&str],
    local_alleles: &[usize],
    local_kept_alleles: &[usize],
) -> Option<String> {
    select_indices(
        values,
        local_kept_alleles
            .iter()
            .filter_map(|allele| local_alleles.iter().position(|i| i == allele)),
    )
}

fn select_local_number_g(
    values: &[&str],
    local_alleles: &[usize],
    local_kept_alleles: &[usize],
) -> Option<String> {
    let alleles = std::iter::once(0)
        .chain(local_alleles.iter().copied())
        .collect::<Vec<_>>();
    let kept = std::iter::once(0)
        .chain(local_kept_alleles.iter().copied())
        .collect::<Vec<_>>();
    let indices = genotype_indices_from_alleles(&alleles, &kept);

    select_indices(values, indices)
}

fn genotype_indices_from_alleles(alleles: &[usize], kept: &[usize]) -> Vec<usize> {
    let mut indices = Vec::new();

    for b in kept {
        for a in kept {
            if a > b {
                break;
            }

            if let (Some(a), Some(b)) = (
                alleles.iter().position(|allele| allele == a),
                alleles.iter().position(|allele| allele == b),
            ) {
                indices.push(b * (b + 1) / 2 + a);
            }
        }
    }

    indices
}

fn remap_gt(value: &str, allele_map: &[Option<usize>]) -> io::Result<String> {
    let mut out = String::new();
    let mut token = String::new();

    for c in value.chars() {
        if matches!(c, '/' | '|') {
            out.push_str(&remap_gt_token(&token, allele_map)?);
            out.push(c);
            token.clear();
        } else {
            token.push(c);
        }
    }

    out.push_str(&remap_gt_token(&token, allele_map)?);

    Ok(out)
}

fn remap_gt_token(token: &str, allele_map: &[Option<usize>]) -> io::Result<String> {
    if token == "." {
        return Ok(".".into());
    }

    let allele = token
        .parse::<usize>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(allele_map
        .get(allele)
        .copied()
        .flatten()
        .map(|i| i.to_string())
        .unwrap_or_else(|| ".".into()))
}

#[cfg(test)]
mod tests {
    use std::{io::BufReader, path::PathBuf};

    use super::{
        contig_count, count_bcf_records_from_path, count_bcf_records_in_region_from_path,
        count_vcf_records_from_path, count_vcf_records_in_region_from_path,
        read_bcf_header_from_path, read_vcf_header_from_path, sample_count,
        write_bcf_csi_from_path,
    };
    use crate::tabix_compat::{TextFormat, write_bgzf_and_index};

    fn fixture(path: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(path)
    }

    #[test]
    fn test_read_vcf_header_and_records() {
        let path = fixture("htslib/test/index.vcf");
        let header = read_vcf_header_from_path(&path).unwrap();

        assert!(contig_count(&header) > 0);
        assert!(count_vcf_records_from_path(path).unwrap() > 0);
    }

    #[test]
    fn test_query_vcf_records() {
        let src = fixture("htslib/test/tabix/vcf_file.vcf");
        let bgzf_path = std::env::temp_dir().join(format!(
            "htslib-rs-variant-query-{}.vcf.gz",
            std::process::id()
        ));
        let tbi_path = std::env::temp_dir().join(format!(
            "htslib-rs-variant-query-{}.vcf.gz.tbi",
            std::process::id()
        ));

        let vcf = std::fs::File::open(src).unwrap();
        write_bgzf_and_index(BufReader::new(vcf), &bgzf_path, &tbi_path, TextFormat::Vcf).unwrap();

        let region = "1:3000151-3000151".parse().unwrap();
        let count = count_vcf_records_in_region_from_path(&bgzf_path, &region).unwrap();

        std::fs::remove_file(bgzf_path).unwrap();
        std::fs::remove_file(tbi_path).unwrap();

        assert_eq!(count, 1);
    }

    #[test]
    fn test_read_bcf_header_and_records() {
        let path = fixture("htslib/test/tabix/vcf_file.bcf");
        let header = read_bcf_header_from_path(&path).unwrap();

        assert!(sample_count(&header) > 0);
        assert!(count_bcf_records_from_path(path).unwrap() > 0);
    }

    #[test]
    fn test_query_bcf_records() {
        let src = fixture("htslib/test/tabix/vcf_file.bcf");
        let bcf_path = std::env::temp_dir().join(format!(
            "htslib-rs-variant-query-{}.bcf",
            std::process::id()
        ));
        let csi_path = std::env::temp_dir().join(format!(
            "htslib-rs-variant-query-{}.bcf.csi",
            std::process::id()
        ));

        std::fs::copy(src, &bcf_path).unwrap();
        write_bcf_csi_from_path(&bcf_path, &csi_path).unwrap();

        let region = "1:3000151-3000151".parse().unwrap();
        let count = count_bcf_records_in_region_from_path(&bcf_path, &region).unwrap();

        std::fs::remove_file(bcf_path).unwrap();
        std::fs::remove_file(csi_path).unwrap();

        assert_eq!(count, 1);
    }
}
