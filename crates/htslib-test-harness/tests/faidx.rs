use std::io::{BufReader, Cursor};

use htslib_rs::{
    bgzf,
    bgzf_compat::{build_gzi, write_all as write_bgzf_all},
    faidx_compat::{
        build_fastq_index, build_index, fastq_line_length, fastq_sequence_len, fastq_sequence_name,
        has_fastq_sequence, has_sequence, line_length, read_index, sequence_len, sequence_name,
        write_fasta_retrieval, write_fastq_as_fasta_retrieval, write_fastq_index,
        write_fastq_retrieval, write_index,
    },
};

const FASTA: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/faidx.fa");
const FASTA_FAI: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/faidx.fa.expected.fai");
const FASTA_EXPECTED: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/faidx.1.expected.fa");
const CE_FASTA: &[u8] = include_bytes!("../../../repos/htslib/test/ce.fa");
const CE_FASTA_FAI: &[u8] = include_bytes!("../../../repos/htslib/test/ce.fa.fai");
const CE_FASTA_EXPECTED: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/ce.1.expected.fa");
const FASTQ: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/fastqs.fq");
const FASTQ_FAI: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/fastqs.fq.expected.fai");
const FASTQ_EXPECTED: &[u8] = include_bytes!("../../../repos/htslib/test/faidx/fastqs.1.expected.fq");
const FASTQ_AS_FASTA_EXPECTED: &[u8] =
    include_bytes!("../../../repos/htslib/test/faidx/fastqs.2.expected.fa");

#[test]
fn builds_fasta_fai_matching_htslib_golden() -> Result<(), Box<dyn std::error::Error>> {
    let index = build_index(BufReader::new(Cursor::new(FASTA)))?;
    let mut actual = Vec::new();
    write_index(&mut actual, &index)?;

    assert_eq!(actual, FASTA_FAI);
    assert!(has_sequence(&index, "foo"));
    assert!(!has_sequence(&index, "absent"));
    assert_eq!(sequence_name(&index, 3), Some(&b"trailingblank3"[..]));
    assert_eq!(sequence_len(&index, "trailingblank1"), Some(33));
    assert_eq!(line_length(&index, "trailingblank2"), Some(24));

    Ok(())
}

#[test]
fn builds_fastq_fai_matching_htslib_golden() -> Result<(), Box<dyn std::error::Error>> {
    let index = build_fastq_index(BufReader::new(Cursor::new(FASTQ)))?;
    let mut actual = Vec::new();
    write_fastq_index(&mut actual, &index)?;

    assert_eq!(actual, FASTQ_FAI);
    assert!(has_fastq_sequence(&index, "SRR014849.203935_3"));
    assert!(!has_fastq_sequence(&index, "absent"));
    assert_eq!(fastq_sequence_name(&index, 0), Some("FAKE0005_1"));
    assert_eq!(fastq_sequence_len(&index, "FSRRS4401CM938_1"), Some(453));
    assert_eq!(fastq_line_length(&index, "SRR014849.203935_3"), Some(144));

    Ok(())
}

#[test]
fn retrieves_fasta_regions_matching_htslib_golden() -> Result<(), Box<dyn std::error::Error>> {
    let index = build_index(BufReader::new(Cursor::new(FASTA)))?;
    let mut reader = Cursor::new(FASTA);
    let mut actual = Vec::new();
    let regions = ["trailingblank2:28-33", "trailingblank3:4-5", "bar:4-5"];

    write_fasta_retrieval(&mut reader, &index, &regions, &mut actual)?;

    assert_eq!(actual, FASTA_EXPECTED);

    Ok(())
}

#[test]
fn retrieves_bgzf_fasta_regions_matching_htslib_golden() -> Result<(), Box<dyn std::error::Error>> {
    let compressed = write_bgzf_all(Vec::new(), CE_FASTA)?;
    let gzi = build_gzi(&mut Cursor::new(&compressed))?;
    let index = read_index(BufReader::new(Cursor::new(CE_FASTA_FAI)))?;
    let mut reader = bgzf::io::IndexedReader::new(Cursor::new(compressed), gzi);
    let mut actual = Vec::new();
    let regions = ["CHROMOSOME_I:5001-5125", "CHROMOSOME_X:101-225"];

    write_fasta_retrieval(&mut reader, &index, &regions, &mut actual)?;

    assert_eq!(actual, CE_FASTA_EXPECTED);

    Ok(())
}

#[test]
fn retrieves_fastq_regions_matching_htslib_golden() -> Result<(), Box<dyn std::error::Error>> {
    let index = build_fastq_index(BufReader::new(Cursor::new(FASTQ)))?;
    let mut reader = Cursor::new(FASTQ);
    let mut actual = Vec::new();
    let regions = [
        "FAKE0006_1:4-12",
        "FSRRS4401BE7HA_1:81-120",
        "FAKE0010_2",
        "SRR014849.50939_3:71-90",
    ];

    write_fastq_retrieval(&mut reader, &index, &regions, &mut actual)?;

    assert_eq!(actual, FASTQ_EXPECTED);

    Ok(())
}

#[test]
fn retrieves_fastq_regions_as_fasta_matching_htslib_golden()
-> Result<(), Box<dyn std::error::Error>> {
    let index = build_fastq_index(BufReader::new(Cursor::new(FASTQ)))?;
    let mut reader = Cursor::new(FASTQ);
    let mut actual = Vec::new();
    let regions = [
        "FAKE0006_1:4-12",
        "FSRRS4401BE7HA_1:81-120",
        "FAKE0010_2",
        "SRR014849.50939_3:71-90",
    ];

    write_fastq_as_fasta_retrieval(&mut reader, &index, &regions, &mut actual)?;

    assert_eq!(actual, FASTQ_AS_FASTA_EXPECTED);

    Ok(())
}
