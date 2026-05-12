use std::{
    error::Error as _,
    path::{Path, PathBuf},
};

use htslib_rs::error::ContextError;
use htslib_rs::format::{Category, Compression, Exact, Format};

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib/test")
        .join(path)
}

fn missing_fixture(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("htslib-rs-missing-{name}-{}", std::process::id()))
}

#[test]
fn detects_htslib_alignment_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        htslib_rs::io::detect_format(fixture("index.sam"))?,
        Format::new(Category::SequenceData, Exact::Sam, Compression::None)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("colons.bam"))?,
        Format::new(Category::SequenceData, Exact::Bam, Compression::Bgzf)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("ce#5b_java.cram"))?,
        Format::new(Category::SequenceData, Exact::Cram, Compression::None)
    );

    Ok(())
}

#[test]
fn detects_htslib_variant_and_index_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        htslib_rs::io::detect_format(fixture("index.vcf"))?,
        Format::new(Category::VariantData, Exact::Vcf, Compression::None)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("index.bam.bai"))?,
        Format::new(Category::IndexFile, Exact::Bai, Compression::None)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("index.bam.csi"))?,
        Format::new(Category::IndexFile, Exact::Csi, Compression::Bgzf)
    );

    Ok(())
}

#[test]
fn detects_htslib_fasta_fastq_and_bed_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        htslib_rs::io::detect_format(fixture("c1.fa"))?,
        Format::new(Category::SequenceData, Exact::Fasta, Compression::None)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("fastq/minimal.fq"))?,
        Format::new(Category::SequenceData, Exact::Fastq, Compression::None)
    );

    assert_eq!(
        htslib_rs::io::detect_format(fixture("tabix/bed_file.bed"))?,
        Format::new(Category::RegionList, Exact::Bed, Compression::None)
    );

    Ok(())
}

#[test]
fn format_detection_errors_include_operation_and_path_context() {
    let path = missing_fixture("format-detection");
    let err = htslib_rs::io::detect_format(&path).unwrap_err();
    let message = err.to_string();

    assert!(message.contains("open input for format detection"));
    assert!(message.contains(&path.display().to_string()));

    let context = err
        .source()
        .and_then(|e| e.downcast_ref::<ContextError>())
        .unwrap();

    assert_eq!(context.operation(), "open input for format detection");
    assert_eq!(context.path(), Some(path.as_path()));
}

#[test]
fn cram_reference_errors_include_index_path_context() {
    let path = missing_fixture("reference.fa");
    let err =
        htslib_rs::alignment_compat::cram_reference_repository_from_fasta_path(&path).unwrap_err();
    let message = err.to_string();
    let mut expected_index = path.as_os_str().to_os_string();
    expected_index.push(".fai");
    let expected_index = PathBuf::from(expected_index);

    assert!(message.contains("read FASTA index"));
    assert!(message.contains(&expected_index.display().to_string()));

    let context = err
        .get_ref()
        .and_then(|e| e.downcast_ref::<ContextError>())
        .unwrap();

    assert_eq!(context.operation(), "read FASTA index");
    assert_eq!(context.path(), Some(expected_index.as_path()));
}
