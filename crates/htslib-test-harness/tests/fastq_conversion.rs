use std::{
    io::{BufReader, Cursor},
    path::Path,
};

use htslib_rs::fastq_compat::{
    FastxToSamOptions, SamToFastxOptions, write_fasta_from_sam, write_fastq_from_sam,
    write_sam_from_fasta, write_sam_from_fastq,
};

fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../repos/htslib/test/fastq")
        .join(name);

    std::fs::read(&path).unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
}

fn assert_fastq_to_sam(
    src: &str,
    expected: &str,
    options: FastxToSamOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut actual = Vec::new();
    write_sam_from_fastq(
        BufReader::new(Cursor::new(fixture(src))),
        &mut actual,
        &options,
    )?;

    assert_eq!(actual, fixture(expected));

    Ok(())
}

fn assert_fasta_to_sam(
    src: &str,
    expected: &str,
    options: FastxToSamOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut actual = Vec::new();
    write_sam_from_fasta(
        BufReader::new(Cursor::new(fixture(src))),
        &mut actual,
        &options,
    )?;

    assert_eq!(actual, fixture(expected));

    Ok(())
}

fn assert_sam_to_fastq(
    src: &str,
    expected: &str,
    options: SamToFastxOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut actual = Vec::new();
    write_fastq_from_sam(
        BufReader::new(Cursor::new(fixture(src))),
        &mut actual,
        &options,
    )?;

    assert_eq!(actual, fixture(expected));

    Ok(())
}

fn assert_sam_to_fasta(
    src: &str,
    expected: &str,
    options: SamToFastxOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut actual = Vec::new();
    write_fasta_from_sam(
        BufReader::new(Cursor::new(fixture(src))),
        &mut actual,
        &options,
    )?;

    assert_eq!(actual, fixture(expected));

    Ok(())
}

#[test]
fn reads_fastq_and_fasta_as_sam() -> Result<(), Box<dyn std::error::Error>> {
    assert_fastq_to_sam("minimal.fq", "minimal.sam", FastxToSamOptions::default())?;
    assert_fasta_to_sam("minimal.fa", "minimal-q.sam", FastxToSamOptions::default())?;
    assert_fastq_to_sam(
        "multiline.fq",
        "multiline.sam",
        FastxToSamOptions::default(),
    )?;
    assert_fasta_to_sam(
        "multiline.fa",
        "multiline-q.sam",
        FastxToSamOptions::default(),
    )?;
    assert_fastq_to_sam(
        "single.fq",
        "single_noaux.sam",
        FastxToSamOptions::default(),
    )?;
    assert_fasta_to_sam(
        "single.fa",
        "single_noaux-q.sam",
        FastxToSamOptions::default(),
    )?;

    Ok(())
}

#[test]
fn reads_aux_and_interleaved_fastq_as_sam() -> Result<(), Box<dyn std::error::Error>> {
    let with_aux = FastxToSamOptions {
        include_aux: true,
        ..Default::default()
    };

    assert_fastq_to_sam("single.fq", "single_aux.sam", with_aux.clone())?;
    assert_fasta_to_sam("single.fa", "single_aux-q.sam", with_aux.clone())?;
    assert_fastq_to_sam("longline.fq", "longline.sam", with_aux.clone())?;
    assert_fastq_to_sam("interleaved.fq", "inter_aux.sam", with_aux.clone())?;
    assert_fasta_to_sam("interleaved.fa", "inter_aux-q.sam", with_aux)?;
    assert_fastq_to_sam(
        "interleaved.fq",
        "inter_noaux.sam",
        FastxToSamOptions::default(),
    )?;
    assert_fasta_to_sam(
        "interleaved.fa",
        "inter_noaux-q.sam",
        FastxToSamOptions::default(),
    )?;

    Ok(())
}

#[test]
fn reads_name2_and_umi_fastq_as_sam() -> Result<(), Box<dyn std::error::Error>> {
    let name2 = FastxToSamOptions {
        name2: true,
        ..Default::default()
    };
    let umi = FastxToSamOptions {
        umi_tag: Some("RX".into()),
        ..Default::default()
    };

    assert_fastq_to_sam("name2.fq", "name2.sam", name2.clone())?;
    assert_fasta_to_sam("name2.fa", "name2-q.sam", name2)?;
    assert_fastq_to_sam("UMI.fq", "UMI.sam", umi)?;

    Ok(())
}

#[test]
fn reads_casava_fastq_and_fasta_as_sam() -> Result<(), Box<dyn std::error::Error>> {
    let casava = FastxToSamOptions {
        casava: true,
        ..Default::default()
    };
    let casava_ox = FastxToSamOptions {
        casava: true,
        barcode_tag: Some("OX".into()),
        ..Default::default()
    };

    assert_fastq_to_sam("interleaved_casava.fq", "inter_casava.sam", casava.clone())?;
    assert_fastq_to_sam(
        "interleaved_casava.fq",
        "inter_casavaOX.sam",
        casava_ox.clone(),
    )?;
    assert_fasta_to_sam(
        "interleaved_casava.fa",
        "inter_casava-q.sam",
        casava.clone(),
    )?;
    assert_fasta_to_sam("interleaved_casava.fa", "inter_casavaOX-q.sam", casava_ox)?;
    assert_fastq_to_sam("filter_casava.fq", "filter_casava.sam", casava.clone())?;
    assert_fasta_to_sam("filter_casava.fa", "filter_casava-q.sam", casava)?;

    Ok(())
}

#[test]
fn writes_sam_as_fastq_and_fasta() -> Result<(), Box<dyn std::error::Error>> {
    assert_sam_to_fastq("minimal.sam", "minimal.fq", SamToFastxOptions::default())?;
    assert_sam_to_fasta("minimal.sam", "minimal.fa", SamToFastxOptions::default())?;

    let with_aux = SamToFastxOptions {
        include_aux: true,
        ..Default::default()
    };
    let with_aux_and_read_number = SamToFastxOptions {
        include_aux: true,
        append_read_number: true,
        ..Default::default()
    };

    assert_sam_to_fastq("single_aux.sam", "single.fq", with_aux.clone())?;
    assert_sam_to_fasta("single_aux.sam", "single.fa", with_aux)?;
    assert_sam_to_fastq(
        "inter_aux.sam",
        "interleaved.fq",
        with_aux_and_read_number.clone(),
    )?;
    assert_sam_to_fasta("inter_aux.sam", "interleaved.fa", with_aux_and_read_number)?;

    Ok(())
}

#[test]
fn writes_sam_as_casava_fastq_and_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let casava = SamToFastxOptions {
        casava: true,
        ..Default::default()
    };
    let casava_ox = SamToFastxOptions {
        casava: true,
        barcode_tag: Some("OX".into()),
        ..Default::default()
    };

    assert_sam_to_fastq("inter_casava.sam", "interleaved_casava.fq", casava.clone())?;
    assert_sam_to_fastq(
        "inter_casavaOX.sam",
        "interleaved_casava.fq",
        casava_ox.clone(),
    )?;
    assert_sam_to_fasta("inter_casava.sam", "interleaved_casava.fa", casava.clone())?;
    assert_sam_to_fasta("inter_casavaOX.sam", "interleaved_casava.fa", casava_ox)?;
    assert_sam_to_fastq("filter_casava.sam", "filter_casava.fq", casava.clone())?;
    assert_sam_to_fasta("filter_casava.sam", "filter_casava.fa", casava)?;

    Ok(())
}

#[test]
fn reads_and_writes_separate_paired_fastq_and_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let with_aux = FastxToSamOptions {
        include_aux: true,
        ..Default::default()
    };
    let with_aux_and_read_number = SamToFastxOptions {
        include_aux: true,
        append_read_number: true,
        ..Default::default()
    };

    assert_fastq_to_sam("r1.fq", "r1.sam", with_aux.clone())?;
    assert_fastq_to_sam("r2.fq", "r2.sam", with_aux.clone())?;
    assert_fasta_to_sam("r1.fa", "r1-q.sam", with_aux.clone())?;
    assert_fasta_to_sam("r2.fa", "r2-q.sam", with_aux)?;

    assert_sam_to_fastq("r1.sam", "r1.fq", with_aux_and_read_number.clone())?;
    assert_sam_to_fastq("r2.sam", "r2.fq", with_aux_and_read_number.clone())?;
    assert_sam_to_fasta("r1.sam", "r1.fa", with_aux_and_read_number.clone())?;
    assert_sam_to_fasta("r2.sam", "r2.fa", with_aux_and_read_number)?;

    Ok(())
}

#[test]
fn writes_sam_as_umi_fastq() -> Result<(), Box<dyn std::error::Error>> {
    let options = SamToFastxOptions {
        append_read_number: true,
        umi_tag: Some("RX".into()),
        ..Default::default()
    };

    assert_sam_to_fastq("UMI.sam", "UMI.fq", options)?;

    Ok(())
}
