use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

use htslib_rs::{
    bgzf_compat::write_all as write_bgzf_all,
    csi::BinningIndex,
    index_compat::{
        bai_reference_sequence_count, build_bai, build_bam_csi_with_min_shift, build_bcf_csi,
        build_bcf_csi_with_min_shift, build_cram_crai, build_sam_bai, build_sam_csi,
        build_sam_csi_with_min_shift, build_vcf_csi, build_vcf_csi_with_min_shift, build_vcf_tbi,
        crai_record_count, csi_reference_sequence_count, read_bai, read_cram_crai, read_csi,
        read_tbi, tbi_reference_sequence_count, write_bai, write_cram_crai, write_csi, write_tbi,
    },
    tabix_compat::{query_csi_records_from_path, query_records_from_path},
    variant_io_compat::query_bcf_records_from_path,
};

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib/test")
        .join(path)
}

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "htslib-rs-index-build-{}-{name}",
        std::process::id()
    ))
}

fn cleanup(paths: &[PathBuf]) {
    for path in paths {
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn builds_csi_for_very_large_reference_and_queries_it() -> Result<(), Box<dyn std::error::Error>> {
    use htslib_rs::bam;
    use htslib_rs::csi::BinningIndex;
    use htslib_rs::sam;
    use htslib_rs::sam::alignment::io::Write as _;

    // ref2 length 541556283 > 2^29 (BAI's limit). The BAM-CSI builder must
    // auto-size depth from the header so the index can address it without
    // a "index out of bounds" panic.
    let sam_text = "@HD\tVN:1.6\tSO:coordinate\n\
                    @SQ\tSN:ref2\tLN:541556283\n\
                    r\t0\tref2\t536880911\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n";
    let mut sam_reader = sam::io::Reader::new(std::io::Cursor::new(sam_text.as_bytes()));
    let header = sam_reader.read_header()?;

    let bam_path = temp_path("large-ref.bam");
    let mut writer = bam::io::Writer::new(File::create(&bam_path)?);
    writer.write_header(&header)?;
    for result in sam_reader.records() {
        writer.write_alignment_record(&header, &result?)?;
    }
    drop(writer);

    let index = build_bam_csi_with_min_shift(&bam_path, 14)?;
    cleanup(&[bam_path]);

    // depth must exceed the old fixed 5 to cover > 2^29.
    assert!(index.depth() > 5);
    assert_eq!(csi_reference_sequence_count(&index), 1);
    Ok(())
}

#[test]
fn builds_bai_for_bam_without_so_coordinate_header() -> Result<(), Box<dyn std::error::Error>> {
    use htslib_rs::bam;
    use htslib_rs::sam;
    use htslib_rs::sam::alignment::io::Write as _;

    // Coordinate-ordered records, but the @HD line carries no SO tag —
    // upstream `samtools index` indexes such BAMs anyway.
    let sam_text = "@HD\tVN:1.6\n\
                    @SQ\tSN:ref0\tLN:1000\n\
                    a\t0\tref0\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
                    b\t0\tref0\t40\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n";
    let mut sam_reader = sam::io::Reader::new(std::io::Cursor::new(sam_text.as_bytes()));
    let header = sam_reader.read_header()?;

    let bam_path = temp_path("no-so.bam");
    let mut writer = bam::io::Writer::new(File::create(&bam_path)?);
    writer.write_header(&header)?;
    for result in sam_reader.records() {
        writer.write_alignment_record(&header, &result?)?;
    }
    drop(writer);

    let index = build_bai(&bam_path)?;
    cleanup(&[bam_path]);

    assert_eq!(bai_reference_sequence_count(&index), 1);
    Ok(())
}

#[test]
fn builds_bai_for_bam_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let expected = read_bai(fixture("range.bam.bai"))?;
    let actual = build_bai(fixture("range.bam"))?;

    assert_eq!(
        bai_reference_sequence_count(&actual),
        bai_reference_sequence_count(&expected)
    );

    let bai_path = temp_path("range.bam.bai");
    write_bai(&bai_path, &actual)?;
    let reread = read_bai(&bai_path)?;
    cleanup(&[bai_path]);

    assert_eq!(
        bai_reference_sequence_count(&reread),
        bai_reference_sequence_count(&actual)
    );

    Ok(())
}

#[test]
fn reads_and_writes_crai_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let index = read_cram_crai(fixture("range.cram.crai"))?;

    let crai_path = temp_path("range.cram.crai");
    write_cram_crai(&crai_path, &index)?;
    let reread = read_cram_crai(&crai_path)?;
    cleanup(&[crai_path]);

    assert_eq!(crai_record_count(&reread), crai_record_count(&index));

    Ok(())
}

#[test]
fn builds_crai_for_cram_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let actual = build_cram_crai(fixture("tlen/a4.cram"))?;

    assert!(crai_record_count(&actual) > 0);

    let crai_path = temp_path("built-a4.cram.crai");
    write_cram_crai(&crai_path, &actual)?;
    let reread = read_cram_crai(&crai_path)?;
    cleanup(&[crai_path]);

    assert_eq!(crai_record_count(&reread), crai_record_count(&actual));

    Ok(())
}

#[test]
fn builds_bai_for_bgzf_sam_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let expected = read_bai(fixture("index.sam.gz.bai"))?;
    let sam_gz_path = temp_path("bai-index.sam.gz");
    let bai_path = temp_path("bai-index.sam.gz.bai");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    File::create(&sam_gz_path)?.write_all(&compressed)?;

    let actual = build_sam_bai(&sam_gz_path)?;
    assert_eq!(
        bai_reference_sequence_count(&actual),
        bai_reference_sequence_count(&expected)
    );

    write_bai(&bai_path, &actual)?;
    let reread = read_bai(&bai_path)?;
    cleanup(&[sam_gz_path, bai_path]);

    assert_eq!(
        bai_reference_sequence_count(&reread),
        bai_reference_sequence_count(&actual)
    );

    Ok(())
}

#[test]
fn builds_csi_for_bam_with_explicit_min_shift() -> Result<(), Box<dyn std::error::Error>> {
    let index = build_bam_csi_with_min_shift(fixture("range.bam"), 10)?;

    assert_eq!(index.min_shift(), 10);
    assert!(csi_reference_sequence_count(&index) > 0);

    let csi_path = temp_path("range.bam.csi");
    write_csi(&csi_path, &index)?;
    let reread = read_csi(&csi_path)?;
    cleanup(&[csi_path]);

    assert_eq!(reread.min_shift(), 10);
    assert_eq!(
        csi_reference_sequence_count(&reread),
        csi_reference_sequence_count(&index)
    );

    Ok(())
}

#[test]
fn builds_csi_for_bgzf_sam_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let expected = read_csi(fixture("index.sam.gz.csi"))?;
    let sam_gz_path = temp_path("index.sam.gz");
    let csi_path = temp_path("index.sam.gz.csi");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    File::create(&sam_gz_path)?.write_all(&compressed)?;

    let actual = build_sam_csi(&sam_gz_path)?;
    assert_eq!(
        csi_reference_sequence_count(&actual),
        csi_reference_sequence_count(&expected)
    );

    write_csi(&csi_path, &actual)?;
    let reread = read_csi(&csi_path)?;
    cleanup(&[sam_gz_path, csi_path]);

    assert_eq!(
        csi_reference_sequence_count(&reread),
        csi_reference_sequence_count(&actual)
    );

    Ok(())
}

#[test]
fn builds_csi_for_bgzf_sam_with_explicit_min_shift() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("sam-min-shift-index.sam.gz");
    let csi_path = temp_path("sam-min-shift-index.sam.gz.csi");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    File::create(&sam_gz_path)?.write_all(&compressed)?;

    let index = build_sam_csi_with_min_shift(&sam_gz_path, 10)?;
    assert_eq!(index.min_shift(), 10);
    assert!(csi_reference_sequence_count(&index) > 0);

    write_csi(&csi_path, &index)?;
    let reread = read_csi(&csi_path)?;
    cleanup(&[sam_gz_path, csi_path]);

    assert_eq!(reread.min_shift(), 10);
    assert_eq!(
        csi_reference_sequence_count(&reread),
        csi_reference_sequence_count(&index)
    );

    Ok(())
}

#[test]
fn builds_tbi_for_vcf_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let vcf = File::open(fixture("index.vcf"))?;
    let (compressed, index) = build_vcf_tbi(BufReader::new(vcf), Vec::new())?;

    assert!(tbi_reference_sequence_count(&index) > 0);

    let bgzf_path = temp_path("tbi-index.vcf.gz");
    let tbi_path = temp_path("tbi-index.vcf.gz.tbi");
    File::create(&bgzf_path)?.write_all(&compressed)?;
    write_tbi(&tbi_path, &index)?;

    let reread = read_tbi(&tbi_path)?;
    let records = query_records_from_path(&bgzf_path, reread, &"1:9999919-9999919".parse()?)?;
    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(records.len(), 1);
    assert!(records[0].starts_with("1\t9999919\t"));

    Ok(())
}

#[test]
fn builds_csi_for_vcf_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let vcf = File::open(fixture("index.vcf"))?;
    let (compressed, index) = build_vcf_csi(BufReader::new(vcf), Vec::new())?;

    assert!(csi_reference_sequence_count(&index) > 0);

    let bgzf_path = temp_path("csi-index.vcf.gz");
    let csi_path = temp_path("csi-index.vcf.gz.csi");
    File::create(&bgzf_path)?.write_all(&compressed)?;
    write_csi(&csi_path, &index)?;

    let reread = read_csi(&csi_path)?;
    let records = query_csi_records_from_path(&bgzf_path, reread, &"1:9999919-9999919".parse()?)?;
    cleanup(&[bgzf_path, csi_path]);

    assert_eq!(records.len(), 1);
    assert!(records[0].starts_with("1\t9999919\t"));

    Ok(())
}

#[test]
fn builds_csi_for_vcf_with_explicit_min_shift() -> Result<(), Box<dyn std::error::Error>> {
    let vcf = File::open(fixture("index.vcf"))?;
    let (compressed, index) = build_vcf_csi_with_min_shift(BufReader::new(vcf), Vec::new(), 10)?;

    assert_eq!(index.min_shift(), 10);
    assert!(csi_reference_sequence_count(&index) > 0);

    let bgzf_path = temp_path("csi-min-shift-index.vcf.gz");
    let csi_path = temp_path("csi-min-shift-index.vcf.gz.csi");
    File::create(&bgzf_path)?.write_all(&compressed)?;
    write_csi(&csi_path, &index)?;

    let reread = read_csi(&csi_path)?;
    let records = query_csi_records_from_path(&bgzf_path, reread, &"1:9999919-9999919".parse()?)?;
    cleanup(&[bgzf_path, csi_path]);

    assert_eq!(records.len(), 1);
    assert!(records[0].starts_with("1\t9999919\t"));

    Ok(())
}

#[test]
fn builds_csi_for_bcf_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("index.bcf");
    let csi_path = temp_path("index.bcf.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    let index = build_bcf_csi(&bcf_path)?;

    assert!(csi_reference_sequence_count(&index) > 0);

    write_csi(&csi_path, &index)?;
    let reread = read_csi(&csi_path)?;

    assert_eq!(
        csi_reference_sequence_count(&reread),
        csi_reference_sequence_count(&index)
    );

    let records = query_bcf_records_from_path(&bcf_path, &"1:3000151-3000151".parse()?)?;
    cleanup(&[bcf_path, csi_path]);

    assert_eq!(records.len(), 1);

    Ok(())
}

#[test]
fn builds_csi_for_bcf_with_explicit_min_shift() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("min-shift-index.bcf");
    let csi_path = temp_path("min-shift-index.bcf.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    let index = build_bcf_csi_with_min_shift(&bcf_path, 10)?;

    assert_eq!(index.min_shift(), 10);
    assert!(csi_reference_sequence_count(&index) > 0);

    write_csi(&csi_path, &index)?;
    let reread = read_csi(&csi_path)?;
    cleanup(&[bcf_path, csi_path]);

    assert_eq!(reread.min_shift(), 10);
    assert_eq!(
        csi_reference_sequence_count(&reread),
        csi_reference_sequence_count(&index)
    );

    Ok(())
}
