use std::path::{Path, PathBuf};

use htslib_rs::alignment_compat::{
    SamAuxValue, SamHeaderAdapter, SamRecordAdapter, add_bam_header_nul_padding,
    apply_existing_baq_from_sam_path, base_modification_pileup_report_from_sam_path,
    base_modification_report_from_sam_path, benchmark_bam_view_from_path,
    benchmark_cram_view_from_path_with_reference, benchmark_sam_view_from_path,
    cap_mapping_qualities_from_sam_path, count_bam_records_from_path, count_sam_records_from_path,
    count_sam_records_matching_filter_from_path, extended_base_modification_report_from_sam_path,
    force_recalculate_baq_from_sam_path, iter_bam_records_from_path,
    iter_cram_records_from_path_with_reference, pileup_report_from_bam_path,
    pileup_report_from_sam_path, query_bam_records_from_path, query_bam_regions_from_path,
    query_cram_records_from_path_with_reference, query_sam_records_from_path,
    read_bam_header_from_path, read_cram_header_from_path, read_sam_header_from_path,
    recalculate_and_apply_baq_from_sam_path, recalculate_baq_from_sam_path,
    recalculate_extended_baq_from_sam_path, reference_sequence_count,
    revert_existing_baq_from_sam_path, sam_aux_get, sam_aux_insert, sam_aux_remove,
    summarize_cram_records_from_path, summarize_cram_records_from_path_with_reference,
    summarize_sam_records_from_path, view_bam_as_sam_text_from_path_with_limit,
    view_bam_regions_as_sam_text_from_path, view_bam_regions_as_sam_text_from_path_with_dedup,
    view_bam_regions_as_sam_text_from_path_with_limit, view_bgzf_sam_text_from_path_with_limit,
    view_cram_as_sam_text_from_path_with_reference_and_limit,
    view_cram_regions_as_sam_text_from_path_with_reference,
    view_cram_regions_as_sam_text_from_path_with_reference_and_limit,
    view_sam_as_fasta_text_from_path_with_limit, view_sam_as_fastq_text_from_path_with_limit,
    view_sam_regions_as_text_from_path, view_sam_text_from_path_with_limit,
    view_sam_text_from_path_with_limit_and_parse_errors, view_sam_text_matching_filter_from_path,
    write_bam_from_path, write_bam_from_sam_path, write_bam_from_sam_path_with_bai,
    write_bam_from_sam_path_with_compression_level, write_bgzf_sam_from_path_with_limit,
    write_cram_from_path_with_reference, write_cram_from_sam_path_with_reference,
    write_cram_from_sam_path_with_reference_and_crai,
};
use htslib_rs::bgzf_compat::write_all as write_bgzf_all;
use htslib_rs::index_compat::{build_bai, build_sam_bai, build_sam_csi, write_bai, write_csi};

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib/test")
        .join(path)
}

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("htslib-rs-alignment-{}-{name}", std::process::id()))
}

fn cleanup(paths: &[PathBuf]) {
    for path in paths {
        let _ = std::fs::remove_file(path);
    }
}

fn strip_sam_header(s: &str) -> String {
    s.lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect()
}

fn sam_text_with_first_records(src: &str, n: usize) -> String {
    let mut expected = src
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    expected.push_str(
        &src.lines()
            .filter(|line| !line.starts_with('@'))
            .take(n)
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );
    expected
}

fn first_sam_record_fields(
    src: impl AsRef<Path>,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(src)?;
    let fields = text
        .lines()
        .find(|line| !line.starts_with('@'))
        .ok_or("missing SAM record")?
        .split('\t')
        .map(String::from)
        .collect();

    Ok(fields)
}

fn sam_record_line_count(src: &str) -> usize {
    src.lines().filter(|line| !line.starts_with('@')).count()
}

#[test]
fn summarizing_cram_without_reference_returns_error() {
    let err = summarize_cram_records_from_path(fixture("range.cram")).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
    assert!(
        err.to_string()
            .contains("missing reference sequence: CHROMOSOME_I")
    );
}

#[test]
fn reads_sam_header_and_records_from_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("xx#minimal.sam");
    let header = read_sam_header_from_path(&path)?;

    assert_eq!(reference_sequence_count(&header), 2);
    assert_eq!(count_sam_records_from_path(path)?, 8);

    Ok(())
}

#[test]
fn adapts_typed_sam_auxiliary_fields_for_record_mutation() {
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    assert_eq!(sam_aux_get(&record, *b"NM"), None);
    assert_eq!(
        sam_aux_insert(&mut record, *b"NM", SamAuxValue::Int32(8)),
        None
    );
    assert_eq!(sam_aux_get(&record, *b"NM"), Some(SamAuxValue::Int32(8)));
    assert_eq!(
        sam_aux_insert(&mut record, *b"NM", SamAuxValue::UInt8(3)),
        Some(SamAuxValue::Int32(8))
    );
    assert_eq!(sam_aux_get(&record, *b"NM"), Some(SamAuxValue::UInt8(3)));

    assert_eq!(
        sam_aux_insert(
            &mut record,
            *b"ML",
            SamAuxValue::UInt8Array(vec![1, 128, 255])
        ),
        None
    );
    assert_eq!(
        sam_aux_get(&record, *b"ML"),
        Some(SamAuxValue::UInt8Array(vec![1, 128, 255]))
    );

    assert_eq!(
        sam_aux_insert(&mut record, *b"RG", SamAuxValue::String("rg0".into())),
        None
    );
    assert_eq!(
        sam_aux_remove(&mut record, *b"RG"),
        Some(SamAuxValue::String("rg0".into()))
    );
    assert_eq!(sam_aux_get(&record, *b"RG"), None);
}

#[test]
fn adapts_sam_header_reference_sequence_mutation() -> Result<(), Box<dyn std::error::Error>> {
    let mut header = htslib_rs::sam::Header::default();

    {
        let mut adapter = SamHeaderAdapter::new(&mut header);

        assert_eq!(adapter.reference_sequence_count(), 0);
        assert_eq!(adapter.reference_sequence_len("sq0"), None);
        assert_eq!(adapter.insert_reference_sequence("sq0", 13)?, None);
        assert_eq!(adapter.reference_sequence_count(), 1);
        assert_eq!(adapter.reference_sequence_len("sq0"), Some(13));
        assert_eq!(adapter.insert_reference_sequence("sq0", 8)?, Some(13));
        assert_eq!(adapter.reference_sequence_len("sq0"), Some(8));
        assert!(adapter.insert_reference_sequence("bad", 0).is_err());
        assert_eq!(adapter.remove_reference_sequence("sq0"), Some(8));
        assert_eq!(adapter.reference_sequence_count(), 0);
    }

    assert!(header.reference_sequences().is_empty());

    Ok(())
}

#[test]
fn adapts_sam_record_core_field_mutation() -> Result<(), Box<dyn std::error::Error>> {
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    {
        let mut adapter = SamRecordAdapter::new(&mut record);

        assert_eq!(adapter.name(), None);
        adapter.set_name(Some("r0"));
        assert_eq!(adapter.name().as_deref(), Some("r0"));

        assert_eq!(adapter.flags(), 0x04);
        adapter.set_flags(0x41);
        assert_eq!(adapter.flags(), 0x41);

        adapter.set_reference_sequence_id(Some(2));
        assert_eq!(adapter.reference_sequence_id(), Some(2));
        adapter.set_alignment_start(Some(13))?;
        assert_eq!(adapter.alignment_start(), Some(13));
        assert!(adapter.set_alignment_start(Some(0)).is_err());

        assert_eq!(adapter.mapping_quality(), 255);
        adapter.set_mapping_quality(42);
        assert_eq!(adapter.mapping_quality(), 42);
        adapter.set_mapping_quality(255);
        assert_eq!(adapter.mapping_quality(), 255);

        adapter.set_mate_reference_sequence_id(Some(3));
        assert_eq!(adapter.mate_reference_sequence_id(), Some(3));
        adapter.set_mate_alignment_start(Some(21))?;
        assert_eq!(adapter.mate_alignment_start(), Some(21));
        assert!(adapter.set_mate_alignment_start(Some(0)).is_err());

        adapter.set_template_length(-8);
        assert_eq!(adapter.template_length(), -8);
        adapter.set_name(None);
        assert_eq!(adapter.name(), None);
    }

    assert_eq!(
        sam_aux_insert(&mut record, *b"NM", SamAuxValue::UInt8(1)),
        None
    );
    assert_eq!(sam_aux_get(&record, *b"NM"), Some(SamAuxValue::UInt8(1)));

    Ok(())
}

#[test]
fn ports_test_view_sam_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let src = fixture("xx#minimal.sam");
    let actual = view_sam_text_from_path_with_limit(&src, Some(2))?;
    let mut expected = std::fs::read_to_string(src)?
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    expected.push_str(
        &std::fs::read_to_string(fixture("xx#minimal.sam"))?
            .lines()
            .filter(|line| !line.starts_with('@'))
            .take(2)
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_sam_ignore_parse_errors() -> Result<(), Box<dyn std::error::Error>> {
    let sam_path = temp_path("ignore-parse-errors.sam");
    let input = concat!(
        "@SQ\tSN:xx\tLN:20\n",
        "a0\t16\txx\t4\t1\t10H\t*\t0\t0\t*\t*\n",
        "malformed\tSAM\tline\n",
        "a1\t16\txx\t4\t1\t5H0M5H\t*\t0\t0\t*\t*\n",
    );
    std::fs::write(&sam_path, input)?;

    assert!(view_sam_text_from_path_with_limit(&sam_path, None).is_err());

    let actual = view_sam_text_from_path_with_limit_and_parse_errors(&sam_path, Some(2), true)?;
    let expected = concat!(
        "@SQ\tSN:xx\tLN:20\n",
        "a0\t16\txx\t4\t1\t10H\t*\t0\t0\t*\t*\n",
        "a1\t16\txx\t4\t1\t5H0M5H\t*\t0\t0\t*\t*\n",
    );

    cleanup(&[sam_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_alignment_benchmark_mode_counts() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        benchmark_sam_view_from_path(fixture("xx#minimal.sam"))?,
        sam_record_line_count(&view_sam_text_from_path_with_limit(
            fixture("xx#minimal.sam"),
            None
        )?)
    );
    assert_eq!(
        benchmark_bam_view_from_path(fixture("range.bam"))?,
        sam_record_line_count(&view_bam_as_sam_text_from_path_with_limit(
            fixture("range.bam"),
            None
        )?)
    );
    assert_eq!(
        benchmark_cram_view_from_path_with_reference(fixture("range.cram"), fixture("ce.fa"))?,
        sam_record_line_count(&view_cram_as_sam_text_from_path_with_reference_and_limit(
            fixture("range.cram"),
            fixture("ce.fa"),
            None
        )?)
    );

    Ok(())
}

#[test]
fn ports_test_realn_existing_baq_apply_and_revert() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        apply_existing_baq_from_sam_path(fixture("realn02_exp.sam"))?,
        std::fs::read_to_string(fixture("realn02_exp-a.sam"))?
    );
    assert_eq!(
        revert_existing_baq_from_sam_path(fixture("realn02_exp-a.sam"))?,
        std::fs::read_to_string(fixture("realn02_exp.sam"))?
    );

    Ok(())
}

#[test]
fn ports_test_realn_existing_baq_noop_states() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        apply_existing_baq_from_sam_path(fixture("realn02_exp-a.sam"))?,
        std::fs::read_to_string(fixture("realn02_exp-a.sam"))?
    );
    assert_eq!(
        revert_existing_baq_from_sam_path(fixture("realn02_exp.sam"))?,
        std::fs::read_to_string(fixture("realn02_exp.sam"))?
    );

    Ok(())
}

#[test]
fn ports_test_realn_non_extended_baq_recalculation() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        recalculate_baq_from_sam_path(fixture("realn01.sam"), fixture("realn01.fa"))?,
        std::fs::read_to_string(fixture("realn01_exp.sam"))?
    );
    assert_eq!(
        recalculate_baq_from_sam_path(fixture("realn02.sam"), fixture("realn02.fa"))?,
        std::fs::read_to_string(fixture("realn02_exp.sam"))?
    );

    Ok(())
}

#[test]
fn ports_test_realn_non_extended_baq_recalculation_apply() -> Result<(), Box<dyn std::error::Error>>
{
    assert_eq!(
        recalculate_and_apply_baq_from_sam_path(fixture("realn01.sam"), fixture("realn01.fa"))?,
        std::fs::read_to_string(fixture("realn01_exp-a.sam"))?
    );
    assert_eq!(
        recalculate_and_apply_baq_from_sam_path(fixture("realn02.sam"), fixture("realn02.fa"))?,
        std::fs::read_to_string(fixture("realn02_exp-a.sam"))?
    );

    Ok(())
}

#[test]
fn ports_test_realn_forced_non_extended_baq_recalculation() -> Result<(), Box<dyn std::error::Error>>
{
    assert_eq!(
        force_recalculate_baq_from_sam_path(fixture("realn02-r.sam"), fixture("realn02.fa"))?,
        std::fs::read_to_string(fixture("realn02_exp.sam"))?
    );

    Ok(())
}

#[test]
fn ports_test_realn_extended_baq_recalculation() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        recalculate_extended_baq_from_sam_path(fixture("realn01.sam"), fixture("realn01.fa"))?,
        std::fs::read_to_string(fixture("realn01_exp-e.sam"))?
    );
    assert_eq!(
        recalculate_extended_baq_from_sam_path(fixture("realn02.sam"), fixture("realn02.fa"))?,
        std::fs::read_to_string(fixture("realn02_exp-e.sam"))?
    );
    assert_eq!(
        recalculate_extended_baq_from_sam_path(fixture("realn03.sam"), fixture("realn03.fa"))?,
        std::fs::read_to_string(fixture("realn03_exp.sam"))?
    );

    Ok(())
}

#[test]
fn ports_sam_cap_mapq_cases() -> Result<(), Box<dyn std::error::Error>> {
    let sam_path = temp_path("cap-mapq.sam");
    let fasta_path = temp_path("cap-mapq.fa");
    let sam = "\
@SQ\tSN:sq0\tLN:8
perfect\t0\tsq0\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII
mismatch\t0\tsq0\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII
softclip\t0\tsq0\t1\t60\t2S6M\t*\t0\t0\tTTACGTAC\tIIIIIIII
";

    std::fs::write(&sam_path, sam)?;
    std::fs::write(&fasta_path, ">sq0\nACGTACGT\n")?;

    assert_eq!(
        cap_mapping_qualities_from_sam_path(&sam_path, &fasta_path, -1)?,
        vec![40, 28, 31]
    );
    assert_eq!(
        cap_mapping_qualities_from_sam_path(&sam_path, &fasta_path, 10)?,
        vec![10, -1, -1]
    );

    cleanup(&[sam_path, fasta_path]);

    Ok(())
}

#[test]
fn ports_test_view_bgzf_sam_view_output() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("view-index.sam.gz");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    std::fs::write(&sam_gz_path, compressed)?;

    let actual = view_bgzf_sam_text_from_path_with_limit(&sam_gz_path, None)?;
    let expected = view_sam_text_from_path_with_limit(fixture("index.sam"), None)?;

    cleanup(&[sam_gz_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_bgzf_sam_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("view-limit-index.sam.gz");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    std::fs::write(&sam_gz_path, compressed)?;

    let actual = view_bgzf_sam_text_from_path_with_limit(&sam_gz_path, Some(2))?;
    let expected = view_sam_text_from_path_with_limit(fixture("index.sam"), Some(2))?;

    cleanup(&[sam_gz_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_generated_large_bam_record_round_trip() -> Result<(), Box<dyn std::error::Error>>
{
    let sam_path = temp_path("large-rec.sam");
    let bam_path = temp_path("large-rec.bam");
    let cigar = "1M1I".repeat(16_000);
    let sequence = "A".repeat(32_000);
    let quality = "Q".repeat(32_000);
    let sam = format!(
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref\tLN:100000
read\t0\tref\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{quality}
"
    );

    std::fs::write(&sam_path, &sam)?;
    std::fs::write(&bam_path, write_bam_from_sam_path(&sam_path, Vec::new())?)?;

    assert_eq!(
        view_bam_as_sam_text_from_path_with_limit(&bam_path, None)?,
        sam
    );

    cleanup(&[sam_path, bam_path]);

    Ok(())
}

#[test]
fn ports_test_view_write_compressed_sam_output() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("write-compressed-sam.sam.gz");
    let encoded = write_bgzf_sam_from_path_with_limit(fixture("ce#1.sam"), Vec::new(), None)?;

    std::fs::write(&sam_gz_path, encoded)?;

    let actual = view_bgzf_sam_text_from_path_with_limit(&sam_gz_path, None)?;
    let expected = view_sam_text_from_path_with_limit(fixture("ce#1.sam"), None)?;

    cleanup(&[sam_gz_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_write_compressed_sam_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("write-compressed-sam-limit.sam.gz");
    let encoded = write_bgzf_sam_from_path_with_limit(fixture("ce#1000.sam"), Vec::new(), Some(3))?;

    std::fs::write(&sam_gz_path, encoded)?;

    let actual = view_bgzf_sam_text_from_path_with_limit(&sam_gz_path, None)?;
    let expected = view_sam_text_from_path_with_limit(fixture("ce#1000.sam"), Some(3))?;

    cleanup(&[sam_gz_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_sam_fastq_output_and_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let actual = view_sam_as_fastq_text_from_path_with_limit(fixture("ce#1.sam"), Some(1))?;
    let fields = first_sam_record_fields(fixture("ce#1.sam"))?;
    let expected = format!("@{}\n{}\n+\n{}\n", fields[0], fields[9], fields[10]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_sam_fasta_output_and_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let actual = view_sam_as_fasta_text_from_path_with_limit(fixture("ce#1.sam"), Some(1))?;
    let fields = first_sam_record_fields(fixture("ce#1.sam"))?;
    let expected = format!(">{}\n{}\n", fields[0], fields[9]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn reads_bam_header_and_records_from_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("range.bam");
    let header = read_bam_header_from_path(&path)?;

    assert!(reference_sequence_count(&header) > 0);
    assert!(count_bam_records_from_path(path)? > 0);

    Ok(())
}

#[test]
fn queries_bam_records_from_htslib_range_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("range.bam");
    let regions = [
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_IV:1500-1500".parse()?,
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_I:1000-1100".parse()?,
    ];

    let records = query_bam_regions_from_path(path, &regions)?;
    let expected_record_count = std::fs::read_to_string(fixture("range.out"))?
        .lines()
        .filter(|line| !line.starts_with('@'))
        .count();

    assert_eq!(records.len(), expected_record_count);

    Ok(())
}

#[test]
fn iterates_bam_records_from_htslib_range_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let region = "CHROMOSOME_II:2980-2980".parse()?;
    let records = query_bam_records_from_path(fixture("range.bam"), &region)?;
    let iter_count = iter_bam_records_from_path(fixture("range.bam"), &region)?.count();

    assert_eq!(iter_count, records.len());
    assert!(iter_count > 0);

    Ok(())
}

#[test]
fn writes_bam_round_trip_records_matching_htslib_fixture() -> Result<(), Box<dyn std::error::Error>>
{
    let bam_path = temp_path("write-round-trip.bam");
    let encoded = write_bam_from_path(fixture("range.bam"), Vec::new())?;
    std::fs::write(&bam_path, encoded)?;

    let actual = view_bam_as_sam_text_from_path_with_limit(&bam_path, None)?;
    let expected = view_bam_as_sam_text_from_path_with_limit(fixture("range.bam"), None)?;

    cleanup(&[bam_path]);

    assert_eq!(strip_sam_header(&actual), strip_sam_header(&expected));

    Ok(())
}

#[test]
fn writes_bam_from_sam_matching_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bam_path = temp_path("minimal-from-sam.bam");
    let encoded = write_bam_from_sam_path(fixture("xx#minimal.sam"), Vec::new())?;
    std::fs::write(&bam_path, encoded)?;

    let actual = view_bam_as_sam_text_from_path_with_limit(&bam_path, None)?;
    let expected = std::fs::read_to_string(fixture("xx#minimal.sam"))?;

    cleanup(&[bam_path]);

    assert_eq!(strip_sam_header(&actual), strip_sam_header(&expected));

    Ok(())
}

#[test]
fn ports_test_view_uncompressed_bam_output() -> Result<(), Box<dyn std::error::Error>> {
    let encoded = write_bam_from_sam_path_with_compression_level(
        fixture("xx#minimal.sam"),
        Vec::new(),
        htslib_rs::bgzf::io::writer::CompressionLevel::NONE,
    )?;
    let records = htslib_rs::bgzf_compat::read_all(&encoded[..])?;

    assert!(records.starts_with(b"BAM\x01"));
    assert_eq!(encoded[18] & 0x06, 0, "expected a stored DEFLATE block");

    Ok(())
}

#[test]
fn ports_test_view_bam_index_output() -> Result<(), Box<dyn std::error::Error>> {
    let bam_path = temp_path("index2-view-write.bam");
    let bai_path = temp_path("index2-view-write.bam.bai");

    write_bam_from_sam_path_with_bai(fixture("index2.sam"), &bam_path, &bai_path)?;

    let region = "1:1000000-1000000".parse()?;
    let records = query_bam_records_from_path(&bam_path, &region)?;

    cleanup(&[bam_path, bai_path]);

    assert_eq!(records.len(), 2);

    Ok(())
}

#[test]
fn ports_test_view_bam_whole_file_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let full = view_bam_as_sam_text_from_path_with_limit(fixture("range.bam"), None)?;
    let actual = view_bam_as_sam_text_from_path_with_limit(fixture("range.bam"), Some(2))?;
    let expected = sam_text_with_first_records(&full, 2);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_padded_bam_header_cases() -> Result<(), Box<dyn std::error::Error>> {
    let source_sam = fixture("ce#1.sam");
    let expected = std::fs::read_to_string(&source_sam)?;
    let encoded = write_bam_from_sam_path(&source_sam, Vec::new())?;

    for extra_nuls in [0, 1, 678] {
        let bam_path = temp_path(&format!("headernul{extra_nuls}.bam"));
        let padded = add_bam_header_nul_padding(&encoded, extra_nuls)?;

        std::fs::write(&bam_path, padded)?;

        let header = read_bam_header_from_path(&bam_path)?;
        let actual = view_bam_as_sam_text_from_path_with_limit(&bam_path, None)?;

        cleanup(&[bam_path]);

        assert!(
            reference_sequence_count(&header) > 0,
            "extra NULs: {extra_nuls}"
        );
        assert_eq!(
            strip_sam_header(&actual),
            strip_sam_header(&expected),
            "extra NULs: {extra_nuls}"
        );
    }

    Ok(())
}

#[test]
fn ports_test_view_bam_region_output() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("range.bam");
    let regions = [
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_IV:1500-1500".parse()?,
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_I:1000-1100".parse()?,
    ];

    let actual = view_bam_regions_as_sam_text_from_path(path, &regions)?;
    let expected = std::fs::read_to_string(fixture("range.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_bam_region_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let regions = [
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_IV:1500-1500".parse()?,
    ];

    let actual =
        view_bam_regions_as_sam_text_from_path_with_limit(fixture("range.bam"), &regions, Some(2))?;
    let expected = sam_text_with_first_records(&std::fs::read_to_string(fixture("range.out"))?, 2);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_bam_multi_region_output() -> Result<(), Box<dyn std::error::Error>> {
    let regions = [
        "CHROMOSOME_I:1122-1122".parse()?,
        "CHROMOSOME_II:1136-1136".parse()?,
        "CHROMOSOME_II:1241-1241".parse()?,
        "CHROMOSOME_II:1267-1267".parse()?,
        "CHROMOSOME_II:1326-1326".parse()?,
        "CHROMOSOME_II:1345-1345".parse()?,
        "CHROMOSOME_II:1353-1353".parse()?,
        "CHROMOSOME_II:1366-1366".parse()?,
        "CHROMOSOME_II:1416-1416".parse()?,
        "CHROMOSOME_II:1459-1459".parse()?,
        "CHROMOSOME_II:1536-1536".parse()?,
    ];

    let actual =
        view_bam_regions_as_sam_text_from_path_with_dedup(fixture("range.bam"), &regions, true)?;
    let expected = std::fs::read_to_string(fixture("range.out2"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_cram_region_output() -> Result<(), Box<dyn std::error::Error>> {
    let regions = [
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_IV:1500-1500".parse()?,
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_I:1000-1100".parse()?,
    ];

    let actual = view_cram_regions_as_sam_text_from_path_with_reference(
        fixture("range.cram"),
        fixture("ce.fa"),
        &regions,
        false,
    )?;
    let expected = std::fs::read_to_string(fixture("range.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_cram_index_output() -> Result<(), Box<dyn std::error::Error>> {
    let cram_path = temp_path("ce-1-view-write.cram");
    let crai_path = temp_path("ce-1-view-write.cram.crai");

    write_cram_from_sam_path_with_reference_and_crai(
        fixture("ce#1.sam"),
        fixture("ce.fa"),
        &cram_path,
        &crai_path,
    )?;

    let region = "CHROMOSOME_I:2-102".parse()?;
    let records =
        query_cram_records_from_path_with_reference(&cram_path, &region, fixture("ce.fa"))?;

    cleanup(&[cram_path, crai_path]);

    assert!(!records.is_empty());

    Ok(())
}

#[test]
fn ports_test_view_cram_whole_file_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let full = view_cram_as_sam_text_from_path_with_reference_and_limit(
        fixture("range.cram"),
        fixture("ce.fa"),
        None,
    )?;
    let actual = view_cram_as_sam_text_from_path_with_reference_and_limit(
        fixture("range.cram"),
        fixture("ce.fa"),
        Some(2),
    )?;
    let expected = sam_text_with_first_records(&full, 2);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_cram_region_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let regions = [
        "CHROMOSOME_II:2980-2980".parse()?,
        "CHROMOSOME_IV:1500-1500".parse()?,
    ];

    let actual = view_cram_regions_as_sam_text_from_path_with_reference_and_limit(
        fixture("range.cram"),
        fixture("ce.fa"),
        &regions,
        Some(2),
    )?;
    let expected = sam_text_with_first_records(&std::fs::read_to_string(fixture("range.out"))?, 2);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_cram_multi_region_output() -> Result<(), Box<dyn std::error::Error>> {
    let regions = [
        "CHROMOSOME_I:1122-1122".parse()?,
        "CHROMOSOME_II:1136-1136".parse()?,
        "CHROMOSOME_II:1241-1241".parse()?,
        "CHROMOSOME_II:1267-1267".parse()?,
        "CHROMOSOME_II:1326-1326".parse()?,
        "CHROMOSOME_II:1345-1345".parse()?,
        "CHROMOSOME_II:1353-1353".parse()?,
        "CHROMOSOME_II:1366-1366".parse()?,
        "CHROMOSOME_II:1416-1416".parse()?,
        "CHROMOSOME_II:1459-1459".parse()?,
        "CHROMOSOME_II:1536-1536".parse()?,
    ];

    let actual = view_cram_regions_as_sam_text_from_path_with_reference(
        fixture("range.cram"),
        fixture("ce.fa"),
        &regions,
        true,
    )?;
    let expected = std::fs::read_to_string(fixture("range.out2"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_pl_bam_view_without_raw_sq_header() -> Result<(), Box<dyn std::error::Error>> {
    let actual = view_bam_as_sam_text_from_path_with_limit(fixture("no_hdr_sq_1.bam"), None)?;
    let expected = std::fs::read_to_string(fixture("no_hdr_sq_1.expected.sam"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_pl_bam_region_view_without_raw_sq_header() -> Result<(), Box<dyn std::error::Error>> {
    let actual = view_bam_regions_as_sam_text_from_path(
        fixture("no_hdr_sq_1.bam"),
        &["CHROMOSOME_I".parse()?],
    )?;
    let expected = std::fs::read_to_string(fixture("no_hdr_sq_1.expected.sam"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_sam_filter_integer_expression_counts() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        count_sam_records_matching_filter_from_path(fixture("ce#1000.sam"), "pos % 23 == 11")?,
        98
    );
    assert_eq!(
        count_sam_records_matching_filter_from_path(
            fixture("ce#1000.sam"),
            "qlen/(flag*mapq+pos)>5"
        )?,
        38
    );
    assert_eq!(
        count_sam_records_matching_filter_from_path(
            fixture("ce#1000.sam"),
            r#"[NM]>=10 || [MD]=~"A.*A.*A""#
        )?,
        72
    );

    Ok(())
}

#[test]
fn ports_sam_filter_function_expression_counts() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        ("ce#5b.sam", "length(seq) != qlen", 1),
        ("ce#1000.sam", "min(qual) >= 20", 42),
        ("ce#1000.sam", "max(qual) <= 20", 2),
        ("ce#1000.sam", "avg(qual) >= 20 && avg(qual) <= 30", 604),
    ];

    for (src, filter, expected) in cases {
        assert_eq!(
            count_sam_records_matching_filter_from_path(fixture(src), filter)?,
            expected,
            "{filter}"
        );
    }

    Ok(())
}

#[test]
fn ports_sam_filter_cigar_metric_output() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        ("realn02.sam", "sclen>=20", "sam_filter/func5.out"),
        ("realn02.sam", "rlen<50", "sam_filter/func6.out"),
        ("realn02.sam", "qlen>100", "sam_filter/func7.out"),
        ("c1#clip.sam", "hclen>=4", "sam_filter/func8.out"),
    ];

    for (src, filter, expected) in cases {
        let actual = view_sam_text_matching_filter_from_path(fixture(src), filter)?;
        let expected = std::fs::read_to_string(fixture(expected))?;

        assert_eq!(strip_sam_header(&actual), expected, "{filter}");
    }

    Ok(())
}

#[test]
fn ports_sam_filter_string_expression_output() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        (
            "ce#1000.sam",
            r#"qname =~ "\.1" && cigar =~ "D""#,
            "sam_filter/string1.out",
        ),
        (
            "ce#5b.sam",
            r#"rname=="CHROMOSOME_II""#,
            "sam_filter/string2.out",
        ),
        (
            "ce#5b.sam",
            r#"rname=~"CHROMOSOME_II""#,
            "sam_filter/string3.out",
        ),
        ("ce#1000.sam", r#"cigar=~"D""#, "sam_filter/string4.out"),
        (
            "ce#1000.sam",
            r#"seq =~ "(AT){2}""#,
            "sam_filter/string5.out",
        ),
        ("xx#rg.sam", r#"library=="x""#, "sam_filter/string6.out"),
        ("xx#rg.sam", r#"library!="x""#, "sam_filter/string7.out"),
    ];

    for (src, filter, expected) in cases {
        let actual = view_sam_text_matching_filter_from_path(fixture(src), filter)?;
        let expected = std::fs::read_to_string(fixture(expected))?;

        assert_eq!(actual, expected, "{filter}");
    }

    Ok(())
}

#[test]
fn queries_bam_with_htslib_index_lookup_fallbacks() -> Result<(), Box<dyn std::error::Error>> {
    let bam_path = temp_path("range.bam");
    let replaced_bai_path = temp_path("range.bai");
    let explicit_bai_path = temp_path("custom-range.bai");

    std::fs::copy(fixture("range.bam"), &bam_path)?;
    let index = build_bai(&bam_path)?;
    write_bai(&replaced_bai_path, &index)?;
    write_bai(&explicit_bai_path, &index)?;

    let region = "CHROMOSOME_II:2980-2980".parse()?;
    let records = query_bam_records_from_path(&bam_path, &region)?;
    let explicit_src = format!(
        "{}##idx##{}",
        bam_path.display(),
        explicit_bai_path.display()
    );
    let explicit_records = query_bam_records_from_path(explicit_src, &region)?;

    cleanup(&[bam_path, replaced_bai_path, explicit_bai_path]);

    assert!(!records.is_empty());
    assert_eq!(explicit_records.len(), records.len());

    Ok(())
}

#[test]
fn ports_test_pl_index2_mapped_unmapped_pair_queries() -> Result<(), Box<dyn std::error::Error>> {
    let bam_path = temp_path("index2.bam");
    let bai_path = temp_path("index2.bam.bai");

    let encoded = write_bam_from_sam_path(fixture("index2.sam"), Vec::new())?;
    std::fs::write(&bam_path, encoded)?;
    let index = build_bai(&bam_path)?;
    write_bai(&bai_path, &index)?;

    for tid in 1..=2 {
        for pos in 1..=2 {
            let region = format!("{tid}:{pos}000000-{pos}000000").parse()?;
            let records = query_bam_records_from_path(&bam_path, &region)?;
            let viewed = view_bam_regions_as_sam_text_from_path(&bam_path, &[region])?;
            let viewed_record_count = viewed.lines().filter(|line| !line.starts_with('@')).count();

            assert_eq!(records.len(), 2, "{tid}:{pos}000000");
            assert_eq!(viewed_record_count, 2, "{tid}:{pos}000000");
        }
    }

    cleanup(&[bam_path, bai_path]);

    Ok(())
}

#[test]
fn queries_bgzf_sam_records_from_htslib_index_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("index.sam.gz");
    let csi_path = temp_path("index.sam.gz.csi");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    std::fs::write(&sam_gz_path, compressed)?;
    let index = build_sam_csi(&sam_gz_path)?;
    write_csi(&csi_path, &index)?;

    let region = "CHROMOSOME_I:999901-1000000".parse()?;
    let records = query_sam_records_from_path(&sam_gz_path, &region)?;

    cleanup(&[sam_gz_path, csi_path]);

    assert!(!records.is_empty());

    Ok(())
}

#[test]
fn queries_bgzf_sam_records_with_bai_index() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("index-bai.sam.gz");
    let bai_path = temp_path("index-bai.sam.gz.bai");
    let sam = std::fs::read(fixture("index.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    std::fs::write(&sam_gz_path, compressed)?;
    let index = build_sam_bai(&sam_gz_path)?;
    write_bai(&bai_path, &index)?;

    let region = "CHROMOSOME_I:999901-1000000".parse()?;
    let records = query_sam_records_from_path(&sam_gz_path, &region)?;

    cleanup(&[sam_gz_path, bai_path]);

    assert!(!records.is_empty());

    Ok(())
}

#[test]
fn ports_test_pl_large_position_sam_iterators() -> Result<(), Box<dyn std::error::Error>> {
    let sam_gz_path = temp_path("longref.sam.gz");
    let csi_path = temp_path("longref.sam.gz.csi");
    let sam = std::fs::read(fixture("longrefs/longref.sam"))?;

    let compressed = write_bgzf_all(Vec::new(), &sam)?;
    std::fs::write(&sam_gz_path, compressed)?;
    let index = build_sam_csi(&sam_gz_path)?;
    write_csi(&csi_path, &index)?;

    let single_region = ["CHROMOSOME_I:10000000000-10000000003".parse()?];
    let actual = view_sam_regions_as_text_from_path(&sam_gz_path, &single_region)?;
    let expected = std::fs::read_to_string(fixture("longrefs/longref_itr.expected.sam"))?;
    assert_eq!(actual, expected);

    let multi_regions = [
        "CHROMOSOME_I:10000000000-10000000003".parse()?,
        "CHROMOSOME_I:10000000100-10000000110".parse()?,
    ];
    let actual = view_sam_regions_as_text_from_path(&sam_gz_path, &multi_regions)?;
    let expected = std::fs::read_to_string(fixture("longrefs/longref_multi.expected.sam"))?;
    assert_eq!(actual, expected);

    cleanup(&[sam_gz_path, csi_path]);

    Ok(())
}

#[test]
fn reads_cram_header_from_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let header = read_cram_header_from_path(fixture("range.cram"))?;

    assert!(reference_sequence_count(&header) > 0);

    Ok(())
}

#[test]
fn queries_cram_records_with_htslib_index_lookup_fallbacks()
-> Result<(), Box<dyn std::error::Error>> {
    let cram_path = temp_path("range.cram");
    let replaced_crai_path = temp_path("range.crai");
    let explicit_crai_path = temp_path("custom-range.crai");

    std::fs::copy(fixture("range.cram"), &cram_path)?;
    std::fs::copy(fixture("range.cram.crai"), &replaced_crai_path)?;
    std::fs::copy(fixture("range.cram.crai"), &explicit_crai_path)?;

    let region = "CHROMOSOME_II:2980-2980".parse()?;
    let records =
        query_cram_records_from_path_with_reference(&cram_path, &region, fixture("ce.fa"))?;
    let explicit_src = format!(
        "{}##idx##{}",
        cram_path.display(),
        explicit_crai_path.display()
    );
    let explicit_records =
        query_cram_records_from_path_with_reference(explicit_src, &region, fixture("ce.fa"))?;

    cleanup(&[cram_path, replaced_crai_path, explicit_crai_path]);

    assert!(!records.is_empty());
    assert_eq!(explicit_records.len(), records.len());

    Ok(())
}

#[test]
fn iterates_cram_records_from_htslib_range_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let region = "CHROMOSOME_II:2980-2980".parse()?;
    let records = query_cram_records_from_path_with_reference(
        fixture("range.cram"),
        &region,
        fixture("ce.fa"),
    )?;
    let iter_count = iter_cram_records_from_path_with_reference(
        fixture("range.cram"),
        &region,
        fixture("ce.fa"),
    )?
    .count();

    assert_eq!(iter_count, records.len());
    assert!(iter_count > 0);

    Ok(())
}

#[test]
fn writes_cram_round_trip_records_matching_htslib_fixture() -> Result<(), Box<dyn std::error::Error>>
{
    let out_path = temp_path("range-roundtrip.cram");

    let encoded =
        write_cram_from_path_with_reference(fixture("range.cram"), fixture("ce.fa"), Vec::new())?;
    std::fs::write(&out_path, encoded)?;

    let expected =
        summarize_cram_records_from_path_with_reference(fixture("range.cram"), fixture("ce.fa"))?;
    let actual = summarize_cram_records_from_path_with_reference(&out_path, fixture("ce.fa"))?;

    cleanup(&[out_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn writes_cram_from_sam_matching_htslib_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        ("ce#1.sam", "ce.fa", "ce-1-from-sam.cram"),
        ("ce#2.sam", "ce.fa", "ce-2-from-sam.cram"),
        ("ce#1000.sam", "ce.fa", "ce-1000-from-sam.cram"),
        ("xx#pair.sam", "xx.fa", "xx-pair-from-sam.cram"),
    ];

    for (sam, reference, temp_name) in cases {
        let out_path = temp_path(temp_name);

        let encoded =
            write_cram_from_sam_path_with_reference(fixture(sam), fixture(reference), Vec::new())?;
        std::fs::write(&out_path, encoded)?;

        let expected = summarize_sam_records_from_path(fixture(sam))?;
        let actual =
            summarize_cram_records_from_path_with_reference(&out_path, fixture(reference))?;

        cleanup(&[out_path]);

        assert_eq!(actual, expected, "case: {sam}");
    }

    Ok(())
}

#[test]
fn ports_tlen_cram_auto_creation_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        "a7", "a7b", "a8", "a8b", "a9", "a9b", "b7", "b7b", "b8", "b8b", "c7", "c7b", "c8", "c8b",
        "d7", "d7b", "d4", "d4b", "d4c", "d4d", "d4e", "d4f", "d5", "d5b", "d5c", "d5d", "d5e",
        "d5f", "a4", "a5",
    ];

    for case in cases {
        let sam = summarize_sam_records_from_path(fixture(format!("tlen/{case}.sam")))?;
        let cram = summarize_cram_records_from_path(fixture(format!("tlen/{case}.cram")))?;

        assert_eq!(cram, sam, "case: {case}");
        assert!(
            cram.iter().all(|record| record.template_length() != 0),
            "case: {case}"
        );
    }

    Ok(())
}

#[test]
fn ports_base_modification_mm_variants() -> Result<(), Box<dyn std::error::Error>> {
    let actual =
        base_modification_report_from_sam_path(fixture("base_mods/MM-variants.sam"), true)?;
    let expected = std::fs::read_to_string(fixture("base_mods/MM-variants.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_base_modification_mm_probability_cases() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        ("MM-explicit.sam", false, "MM-explicit.out"),
        ("MM-explicit.sam", true, "MM-explicit-f.out"),
        ("MM-double.sam", false, "MM-double.out"),
        ("MM-chebi.sam", false, "MM-chebi.out"),
        ("MM-multi.sam", false, "MM-multi.out"),
        ("MM-not-all-modded.sam", false, "MM-not-all-modded.out"),
    ];

    for (src, report_unchecked, expected) in cases {
        let actual = base_modification_report_from_sam_path(
            fixture(format!("base_mods/{src}")),
            report_unchecked,
        )?;
        let expected = std::fs::read_to_string(fixture(format!("base_mods/{expected}")))?;

        assert_eq!(
            actual, expected,
            "case: {src}, report_unchecked: {report_unchecked}"
        );
    }

    Ok(())
}

#[test]
fn ports_base_modification_extended_metadata_case() -> Result<(), Box<dyn std::error::Error>> {
    let actual =
        extended_base_modification_report_from_sam_path(fixture("base_mods/MM-explicit.sam"))?;
    let expected = std::fs::read_to_string(fixture("base_mods/MM-explicit-x.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_base_modification_pileup_cases() -> Result<(), Box<dyn std::error::Error>> {
    for (src, expected) in [
        ("MM-pileup.sam", "MM-pileup.out"),
        ("MM-pileup2.sam", "MM-pileup2.out"),
        ("MM-MNp.sam", "MM-pileup.out"),
    ] {
        let actual =
            base_modification_pileup_report_from_sam_path(fixture(format!("base_mods/{src}")))?;
        let expected = std::fs::read_to_string(fixture(format!("base_mods/{expected}")))?;

        assert_eq!(actual, expected, "case: {src}");
    }

    Ok(())
}

#[test]
fn ports_mpileup_sam_cigar_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    for (src, expected) in [
        ("mp_D.sam", "mp_D.out"),
        ("mp_I.sam", "mp_I.out"),
        ("mp_DI.sam", "mp_DI.out"),
        ("mp_ID.sam", "mp_ID.out"),
        ("mp_N.sam", "mp_N.out"),
        ("mp_N2.sam", "mp_N2.out"),
        ("mp_P.sam", "mp_P.out"),
        ("c1#pad1.sam", "c1#pad1.out"),
        ("c1#pad2.sam", "c1#pad2.out"),
        ("c1#pad3.sam", "c1#pad3.out"),
    ] {
        let actual = pileup_report_from_sam_path(fixture(format!("mpileup/{src}")))?;
        let expected = std::fs::read_to_string(fixture(format!("mpileup/{expected}")))?;

        assert_eq!(actual, expected, "case: {src}");
    }

    Ok(())
}

#[test]
fn ports_mpileup_overlap_removal_fixtures() -> Result<(), Box<dyn std::error::Error>> {
    for (src, expected) in [
        ("mp_overlap1.sam", "mp_overlap1.out"),
        ("mp_overlap2.sam", "mp_overlap2.out"),
    ] {
        let actual = pileup_report_from_sam_path(fixture(format!("mpileup/{src}")))?;
        let expected = std::fs::read_to_string(fixture(format!("mpileup/{expected}")))?;

        assert_eq!(actual, expected, "case: {src}");
    }

    Ok(())
}

#[test]
fn ports_mpileup_small_bam_edge_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let actual = pileup_report_from_bam_path(fixture("mpileup/small.bam"))?;
    let expected = std::fs::read_to_string(fixture("mpileup/small.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn rejects_invalid_base_modification_mn_and_bounds() {
    for src in [
        "MM-MNf1.sam",
        "MM-MNf2.sam",
        "MM-bounds+.sam",
        "MM-bounds-.sam",
    ] {
        assert!(
            base_modification_report_from_sam_path(fixture(format!("base_mods/{src}")), false)
                .is_err(),
            "case: {src}"
        );
        assert!(
            base_modification_pileup_report_from_sam_path(fixture(format!("base_mods/{src}")))
                .is_err(),
            "pileup case: {src}"
        );
    }
}
