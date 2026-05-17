use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use htslib_rs::tabix_compat::{TextFormat, write_bgzf_and_index};
use htslib_rs::variant_io_compat::{
    BcfFloatValue, HeaderNumber, VcfHeaderId, VcfRecordAdapter, benchmark_bcf_view_from_path,
    benchmark_vcf_view_from_path, build_bcf_csi_from_path, contig_count,
    count_bcf_records_from_path, count_vcf_records_from_path, filter_vcf_text_by_region_from_path,
    filter_vcf_text_by_target_from_path, format_header_number, header_has_contig,
    header_has_filter, header_has_format, header_has_info, header_has_other_record_id,
    header_has_unstructured_other_record, htslib_variant_span_from_vcf_line,
    htslib_variant_spans_from_vcf, info_header_number, iter_bcf_records_from_path,
    iter_vcf_records_from_path, normalize_bcf_float_vector, parse_numbered_header_record,
    query_bcf_records_by_reference_id_from_path, query_bcf_records_from_path,
    query_vcf_records_from_path, read_bcf_header_from_path, read_bcf_sweep_from_path,
    read_vcf_header_from_path, read_vcf_sweep_from_path, remove_header_contig,
    remove_header_filter, remove_header_format, remove_header_info, remove_header_other_record_id,
    remove_header_other_records, remove_vcf_allele_set_from_line, sample_count,
    synced_bcf_output_no_index_from_paths, synced_vcf_output_no_index_from_paths,
    synced_vcf_summary_no_index_from_paths, synced_vcf_summary_no_index_from_readers,
    test_vcf_api_record_serialization_text, translated_bcf_record_fixture_vcf_text,
    update_vcf_line_alleles, update_vcf_line_format_i32, update_vcf_line_info_i32,
    vcf_format_i32_values_from_line, vcf_header_has_id, vcf_header_remove_id,
    vcf_info_float_values_from_line, vcf_info_i32_values_from_line, vcf_open_mode_suffix,
    view_bcf_as_vcf_text_from_path_with_limit, view_bcf_regions_as_vcf_text_from_path,
    view_bcf_regions_as_vcf_text_from_path_with_limit, view_bcf_targets_as_vcf_text_from_path,
    view_vcf_regions_as_text_from_path, view_vcf_regions_as_text_from_path_with_limit,
    view_vcf_text_from_path_with_limit, write_bcf_csi_from_path, write_bcf_from_vcf_path,
    write_bcf_from_vcf_path_with_csi, write_vcf_bgzf_from_path_with_tbi, write_vcf_from_path,
    write_vcf_header_from_path,
};
use htslib_rs::{index_compat::write_csi, tabix_compat::write_bgzf_and_csi};

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib/test")
        .join(path)
}

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("htslib-rs-variant-{}-{name}", std::process::id()))
}

fn cleanup(paths: &[PathBuf]) {
    for path in paths {
        let _ = std::fs::remove_file(path);
    }
}

fn parse_header(lines: &[&str]) -> htslib_rs::vcf::Header {
    let mut raw = String::from("##fileformat=VCFv4.3\n");

    for line in lines {
        raw.push_str(line);
        raw.push('\n');
    }

    raw.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    raw.parse().expect("header should parse")
}

fn vcf_text_with_first_records(src: &str, n: usize) -> String {
    let mut expected = src
        .lines()
        .take_while(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    expected.push_str(
        &src.lines()
            .filter(|line| !line.starts_with('#'))
            .take(n)
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );
    expected
}

fn vcf_record_lines(src: &str) -> Vec<&str> {
    src.lines().filter(|line| !line.starts_with('#')).collect()
}

fn vcf_record_line_count(src: &str) -> usize {
    vcf_record_lines(src).len()
}

fn vcf_miniview_filtered_text(src: &str) -> String {
    const ERASE_TAGS: &[&str] = &[
        "IMF=", "DP=", "IDV=", "IMP=", "IS=", "VDB=", "SGB=", "MQB=", "BQB=", "RPB=", "MQ0F=",
        "MQSB=",
    ];

    let mut out = String::new();
    for line in src.lines() {
        if line.starts_with("##") {
            continue;
        }

        let mut line = line.to_string();
        if !line.starts_with('#') {
            for tag in ERASE_TAGS {
                erase_vcf_miniview_tag(&mut line, tag);
            }
        }
        out.push_str(&line);
        out.push('\n');
    }
    out
}

fn erase_vcf_miniview_tag(line: &mut String, tag: &str) {
    let Some(mut begin) = line.get(1..).and_then(|s| s.find(tag).map(|i| i + 1)) else {
        return;
    };

    let bytes = line.as_bytes();
    let mut end = begin;
    while end < bytes.len() && bytes[end] != b'\t' && bytes[end] != b';' {
        end += 1;
    }
    if begin > 0 && bytes[begin - 1] == b';' {
        begin -= 1;
    }
    line.replace_range(begin..end, "");
}

fn bcf_sr_weird_chromosome_cases() -> [(&'static str, &'static str); 12] {
    [
        ("1", "bcf-sr/weird-chr-names.1.out"),
        ("1:1-2", "bcf-sr/weird-chr-names.1.out"),
        ("1:1,1:2", "bcf-sr/weird-chr-names.1.out"),
        ("1:1-1", "bcf-sr/weird-chr-names.2.out"),
        ("{1:1}", "bcf-sr/weird-chr-names.3.out"),
        ("{1:1}:1-2", "bcf-sr/weird-chr-names.3.out"),
        ("{1:1}:1,{1:1}:2", "bcf-sr/weird-chr-names.3.out"),
        ("{1:1}:1-1", "bcf-sr/weird-chr-names.4.out"),
        ("{1:1-1}", "bcf-sr/weird-chr-names.5.out"),
        ("{1:1-1}:1-2", "bcf-sr/weird-chr-names.5.out"),
        ("{1:1-1}:1,{1:1-1}:2", "bcf-sr/weird-chr-names.5.out"),
        ("{1:1-1}:1-1", "bcf-sr/weird-chr-names.6.out"),
    ]
}

#[test]
fn reads_vcf_header_and_records_from_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("index.vcf");
    let header = read_vcf_header_from_path(&path)?;

    assert!(contig_count(&header) > 0);
    assert!(count_vcf_records_from_path(path)? > 0);

    Ok(())
}

#[test]
fn queries_vcf_records_from_htslib_tabix_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = std::env::temp_dir().join(format!(
        "htslib-rs-variant-vcf-{}.vcf.gz",
        std::process::id()
    ));
    let tbi_path = std::env::temp_dir().join(format!(
        "htslib-rs-variant-vcf-{}.vcf.gz.tbi",
        std::process::id()
    ));

    let vcf = std::fs::File::open(fixture("tabix/vcf_file.vcf"))?;
    write_bgzf_and_index(BufReader::new(vcf), &bgzf_path, &tbi_path, TextFormat::Vcf)?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_vcf_records_from_path(&bgzf_path, &region)?;
    let expected_record_count = std::fs::read_to_string(fixture("tabix/vcf_file.1.3000151.out"))?
        .lines()
        .count();

    std::fs::remove_file(bgzf_path)?;
    std::fs::remove_file(tbi_path)?;

    assert_eq!(records.len(), expected_record_count);

    Ok(())
}

#[test]
fn iterates_vcf_records_from_htslib_tabix_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("iter.vcf.gz");
    let tbi_path = temp_path("iter.vcf.gz.tbi");

    let vcf = std::fs::File::open(fixture("tabix/vcf_file.vcf"))?;
    write_bgzf_and_index(BufReader::new(vcf), &bgzf_path, &tbi_path, TextFormat::Vcf)?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_vcf_records_from_path(&bgzf_path, &region)?;
    let iter_count = iter_vcf_records_from_path(&bgzf_path, &region)?.count();

    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(iter_count, records.len());
    assert!(iter_count > 0);

    Ok(())
}

#[test]
fn queries_vcf_with_htslib_index_lookup_fallbacks() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("lookup.vcf.gz");
    let replaced_tbi_path = temp_path("lookup.vcf.tbi");
    let explicit_csi_path = temp_path("vcf-custom-lookup.csi");

    let vcf = std::fs::File::open(fixture("tabix/vcf_file.vcf"))?;
    write_bgzf_and_index(
        BufReader::new(vcf),
        &bgzf_path,
        &replaced_tbi_path,
        TextFormat::Vcf,
    )?;

    let vcf = std::fs::File::open(fixture("tabix/vcf_file.vcf"))?;
    write_bgzf_and_csi(
        BufReader::new(vcf),
        &bgzf_path,
        &explicit_csi_path,
        TextFormat::Vcf,
    )?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_vcf_records_from_path(&bgzf_path, &region)?;
    let explicit_src = format!(
        "{}##idx##{}",
        bgzf_path.display(),
        explicit_csi_path.display()
    );
    let explicit_records = query_vcf_records_from_path(explicit_src, &region)?;

    cleanup(&[bgzf_path, replaced_tbi_path, explicit_csi_path]);

    assert_eq!(records.len(), 1);
    assert_eq!(explicit_records.len(), 1);

    Ok(())
}

#[test]
fn writes_vcf_text_matching_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("tabix/vcf_file.vcf");
    let actual = write_vcf_from_path(&path, Vec::new())?;
    let expected = std::fs::read(path)?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_pl_vcf_various_canonical_output() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        ("formatcols.vcf", "formatcols.vcf"),
        ("noroundtrip.vcf", "noroundtrip-out.vcf"),
        ("formatmissing.vcf", "formatmissing-out.vcf"),
        ("vcf_meta_meta.vcf", "vcf_meta_meta.vcf"),
        ("vcf44_1.vcf", "vcf44_1.expected"),
    ];

    for (input, expected) in cases {
        let actual = write_vcf_from_path(fixture(input), Vec::new()).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{input}: {e}"))
        })?;
        let expected = std::fs::read(fixture(expected))?;

        assert_eq!(actual, expected, "case: {input}");
    }

    Ok(())
}

#[test]
fn ports_test_pl_vcf_various_header_canonicalization() -> Result<(), Box<dyn std::error::Error>> {
    let actual = write_vcf_header_from_path(fixture("test-vcf-hdr-in.vcf"), Vec::new())?;
    let expected = std::fs::read(fixture("test-vcf-hdr.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_vcf_indexed_region_header_output() -> Result<(), Box<dyn std::error::Error>> {
    let region = "chr22:1-2".parse()?;
    let actual = view_vcf_regions_as_text_from_path(fixture("modhdr.vcf.gz"), &[region])?;
    let expected = std::fs::read_to_string(fixture("modhdr.expected.vcf"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_view_vcf_indexed_region_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("region-limit.vcf.gz");
    let tbi_path = temp_path("region-limit.vcf.gz.tbi");
    let vcf = File::open(fixture("tabix/vcf_file.vcf"))?;

    write_bgzf_and_index(BufReader::new(vcf), &bgzf_path, &tbi_path, TextFormat::Vcf)?;

    let regions = ["1:3000150-3258501".parse()?];
    let actual = view_vcf_regions_as_text_from_path_with_limit(&bgzf_path, &regions, Some(2))?;
    let source = std::fs::read_to_string(fixture("tabix/vcf_file.vcf"))?;
    let expected_records = vcf_record_lines(&source)
        .into_iter()
        .take(2)
        .collect::<Vec<_>>();

    cleanup(&[bgzf_path, tbi_path]);

    assert!(actual.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    assert_eq!(vcf_record_lines(&actual), expected_records);

    Ok(())
}

#[test]
fn reads_bcf_header_and_records_from_htslib_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let path = fixture("tabix/vcf_file.bcf");
    let header = read_bcf_header_from_path(&path)?;

    assert!(sample_count(&header) > 0);
    assert!(count_bcf_records_from_path(path)? > 0);

    Ok(())
}

#[test]
fn ports_test_view_vcf_compressed_index_output() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("view-write-index.vcf.gz");
    let tbi_path = temp_path("view-write-index.vcf.gz.tbi");

    write_vcf_bgzf_from_path_with_tbi(fixture("tabix/vcf_file.vcf"), &bgzf_path, &tbi_path)?;

    let actual = view_vcf_regions_as_text_from_path(&bgzf_path, &["1:3000150-3258501".parse()?])?;
    let expected =
        filter_vcf_text_by_region_from_path(fixture("tabix/vcf_file.vcf"), "1:3000150-3258501")?;

    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(vcf_record_lines(&actual), vcf_record_lines(&expected));

    Ok(())
}

#[test]
fn ports_test_view_bcf_to_vcf_output() -> Result<(), Box<dyn std::error::Error>> {
    let actual = view_bcf_as_vcf_text_from_path_with_limit(fixture("tabix/vcf_file.bcf"), None)?;
    let expected = std::fs::read_to_string(fixture("tabix/vcf_file.vcf"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_vcf_miniview_filtered_bcf_output() -> Result<(), Box<dyn std::error::Error>> {
    let vcf_path = temp_path("miniview.vcf");
    let bcf_path = temp_path("miniview.bcf");
    let vcf = concat!(
        "##fileformat=VCFv4.3\n",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
        "##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Bias\">\n",
        "##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Bias\">\n",
        "##INFO=<ID=KEEP,Number=1,Type=Integer,Description=\"Keep\">\n",
        "##contig=<ID=1,length=100>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "1\t10\t.\tA\tC\t.\tPASS\tDP=9;KEEP=1;VDB=0.5;MQSB=0.7\n",
        "1\t11\t.\tG\tT\t.\tPASS\tKEEP=2;DP=10\n",
    );
    std::fs::write(&vcf_path, vcf)?;
    std::fs::write(&bcf_path, write_bcf_from_vcf_path(&vcf_path, Vec::new())?)?;

    let viewed = view_bcf_as_vcf_text_from_path_with_limit(&bcf_path, None)?;
    let filtered = vcf_miniview_filtered_text(&viewed);

    cleanup(&[vcf_path, bcf_path]);

    assert!(filtered.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(!filtered.contains("##INFO"));
    assert!(!filtered.contains("DP="));
    assert!(!filtered.contains("VDB="));
    assert!(!filtered.contains("MQSB="));
    assert!(filtered.contains("KEEP=1"));
    assert!(filtered.contains("KEEP=2"));

    Ok(())
}

#[test]
fn ports_test_view_bcf_indexed_region_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("region-limit.bcf");
    let csi_path = temp_path("region-limit.bcf.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    let actual =
        view_bcf_regions_as_vcf_text_from_path_with_limit(&bcf_path, "1:3000150-3258501", Some(2))?;
    let source = std::fs::read_to_string(fixture("tabix/vcf_file.vcf"))?;
    let expected_records = vcf_record_lines(&source)
        .into_iter()
        .take(2)
        .collect::<Vec<_>>();

    cleanup(&[bcf_path, csi_path]);

    assert!(actual.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    assert_eq!(vcf_record_lines(&actual), expected_records);

    Ok(())
}

#[test]
fn ports_test_view_bcf_index_output() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("view-write-index.bcf");
    let csi_path = temp_path("view-write-index.bcf.csi");

    write_bcf_from_vcf_path_with_csi(fixture("tabix/vcf_file.vcf"), &bcf_path, &csi_path)?;

    let actual = view_bcf_regions_as_vcf_text_from_path(&bcf_path, "1:3000150-3258501")?;
    let expected =
        filter_vcf_text_by_region_from_path(fixture("tabix/vcf_file.vcf"), "1:3000150-3258501")?;

    cleanup(&[bcf_path, csi_path]);

    assert_eq!(vcf_record_lines(&actual), vcf_record_lines(&expected));

    Ok(())
}

#[test]
fn ports_test_view_variant_record_limit() -> Result<(), Box<dyn std::error::Error>> {
    let src = std::fs::read_to_string(fixture("tabix/vcf_file.vcf"))?;
    let expected = vcf_text_with_first_records(&src, 2);

    assert_eq!(
        view_vcf_text_from_path_with_limit(fixture("tabix/vcf_file.vcf"), Some(2))?,
        expected
    );
    assert_eq!(
        view_bcf_as_vcf_text_from_path_with_limit(fixture("tabix/vcf_file.bcf"), Some(2))?,
        expected
    );

    Ok(())
}

#[test]
fn ports_test_view_variant_benchmark_mode_counts() -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(
        benchmark_vcf_view_from_path(fixture("tabix/vcf_file.vcf"))?,
        vcf_record_line_count(&view_vcf_text_from_path_with_limit(
            fixture("tabix/vcf_file.vcf"),
            None
        )?)
    );
    assert_eq!(
        benchmark_bcf_view_from_path(fixture("tabix/vcf_file.bcf"))?,
        vcf_record_line_count(&view_bcf_as_vcf_text_from_path_with_limit(
            fixture("tabix/vcf_file.bcf"),
            None
        )?)
    );

    Ok(())
}

#[test]
fn ports_vcf_sweep_forward_backward_checksums() -> Result<(), Box<dyn std::error::Error>> {
    let mut sweep = read_vcf_sweep_from_path(fixture("index.vcf"))?;

    assert_eq!(sample_count(sweep.header()), 1);
    assert_eq!(sweep.len(), 621);

    let mut position_sum = 0_u64;
    while let Some(record) = sweep.forward() {
        position_sum += u64::try_from(usize::from(record.position()))?;
    }
    assert_eq!(position_sum, 3_638_024_113);

    let mut reverse_position_sum = 0_u64;
    while let Some(record) = sweep.backward() {
        reverse_position_sum += u64::try_from(usize::from(record.position()))?;
    }
    assert_eq!(reverse_position_sum, position_sum);

    let mut pl_sum = 0_i64;
    while let Some(record) = sweep.forward() {
        pl_sum += record.pl_sum();
    }
    assert_eq!(pl_sum, 42_014);

    let mut reverse_pl_sum = 0_i64;
    while let Some(record) = sweep.backward() {
        reverse_pl_sum += record.pl_sum();
    }
    assert_eq!(reverse_pl_sum, pl_sum);

    Ok(())
}

#[test]
fn ports_bcf_sweep_forward_backward_positions() -> Result<(), Box<dyn std::error::Error>> {
    let mut sweep = read_bcf_sweep_from_path(fixture("tabix/vcf_file.bcf"))?;

    assert_eq!(sample_count(sweep.header()), 2);
    assert_eq!(sweep.len(), 15);

    let mut positions = Vec::new();
    while let Some(record) = sweep.forward() {
        positions.push(usize::from(record.position()));
        assert_eq!(record.pl_sum(), 0);
    }

    assert_eq!(
        positions,
        [
            3_000_150, 3_000_151, 3_062_915, 3_062_915, 3_106_154, 3_106_154, 3_157_410, 3_162_006,
            3_177_144, 3_177_144, 3_184_885, 3_199_812, 3_212_016, 3_258_448, 3_258_501
        ]
    );

    let mut reverse_positions = Vec::new();
    while let Some(record) = sweep.backward() {
        reverse_positions.push(usize::from(record.position()));
    }

    assert_eq!(
        reverse_positions,
        positions.into_iter().rev().collect::<Vec<_>>()
    );

    Ok(())
}

#[test]
fn ports_bcf_sr_no_index_summary() -> Result<(), Box<dyn std::error::Error>> {
    let paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.b.vcf"),
        fixture("bcf-sr/merge.noidx.c.vcf"),
    ];
    let actual = synced_vcf_summary_no_index_from_paths(&paths)?;
    let expected = std::fs::read_to_string(fixture("bcf-sr/merge.noidx.abc.expected.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_bcf_sr_reader_addition_summary() -> Result<(), Box<dyn std::error::Error>> {
    let paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.b.vcf"),
        fixture("bcf-sr/merge.noidx.c.vcf"),
    ];
    let expected = synced_vcf_summary_no_index_from_paths(&paths)?;
    let readers = paths
        .iter()
        .map(File::open)
        .map(|result| result.map(BufReader::new))
        .collect::<Result<Vec<_>, _>>()?;

    let actual = synced_vcf_summary_no_index_from_readers(readers)?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_bcf_sr_no_index_vcf_output() -> Result<(), Box<dyn std::error::Error>> {
    let paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.b.vcf"),
        fixture("bcf-sr/merge.noidx.c.vcf"),
    ];
    let actual = synced_vcf_output_no_index_from_paths(&paths)?;
    let expected = concat!(
        "##fileformat=VCFv4.3\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##contig=<ID=1,assembly=b37,length=249250621>\n",
        "##contig=<ID=2,assembly=b37,length=249250621>\n",
        "##contig=<ID=3,assembly=b37,length=198022430>\n",
        "##contig=<ID=4,assembly=b37,length=191154276>\n",
        "##reference=file:///lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n",
        "1\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/0\n",
        "1\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "1\t3000150\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "1\t3000151\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "1\t3000151\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "2\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/0\n",
        "2\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "2\t3000150\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "2\t3000151\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "2\t3000151\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "3\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/0\n",
        "3\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "3\t3000150\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "3\t3000151\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "3\t3000151\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "4\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/0\n",
        "4\t3000150\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "4\t3000150\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
        "4\t3000151\t.\tC\tT\t.\t.\t.\tGT\t0/1\n",
        "4\t3000151\t.\tC\tT\t.\t.\t.\tGT\t1/1\n",
    );

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_bcf_sr_no_index_bcf_output() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("bcf-sr-no-index-output.bcf");
    let paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.b.vcf"),
        fixture("bcf-sr/merge.noidx.c.vcf"),
    ];

    let encoded = synced_bcf_output_no_index_from_paths(&paths)?;
    std::fs::write(&bcf_path, encoded)?;

    let actual = view_bcf_as_vcf_text_from_path_with_limit(&bcf_path, None)?;
    let expected = synced_vcf_output_no_index_from_paths(&paths)?;

    cleanup(&[bcf_path]);

    assert!(actual.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"));
    assert_eq!(vcf_record_lines(&actual), vcf_record_lines(&expected));

    Ok(())
}

#[test]
fn ports_bcf_sr_no_index_order_errors() {
    let header_order_paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.hdr_order.vcf"),
    ];
    let record_order_paths = [
        fixture("bcf-sr/merge.noidx.a.vcf"),
        fixture("bcf-sr/merge.noidx.rec_order.vcf"),
    ];
    let reverse_record_order_paths = [
        fixture("bcf-sr/merge.noidx.rec_order.vcf"),
        fixture("bcf-sr/merge.noidx.a.vcf"),
    ];

    assert!(synced_vcf_summary_no_index_from_paths(&header_order_paths).is_err());
    assert!(synced_vcf_summary_no_index_from_paths(&record_order_paths).is_err());
    assert!(synced_vcf_summary_no_index_from_paths(&reverse_record_order_paths).is_err());
}

#[test]
fn ports_bcf_sr_weird_chromosome_region_queries() -> Result<(), Box<dyn std::error::Error>> {
    for (region, expected_path) in bcf_sr_weird_chromosome_cases() {
        let actual =
            filter_vcf_text_by_region_from_path(fixture("bcf-sr/weird-chr-names.vcf"), region)?;
        let expected = std::fs::read_to_string(fixture(expected_path))?;
        assert_eq!(actual, expected, "region: {region}");
    }

    assert!(
        filter_vcf_text_by_region_from_path(fixture("bcf-sr/weird-chr-names.vcf"), "{1:1-1}-2")
            .is_err()
    );

    Ok(())
}

#[test]
fn ports_bcf_sr_weird_chromosome_target_queries() -> Result<(), Box<dyn std::error::Error>> {
    for (target, expected_path) in bcf_sr_weird_chromosome_cases() {
        let actual =
            filter_vcf_text_by_target_from_path(fixture("bcf-sr/weird-chr-names.vcf"), target)?;
        let expected = std::fs::read_to_string(fixture(expected_path))?;
        assert_eq!(actual, expected, "target: {target}");
    }

    assert!(
        filter_vcf_text_by_target_from_path(fixture("bcf-sr/weird-chr-names.vcf"), "{1:1-1}-2")
            .is_err()
    );

    Ok(())
}

#[test]
fn ports_bcf_sr_indexed_bcf_weird_chromosome_region_queries()
-> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("weird-chr-names.bcf");
    let csi_path = temp_path("weird-chr-names.bcf.csi");

    let encoded = write_bcf_from_vcf_path(fixture("bcf-sr/weird-chr-names.vcf"), Vec::new())?;
    std::fs::write(&bcf_path, encoded)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    for (region, expected_path) in bcf_sr_weird_chromosome_cases() {
        let actual = view_bcf_regions_as_vcf_text_from_path(&bcf_path, region)?;
        let expected = std::fs::read_to_string(fixture(expected_path))?;
        assert_eq!(actual, expected, "region: {region}");
    }

    assert!(view_bcf_regions_as_vcf_text_from_path(&bcf_path, "{1:1-1}-2").is_err());

    cleanup(&[bcf_path, csi_path]);

    Ok(())
}

#[test]
fn ports_bcf_sr_indexed_bcf_weird_chromosome_target_queries()
-> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("weird-chr-names-targets.bcf");
    let csi_path = temp_path("weird-chr-names-targets.bcf.csi");

    let encoded = write_bcf_from_vcf_path(fixture("bcf-sr/weird-chr-names.vcf"), Vec::new())?;
    std::fs::write(&bcf_path, encoded)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    for (target, expected_path) in bcf_sr_weird_chromosome_cases() {
        let actual = view_bcf_targets_as_vcf_text_from_path(&bcf_path, target)?;
        let expected = std::fs::read_to_string(fixture(expected_path))?;
        assert_eq!(actual, expected, "target: {target}");
    }

    assert!(view_bcf_targets_as_vcf_text_from_path(&bcf_path, "{1:1-1}-2").is_err());

    cleanup(&[bcf_path, csi_path]);

    Ok(())
}

#[test]
fn writes_bcf_round_trip_records_matching_htslib_vcf_fixture()
-> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("write-round-trip.bcf");
    let encoded = write_bcf_from_vcf_path(fixture("tabix/vcf_file.vcf"), Vec::new())?;
    std::fs::write(&bcf_path, encoded)?;

    let actual = view_bcf_as_vcf_text_from_path_with_limit(&bcf_path, None)?;
    let expected = std::fs::read_to_string(fixture("tabix/vcf_file.vcf"))?;

    cleanup(&[bcf_path]);

    assert_eq!(vcf_record_lines(&actual), vcf_record_lines(&expected));

    Ok(())
}

#[test]
fn ports_bcf_translate_synthetic_header_and_record() -> Result<(), Box<dyn std::error::Error>> {
    let actual = translated_bcf_record_fixture_vcf_text()?;
    let expected = std::fs::read_to_string(fixture("test-bcf-translate.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_vcf_header_get_and_remove_cases() {
    let mut header = parse_header(&[
        r#"##FILTER=<ID=Flt,Description="Unused FILTER">"#,
        r#"##INFO=<ID=UI,Number=1,Type=Integer,Description="Unused INFO">"#,
        r#"##FORMAT=<ID=UF,Number=1,Type=Integer,Description="Unused FORMAT">"#,
        "##contig=<ID=Unused,length=1>",
    ]);

    assert!(header_has_filter(&header, "Flt"));
    assert!(remove_header_filter(&mut header, "Flt"));
    assert!(!header_has_filter(&header, "Flt"));

    assert!(header_has_info(&header, "UI"));
    assert!(remove_header_info(&mut header, "UI"));
    assert!(!header_has_info(&header, "UI"));

    assert!(header_has_format(&header, "UF"));
    assert!(remove_header_format(&mut header, "UF"));
    assert!(!header_has_format(&header, "UF"));

    assert!(header_has_contig(&header, "Unused"));
    assert!(remove_header_contig(&mut header, "Unused"));
    assert!(!header_has_contig(&header, "Unused"));

    let mut header = parse_header(&[r#"##unused=<ID=BB,Description="Unused generic with ID">"#]);
    assert!(header_has_other_record_id(&header, "unused", "BB"));
    assert!(remove_header_other_record_id(&mut header, "unused", "BB"));
    assert!(!header_has_other_record_id(&header, "unused", "BB"));

    let mut header = parse_header(&["##unused=unformatted text 1", "##unused=unformatted text 2"]);
    assert!(header_has_unstructured_other_record(
        &header,
        "unused",
        "unformatted text 1"
    ));
    assert!(remove_header_other_records(&mut header, "unused"));
    assert!(!header_has_unstructured_other_record(
        &header,
        "unused",
        "unformatted text 1"
    ));
}

#[test]
fn adapts_vcf_header_ids() {
    let mut header = parse_header(&[
        r#"##FILTER=<ID=Flt,Description="Unused FILTER">"#,
        r#"##INFO=<ID=UI,Number=1,Type=Integer,Description="Unused INFO">"#,
        r#"##FORMAT=<ID=UF,Number=1,Type=Integer,Description="Unused FORMAT">"#,
        "##contig=<ID=Unused,length=1>",
        r#"##unused=<ID=BB,Description="Unused generic with ID">"#,
    ]);
    let ids = [
        VcfHeaderId::Filter("Flt".into()),
        VcfHeaderId::Info("UI".into()),
        VcfHeaderId::Format("UF".into()),
        VcfHeaderId::Contig("Unused".into()),
        VcfHeaderId::Other {
            key: "unused".into(),
            id: "BB".into(),
        },
    ];

    for id in ids {
        assert!(vcf_header_has_id(&header, &id));
        assert!(vcf_header_remove_id(&mut header, &id));
        assert!(!vcf_header_has_id(&header, &id));
    }
}

#[test]
fn ports_vcf_open_mode_extension_cases() {
    assert_eq!(vcf_open_mode_suffix("mode1.bcf").unwrap(), "b");
    assert_eq!(vcf_open_mode_suffix("mode1.vcf").unwrap(), "");
    assert_eq!(vcf_open_mode_suffix("mode1.vcf.gz").unwrap(), "z");
    assert_eq!(vcf_open_mode_suffix("mode1.vcf.bgz").unwrap(), "z");

    assert!(vcf_open_mode_suffix("mode1.xcf").is_err());
    assert!(vcf_open_mode_suffix("mode1.vcf.gbz").is_err());
    assert!(vcf_open_mode_suffix("mode1.bvcf.bgz").is_err());
}

#[test]
fn ports_vcf_header_number_classification_cases() {
    let supported_header = parse_header(&[
        r#"##INFO=<ID=FIXED_1_INFO,Number=1,Type=Integer,Description="Fixed number 1">"#,
        r#"##INFO=<ID=FIXED_4_INFO,Number=4,Type=Float,Description="Fixed number 4">"#,
        r#"##INFO=<ID=VL_DOT_INFO,Number=.,Type=Integer,Description="Variable number">"#,
        r#"##INFO=<ID=VL_A_INFO,Number=A,Type=Integer,Description="One value for each ALT allele">"#,
        r#"##INFO=<ID=VL_G_INFO,Number=G,Type=Integer,Description="One value for each possible genotype">"#,
        r#"##INFO=<ID=VL_R_INFO,Number=R,Type=Integer,Description="One value for each allele including REF">"#,
        r#"##FORMAT=<ID=FIXED_1_FMT,Number=1,Type=String,Description="Fixed number 1">"#,
        r#"##FORMAT=<ID=FIXED_4_FMT,Number=4,Type=String,Description="Fixed number 4">"#,
        r#"##FORMAT=<ID=VL_DOT_FMT,Number=.,Type=String,Description="Variable number">"#,
        r#"##FORMAT=<ID=VL_A_FMT,Number=A,Type=Integer,Description="One value for each ALT allele">"#,
        r#"##FORMAT=<ID=VL_G_FMT,Number=G,Type=Integer,Description="One value for each possible genotype">"#,
        r#"##FORMAT=<ID=VL_R_FMT,Number=R,Type=Integer,Description="One value for each allele including REF">"#,
    ]);
    let extended_lines = [
        r#"##FORMAT=<ID=VL_P_FMT,Number=P,Type=String,Description="One value for each allele value defined in GT">"#,
        r#"##FORMAT=<ID=VL_LA_FMT,Number=LA,Type=Integer,Description="One value for each local ALT allele">"#,
        r#"##FORMAT=<ID=VL_LG_FMT,Number=LG,Type=Integer,Description="One value for each local genotype">"#,
        r#"##FORMAT=<ID=VL_LR_FMT,Number=LR,Type=Integer,Description="One value for each local allele including REF">"#,
        r#"##FORMAT=<ID=VL_M_FMT,Number=M,Type=Integer,Description="One value for each posible base modification of the given type">"#,
    ];

    assert_eq!(
        info_header_number(&supported_header, "FIXED_1_INFO"),
        Some(HeaderNumber::Fixed(1))
    );
    assert_eq!(
        info_header_number(&supported_header, "FIXED_4_INFO"),
        Some(HeaderNumber::Fixed(4))
    );
    assert_eq!(
        info_header_number(&supported_header, "VL_DOT_INFO"),
        Some(HeaderNumber::Variable)
    );
    assert_eq!(
        info_header_number(&supported_header, "VL_A_INFO"),
        Some(HeaderNumber::AlternateBases)
    );
    assert_eq!(
        info_header_number(&supported_header, "VL_G_INFO"),
        Some(HeaderNumber::Genotypes)
    );
    assert_eq!(
        info_header_number(&supported_header, "VL_R_INFO"),
        Some(HeaderNumber::ReferenceAlternateBases)
    );

    assert_eq!(
        format_header_number(&supported_header, "FIXED_1_FMT"),
        Some(HeaderNumber::Fixed(1))
    );
    assert_eq!(
        format_header_number(&supported_header, "FIXED_4_FMT"),
        Some(HeaderNumber::Fixed(4))
    );
    assert_eq!(
        format_header_number(&supported_header, "VL_DOT_FMT"),
        Some(HeaderNumber::Variable)
    );
    assert_eq!(
        format_header_number(&supported_header, "VL_A_FMT"),
        Some(HeaderNumber::AlternateBases)
    );
    assert_eq!(
        format_header_number(&supported_header, "VL_G_FMT"),
        Some(HeaderNumber::Genotypes)
    );
    assert_eq!(
        format_header_number(&supported_header, "VL_R_FMT"),
        Some(HeaderNumber::ReferenceAlternateBases)
    );

    let expected = [
        ("VL_P_FMT", HeaderNumber::Ploidy),
        ("VL_LA_FMT", HeaderNumber::LocalAlternateBases),
        ("VL_LG_FMT", HeaderNumber::LocalGenotypes),
        ("VL_LR_FMT", HeaderNumber::LocalReferenceAlternateBases),
        ("VL_M_FMT", HeaderNumber::BaseModifications),
    ];

    for (line, (expected_id, expected_number)) in extended_lines.iter().zip(expected) {
        let (_, id, number) = parse_numbered_header_record(line)
            .unwrap()
            .expect("FORMAT header record");

        assert_eq!(id, expected_id);
        assert_eq!(number, expected_number);
    }
}

#[test]
fn ports_vcf_record_rlen_table_cases() -> Result<(), Box<dyn std::error::Error>> {
    let body = concat!(
        "##reference=file://tmp\n",
        r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
        "\n",
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="end">"#,
        "\n",
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">"#,
        "\n",
        r#"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="svlen">"#,
        "\n",
        r#"##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number">"#,
        "\n",
        r#"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="svclaim">"#,
        "\n",
        r#"##FORMAT=<ID=LEN,Number=1,Type=Integer,Description="fmt len">"#,
        "\n",
        "##contig=<ID=1,Length=40>\n",
        r#"##ALT=<ID=INS,Description="INS">"#,
        "\n",
        r#"##ALT=<ID=DEL,Description="DEL">"#,
        "\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n",
        "1\t4310\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0|1\n",
        "1\t4311\t.\tC\tCT\t.\t.\t.\tGT\t0/0\t1/0\n",
        "1\t4312\t.\tTC\tT\t213.73\t.\t.\tGT\t0/1\t0|0\n",
        "1\t4314\t.\tG\t<INS>\t213.73\t.\tSVLEN=10;SVCLAIM=J\tGT\t0/1\t0|0\n",
        "1\t4315\t.\tG\t<DEL>\t213.73\t.\tSVLEN=-10;SVCLAIM=D\tGT\t0/1\t0|0\n",
        "1\t4326\t.\tG\t<INS>\t213.73\t.\tEND=4326;SVLEN=10;SVCLAIM=J\tGT\t0/1\t0|0\n",
        "1\t4327\t.\tG\t<DEL>\t213.73\t.\tEND=4337;SVLEN=-10;SVCLAIM=J\tGT\t0/1\t0|0\n",
        "1\t4338\t.\tG\t<*>\t213.73\t.\tEND=4342;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:8\n",
        "1\t4353\t.\tG\t<*>\t213.73\t.\tEND=4357;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n",
        "1\t4363\t.\tG\t<*>\t213.73\t.\tEND=4367;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n",
        "1\t4370\t.\tG\t<INS>,<*>\t213.73\t.\tEND=4371;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n",
        "1\t4378\t.\tG\t<DEL>,<INS>,<*>\t213.73\t.\tEND=4379;SVLEN=3,5,.;SVCLAIM=D,J,.\tGT:LEN\t0/1:7\t0|0:.\n",
        "1\t4385\t.\tG\tT,<DEL>\t213.73\t.\tEND=4387;SVLEN=.,180\tGT\t0/1\t0|0\n",
        "1\t4585\t.\tG\tT,<DEL:ME>\t213.73\t.\tEND=4587;SVLEN=.,180\tGT\t0/1\t0|0\n",
        "1\t4685\t.\tG\t<DUP>,<DUP>\t213.73\t.\tEND=4687;SVLEN=10,10\tGT\t0/1\t0|0\n",
        "1\t4705\t.\tG\t<CNV>\t213.73\t.\tEND=4707;SVLEN=11;CN=2\tGT\t0/1\t0|0\n",
        "1\t4725\t.\tG\t<CNV:TR>\t213.73\t.\tEND=4727;SVLEN=12;CN=1.5\tGT\t0/1\t0|0\n",
        "1\t4745\t.\tG\t<INV>\t213.73\t.\tEND=4747;SVLEN=10\tGT\t0/1\t0|0\n",
        "1\t4885\t.\tG\tT,<*>\t213.73\t.\tEND=4887\tGT:LEN\t0/1:190\t0|0:.\n",
        "1\t5885\t.\tG\tT\t213.73\t.\tEND=5887;SVLEN=8;SVCLAIM=.\tGT:LEN\t0/1:.\t0|0:10\n",
    );
    let expected = [
        1, 1, 2, 1, 11, 1, 11, 8, 7, 7, 7, 7, 181, 181, 11, 12, 13, 11, 190, 3,
    ];

    for version in ["VCFv4.3", "VCFv4.4", "VCFv4.5"] {
        let vcf = format!("##fileformat={version}\n{body}");
        let spans = htslib_variant_spans_from_vcf(std::io::Cursor::new(vcf))?;

        assert_eq!(spans, expected, "version: {version}");
    }

    Ok(())
}

#[test]
fn ports_vcf_record_rlen_update_cases() -> Result<(), Box<dyn std::error::Error>> {
    let mut record = "1\t4310\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0|1".to_string();

    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 1);

    record = update_vcf_line_alleles(&record, &["G", "AT"])?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 1);

    record = update_vcf_line_alleles(&record, &["GC", "A"])?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 2);

    record = update_vcf_line_alleles(&record, &["G", "<*>"])?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 1);

    record = update_vcf_line_info_i32(&record, "END", Some(&[4323]))?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 14);

    record = update_vcf_line_format_i32(&record, "LEN", Some(&[1, 15]))?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 15);

    record = update_vcf_line_info_i32(&record, "END", None)?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 15);

    record = update_vcf_line_format_i32(&record, "LEN", None)?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 1);

    record = update_vcf_line_alleles(&record, &["G", "T", "<DEL>"])?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 1);

    record = update_vcf_line_info_i32(&record, "SVLEN", Some(&[0, -5]))?;
    assert_eq!(htslib_variant_span_from_vcf_line(&record)?, 6);

    let copied_record = record.clone();
    assert_eq!(htslib_variant_span_from_vcf_line(&copied_record)?, 6);

    Ok(())
}

#[test]
fn adapts_vcf_record_typed_values_genotypes_and_mutation() -> Result<(), Box<dyn std::error::Error>>
{
    let mut record = VcfRecordAdapter::new(
        "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;NEG=-127;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.",
    )?;

    assert_eq!(record.info_i32_values("NEG")?, Some(vec![Some(-127)]));
    assert_eq!(
        record.format_i32_values("GQ")?,
        Some(vec![vec![Some(48)], vec![Some(48)], vec![Some(43)]])
    );
    assert_eq!(
        record.genotypes()?,
        Some(vec!["0|0".into(), "1|0".into(), "1/1".into()])
    );
    assert_eq!(record.htslib_span()?, 1);

    record.set_alleles(&["GC", "A"])?;
    assert_eq!(record.htslib_span()?, 2);
    record.set_info_i32("END", Some(&[14375]))?;
    assert_eq!(record.htslib_span()?, 6);
    record.set_format_i32("DP", Some(&[10, 11, 12]))?;
    assert_eq!(
        record.format_i32_values("DP")?,
        Some(vec![vec![Some(10)], vec![Some(11)], vec![Some(12)]])
    );

    record.remove_alleles(&[1])?;
    assert_eq!(record.as_str().split('\t').nth(4), Some("."));

    Ok(())
}

#[test]
fn ports_vcf_invalid_end_tag_rlen_cases() -> Result<(), Box<dyn std::error::Error>> {
    let vcf = concat!(
        "##fileformat=VCFv4.1\n",
        "##contig=<ID=X,length=155270560>\n",
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">"#,
        "\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "X\t86470037\trs59780433a\tTTTCA\tTGGTT,T\t.\t.\tEND=85725113\n",
        "X\t86470038\trs59780433b\tT\tTGGTT,T\t.\t.\tEND=86470047\n",
    );

    assert_eq!(
        htslib_variant_spans_from_vcf(std::io::Cursor::new(vcf))?,
        [5, 10]
    );

    Ok(())
}

#[test]
fn ports_vcf_get_info_typed_values() -> Result<(), Box<dyn std::error::Error>> {
    let records = [
        (
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;NEG=-127;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.",
            vec![Some(0.5)],
            vec![Some(-127)],
        ),
        (
            "20\t1110696\t.\tA\tG,T\t67\t.\tNS=2;DP=10;NEG=-128;AF=0.333,.;AA=T;DB\tGT\t2\t1\t./.",
            vec![Some(0.333), None],
            vec![Some(-128)],
        ),
    ];

    for (record, expected_af, expected_neg) in records {
        assert_eq!(
            vcf_info_float_values_from_line(record, "AF")?,
            Some(expected_af)
        );
        assert_eq!(
            vcf_info_i32_values_from_line(record, "NEG")?,
            Some(expected_neg)
        );
    }

    Ok(())
}

#[test]
fn ports_vcf_get_format_integer_values() -> Result<(), Box<dyn std::error::Error>> {
    let record = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;NEG=-127;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";

    assert_eq!(
        vcf_format_i32_values_from_line(record, "GQ")?,
        Some(vec![vec![Some(48)], vec![Some(48)], vec![Some(43)]])
    );
    assert_eq!(
        vcf_format_i32_values_from_line(record, "DP")?,
        Some(vec![vec![Some(1)], vec![Some(8)], vec![Some(5)]])
    );
    assert_eq!(
        vcf_format_i32_values_from_line(record, "HQ")?,
        Some(vec![
            vec![Some(51), Some(51)],
            vec![Some(51), Some(51)],
            vec![None, None],
        ])
    );

    Ok(())
}

#[test]
fn ports_vcf_api_record_serialization_output() -> Result<(), Box<dyn std::error::Error>> {
    let actual = test_vcf_api_record_serialization_text();
    let expected = std::fs::read_to_string(fixture("test-vcf-api.out"))?;

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_bcf_get_format_float_vector_end_values() {
    let values = normalize_bcf_float_vector(&[
        BcfFloatValue::Missing,
        BcfFloatValue::Value(47.11),
        BcfFloatValue::VectorEnd,
        BcfFloatValue::Value(-1.2e-13),
    ]);

    assert_eq!(
        values,
        [
            BcfFloatValue::Missing,
            BcfFloatValue::Value(47.11),
            BcfFloatValue::VectorEnd,
            BcfFloatValue::VectorEnd,
        ]
    );
}

#[test]
fn ports_bcf_remove_allele_set_cases() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        (
            "5\t110285\t.\tT\tC,<*>\t.\tPASS\tAC=1,0;AD=6,5,0;AF=0.99,0.01;VL_A_STR_INFO=alt_c,alt_nonref;VL_R_STR_INFO=ref,alt_c,alt_nonref\tGT:AD:EC:PL:VL_A_STR_FMT:VL_G_STR_FMT:VL_R_STR_FMT\t.:.:.:.:.:.:.\t0/1:6,5,0:4,0:114,0,15,35,73,113:alt_c,alt_nonref:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_nonref\t.:.:.:.:.:.:.",
            &[2][..],
            "5\t110285\t.\tT\tC\t.\tPASS\tAC=1;AD=6,5;AF=0.99;VL_A_STR_INFO=alt_c;VL_R_STR_INFO=ref,alt_c\tGT:AD:EC:PL:VL_A_STR_FMT:VL_G_STR_FMT:VL_R_STR_FMT\t.:.:.:.:.:.:.\t0/1:6,5:4:114,0,15:alt_c:gt_00,gt_01,gt_11:ref,alt_c\t.:.:.:.:.:.:.",
        ),
        (
            "5\t110290\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
            &[2][..],
            "5\t110290\t.\tT\tC\t.\tPASS\tAC=90;AD=6,5;AF=0.009;VL_A_STR_INFO=alt_c;VL_R_STR_INFO=ref,alt_c\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1:3,2:44:114,0,15:alt_c:gt_00,gt_01,gt_11:ref,alt_c\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
        ),
        (
            "5\t110350\t.\tT\t<INS>,<INS>\t.\tPASS\tIMPRECISE;SVLEN=100,200;CIEND=-50,50,-25,25;CIPOS=-10,10,-20,20\tGT\t0/1\t0/1\t0/1",
            &[2][..],
            "5\t110350\t.\tT\t<INS>\t.\tPASS\tIMPRECISE;SVLEN=100;CIEND=-50,50;CIPOS=-10,10\tGT\t0/1\t0/1\t0/1",
        ),
        (
            "5\t110500\t.\tT\t<CNV>,<CNV>\t.\tPASS\tIMPRECISE;SVLEN=50,100;CILEN=0,25,-25,25;CN=2,4;CICN=-0.5,1,-1.5,1.5\tGT\t0/1\t0/1\t0/1",
            &[2][..],
            "5\t110500\t.\tT\t<CNV>\t.\tPASS\tIMPRECISE;SVLEN=50;CILEN=0,25;CN=2;CICN=-0.5,1\tGT\t0/1\t0/1\t0/1",
        ),
        (
            "5\t110700\t.\tA\t<INS:ME>,<INS:ME>\t.\tPASS\tMEINFO=AluY,1,260,+,FLAM_C,1,110,-;METRANS=1,94820,95080,+,1,129678,129788,-\tGT\t0/1\t0/1\t0/1",
            &[2][..],
            "5\t110700\t.\tA\t<INS:ME>\t.\tPASS\tMEINFO=AluY,1,260,+;METRANS=1,94820,95080,+\tGT\t0/1\t0/1\t0/1",
        ),
        (
            "5\t112000\t.\tC\t<CNV:TR>,<CNV:TR>\t.\tPASS\tRN=2,1;RUS=CAG,TTG,CA;RUL=3,3,2;RB=12,6,6;RUC=4,2,3;RUB=3,3,3,3,3,3,2,2,2;SVLEN=18,6",
            &[2][..],
            "5\t112000\t.\tC\t<CNV:TR>\t.\tPASS\tRN=2;RUS=CAG,TTG;RUL=3,3;RB=12,6;RUC=4,2;RUB=3,3,3,3,3,3;SVLEN=18",
        ),
        (
            "5\t113000\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
            &[1][..],
            "5\t113000\t.\tT\tA\t.\tPASS\tAC=1;AD=6,6;AF=0.0001;VL_A_STR_INFO=alt_a;VL_R_STR_INFO=ref,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/.:1:3,0:27:114,35,113:alt_a:gt_00,gt_02,gt_22:ref,alt_a\t./.:.:0:.:110:.:gt_00:ref",
        ),
        (
            "5\t114000\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
            &[1, 2][..],
            "5\t114000\t.\tT\t.\t.\tPASS\tAD=6;VL_R_STR_INFO=ref\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/.:.:3:.:114:.:gt_00:ref\t./.:.:0:.:110:.:gt_00:ref",
        ),
        (
            "5\t115000\t.\tC\t<CNV:TR>,<CNV:TR>\t.\tPASS\tRN=2,1;RUS=CAG,TTG,CA;RUL=3,3,2;RB=12,6,6;RUC=4,2,3;RUB=3,3,3,3,3,3,2,2,2;SVLEN=18,6",
            &[1, 2][..],
            "5\t115000\t.\tC\t.\t.\tPASS\t.",
        ),
    ];

    for (input, remove, expected) in cases {
        let actual = remove_vcf_allele_set_from_line(input, remove)?;
        assert_eq!(actual, expected, "input: {input}");
    }

    Ok(())
}

#[test]
fn queries_bcf_records_from_htslib_tabix_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path =
        std::env::temp_dir().join(format!("htslib-rs-variant-bcf-{}.bcf", std::process::id()));
    let csi_path = std::env::temp_dir().join(format!(
        "htslib-rs-variant-bcf-{}.bcf.csi",
        std::process::id()
    ));

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_bcf_records_from_path(&bcf_path, &region)?;
    let expected_record_count = std::fs::read_to_string(fixture("tabix/vcf_file.1.3000151.out"))?
        .lines()
        .count();

    std::fs::remove_file(bcf_path)?;
    std::fs::remove_file(csi_path)?;

    assert_eq!(records.len(), expected_record_count);

    Ok(())
}

#[test]
fn iterates_bcf_records_from_htslib_tabix_fixture() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("iter.bcf");
    let csi_path = temp_path("iter.bcf.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_bcf_records_from_path(&bcf_path, &region)?;
    let iter_count = iter_bcf_records_from_path(&bcf_path, &region)?.count();

    cleanup(&[bcf_path, csi_path]);

    assert_eq!(iter_count, records.len());
    assert!(iter_count > 0);

    Ok(())
}

#[test]
fn ports_test_vcf_api_bcf_iterator_creation() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("iterator.bcf");
    let csi_path = temp_path("iterator.bcf.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    write_bcf_csi_from_path(&bcf_path, &csi_path)?;

    let region = "1:3000151-3000151".parse()?;
    let by_name = query_bcf_records_from_path(&bcf_path, &region)?;
    let by_id = query_bcf_records_by_reference_id_from_path(&bcf_path, 0, 3_000_150, 3_000_151)?;

    cleanup(&[bcf_path, csi_path]);

    assert_eq!(by_id.len(), by_name.len());
    assert_eq!(by_id.len(), 1);

    Ok(())
}

#[test]
fn queries_bcf_with_htslib_index_lookup_fallbacks() -> Result<(), Box<dyn std::error::Error>> {
    let bcf_path = temp_path("lookup.bcf");
    let replaced_csi_path = temp_path("lookup.csi");
    let explicit_csi_path = temp_path("bcf-custom-lookup.csi");

    std::fs::copy(fixture("tabix/vcf_file.bcf"), &bcf_path)?;
    let index = build_bcf_csi_from_path(&bcf_path)?;
    write_csi(&replaced_csi_path, &index)?;
    write_csi(&explicit_csi_path, &index)?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_bcf_records_from_path(&bcf_path, &region)?;
    let explicit_src = format!(
        "{}##idx##{}",
        bcf_path.display(),
        explicit_csi_path.display()
    );
    let explicit_records = query_bcf_records_from_path(explicit_src, &region)?;

    cleanup(&[bcf_path, replaced_csi_path, explicit_csi_path]);

    assert_eq!(records.len(), 1);
    assert_eq!(explicit_records.len(), 1);

    Ok(())
}
