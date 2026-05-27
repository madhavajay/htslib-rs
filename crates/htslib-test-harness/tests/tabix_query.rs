use std::{
    fs,
    io::BufReader,
    path::{Path, PathBuf},
};

use htslib_rs::tabix_compat::{
    TextFormat, query_csi_records_from_path, query_records_from_path,
    query_records_from_path_separate_regions, query_vcf_csi_records_from_path,
    query_vcf_records_from_associated_path, write_bed_bgzf_and_index, write_bgzf_and_csi,
    write_bgzf_and_index,
};

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("repos/htslib/test/tabix")
        .join(path)
}

fn htslib_fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("repos/htslib/test")
        .join(path)
}

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("htslib-rs-{name}-{}", std::process::id()))
}

fn cleanup(paths: &[PathBuf]) {
    for path in paths {
        let _ = fs::remove_file(path);
    }
}

#[test]
fn queries_bed_records_with_tabix_index() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("bed-file.bed.gz");
    let tbi_path = temp_path("bed-file.bed.gz.tbi");

    let bed = fs::File::open(fixture("bed_file.bed"))?;
    write_bed_bgzf_and_index(BufReader::new(bed), &bgzf_path, &tbi_path)?;

    let index = htslib_rs::tabix::fs::read(&tbi_path)?;
    let region = "Y:100200-100200".parse()?;
    let records = query_records_from_path(&bgzf_path, index, &region)?;

    let expected = fs::read_to_string(fixture("bed_file.Y.100200.out"))?;
    let actual = records.join("\n") + "\n";

    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn queries_bed_records_with_separate_region_headers() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("bed-file-separate.bed.gz");
    let tbi_path = temp_path("bed-file-separate.bed.gz.tbi");

    let bed = fs::File::open(fixture("bed_file.bed"))?;
    write_bed_bgzf_and_index(BufReader::new(bed), &bgzf_path, &tbi_path)?;

    let index = htslib_rs::tabix::fs::read(&tbi_path)?;
    let regions = [
        ("X:1100-1400", "X:1100-1400".parse()?),
        ("Y:100000-100550", "Y:100000-100550".parse()?),
        ("Z:100000-100005", "Z:100000-100005".parse()?),
    ];
    let actual = query_records_from_path_separate_regions(&bgzf_path, index, &regions)?;
    let expected = fs::read_to_string(fixture("bed_file.separate.out"))?;

    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn queries_vcf_records_with_tabix_index() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("vcf-file.vcf.gz");
    let tbi_path = temp_path("vcf-file.vcf.gz.tbi");

    let vcf = fs::File::open(fixture("vcf_file.vcf"))?;
    write_bgzf_and_index(BufReader::new(vcf), &bgzf_path, &tbi_path, TextFormat::Vcf)?;

    let index = htslib_rs::tabix::fs::read(&tbi_path)?;

    for (region, expected_name) in [
        ("1:3000151-3000151", "vcf_file.1.3000151.out"),
        ("2:3199812-3199812", "vcf_file.2.3199812.out"),
    ] {
        let records = query_records_from_path(&bgzf_path, index.clone(), &region.parse()?)?;
        let expected = fs::read_to_string(fixture(expected_name))?;
        let actual = records.join("\n") + "\n";

        assert_eq!(actual, expected);
    }

    cleanup(&[bgzf_path, tbi_path]);

    Ok(())
}

#[test]
fn queries_gff_records_with_tabix_index() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("gff-file.gff.gz");
    let tbi_path = temp_path("gff-file.gff.gz.tbi");

    let gff = fs::File::open(fixture("gff_file.gff"))?;
    write_bgzf_and_index(BufReader::new(gff), &bgzf_path, &tbi_path, TextFormat::Gff)?;

    let index = htslib_rs::tabix::fs::read(&tbi_path)?;
    let region = "X:2934832-2935190".parse()?;
    let records = query_records_from_path(&bgzf_path, index, &region)?;

    let expected = fs::read_to_string(fixture("gff_file.X.2934832.2935190.out"))?;
    let actual = records.join("\n") + "\n";

    cleanup(&[bgzf_path, tbi_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn queries_large_vcf_records_with_csi_index() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("large-chr.vcf.gz");
    let csi_path = temp_path("large-chr.vcf.gz.csi");

    let vcf = fs::File::open(fixture("large_chr.vcf"))?;
    write_bgzf_and_csi(BufReader::new(vcf), &bgzf_path, &csi_path, TextFormat::Vcf)?;

    let index = htslib_rs::csi::fs::read(&csi_path)?;
    let region = "chr20:1-2147483647".parse()?;
    let records = query_csi_records_from_path(&bgzf_path, index, &region)?;

    let expected = fs::read_to_string(fixture("large_chr.20.1.2147483647.out"))?;
    let actual = records.join("\n") + "\n";

    cleanup(&[bgzf_path, csi_path]);

    assert_eq!(actual, expected);

    Ok(())
}

#[test]
fn ports_test_pl_large_reference_vcf_csi_queries() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("longref-index.vcf.gz");
    let csi_path = std::env::temp_dir().join(format!(
        "htslib-rs-longref-index-{}.vcf.gz.csi",
        std::process::id()
    ));

    let vcf = fs::File::open(htslib_fixture("longrefs/index.vcf"))?;
    write_bgzf_and_csi(BufReader::new(vcf), &bgzf_path, &csi_path, TextFormat::Vcf)?;

    let index = htslib_rs::csi::fs::read(&csi_path)?;

    for (region, expected_path) in [
        ("1:10010000100-10010000105", "longrefs/index.expected1.vcf"),
        ("1:10010000120-10010000130", "longrefs/index.expected2.vcf"),
    ] {
        let records = query_vcf_csi_records_from_path(&bgzf_path, index.clone(), &region.parse()?)?;
        let expected = fs::read_to_string(htslib_fixture(expected_path))?;
        let actual = records.join("\n") + "\n";

        assert_eq!(actual, expected, "region: {region}");
    }

    let explicit_src = format!("{}##idx##{}", bgzf_path.display(), csi_path.display());
    let associated_records = query_vcf_records_from_associated_path(
        explicit_src,
        &"1:10010000120-10010000130".parse()?,
    )?;
    let expected = fs::read_to_string(htslib_fixture("longrefs/index.expected2.vcf"))?;

    assert_eq!(associated_records.join("\n") + "\n", expected);

    cleanup(&[bgzf_path, csi_path]);

    Ok(())
}

#[test]
fn queries_vcf_text_with_htslib_index_lookup_fallbacks() -> Result<(), Box<dyn std::error::Error>> {
    let bgzf_path = temp_path("associated-vcf-file.vcf.gz");
    let mut replaced_tbi_path = bgzf_path.clone();
    replaced_tbi_path.set_extension("tbi");

    let vcf = fs::File::open(fixture("vcf_file.vcf"))?;
    write_bgzf_and_index(
        BufReader::new(vcf),
        &bgzf_path,
        &replaced_tbi_path,
        TextFormat::Vcf,
    )?;

    let region = "1:3000151-3000151".parse()?;
    let records = query_vcf_records_from_associated_path(&bgzf_path, &region)?;

    let expected = fs::read_to_string(fixture("vcf_file.1.3000151.out"))?;
    let actual = records.join("\n") + "\n";

    cleanup(&[bgzf_path, replaced_tbi_path]);

    assert_eq!(actual, expected);

    Ok(())
}
