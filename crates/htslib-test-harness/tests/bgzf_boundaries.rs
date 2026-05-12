use std::path::{Path, PathBuf};

use htslib_rs::alignment_compat::count_bam_records_from_path;

fn fixture(path: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib/test")
        .join(path)
}

#[test]
fn reads_bam_records_split_across_bgzf_blocks() -> Result<(), Box<dyn std::error::Error>> {
    for path in [
        "bgzf_boundaries/bgzf_boundaries1.bam",
        "bgzf_boundaries/bgzf_boundaries2.bam",
        "bgzf_boundaries/bgzf_boundaries3.bam",
    ] {
        assert_eq!(count_bam_records_from_path(fixture(path))?, 1, "{path}");
    }

    Ok(())
}
