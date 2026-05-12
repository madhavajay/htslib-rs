use std::path::{Path, PathBuf};

use htslib_rs::log_compat::{
    LogLevel, check_htslib_log_message_style, hts_get_log_level, hts_log_enabled,
    hts_set_log_level, hts_verbose, set_hts_verbose,
};

fn htslib_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("htslib")
}

#[test]
fn ports_test_logging_pl_log_message_style() -> Result<(), Box<dyn std::error::Error>> {
    let report = check_htslib_log_message_style(htslib_root())?;

    assert!(report.files_scanned > 0);
    assert!(report.messages_checked > 0);
    assert!(report.issues.is_empty());

    Ok(())
}

#[test]
fn ports_hts_log_level_state() {
    let original = hts_verbose();

    hts_set_log_level(LogLevel::Warning);
    assert_eq!(hts_get_log_level(), LogLevel::Warning);
    assert_eq!(hts_verbose(), 3);
    assert!(hts_log_enabled(LogLevel::Error));
    assert!(hts_log_enabled(LogLevel::Warning));
    assert!(!hts_log_enabled(LogLevel::Info));

    hts_set_log_level(LogLevel::Off);
    assert!(!hts_log_enabled(LogLevel::Error));

    set_hts_verbose(original);
}
