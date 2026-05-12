//! HTSlib-compatible logging helpers and source-style validation.

use std::{
    fs, io,
    path::{Path, PathBuf},
    sync::atomic::{AtomicI32, Ordering},
};

use regex::Regex;

static HTS_VERBOSE: AtomicI32 = AtomicI32::new(LogLevel::Warning as i32);

/// HTSlib log levels.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
#[repr(i32)]
pub enum LogLevel {
    /// All logging disabled.
    Off = 0,
    /// Log errors only.
    Error = 1,
    /// Log errors and warnings.
    Warning = 3,
    /// Log normal but significant events.
    Info = 4,
    /// Log debug events.
    Debug = 5,
    /// Log all events.
    Trace = 6,
}

impl LogLevel {
    /// Returns the HTSlib integer value for this level.
    pub const fn as_i32(self) -> i32 {
        self as i32
    }
}

impl TryFrom<i32> for LogLevel {
    type Error = i32;

    fn try_from(value: i32) -> Result<Self, i32> {
        match value {
            0 => Ok(Self::Off),
            1 => Ok(Self::Error),
            3 => Ok(Self::Warning),
            4 => Ok(Self::Info),
            5 => Ok(Self::Debug),
            6 => Ok(Self::Trace),
            _ => Err(value),
        }
    }
}

/// Sets the selected HTSlib-compatible log level.
pub fn hts_set_log_level(level: LogLevel) {
    HTS_VERBOSE.store(level.as_i32(), Ordering::SeqCst);
}

/// Gets the selected HTSlib-compatible log level.
pub fn hts_get_log_level() -> LogLevel {
    LogLevel::try_from(HTS_VERBOSE.load(Ordering::SeqCst))
        .expect("stored log level should be set through hts_set_log_level")
}

/// Returns the current `hts_verbose` integer value.
pub fn hts_verbose() -> i32 {
    HTS_VERBOSE.load(Ordering::SeqCst)
}

/// Sets the current `hts_verbose` integer value.
pub fn set_hts_verbose(value: i32) {
    HTS_VERBOSE.store(value, Ordering::SeqCst);
}

/// Returns whether an event with `severity` would be logged at the current level.
pub fn hts_log_enabled(severity: LogLevel) -> bool {
    let level = HTS_VERBOSE.load(Ordering::SeqCst);

    level != LogLevel::Off.as_i32() && severity.as_i32() <= level
}

/// A source location whose HTSlib log message does not match project style.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LogMessageIssue {
    /// Source file containing the message.
    pub path: PathBuf,
    /// One-based source line number.
    pub line: usize,
    /// The quoted log message literal.
    pub message: String,
    /// The style violation found for the literal.
    pub kind: LogMessageIssueKind,
}

/// A kind of HTSlib log-message style violation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LogMessageIssueKind {
    /// The literal does not begin with a capital letter, punctuation, digit, or `%s`.
    InvalidStart,
    /// The literal ends with an escaped newline.
    TrailingNewline,
    /// The literal ends with a full stop.
    TrailingFullStop,
}

/// Summary of a source-tree log-message style check.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LogMessageStyleReport {
    /// Number of source files scanned.
    pub files_scanned: usize,
    /// Number of HTSlib log message literals checked.
    pub messages_checked: usize,
    /// Style violations found while scanning.
    pub issues: Vec<LogMessageIssue>,
}

/// Checks C source files under an HTSlib source tree like `test-logging.pl`.
///
/// This scans `*.c` files in the source root and its `cram` subdirectory,
/// extracts single-line `hts_log_*("...")` message literals, and applies the
/// same user-facing style rules as HTSlib's `test/test-logging.pl`.
pub fn check_htslib_log_message_style<P>(src_root: P) -> io::Result<LogMessageStyleReport>
where
    P: AsRef<Path>,
{
    let src_root = src_root.as_ref();
    let pattern = Regex::new(r#"hts_log_\w+\s*\(\s*("[^"]*")"#)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let mut report = LogMessageStyleReport::default();

    check_source_dir(src_root, &pattern, &mut report)?;
    check_source_dir(src_root.join("cram"), &pattern, &mut report)?;

    Ok(report)
}

fn check_source_dir(
    dir: impl AsRef<Path>,
    pattern: &Regex,
    report: &mut LogMessageStyleReport,
) -> io::Result<()> {
    let dir = dir.as_ref();

    if !dir.exists() {
        return Ok(());
    }

    let mut paths = fs::read_dir(dir)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<io::Result<Vec<_>>>()?;
    paths.sort();

    for path in paths {
        if path.extension().and_then(|ext| ext.to_str()) != Some("c") {
            continue;
        }

        check_source_file(&path, pattern, report)?;
    }

    Ok(())
}

fn check_source_file(
    path: &Path,
    pattern: &Regex,
    report: &mut LogMessageStyleReport,
) -> io::Result<()> {
    let src = fs::read_to_string(path)?;

    report.files_scanned += 1;

    for (i, line) in src.lines().enumerate() {
        let Some(captures) = pattern.captures(line) else {
            continue;
        };

        if line.trim_end().ends_with(r#"\n""#) {
            continue;
        }

        let message = captures[1].to_string();
        check_log_message(&message, path, i + 1, report);
    }

    Ok(())
}

fn check_log_message(message: &str, path: &Path, line: usize, report: &mut LogMessageStyleReport) {
    report.messages_checked += 1;

    if !has_valid_log_message_start(message) {
        report.issues.push(LogMessageIssue {
            path: path.into(),
            line,
            message: message.into(),
            kind: LogMessageIssueKind::InvalidStart,
        });
    }

    if message.ends_with(r#"\n""#) {
        report.issues.push(LogMessageIssue {
            path: path.into(),
            line,
            message: message.into(),
            kind: LogMessageIssueKind::TrailingNewline,
        });
    }

    if message.ends_with(r#".""#) {
        report.issues.push(LogMessageIssue {
            path: path.into(),
            line,
            message: message.into(),
            kind: LogMessageIssueKind::TrailingFullStop,
        });
    }
}

fn has_valid_log_message_start(message: &str) -> bool {
    let Some(body) = message.strip_prefix('"') else {
        return false;
    };

    body.starts_with("%s")
        || body
            .as_bytes()
            .first()
            .is_some_and(|b| b.is_ascii_uppercase() || (b'!'..=b'@').contains(b))
}

#[cfg(test)]
mod tests {
    use super::{
        LogLevel, LogMessageIssueKind, LogMessageStyleReport, check_log_message, hts_get_log_level,
        hts_log_enabled, hts_set_log_level, hts_verbose, set_hts_verbose,
    };
    use std::path::Path;

    #[test]
    fn manages_htslib_log_level() {
        let original = hts_verbose();

        hts_set_log_level(LogLevel::Error);
        assert_eq!(hts_get_log_level(), LogLevel::Error);
        assert_eq!(hts_verbose(), 1);
        assert!(hts_log_enabled(LogLevel::Error));
        assert!(!hts_log_enabled(LogLevel::Warning));

        hts_set_log_level(LogLevel::Trace);
        assert!(hts_log_enabled(LogLevel::Debug));
        assert!(hts_log_enabled(LogLevel::Trace));

        set_hts_verbose(original);
    }

    #[test]
    fn accepts_htslib_log_message_style() {
        let mut report = LogMessageStyleReport::default();

        check_log_message(r#""Could not open file""#, Path::new("x.c"), 1, &mut report);
        check_log_message(r#""%s failed""#, Path::new("x.c"), 2, &mut report);
        check_log_message(r#""2 records skipped""#, Path::new("x.c"), 3, &mut report);

        assert_eq!(report.messages_checked, 3);
        assert!(report.issues.is_empty());
    }

    #[test]
    fn reports_htslib_log_message_style_issues() {
        let mut report = LogMessageStyleReport::default();

        check_log_message(r#""lowercase.""#, Path::new("x.c"), 7, &mut report);

        assert_eq!(report.messages_checked, 1);
        assert_eq!(report.issues.len(), 2);
        assert_eq!(report.issues[0].kind, LogMessageIssueKind::InvalidStart);
        assert_eq!(report.issues[1].kind, LogMessageIssueKind::TrailingFullStop);
    }
}
