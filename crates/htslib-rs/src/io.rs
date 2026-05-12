//! Local I/O helpers.

use std::path::Path;

use crate::format::{self, DetectError, Format};

/// Detects the format of a local file.
pub fn detect_format<P>(src: P) -> Result<Format, DetectError>
where
    P: AsRef<Path>,
{
    format::detect_path(src)
}
