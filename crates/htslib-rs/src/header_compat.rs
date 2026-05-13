//! HTSlib-shaped helpers for VCF/BCF header mutation.
//!
//! noodles' [`vcf::Header`] exposes typed APIs for adding `##key=value`
//! records (see `Header::insert`), but the typed-key surface requires callers
//! to produce a parsed `record::key::Other` and a `record::Value`. HTSlib /
//! bcftools work in terms of literal `##key=value` strings — the consumer
//! already has the rendered text in hand. This module bridges that.
//!
//! # HTSlib analogues
//!
//! The functions here cover the HTSlib `bcf_hdr_append(...)` C entry point,
//! which accepts a literal `##key=value` line (or just `key=value`) and
//! routes it to the appropriate typed or unstructured slot in the header.
//! The shape was not previously ported because noodles already covers the
//! typed records (`INFO`, `FILTER`, `FORMAT`, `ALT`, `contig`) via its
//! header-record types — only the unstructured `Other` form needed a
//! string-keyed bridge.
//!
//! # Upstream consumers
//!
//! - **bcftools-rs** uses [`append_other_record`] from
//!   `bcftools_rs::commands::view` (and the same path will fan out to every
//!   subcommand that emits `##bcftools_<cmd>Version` / `##bcftools_<cmd>Command`
//!   lines, mirroring upstream `bcf_hdr_append_version` in
//!   `bcftools/vcfmerge.c:3362`).
//! - The line-form [`append_line`] is the direct analogue of HTSlib's
//!   `bcf_hdr_append` and is a building block for arbitrary callers that
//!   already hold rendered VCF header text (e.g. `bcftools annotate
//!   --header-lines`).
//!
//! # Implementation
//!
//! Backed by noodles' [`vcf::Header::insert`] for the storage; this module
//! adds only the string-keyed surface and the `##key=value` line splitter.
//! No format-level logic is reimplemented.

use std::io;

use crate::vcf;

/// Append a single `##<key>=<value>` line to a VCF/BCF header.
///
/// `key` must be a valid VCF other-key identifier (HTSlib accepts any
/// printable token; noodles validates against
/// `record::key::other::ParseError`). The value is recorded as an
/// unstructured string entry, matching how HTSlib's `bcf_hdr_append` stores
/// unstructured records.
///
/// Returns an `io::Error` with `InvalidData` for either a malformed key or
/// a noodles `AddError` (e.g. mixing structured and unstructured values
/// under the same key).
///
/// # Examples
///
/// ```
/// use htslib_rs::{header_compat, vcf};
///
/// let mut header = vcf::Header::default();
/// header_compat::append_other_record(&mut header, "fileDate", "20260513")?;
/// header_compat::append_other_record(
///     &mut header,
///     "bcftools_viewVersion",
///     "1.23.1+htslib-0.1.0",
/// )?;
///
/// let collection = header.get("fileDate").expect("fileDate present");
/// assert!(matches!(
///     collection,
///     vcf::header::record::value::Collection::Unstructured(_),
/// ));
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn append_other_record(header: &mut vcf::Header, key: &str, value: &str) -> io::Result<()> {
    let parsed_key = key
        .parse::<vcf::header::record::key::Other>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    let parsed_value = vcf::header::record::Value::from(value);
    header
        .insert(parsed_key, parsed_value)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    Ok(())
}

/// Append many `##<key>=<value>` lines in one call. Stops at the first error.
pub fn append_other_records<I, K, V>(header: &mut vcf::Header, lines: I) -> io::Result<()>
where
    I: IntoIterator<Item = (K, V)>,
    K: AsRef<str>,
    V: AsRef<str>,
{
    for (k, v) in lines {
        append_other_record(header, k.as_ref(), v.as_ref())?;
    }
    Ok(())
}

/// Append a literal VCF header line, accepting either `"##key=value"` or
/// `"key=value"` (HTSlib's `bcf_hdr_append` accepts both). A trailing newline
/// is allowed and stripped.
///
/// Direct analogue of HTSlib's `bcf_hdr_append(hdr, line)`. Routes via
/// [`append_other_record`] for the unstructured cases; typed records (`INFO`,
/// `FILTER`, `FORMAT`, `ALT`, `contig`) are NOT yet routed through here and
/// returning an `InvalidData` error — those should go through noodles' typed
/// header-record API. (Tracked: extend to dispatch into typed slots when
/// bcftools-rs's `annotate --header-lines` lands.)
///
/// Returns `InvalidData` if the line does not contain `=` after the leading
/// `##`.
pub fn append_line(header: &mut vcf::Header, line: &str) -> io::Result<()> {
    let body = line.trim_end_matches(['\r', '\n']);
    let body = body.strip_prefix("##").unwrap_or(body);
    let (key, value) = body.split_once('=').ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("not a `##key=value` header line: {body:?}"),
        )
    })?;
    append_other_record(header, key, value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn appends_single_string_record() {
        let mut header = vcf::Header::default();
        append_other_record(&mut header, "fileDate", "20260513").unwrap();

        match header.get("fileDate") {
            Some(vcf::header::record::value::Collection::Unstructured(values)) => {
                assert_eq!(values.as_slice(), &[String::from("20260513")]);
            }
            other => panic!("unexpected collection: {other:?}"),
        }
    }

    #[test]
    fn appends_multiple_records_under_distinct_keys() {
        let mut header = vcf::Header::default();
        append_other_records(
            &mut header,
            [
                ("bcftools_viewVersion", "1.23.1+htslib-0.1.0"),
                ("bcftools_viewCommand", "view -Oz in.vcf.gz; Date=now"),
            ],
        )
        .unwrap();

        assert!(header.get("bcftools_viewVersion").is_some());
        assert!(header.get("bcftools_viewCommand").is_some());
    }

    #[test]
    fn appending_same_key_twice_extends_collection() {
        let mut header = vcf::Header::default();
        append_other_record(&mut header, "fileDate", "20240101").unwrap();
        append_other_record(&mut header, "fileDate", "20260513").unwrap();

        match header.get("fileDate") {
            Some(vcf::header::record::value::Collection::Unstructured(values)) => {
                assert_eq!(
                    values.as_slice(),
                    &[String::from("20240101"), String::from("20260513")],
                );
            }
            other => panic!("unexpected collection: {other:?}"),
        }
    }

    #[test]
    fn rejects_reserved_key() {
        let mut header = vcf::Header::default();
        // `INFO`, `FILTER`, `FORMAT`, `ALT`, and `contig` are typed records,
        // not free-form `##key=value` lines. noodles' `Other` key parser
        // rejects them so they can't bypass the typed surface.
        for reserved in ["INFO", "FILTER", "FORMAT", "contig"] {
            let err = append_other_record(&mut header, reserved, "value");
            assert!(err.is_err(), "expected error for reserved key {reserved}");
        }
    }

    #[test]
    fn append_line_accepts_double_hash_form() {
        let mut header = vcf::Header::default();
        append_line(&mut header, "##fileDate=20260513").unwrap();
        assert!(header.get("fileDate").is_some());
    }

    #[test]
    fn append_line_accepts_bare_form() {
        let mut header = vcf::Header::default();
        append_line(&mut header, "fileDate=20260513").unwrap();
        assert!(header.get("fileDate").is_some());
    }

    #[test]
    fn append_line_strips_trailing_newline() {
        let mut header = vcf::Header::default();
        append_line(&mut header, "##fileDate=20260513\n").unwrap();
        match header.get("fileDate") {
            Some(vcf::header::record::value::Collection::Unstructured(values)) => {
                assert_eq!(values.as_slice(), &[String::from("20260513")]);
            }
            other => panic!("unexpected: {other:?}"),
        }
    }

    #[test]
    fn append_line_rejects_lines_without_equals() {
        let mut header = vcf::Header::default();
        let err = append_line(&mut header, "##fileDate");
        assert!(err.is_err());
    }

    #[test]
    fn round_trips_through_vcf_text() {
        let mut header = vcf::Header::default();
        append_other_record(&mut header, "fileDate", "20260513").unwrap();
        append_other_record(
            &mut header,
            "bcftools_viewCommand",
            "view -Oz in.vcf; Date=now",
        )
        .unwrap();

        let mut buf = Vec::new();
        crate::vcf::io::Writer::new(&mut buf)
            .write_header(&header)
            .unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("##fileDate=20260513"));
        assert!(text.contains("##bcftools_viewCommand=view -Oz in.vcf; Date=now"));
    }
}
