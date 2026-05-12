//! HTSlib-compatible in-memory region index helpers.

use std::collections::BTreeMap;

/// Maximum coordinate accepted by HTSlib's region index.
pub const REGIDX_MAX: i64 = 1_i64 << 35;

/// Supported HTSlib region-index line parsers.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Parser {
    /// Whitespace-separated `CHROM FROM TO`, 0-based, right-open.
    Bed,
    /// Whitespace-separated `CHROM POS [TO]`, 1-based, inclusive.
    Tab,
    /// `CHROM`, `CHROM:POS`, `CHROM:FROM-TO`, or `CHROM:FROM-`, 1-based, inclusive.
    Region,
}

/// A parsed or indexed region.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegionRecord {
    /// Reference sequence name.
    pub seq: String,
    /// 0-based, inclusive start.
    pub start: i64,
    /// 0-based, inclusive end.
    pub end: i64,
    /// Optional payload. For custom/tabular tests this is the fourth field.
    pub payload: Option<String>,
}

impl RegionRecord {
    fn new(seq: &str, start: i64, end: i64, payload: Option<String>) -> Self {
        Self {
            seq: seq.to_string(),
            start: start.clamp(0, REGIDX_MAX),
            end: end.clamp(0, REGIDX_MAX),
            payload,
        }
    }

    fn overlaps(&self, start: i64, end: i64) -> bool {
        self.end >= start && self.start <= end
    }
}

/// Region-index parse error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Blank or comment line that should be skipped.
    Skip,
    /// Malformed input.
    Invalid,
    /// Coordinate is zero where HTSlib expects a 1-based coordinate.
    ZeroCoordinate,
}

/// An in-memory region index.
#[derive(Clone, Debug, Default)]
pub struct RegionIndex {
    records: BTreeMap<String, Vec<RegionRecord>>,
}

impl RegionIndex {
    /// Creates an empty index.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parses and inserts a line using the given parser.
    pub fn insert_line(&mut self, line: &str, parser: Parser) -> Result<bool, ParseError> {
        match parse_line(line, parser)? {
            Some(record) => {
                self.push(record);
                Ok(true)
            }
            None => Ok(false),
        }
    }

    /// Inserts a pre-parsed record.
    pub fn push(&mut self, record: RegionRecord) {
        let records = self.records.entry(record.seq.clone()).or_default();
        let pos = records
            .binary_search_by_key(&(record.start, record.end), |record| {
                (record.start, record.end)
            })
            .unwrap_or_else(|pos| pos);
        records.insert(pos, record);
    }

    /// Returns all records overlapping `seq:start-end`, where coordinates are 0-based inclusive.
    pub fn overlaps<'a>(&'a self, seq: &str, start: i64, end: i64) -> Vec<&'a RegionRecord> {
        self.records
            .get(seq)
            .map(|records| {
                records
                    .iter()
                    .skip_while(|record| record.end < start)
                    .take_while(|record| record.start <= end)
                    .filter(|record| record.overlaps(start, end))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Returns whether any record overlaps `seq:start-end`.
    pub fn has_overlap(&self, seq: &str, start: i64, end: i64) -> bool {
        self.records
            .get(seq)
            .is_some_and(|records| records.iter().any(|record| record.overlaps(start, end)))
    }

    /// Returns all indexed records in sequence-name order and then coordinate order.
    pub fn iter(&self) -> impl Iterator<Item = &RegionRecord> {
        self.records.values().flat_map(|records| records.iter())
    }
}

/// Parses a line using one of HTSlib's built-in region-index parsers.
pub fn parse_line(line: &str, parser: Parser) -> Result<Option<RegionRecord>, ParseError> {
    match parser {
        Parser::Bed => parse_bed(line),
        Parser::Tab => parse_tab(line),
        Parser::Region => parse_region_line(line),
    }
}

/// Parses a BED-style line: `CHROM FROM TO`, 0-based, right-open.
pub fn parse_bed(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    let Some((seq, fields)) = split_record(line) else {
        return Ok(None);
    };

    if fields.is_empty() {
        return Ok(Some(RegionRecord::new(seq, 0, REGIDX_MAX, None)));
    }

    let start = parse_i64(fields[0])?;
    let end = fields
        .get(1)
        .ok_or(ParseError::Invalid)
        .and_then(|s| parse_i64(s))?
        - 1;

    Ok(Some(RegionRecord::new(seq, start, end, payload(&fields))))
}

/// Parses a tabular line: `CHROM POS [TO]`, 1-based, inclusive.
pub fn parse_tab(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    let Some((seq, fields)) = split_record(line) else {
        return Ok(None);
    };

    if fields.is_empty() {
        return Ok(Some(RegionRecord::new(seq, 0, REGIDX_MAX, None)));
    }

    let pos = parse_i64(fields[0])?;
    if pos == 0 {
        return Err(ParseError::ZeroCoordinate);
    }

    let start = pos - 1;
    let end = match fields.get(1).map(|s| parse_i64(s)) {
        Some(Ok(0)) => return Err(ParseError::ZeroCoordinate),
        Some(Ok(value)) => value - 1,
        Some(Err(_)) | None => start,
    };

    Ok(Some(RegionRecord::new(seq, start, end, payload(&fields))))
}

/// Parses a region line: `CHROM`, `CHROM:POS`, `CHROM:FROM-TO`, or `CHROM:FROM-`.
pub fn parse_region_line(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    let line = line.trim_start();

    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }

    let token = line.split_whitespace().next().ok_or(ParseError::Invalid)?;
    let Some((seq, range)) = token.split_once(':') else {
        return Ok(Some(RegionRecord::new(token, 0, REGIDX_MAX, None)));
    };

    let (from, to) = range.split_once('-').unwrap_or((range, ""));
    let pos = parse_i64(from)?;
    if pos == 0 {
        return Err(ParseError::ZeroCoordinate);
    }

    let start = pos - 1;
    let end = if range.ends_with('-') {
        REGIDX_MAX
    } else if to.is_empty() {
        start
    } else {
        let value = parse_i64(to)?;
        if value == 0 {
            return Err(ParseError::ZeroCoordinate);
        }
        value - 1
    };

    Ok(Some(RegionRecord::new(seq, start, end, None)))
}

fn split_record(line: &str) -> Option<(&str, Vec<&str>)> {
    let line = line.trim_start();

    if line.is_empty() || line.starts_with('#') {
        return None;
    }

    let mut fields = line.split_whitespace();
    let seq = fields.next()?;

    Some((seq, fields.collect()))
}

fn payload(fields: &[&str]) -> Option<String> {
    fields.get(2).map(|s| (*s).to_string())
}

fn parse_i64(s: &str) -> Result<i64, ParseError> {
    let value = s
        .bytes()
        .take_while(u8::is_ascii_digit)
        .fold((0_i64, 0_usize), |(value, len), b| {
            (value.saturating_mul(10) + i64::from(b - b'0'), len + 1)
        });

    if value.1 == 0 {
        Err(ParseError::Invalid)
    } else {
        Ok(value.0)
    }
}

#[cfg(test)]
mod tests {
    use super::{Parser, RegionIndex, parse_bed, parse_region_line, parse_tab};

    #[test]
    fn test_parse_builtins() {
        let record = parse_bed("sq0\t8\t13\tpayload").unwrap().unwrap();
        assert_eq!(
            (record.seq.as_str(), record.start, record.end),
            ("sq0", 8, 12)
        );
        assert_eq!(record.payload.as_deref(), Some("payload"));

        let record = parse_tab("sq0\t9\t13\tpayload").unwrap().unwrap();
        assert_eq!(
            (record.seq.as_str(), record.start, record.end),
            ("sq0", 8, 12)
        );
        assert_eq!(record.payload.as_deref(), Some("payload"));

        let record = parse_region_line("sq0:9-13").unwrap().unwrap();
        assert_eq!(
            (record.seq.as_str(), record.start, record.end),
            ("sq0", 8, 12)
        );
    }

    #[test]
    fn test_index_overlap() {
        let mut index = RegionIndex::new();
        index.insert_line("1\t10\t10", Parser::Tab).unwrap();
        index.insert_line("1\t11\t11", Parser::Tab).unwrap();

        assert!(!index.has_overlap("1", 8, 8));
        assert_eq!(index.overlaps("1", 9, 9).len(), 1);
        assert_eq!(index.overlaps("1", 9, 10).len(), 2);
    }
}
