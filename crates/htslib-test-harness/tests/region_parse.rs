use htslib_rs::region::{HTS_POS_MAX, ParseError, ParseFlags, ParsedRegion, parse_region};

fn name_to_id(name: &str) -> Option<i32> {
    [
        "chr1",
        "chr1:100",
        "chr1:100-200",
        "chr2:100-200",
        "chr3",
        "chr1,chr3",
    ]
    .iter()
    .position(|&reference_sequence_name| reference_sequence_name == name)
    .map(|i| i as i32)
}

fn parse(input: &str, flags: ParseFlags) -> Result<ParsedRegion<'_>, ParseError> {
    parse_region(input, name_to_id, flags)
}

fn assert_region(input: &str, flags: ParseFlags, rest: &str, tid: i32, start: i64, end: i64) {
    assert_eq!(
        parse(input, flags),
        Ok(ParsedRegion {
            tid,
            start,
            end,
            rest
        })
    );
}

#[test]
fn parses_range_extensions() {
    assert_region("chr1", ParseFlags::default(), "", 0, 0, HTS_POS_MAX);
    assert_region("chr1:50", ParseFlags::default(), "", 0, 49, HTS_POS_MAX);
    assert_region("chr1:50", ParseFlags::ONE_COORD, "", 0, 49, 50);
    assert_region("chr1:50-100", ParseFlags::default(), "", 0, 49, 100);
    assert_region("chr1:50-", ParseFlags::default(), "", 0, 49, HTS_POS_MAX);
    assert_region("chr1:-50", ParseFlags::default(), "", 0, 0, 50);
}

#[test]
fn parses_quoted_and_ambiguous_regions() {
    assert_eq!(
        parse("chr1:100-200", ParseFlags::default()),
        Err(ParseError::Ambiguous)
    );

    assert_region("{chr1}:100-200", ParseFlags::default(), "", 0, 99, 200);
    assert_region(
        "{chr1:100-200}",
        ParseFlags::default(),
        "",
        2,
        0,
        HTS_POS_MAX,
    );
    assert_region(
        "{chr1:100-200}:100-200",
        ParseFlags::default(),
        "",
        2,
        99,
        200,
    );
    assert_region(
        "{chr2:100-200}:100-200",
        ParseFlags::default(),
        "",
        3,
        99,
        200,
    );
    assert_region(
        "chr2:100-200:100-200",
        ParseFlags::default(),
        "",
        3,
        99,
        200,
    );
    assert_region("chr2:100-200", ParseFlags::default(), "", 3, 0, HTS_POS_MAX);
}

#[test]
fn parses_numeric_forms() {
    assert_region("chr3", ParseFlags::default(), "", 4, 0, HTS_POS_MAX);
    assert_region("chr3:", ParseFlags::default(), "", 4, 0, HTS_POS_MAX);
    assert_region("chr3:1000-1500", ParseFlags::default(), "", 4, 999, 1500);
    assert_region("chr3:1,000-1,500", ParseFlags::default(), "", 4, 999, 1500);
    assert_region("chr3:1k-1.5K", ParseFlags::default(), "", 4, 999, 1500);
    assert_region("chr3:1e3-1.5e3", ParseFlags::default(), "", 4, 999, 1500);
    assert_region("chr3:1e3-15e2", ParseFlags::default(), "", 4, 999, 1500);
}

#[test]
fn parses_list_mode() {
    assert_region("chr1,chr3", ParseFlags::LIST, "chr3", 0, 0, HTS_POS_MAX);
    assert_eq!(
        parse("chr1:100-200,chr3", ParseFlags::LIST),
        Err(ParseError::Ambiguous)
    );
    assert_region("{chr1,chr3}", ParseFlags::LIST, "", 5, 0, HTS_POS_MAX);
    assert_region(
        "{chr1,chr3},chr1",
        ParseFlags::LIST,
        "chr1",
        5,
        0,
        HTS_POS_MAX,
    );
    assert_region(
        "chr3:1,000-1,500",
        ParseFlags::LIST | ParseFlags::ONE_COORD,
        "000-1,500",
        4,
        0,
        1,
    );
}

#[test]
fn rejects_invalid_regions() {
    assert_eq!(
        parse("chr2", ParseFlags::default()),
        Err(ParseError::UnknownReference)
    );
    assert_eq!(
        parse("chr1,", ParseFlags::default()),
        Err(ParseError::UnknownReference)
    );
    assert_eq!(
        parse("{chr1", ParseFlags::default()),
        Err(ParseError::MismatchedBraces)
    );

    assert_region("chr1:10-10", ParseFlags::default(), "", 0, 9, 10);
    assert_eq!(
        parse("chr1:10-9", ParseFlags::default()),
        Err(ParseError::InvalidCoordinates)
    );
    assert_eq!(
        parse("chr1:x", ParseFlags::default()),
        Err(ParseError::InvalidCoordinates)
    );
    assert_eq!(
        parse("chr1:1-y", ParseFlags::default()),
        Err(ParseError::InvalidCoordinates)
    );
    assert_eq!(
        parse("chr1:1,chr3", ParseFlags::default()),
        Err(ParseError::InvalidCoordinates)
    );
}
