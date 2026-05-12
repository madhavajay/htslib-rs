//! HTSlib-compatible region parsing.

use std::ops::{BitOr, BitOrAssign};

/// The maximum HTSlib position value.
pub const HTS_POS_MAX: i64 = ((i32::MAX as i64) << 32) | i32::MAX as i64;

/// Region parser flags.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ParseFlags(u8);

impl ParseFlags {
    /// Ignore commas in numbers.
    pub const THOUSANDS_SEP: Self = Self(1);
    /// Treat `chr:pos` as the single-base region `chr:pos-pos`.
    pub const ONE_COORD: Self = Self(2);
    /// Parse a comma-separated list of regions.
    pub const LIST: Self = Self(4);

    /// Returns true if all flags in `other` are set.
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }

    fn insert(&mut self, other: Self) {
        self.0 |= other.0;
    }

    fn remove(&mut self, other: Self) {
        self.0 &= !other.0;
    }
}

impl BitOr for ParseFlags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for ParseFlags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.insert(rhs);
    }
}

/// A parsed region.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParsedRegion<'a> {
    /// Reference sequence ID.
    pub tid: i32,
    /// 0-based, inclusive start.
    pub start: i64,
    /// 0-based, exclusive end.
    pub end: i64,
    /// Remaining input after the parsed region.
    pub rest: &'a str,
}

/// Region parse error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is malformed.
    Invalid,
    /// Braces are mismatched.
    MismatchedBraces,
    /// The reference name is ambiguous with a ranged expression.
    Ambiguous,
    /// The reference name does not exist.
    UnknownReference,
    /// The coordinates are invalid.
    InvalidCoordinates,
}

/// Parses an HTSlib-style region using a reference-name lookup callback.
pub fn parse_region<'a, F>(
    input: &'a str,
    mut name_to_id: F,
    mut flags: ParseFlags,
) -> Result<ParsedRegion<'a>, ParseError>
where
    F: FnMut(&str) -> Option<i32>,
{
    if flags.contains(ParseFlags::LIST) {
        flags.remove(ParseFlags::THOUSANDS_SEP);
    } else {
        flags.insert(ParseFlags::THOUSANDS_SEP);
    }

    let bytes = input.as_bytes();
    let len = bytes.len();

    if len == 0 {
        return Err(ParseError::UnknownReference);
    }

    let mut quoted = false;
    let mut name_start = 0;
    let mut item_end = len;
    let mut rest_start = len;
    let colon;

    if bytes[0] == b'{' {
        let close = bytes
            .iter()
            .position(|&b| b == b'}')
            .ok_or(ParseError::MismatchedBraces)?;

        quoted = true;
        name_start = 1;

        colon = if bytes.get(close + 1) == Some(&b':') {
            Some(close + 1)
        } else {
            None
        };

        if flags.contains(ParseFlags::LIST)
            && let Some(comma) = bytes[close..].iter().position(|&b| b == b',')
        {
            item_end = close + comma;
            rest_start = item_end + 1;
        }
    } else {
        if flags.contains(ParseFlags::LIST)
            && let Some(comma) = bytes.iter().position(|&b| b == b',')
        {
            item_end = comma;
            rest_start = comma + 1;
        }

        colon = bytes[..item_end].iter().rposition(|&b| b == b':');
    }

    if colon.is_none() {
        let name_end = if quoted {
            bytes[name_start..item_end]
                .iter()
                .position(|&b| b == b'}')
                .map(|i| name_start + i)
                .unwrap_or(item_end.saturating_sub(1))
        } else {
            item_end
        };
        let name = &input[name_start..name_end];
        let tid = name_to_id(name).ok_or(ParseError::UnknownReference)?;

        return Ok(ParsedRegion {
            tid,
            start: 0,
            end: HTS_POS_MAX,
            rest: &input[rest_start..],
        });
    }

    let colon = colon.expect("checked above");

    if !quoted {
        let whole_name = &input[..item_end];

        if let Some(tid) = name_to_id(whole_name) {
            let prefix = &input[..colon];

            if name_to_id(prefix).is_some() {
                return Err(ParseError::Ambiguous);
            }

            return Ok(ParsedRegion {
                tid,
                start: 0,
                end: HTS_POS_MAX,
                rest: &input[rest_start..],
            });
        }
    }

    let name_end = if quoted { colon - 1 } else { colon };
    let name = &input[name_start..name_end];
    let tid = name_to_id(name).ok_or(ParseError::UnknownReference)?;

    let coord_start = colon + 1;
    let (n, mut offset, digits) = parse_decimal(&input[coord_start..], flags);
    offset += coord_start;

    let mut start = n - 1;
    let mut end;
    let next = bytes.get(offset).copied();

    if start < 0 {
        if start != -1 && next == Some(b'-') && coord_start < len {
            return Err(ParseError::InvalidCoordinates);
        }

        if matches!(next, Some(b'0'..=b'9') | Some(b',') | None) {
            end = if start == -1 {
                HTS_POS_MAX
            } else {
                -(start + 1)
            };
            start = 0;
            return finish_region(input, rest_start, tid, start, end);
        } else if start < -1 {
            return Err(ParseError::InvalidCoordinates);
        }
    }

    if !digits && !matches!(next, Some(b'-') | Some(b',') | None) {
        return Err(ParseError::InvalidCoordinates);
    }

    if next.is_none() || (flags.contains(ParseFlags::LIST) && next == Some(b',')) {
        end = if flags.contains(ParseFlags::ONE_COORD) {
            start + 1
        } else {
            HTS_POS_MAX
        };
    } else if next == Some(b'-') {
        let (n, mut end_offset, _) = parse_decimal(&input[offset + 1..], flags);
        end_offset += offset + 1;

        if !matches!(bytes.get(end_offset), None | Some(b',')) {
            return Err(ParseError::InvalidCoordinates);
        }

        end = n;
    } else {
        return Err(ParseError::InvalidCoordinates);
    }

    if end == 0 {
        end = HTS_POS_MAX;
    }

    finish_region(input, rest_start, tid, start, end)
}

fn finish_region(
    input: &str,
    rest_start: usize,
    tid: i32,
    start: i64,
    end: i64,
) -> Result<ParsedRegion<'_>, ParseError> {
    if start >= end {
        return Err(ParseError::InvalidCoordinates);
    }

    Ok(ParsedRegion {
        tid,
        start,
        end,
        rest: &input[rest_start..],
    })
}

fn parse_decimal(s: &str, flags: ParseFlags) -> (i64, usize, bool) {
    let bytes = s.as_bytes();
    let mut i = 0;

    while bytes.get(i).is_some_and(u8::is_ascii_whitespace) {
        i += 1;
    }

    let mut sign = 1i64;

    if bytes.get(i) == Some(&b'+') {
        i += 1;
    } else if bytes.get(i) == Some(&b'-') {
        sign = -1;
        i += 1;
    }

    let mut digits = 0;
    let mut decimals = 0;
    let mut value = 0i64;

    while let Some(&b) = bytes.get(i) {
        if b.is_ascii_digit() {
            digits += 1;
            value = value.saturating_mul(10).saturating_add(i64::from(b - b'0'));
            i += 1;
        } else if b == b',' && flags.contains(ParseFlags::THOUSANDS_SEP) {
            i += 1;
        } else {
            break;
        }
    }

    if bytes.get(i) == Some(&b'.') {
        i += 1;

        while let Some(&b) = bytes.get(i) {
            if b.is_ascii_digit() {
                digits += 1;
                decimals += 1;
                value = value.saturating_mul(10).saturating_add(i64::from(b - b'0'));
                i += 1;
            } else {
                break;
            }
        }
    }

    let mut exponent = 0i32;

    match bytes.get(i).copied() {
        Some(b'e' | b'E') => {
            i += 1;

            let mut exponent_sign = 1i32;

            if bytes.get(i) == Some(&b'+') {
                i += 1;
            } else if bytes.get(i) == Some(&b'-') {
                exponent_sign = -1;
                i += 1;
            }

            while let Some(&b) = bytes.get(i) {
                if b.is_ascii_digit() {
                    exponent = exponent
                        .saturating_mul(10)
                        .saturating_add(i32::from(b - b'0'));
                    i += 1;
                } else {
                    break;
                }
            }

            exponent *= exponent_sign;
        }
        Some(b'k' | b'K') => {
            exponent += 3;
            i += 1;
        }
        Some(b'm' | b'M') => {
            exponent += 6;
            i += 1;
        }
        Some(b'g' | b'G') => {
            exponent += 9;
            i += 1;
        }
        _ => {}
    }

    exponent -= decimals;

    while exponent > 0 {
        value = value.saturating_mul(10);
        exponent -= 1;
    }

    while exponent < 0 {
        value /= 10;
        exponent += 1;
    }

    (sign * value, if digits > 0 { i } else { 0 }, digits > 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_decimal() {
        let flags = ParseFlags::THOUSANDS_SEP;

        assert_eq!(parse_decimal("1,000", flags), (1000, 5, true));
        assert_eq!(parse_decimal("1k", flags), (1000, 2, true));
        assert_eq!(parse_decimal("1.5K", flags), (1500, 4, true));
        assert_eq!(parse_decimal("1e3", flags), (1000, 3, true));
        assert_eq!(parse_decimal("15e2", flags), (1500, 4, true));
        assert_eq!(parse_decimal("x", flags), (0, 0, false));
    }
}
