//! HTSlib-compatible text utilities.

/// Result of an integer conversion.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IntParse<T> {
    /// Parsed value, clamped on overflow.
    pub value: T,
    /// Byte index after the parsed number.
    pub end: usize,
    /// Whether overflow occurred.
    pub failed: bool,
}

/// Converts a string to an unsigned integer using HTSlib clamping semantics.
pub fn str_to_uint(input: &str, bits: u32) -> IntParse<u64> {
    let bytes = input.as_bytes();
    let mut i = 0;

    if bytes.get(i) == Some(&b'+') {
        i += 1;
    }

    let limit = if bits >= 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    let (value, end, failed) = parse_unsigned_digits(bytes, i, limit);

    IntParse { value, end, failed }
}

/// Converts a string to a signed integer using HTSlib clamping semantics.
pub fn str_to_int(input: &str, bits: u32) -> IntParse<i64> {
    let bytes = input.as_bytes();
    let mut i = 0;
    let mut negative = false;

    if bytes.get(i) == Some(&b'-') {
        negative = true;
        i += 1;
    } else if bytes.get(i) == Some(&b'+') {
        i += 1;
    }

    let magnitude_limit = if negative {
        1u64 << (bits - 1)
    } else {
        (1u64 << (bits - 1)) - 1
    };

    let (magnitude, end, failed) = parse_unsigned_digits(bytes, i, magnitude_limit);
    let value = if negative {
        if magnitude == (1u64 << 63) {
            i64::MIN
        } else {
            -(magnitude as i64)
        }
    } else {
        magnitude as i64
    };

    IntParse { value, end, failed }
}

fn parse_unsigned_digits(bytes: &[u8], start: usize, limit: u64) -> (u64, usize, bool) {
    let mut value = 0u64;
    let mut i = start;
    let mut failed = false;

    while let Some(&b) = bytes.get(i) {
        if !b.is_ascii_digit() {
            break;
        }

        let digit = u64::from(b - b'0');

        if value < limit / 10 || (value == limit / 10 && digit <= limit % 10) {
            value = value * 10 + digit;
            i += 1;
        } else {
            while bytes.get(i).is_some_and(u8::is_ascii_digit) {
                i += 1;
            }

            value = limit;
            failed = true;
            break;
        }
    }

    (value, i, failed)
}

/// Escapes and truncates possibly malicious text using HTSlib `hts_strprint`
/// semantics.
///
/// `buflen` is the C destination buffer length, including the trailing NUL.
pub fn strprint(input: &[u8], buflen: usize, quote: Option<u8>) -> String {
    let mut out = Vec::with_capacity(buflen.saturating_sub(1));
    let quote_len = usize::from(quote.is_some());

    if let Some(q) = quote {
        out.push(q);
    }

    for &b in input {
        let escaped = escaped_byte(b, quote);

        if out.len() + escaped.len() + quote_len >= buflen {
            while out.len() + 3 + quote_len >= buflen {
                out.pop();
            }

            if let Some(q) = quote {
                out.push(q);
            }

            out.extend_from_slice(b"...");
            return String::from_utf8(out).expect("escaped text is valid UTF-8");
        }

        out.extend_from_slice(&escaped);
    }

    if let Some(q) = quote {
        out.push(q);
    }

    String::from_utf8(out).expect("escaped text is valid UTF-8")
}

fn escaped_byte(b: u8, quote: Option<u8>) -> Vec<u8> {
    match b {
        b'\n' => b"\\n".to_vec(),
        b'\r' => b"\\r".to_vec(),
        b'\t' => b"\\t".to_vec(),
        b'\0' => b"\\0".to_vec(),
        b'\\' => b"\\\\".to_vec(),
        b if Some(b) == quote => vec![b'\\', b],
        b if b.is_ascii_graphic() || b == b' ' => vec![b],
        b => format!("\\x{b:02X}").into_bytes(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_str_to_uint() {
        assert_eq!(
            str_to_uint("255#", 8),
            IntParse {
                value: 255,
                end: 3,
                failed: false
            }
        );
        assert_eq!(
            str_to_uint("256#", 8),
            IntParse {
                value: 255,
                end: 3,
                failed: true
            }
        );
    }

    #[test]
    fn test_str_to_int() {
        assert_eq!(
            str_to_int("127#", 8),
            IntParse {
                value: 127,
                end: 3,
                failed: false
            }
        );
        assert_eq!(
            str_to_int("-128#", 8),
            IntParse {
                value: -128,
                end: 4,
                failed: false
            }
        );
        assert_eq!(
            str_to_int("-129#", 8),
            IntParse {
                value: -128,
                end: 4,
                failed: true
            }
        );
    }

    #[test]
    fn test_strprint() {
        assert_eq!(strprint(b"tab\twxyz", 10, None), "tab\\twxyz");
        assert_eq!(strprint(b"\xab", 5, None), "\\xAB");
        assert_eq!(strprint(b"chr10", 7, Some(b'\'')), "'c'...");
    }
}
