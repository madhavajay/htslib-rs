//! HTSlib-compatible local hFILE helpers.

use std::io;

/// Appends an extension or replaces one extension on a path or URL.
///
/// This follows HTSlib `haddextension` behavior: local paths are adjusted at
/// the end of the filename, while URLs are adjusted before a query or fragment
/// suffix. S3 URLs preserve `#` as part of the path.
pub fn add_extension(filename: &str, replace: bool, extension: &str) -> String {
    let trailing_start = trailing_start(filename);
    let end = if replace {
        strip_extension(filename, trailing_start)
    } else {
        trailing_start
    };

    let mut s = String::with_capacity(filename.len() + extension.len());
    s.push_str(&filename[..end]);
    s.push_str(extension);
    s.push_str(&filename[trailing_start..]);
    s
}

/// Decodes the `data:` URLs supported by HTSlib's hFILE data plugin.
pub fn decode_data_url(src: &str) -> io::Result<Vec<u8>> {
    let payload = src.strip_prefix("data:").ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "data URL must start with data:",
        )
    })?;
    let (metadata, data) = payload
        .split_once(',')
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing data URL comma"))?;

    if metadata
        .split(';')
        .any(|part| part.eq_ignore_ascii_case("base64"))
    {
        decode_base64(data)
    } else {
        percent_decode(data)
    }
}

fn trailing_start(filename: &str) -> usize {
    if is_url(filename) {
        if filename.starts_with("s3://")
            || filename.starts_with("s3+http://")
            || filename.starts_with("s3+https://")
        {
            filename.find('?').unwrap_or(filename.len())
        } else {
            filename.find(['?', '#']).unwrap_or(filename.len())
        }
    } else {
        filename.len()
    }
}

fn strip_extension(filename: &str, limit: usize) -> usize {
    let prefix = &filename[..limit];
    let last_slash = prefix.rfind('/').map(|i| i + 1).unwrap_or(0);

    prefix[last_slash..]
        .rfind('.')
        .map(|i| last_slash + i)
        .unwrap_or(limit)
}

fn is_url(s: &str) -> bool {
    let Some(i) = s.find(':') else {
        return false;
    };

    if i <= 1 || i >= 12 {
        return false;
    }

    s[..i]
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'+' | b'-' | b'.'))
}

fn percent_decode(src: &str) -> io::Result<Vec<u8>> {
    let mut decoded = Vec::with_capacity(src.len());
    let mut bytes = src.as_bytes().iter().copied();

    while let Some(b) = bytes.next() {
        if b != b'%' {
            decoded.push(b);
            continue;
        }

        let hi = bytes.next().ok_or_else(invalid_percent_encoding)?;
        let lo = bytes.next().ok_or_else(invalid_percent_encoding)?;
        let hi = hex_value(hi).ok_or_else(invalid_percent_encoding)?;
        let lo = hex_value(lo).ok_or_else(invalid_percent_encoding)?;

        decoded.push((hi << 4) | lo);
    }

    Ok(decoded)
}

fn invalid_percent_encoding() -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, "invalid percent encoding")
}

fn hex_value(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'A'..=b'F' => Some(b - b'A' + 10),
        b'a'..=b'f' => Some(b - b'a' + 10),
        _ => None,
    }
}

fn decode_base64(src: &str) -> io::Result<Vec<u8>> {
    let mut decoded = Vec::with_capacity(src.len() / 4 * 3);
    let mut quantum = [0; 4];
    let mut n = 0;

    for b in src.bytes().filter(|b| !b.is_ascii_whitespace()) {
        quantum[n] = match b {
            b'A'..=b'Z' => b - b'A',
            b'a'..=b'z' => b - b'a' + 26,
            b'0'..=b'9' => b - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => 64,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid base64",
                ));
            }
        };
        n += 1;

        if n == 4 {
            push_base64_quantum(&mut decoded, quantum)?;
            n = 0;
        }
    }

    if n == 0 {
        Ok(decoded)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "truncated base64",
        ))
    }
}

fn push_base64_quantum(dst: &mut Vec<u8>, quantum: [u8; 4]) -> io::Result<()> {
    if quantum[0] == 64 || quantum[1] == 64 || (quantum[2] == 64 && quantum[3] != 64) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid base64 padding",
        ));
    }

    dst.push((quantum[0] << 2) | (quantum[1] >> 4));

    if quantum[2] != 64 {
        dst.push((quantum[1] << 4) | (quantum[2] >> 2));
    }

    if quantum[3] != 64 {
        dst.push((quantum[2] << 6) | quantum[3]);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{add_extension, decode_data_url};

    #[test]
    fn test_add_extension() {
        assert_eq!(
            add_extension("foo/bar.bam", false, ".bai"),
            "foo/bar.bam.bai"
        );
        assert_eq!(add_extension("foo/bar.bam", true, ".bai"), "foo/bar.bai");
        assert_eq!(
            add_extension("foo.bar/baz", true, ".bai"),
            "foo.bar/baz.bai"
        );
        assert_eq!(
            add_extension("foo#bar.bam", false, ".bai"),
            "foo#bar.bam.bai"
        );
        assert_eq!(add_extension(".bam", true, ".bai"), ".bai");
        assert_eq!(add_extension("foo", true, ".csi"), "foo.csi");
        assert_eq!(
            add_extension("http://host/bar.cram?a&b&c", false, ".crai"),
            "http://host/bar.cram.crai?a&b&c"
        );
        assert_eq!(
            add_extension("http://host/bar.cram#frag", true, ".crai"),
            "http://host/bar.crai#frag"
        );
    }

    #[test]
    fn test_decode_data_url() {
        assert_eq!(
            decode_data_url("data:,hello, world!%0A").unwrap(),
            b"hello, world!\n"
        );
        assert_eq!(decode_data_url("data:,").unwrap(), b"");
        assert_eq!(decode_data_url("data:;base64,SGVsbG8=").unwrap(), b"Hello");
    }
}
