//! Rust replacements for selected HTSlib kstring behavior.

/// A growable byte string with HTSlib-like insertion and append helpers.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct KString {
    buf: Vec<u8>,
}

impl KString {
    /// Creates an empty string.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the string length.
    pub fn len(&self) -> usize {
        self.buf.len()
    }

    /// Returns the current allocation capacity.
    pub fn capacity(&self) -> usize {
        self.buf.capacity()
    }

    /// Returns whether the string is empty.
    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    /// Clears the string.
    pub fn clear(&mut self) {
        self.buf.clear();
    }

    /// Reserves capacity for at least `size` bytes.
    pub fn resize_capacity(&mut self, size: usize) {
        if self.buf.capacity() < size {
            self.buf.reserve(size - self.buf.capacity());
        }
    }

    /// Reserves capacity for `expansion` additional bytes.
    pub fn expand_capacity(&mut self, expansion: usize) {
        self.resize_capacity(self.buf.len().saturating_add(expansion));
    }

    /// Returns the bytes.
    pub fn as_bytes(&self) -> &[u8] {
        &self.buf
    }

    /// Appends bytes.
    pub fn push_bytes(&mut self, bytes: &[u8]) {
        self.buf.extend_from_slice(bytes);
    }

    /// Appends a byte.
    pub fn push_byte(&mut self, byte: u8) {
        self.buf.push(byte);
    }

    /// Appends an unsigned integer as decimal ASCII.
    pub fn push_u32(&mut self, value: u32) {
        self.push_bytes(value.to_string().as_bytes());
    }

    /// Appends a signed 32-bit integer as decimal ASCII.
    pub fn push_i32(&mut self, value: i32) {
        self.push_bytes(value.to_string().as_bytes());
    }

    /// Appends a signed 64-bit integer as decimal ASCII.
    pub fn push_i64(&mut self, value: i64) {
        self.push_bytes(value.to_string().as_bytes());
    }

    /// Inserts a byte at `position`.
    pub fn insert_byte(&mut self, position: usize, byte: u8) -> Result<(), InsertError> {
        if position > self.buf.len() {
            return Err(InsertError::OutOfBounds);
        }

        self.buf.insert(position, byte);
        Ok(())
    }

    /// Inserts bytes at `position`.
    pub fn insert_bytes(&mut self, position: usize, bytes: &[u8]) -> Result<(), InsertError> {
        if position > self.buf.len() {
            return Err(InsertError::OutOfBounds);
        }

        if !bytes.is_empty() {
            self.buf.splice(position..position, bytes.iter().copied());
        }

        Ok(())
    }
}

/// HTSlib-style initializer.
pub fn ks_initialize(s: &mut KString) {
    s.clear();
}

/// HTSlib-style capacity resize.
pub fn ks_resize(s: &mut KString, size: usize) {
    s.resize_capacity(size);
}

/// HTSlib-style capacity expansion.
pub fn ks_expand(s: &mut KString, expansion: usize) {
    s.expand_capacity(expansion);
}

/// HTSlib-style buffer accessor.
pub fn ks_str(s: &KString) -> &[u8] {
    s.as_bytes()
}

/// HTSlib-style non-null buffer accessor.
pub fn ks_c_str(s: &KString) -> &[u8] {
    s.as_bytes()
}

/// HTSlib-style length accessor.
pub fn ks_len(s: &KString) -> usize {
    s.len()
}

/// HTSlib-style clear operation.
pub fn ks_clear(s: &mut KString) -> &mut KString {
    s.clear();
    s
}

/// HTSlib-style release operation.
pub fn ks_release(s: &mut KString) -> Vec<u8> {
    std::mem::take(&mut s.buf)
}

/// HTSlib-style free operation.
pub fn ks_free(s: &mut KString) {
    s.buf = Vec::new();
}

/// HTSlib-style append bytes with trailing-NUL semantics represented by bytes only.
pub fn kputsn(p: &[u8], s: &mut KString) -> usize {
    s.push_bytes(p);
    p.len()
}

/// HTSlib-style append string.
pub fn kputs(p: &str, s: &mut KString) -> usize {
    kputsn(p.as_bytes(), s)
}

/// HTSlib-style append byte.
pub fn kputc(c: u8, s: &mut KString) -> u8 {
    s.push_byte(c);
    c
}

/// HTSlib-style append byte without trailing-NUL semantics.
pub fn kputc_(c: u8, s: &mut KString) -> usize {
    s.push_byte(c);
    1
}

/// HTSlib-style append raw bytes without trailing-NUL semantics.
pub fn kputsn_(p: &[u8], s: &mut KString) -> usize {
    s.push_bytes(p);
    p.len()
}

/// HTSlib-style append unsigned decimal integer.
pub fn kputuw(value: u32, s: &mut KString) -> usize {
    let before = s.len();
    s.push_u32(value);
    s.len() - before
}

/// HTSlib-style append signed 32-bit decimal integer.
pub fn kputw(value: i32, s: &mut KString) -> usize {
    let before = s.len();
    s.push_i32(value);
    s.len() - before
}

/// HTSlib-style append signed 64-bit decimal integer.
pub fn kputll(value: i64, s: &mut KString) -> usize {
    let before = s.len();
    s.push_i64(value);
    s.len() - before
}

/// HTSlib-style append signed long decimal integer.
pub fn kputl(value: isize, s: &mut KString) -> usize {
    kputll(value as i64, s)
}

/// HTSlib-style character insertion.
pub fn kinsert_char(c: u8, pos: usize, s: &mut KString) -> Result<(), InsertError> {
    s.insert_byte(pos, c)
}

/// HTSlib-style string insertion.
pub fn kinsert_str(str_: &str, pos: usize, s: &mut KString) -> Result<(), InsertError> {
    s.insert_bytes(pos, str_.as_bytes())
}

/// Insertion error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InsertError {
    /// The insertion position is outside the current string.
    OutOfBounds,
}

/// Rounds up to the next power of two, matching HTSlib's `kroundup_size_t`.
pub fn roundup_size_t(value: usize) -> usize {
    crate::kroundup::roundup_size_t(value)
}

/// Rounds up to the next signed 32-bit power-of-two boundary.
///
/// Values above `1 << 30` saturate to `i32::MAX`, matching the upstream test's
/// signed-overflow avoidance expectation.
pub fn roundup_i32(value: i32) -> u32 {
    crate::kroundup::roundup_i32(value)
}

/// Finds `needle` in `haystack`, including NUL bytes.
pub fn memmem(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }

    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

/// HTSlib-style memory search alias.
pub fn kmemmem(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    memmem(haystack, needle)
}

/// Finds `needle` in a UTF-8 string.
pub fn strstr(haystack: &str, needle: &str) -> Option<usize> {
    haystack.find(needle)
}

/// HTSlib-style string search alias.
pub fn kstrstr(haystack: &str, needle: &str) -> Option<usize> {
    strstr(haystack, needle)
}

/// Finds `needle` in at most the first `limit` bytes of a C-style string.
pub fn strnstr(haystack: &[u8], needle: &[u8], limit: usize) -> Option<usize> {
    let capped = &haystack[..haystack.len().min(limit)];
    let c_string = match capped.iter().position(|&b| b == 0) {
        Some(nul) => &capped[..nul],
        None => capped,
    };

    memmem(c_string, needle)
}

/// HTSlib-style bounded string search alias.
pub fn kstrnstr(haystack: &[u8], needle: &[u8], limit: usize) -> Option<usize> {
    strnstr(haystack, needle, limit)
}

/// Reads a single HTSlib-style line from a chunk provider into `dst`.
///
/// Returns `true` when a line or EOF-terminated partial line was read, and
/// `false` when EOF was reached before reading any bytes.
pub fn getline_from_chunks<I>(dst: &mut KString, chunks: &mut I) -> bool
where
    I: Iterator,
    I::Item: AsRef<[u8]>,
{
    for chunk in chunks.by_ref() {
        let chunk = chunk.as_ref();

        if let Some(newline) = chunk.iter().position(|&b| b == b'\n') {
            let line = &chunk[..newline];
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            dst.push_bytes(line);
            return true;
        }

        dst.push_bytes(chunk);
    }

    !dst.is_empty()
}

/// HTSlib-style line-read alias over Rust chunk iterators.
pub fn kgetline_from_chunks<I>(dst: &mut KString, chunks: &mut I) -> bool
where
    I: Iterator,
    I::Item: AsRef<[u8]>,
{
    getline_from_chunks(dst, chunks)
}

#[cfg(test)]
mod tests {
    use super::{
        KString, kgetline_from_chunks, kinsert_char, kinsert_str, kmemmem, kputc, kputc_, kputl,
        kputll, kputs, kputsn, kputsn_, kputuw, kputw, ks_c_str, ks_clear, ks_expand, ks_free,
        ks_initialize, ks_len, ks_release, ks_resize, ks_str, kstrnstr, kstrstr, memmem,
        roundup_i32, roundup_size_t, strnstr, strstr,
    };

    #[test]
    fn test_roundup() {
        assert_eq!(roundup_size_t(0), 0);
        assert_eq!(roundup_size_t(3), 4);
        assert_eq!(roundup_i32((1 << 30) + 1), i32::MAX as u32);
    }

    #[test]
    fn test_kstring() {
        let mut s = KString::new();

        s.push_i32(-12);
        assert_eq!(s.as_bytes(), b"-12");

        s.insert_bytes(0, b"v=").unwrap();
        assert_eq!(s.as_bytes(), b"v=-12");
    }

    #[test]
    fn test_searches() {
        assert_eq!(memmem(b"a\0bc", b"\0b"), Some(1));
        assert_eq!(strstr("foo", "oo"), Some(1));
        assert_eq!(strnstr(b"foo\0bar", b"bar", 7), None);
    }

    #[test]
    fn test_c_shaped_aliases() {
        let mut s = KString::new();
        ks_resize(&mut s, 16);
        assert!(s.capacity() >= 16);
        ks_expand(&mut s, 4);

        assert_eq!(kputs("a", &mut s), 1);
        assert_eq!(kputsn(b"bc", &mut s), 2);
        assert_eq!(kputc(b'd', &mut s), b'd');
        assert_eq!(kputc_(b'e', &mut s), 1);
        assert_eq!(kputsn_(b"fg", &mut s), 2);
        assert_eq!(ks_str(&s), b"abcdefg");
        assert_eq!(ks_c_str(&s), b"abcdefg");
        assert_eq!(ks_len(&s), 7);

        assert_eq!(kputuw(12, &mut s), 2);
        assert_eq!(kputw(-3, &mut s), 2);
        assert_eq!(kputll(i64::MIN, &mut s), i64::MIN.to_string().len());
        assert_eq!(kputl(-4, &mut s), 2);

        kinsert_char(b'X', 0, &mut s).unwrap();
        kinsert_str("YZ", 1, &mut s).unwrap();
        assert!(s.as_bytes().starts_with(b"XYZ"));

        assert_eq!(kmemmem(b"a\0bc", b"\0b"), Some(1));
        assert_eq!(kstrstr("foo", "oo"), Some(1));
        assert_eq!(kstrnstr(b"foo\0bar", b"bar", 7), None);

        let released = ks_release(&mut s);
        assert!(!released.is_empty());
        assert!(s.is_empty());

        kputs("line\n", &mut s);
        ks_clear(&mut s);
        assert!(s.is_empty());
        kputs("free", &mut s);
        ks_free(&mut s);
        assert!(s.is_empty());

        kputs("reset", &mut s);
        ks_initialize(&mut s);
        assert!(s.is_empty());

        let mut chunks = [b"abc\n".as_slice()].into_iter();
        assert!(kgetline_from_chunks(&mut s, &mut chunks));
        assert_eq!(s.as_bytes(), b"abc");
    }
}
