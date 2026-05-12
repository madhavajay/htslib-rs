//! HTSlib-compatible little-endian conversion helpers.

/// Reads a little-endian `u8`.
pub fn le_to_u8(buf: &[u8]) -> u8 {
    buf[0]
}

/// Reads a little-endian `u16`.
pub fn le_to_u16(buf: &[u8]) -> u16 {
    u16::from_le_bytes([buf[0], buf[1]])
}

/// Reads a little-endian `u32`.
pub fn le_to_u32(buf: &[u8]) -> u32 {
    u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]])
}

/// Reads a little-endian `u64`.
pub fn le_to_u64(buf: &[u8]) -> u64 {
    u64::from_le_bytes([
        buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
    ])
}

/// Reads a two's-complement little-endian `i8`.
pub fn le_to_i8(buf: &[u8]) -> i8 {
    i8::from_le_bytes([buf[0]])
}

/// Reads a two's-complement little-endian `i16`.
pub fn le_to_i16(buf: &[u8]) -> i16 {
    i16::from_le_bytes([buf[0], buf[1]])
}

/// Reads a two's-complement little-endian `i32`.
pub fn le_to_i32(buf: &[u8]) -> i32 {
    i32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]])
}

/// Reads a two's-complement little-endian `i64`.
pub fn le_to_i64(buf: &[u8]) -> i64 {
    i64::from_le_bytes([
        buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
    ])
}

/// Reads an IEEE 754 little-endian `f32`.
pub fn le_to_f32(buf: &[u8]) -> f32 {
    f32::from_bits(le_to_u32(buf))
}

/// Reads an IEEE 754 little-endian `f32`.
pub fn le_to_float(buf: &[u8]) -> f32 {
    le_to_f32(buf)
}

/// Reads an IEEE 754 little-endian `f64`.
pub fn le_to_f64(buf: &[u8]) -> f64 {
    f64::from_bits(le_to_u64(buf))
}

/// Reads an IEEE 754 little-endian `f64`.
pub fn le_to_double(buf: &[u8]) -> f64 {
    le_to_f64(buf)
}

/// Writes a little-endian `u16`.
pub fn u16_to_le(value: u16, buf: &mut [u8]) {
    buf[..2].copy_from_slice(&value.to_le_bytes());
}

/// Writes a little-endian `u32`.
pub fn u32_to_le(value: u32, buf: &mut [u8]) {
    buf[..4].copy_from_slice(&value.to_le_bytes());
}

/// Writes a little-endian `u64`.
pub fn u64_to_le(value: u64, buf: &mut [u8]) {
    buf[..8].copy_from_slice(&value.to_le_bytes());
}

/// Writes a two's-complement little-endian `i16`.
pub fn i16_to_le(value: i16, buf: &mut [u8]) {
    buf[..2].copy_from_slice(&value.to_le_bytes());
}

/// Writes a two's-complement little-endian `i32`.
pub fn i32_to_le(value: i32, buf: &mut [u8]) {
    buf[..4].copy_from_slice(&value.to_le_bytes());
}

/// Writes a two's-complement little-endian `i64`.
pub fn i64_to_le(value: i64, buf: &mut [u8]) {
    buf[..8].copy_from_slice(&value.to_le_bytes());
}

/// Writes an IEEE 754 little-endian `f32`.
pub fn f32_to_le(value: f32, buf: &mut [u8]) {
    u32_to_le(value.to_bits(), buf);
}

/// Writes an IEEE 754 little-endian `f32`.
pub fn float_to_le(value: f32, buf: &mut [u8]) {
    f32_to_le(value, buf);
}

/// Writes an IEEE 754 little-endian `f64`.
pub fn f64_to_le(value: f64, buf: &mut [u8]) {
    u64_to_le(value.to_bits(), buf);
}

/// Writes an IEEE 754 little-endian `f64`.
pub fn double_to_le(value: f64, buf: &mut [u8]) {
    f64_to_le(value, buf);
}

#[cfg(test)]
mod tests {
    use super::{i16_to_le, le_to_i16, le_to_u16, u16_to_le};

    #[test]
    fn test_16_bit_helpers() {
        assert_eq!(le_to_u16(&[0xff, 0xff]), u16::MAX);
        assert_eq!(le_to_i16(&[0x00, 0x80]), i16::MIN);

        let mut buf = [0; 3];
        u16_to_le(0x7fff, &mut buf[1..]);
        assert_eq!(&buf[1..3], &[0xff, 0x7f]);

        i16_to_le(-1, &mut buf[1..]);
        assert_eq!(&buf[1..3], &[0xff, 0xff]);
    }
}
