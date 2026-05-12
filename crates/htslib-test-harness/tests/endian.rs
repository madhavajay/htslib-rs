use htslib_rs::endian::{
    double_to_le, f32_to_le, f64_to_le, float_to_le, i16_to_le, i32_to_le, i64_to_le, le_to_double,
    le_to_f32, le_to_f64, le_to_float, le_to_i8, le_to_i16, le_to_i32, le_to_i64, le_to_u8,
    le_to_u16, le_to_u32, le_to_u64, u16_to_le, u32_to_le, u64_to_le,
};

#[test]
fn ports_hts_endian_8_bit_cases() {
    let cases = [
        ([0x00], 0, 0),
        ([0x01], 1, 1),
        ([0x7f], 127, 127),
        ([0x80], -128, 128),
        ([0xff], -1, 255),
    ];

    for (bytes, signed, unsigned) in cases {
        assert_eq!(le_to_u8(&bytes), unsigned);
        assert_eq!(le_to_i8(&bytes), signed);

        let mut unaligned = [0; 2];
        unaligned[1] = bytes[0];
        assert_eq!(le_to_u8(&unaligned[1..]), unsigned);
        assert_eq!(le_to_i8(&unaligned[1..]), signed);
    }
}

#[test]
fn ports_hts_endian_16_bit_cases() {
    let cases = [
        ([0x00, 0x00], 0, 0),
        ([0x01, 0x00], 1, 1),
        ([0x00, 0x01], 256, 256),
        ([0xff, 0x7f], 32767, 32767),
        ([0x00, 0x80], -32768, 32768),
        ([0xff, 0xff], -1, 65535),
    ];

    for (bytes, signed, unsigned) in cases {
        assert_eq!(le_to_u16(&bytes), unsigned);
        assert_eq!(le_to_i16(&bytes), signed);

        let mut unaligned = [0; 3];
        unaligned[1..].copy_from_slice(&bytes);
        assert_eq!(le_to_u16(&unaligned[1..]), unsigned);
        assert_eq!(le_to_i16(&unaligned[1..]), signed);

        let mut buf = [0; 3];
        u16_to_le(unsigned, &mut buf);
        assert_eq!(&buf[..2], &bytes);
        i16_to_le(signed, &mut buf);
        assert_eq!(&buf[..2], &bytes);

        u16_to_le(unsigned, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
        i16_to_le(signed, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
    }
}

#[test]
fn ports_hts_endian_32_bit_cases() {
    let cases = [
        ([0x00, 0x00, 0x00, 0x00], 0, 0),
        ([0x01, 0x00, 0x00, 0x00], 1, 1),
        ([0x00, 0x01, 0x00, 0x00], 256, 256),
        ([0x00, 0x00, 0x01, 0x00], 65536, 65536),
        ([0xff, 0xff, 0xff, 0x7f], i32::MAX, i32::MAX as u32),
        ([0x00, 0x00, 0x00, 0x80], i32::MIN, 2147483648),
        ([0xff, 0xff, 0xff, 0xff], -1, u32::MAX),
    ];

    for (bytes, signed, unsigned) in cases {
        assert_eq!(le_to_u32(&bytes), unsigned);
        assert_eq!(le_to_i32(&bytes), signed);

        let mut unaligned = [0; 5];
        unaligned[1..].copy_from_slice(&bytes);
        assert_eq!(le_to_u32(&unaligned[1..]), unsigned);
        assert_eq!(le_to_i32(&unaligned[1..]), signed);

        let mut buf = [0; 5];
        u32_to_le(unsigned, &mut buf);
        assert_eq!(&buf[..4], &bytes);
        i32_to_le(signed, &mut buf);
        assert_eq!(&buf[..4], &bytes);

        u32_to_le(unsigned, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
        i32_to_le(signed, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
    }
}

#[test]
fn ports_hts_endian_64_bit_cases() {
    let cases = [
        ([0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], 0, 0),
        ([0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], 1, 1),
        ([0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], 256, 256),
        (
            [0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            65536,
            65536,
        ),
        (
            [0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00],
            4294967296,
            4294967296,
        ),
        (
            [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f],
            i64::MAX,
            i64::MAX as u64,
        ),
        (
            [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80],
            i64::MIN,
            9223372036854775808,
        ),
        (
            [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
            -1,
            u64::MAX,
        ),
    ];

    for (bytes, signed, unsigned) in cases {
        assert_eq!(le_to_u64(&bytes), unsigned);
        assert_eq!(le_to_i64(&bytes), signed);

        let mut unaligned = [0; 9];
        unaligned[1..].copy_from_slice(&bytes);
        assert_eq!(le_to_u64(&unaligned[1..]), unsigned);
        assert_eq!(le_to_i64(&unaligned[1..]), signed);

        let mut buf = [0; 9];
        u64_to_le(unsigned, &mut buf);
        assert_eq!(&buf[..8], &bytes);
        i64_to_le(signed, &mut buf);
        assert_eq!(&buf[..8], &bytes);

        u64_to_le(unsigned, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
        i64_to_le(signed, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
    }
}

#[test]
fn ports_hts_endian_float_cases() {
    let cases = [
        ([0x00, 0x00, 0x00, 0x00], 0.0f32),
        ([0x00, 0x00, 0x80, 0x3f], 1.0),
        ([0x00, 0x00, 0x80, 0xbf], -1.0),
        ([0x00, 0x00, 0x20, 0x41], 10.0),
        ([0xd0, 0x0f, 0x49, 0x40], f32::from_bits(0x40490fd0)),
        ([0xa8, 0x0a, 0xff, 0x66], 6.022e23),
        ([0xcd, 0x84, 0x03, 0x13], 1.66e-27),
    ];

    for (bytes, value) in cases {
        assert_eq!(le_to_f32(&bytes), value);
        assert_eq!(le_to_float(&bytes), value);

        let mut unaligned = [0; 5];
        unaligned[1..].copy_from_slice(&bytes);
        assert_eq!(le_to_f32(&unaligned[1..]), value);
        assert_eq!(le_to_float(&unaligned[1..]), value);

        let mut buf = [0; 5];
        f32_to_le(value, &mut buf);
        assert_eq!(&buf[..4], &bytes);
        float_to_le(value, &mut buf);
        assert_eq!(&buf[..4], &bytes);

        f32_to_le(value, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
        float_to_le(value, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
    }
}

#[test]
fn ports_hts_endian_double_cases() {
    let cases = [
        ([0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00], 0.0f64),
        ([0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f], 1.0),
        ([0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf], -1.0),
        ([0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x40], 10.0),
        (
            [0x18, 0x2d, 0x44, 0x54, 0xfb, 0x21, 0x09, 0x40],
            std::f64::consts::PI,
        ),
        (
            [0x2b, 0x08, 0x0c, 0xd3, 0x85, 0xe1, 0xdf, 0x44],
            6.022140858e23,
        ),
        (
            [0x55, 0xfa, 0x81, 0x74, 0xf7, 0x71, 0x60, 0x3a],
            1.66053904e-27,
        ),
    ];

    for (bytes, value) in cases {
        assert_eq!(le_to_f64(&bytes), value);
        assert_eq!(le_to_double(&bytes), value);

        let mut unaligned = [0; 9];
        unaligned[1..].copy_from_slice(&bytes);
        assert_eq!(le_to_f64(&unaligned[1..]), value);
        assert_eq!(le_to_double(&unaligned[1..]), value);

        let mut buf = [0; 9];
        f64_to_le(value, &mut buf);
        assert_eq!(&buf[..8], &bytes);
        double_to_le(value, &mut buf);
        assert_eq!(&buf[..8], &bytes);

        f64_to_le(value, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
        double_to_le(value, &mut buf[1..]);
        assert_eq!(&buf[1..], &bytes);
    }
}
