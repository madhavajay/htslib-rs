use htslib_rs::kstring::{
    KString, getline_from_chunks, kgetline_from_chunks, kinsert_char, kinsert_str, kmemmem, kputc,
    kputc_, kputl, kputll, kputs, kputsn, kputsn_, kputuw, kputw, ks_c_str, ks_clear, ks_expand,
    ks_free, ks_initialize, ks_len, ks_release, ks_resize, ks_str, kstrnstr, kstrstr, memmem,
    roundup_i32, roundup_size_t, strnstr, strstr,
};

fn check_u32_range(start: u32, end: u32) {
    for value in start..=end {
        let mut s = KString::new();
        s.push_u32(value);
        assert_eq!(s.as_bytes(), value.to_string().as_bytes());
    }
}

fn check_i32_range(start: i32, end: i32) {
    for value in start..=end {
        let mut s = KString::new();
        s.push_i32(value);
        assert_eq!(s.as_bytes(), value.to_string().as_bytes());
    }
}

fn check_i64_range(start: i64, end: i64) {
    for value in start..=end {
        let mut s = KString::new();
        s.push_i64(value);
        assert_eq!(s.as_bytes(), value.to_string().as_bytes());
    }
}

#[test]
fn ports_test_kstring_roundup_cases() {
    assert_eq!(roundup_size_t(0), 0);

    for exp in 0..usize::BITS {
        let expected = 1usize << exp;
        let deltas: &[isize] = if exp > 1 { &[-1, 0, 1] } else { &[0] };

        for &delta in deltas {
            let input = expected.wrapping_add_signed(delta);
            let actual = roundup_size_t(input);

            if delta <= 0 {
                assert_eq!(actual, expected, "roundup_size_t({input:#x})");
            } else {
                let expected = expected.saturating_mul(2);
                assert_eq!(actual, expected, "roundup_size_t({input:#x})");
            }
        }
    }

    for exp in 0..31 {
        let expected = 1u32 << exp;
        let deltas: &[i32] = if exp > 1 { &[-1, 0, 1] } else { &[0] };

        for &delta in deltas {
            let input = expected as i32 + delta;
            let actual = roundup_i32(input);

            if delta <= 0 {
                assert_eq!(actual, expected, "roundup_i32({input})");
            } else {
                let expected = if exp < 30 {
                    expected * 2
                } else {
                    ((expected - 1) << 1) | 1
                };
                assert_eq!(actual, expected, "roundup_i32({input})");
            }
        }
    }
}

#[test]
fn ports_test_kstring_integer_format_cases() {
    let mut value = 0u32;
    loop {
        let start = value.saturating_sub(5);
        let end = value.saturating_add(5);
        check_u32_range(start, end);

        if value == 0 {
            value = 1;
        } else if value > u32::MAX / 10 {
            break;
        } else {
            value *= 10;
        }
    }

    check_u32_range(u32::MAX - 5, u32::MAX);

    let mut value = 1i32;
    while value <= i32::MAX / 10 {
        check_i32_range((value - 5).max(0), value + 5);
        value *= 10;
    }

    let mut value = -1i32;
    while value >= i32::MIN / 10 {
        check_i32_range(value - 5, (value + 5).min(0));
        value *= 10;
    }

    check_i32_range(i32::MAX - 5, i32::MAX);
    check_i32_range(i32::MIN, i32::MIN + 5);

    let mut value = 1i64;
    while value <= (i64::MAX - 5) / 10 {
        check_i64_range(if value >= 5 { value - 5 } else { value }, value);
        value *= 10;
    }

    let mut value = 1i64;
    while value <= (i64::MAX - 5) / 10 {
        let negative = -value;
        check_i64_range(negative, negative);
        value *= 10;
    }

    check_i64_range(i64::MAX - 5, i64::MAX);
    check_i64_range(i64::MIN, i64::MIN + 5);
}

#[test]
fn ports_test_kgetline_cases() {
    let chunks = [
        b"ABCD".as_slice(),
        b"\n",
        b"\n",
        b"ABCD",
        b"\r\n",
        b"\r\n",
        b"ABCD",
    ];
    let mut chunks = chunks.into_iter();
    let mut s = KString::new();

    s.push_bytes(b"_");
    assert!(getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"_ABCD");

    s.clear();
    assert!(getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"");

    s.clear();
    assert!(getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"ABCD");

    s.clear();
    assert!(getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"");

    s.clear();
    assert!(getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"ABCD");

    s.clear();
    assert!(!getline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"");
}

#[test]
fn ports_test_kstring_insertion_cases() {
    let expected = [
        None,
        Some("X0123"),
        Some("0X123"),
        Some("01X23"),
        Some("012X3"),
        Some("0123X"),
        None,
    ];

    for (i, expected) in (-1..6).zip(expected) {
        let mut s = KString::new();
        s.push_bytes(b"0123");

        let result = usize::try_from(i)
            .ok()
            .map(|i| s.insert_byte(i, b'X'))
            .unwrap_or(Err(htslib_rs::kstring::InsertError::OutOfBounds));

        match expected {
            Some(expected) => {
                result.unwrap();
                assert_eq!(s.as_bytes(), expected.as_bytes());
            }
            None => assert!(result.is_err()),
        }
    }

    let expected = [
        None,
        Some("XYZ0123"),
        Some("0XYZ123"),
        Some("01XYZ23"),
        Some("012XYZ3"),
        Some("0123XYZ"),
        None,
    ];

    for (i, expected) in (-1..6).zip(expected) {
        let mut s = KString::new();
        s.push_bytes(b"0123");

        let result = usize::try_from(i)
            .ok()
            .map(|i| s.insert_bytes(i, b"XYZ"))
            .unwrap_or(Err(htslib_rs::kstring::InsertError::OutOfBounds));

        match expected {
            Some(expected) => {
                result.unwrap();
                assert_eq!(s.as_bytes(), expected.as_bytes());
            }
            None => assert!(result.is_err()),
        }
    }
}

#[test]
fn ports_test_kstring_search_cases() {
    let mem_cases = [
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"f\0\0".as_slice(),
            Some(0),
        ),
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"\0\0f".as_slice(),
            Some(1),
        ),
        (
            b"\0\0f\0\0f\0\0fbar\0\0f\0\0f".as_slice(),
            b"\0\0f".as_slice(),
            Some(0),
        ),
        (
            b"\0\0f\0\0f\0\0fbar\0\0f\0\0f".as_slice(),
            b"f\0\0".as_slice(),
            Some(2),
        ),
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"bar".as_slice(),
            Some(9),
        ),
        (
            b"f\0\0f\0\0f\0\0baz\0\0f\0\0f".as_slice(),
            b"bar".as_slice(),
            None,
        ),
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"".as_slice(),
            Some(0),
        ),
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"\0\0b".as_slice(),
            Some(7),
        ),
        (
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            b"r\0\0".as_slice(),
            Some(11),
        ),
        (
            b"bar".as_slice(),
            b"f\0\0f\0\0f\0\0bar\0\0f\0\0f".as_slice(),
            None,
        ),
        (b"".as_slice(), b"bar".as_slice(), None),
        (b"".as_slice(), b"".as_slice(), Some(0)),
    ];

    for (haystack, needle, expected) in mem_cases {
        assert_eq!(memmem(haystack, needle), expected);
    }

    let strstr_cases = [
        ("foofoofoobaroofoof", "bar", Some(9)),
        ("foofoofoobazoofoof", "bar", None),
        ("foofoofoobaroofoof", "", Some(0)),
        ("foofoofoobaroofoof", "oob", Some(7)),
        ("foofoofoobaroofoof", "roo", Some(11)),
        ("bar", "foofoofoobaroofoof", None),
        ("", "bar", None),
        ("", "", Some(0)),
    ];

    for (haystack, needle, expected) in strstr_cases {
        assert_eq!(strstr(haystack, needle), expected);
    }

    let strnstr_cases = [
        (
            b"foofoofoobaroofoof".as_slice(),
            b"bar".as_slice(),
            18,
            Some(9),
        ),
        (
            b"foofoofoobazoofoof".as_slice(),
            b"bar".as_slice(),
            18,
            None,
        ),
        (b"foofoofoobaroofoof".as_slice(), b"bar".as_slice(), 9, None),
        (
            b"foofoofoobaroofoof".as_slice(),
            b"".as_slice(),
            18,
            Some(0),
        ),
        (
            b"bar".as_slice(),
            b"foofoofoobaroofoof".as_slice(),
            18,
            None,
        ),
        (
            b"foofoof\0obaroofoof".as_slice(),
            b"bar".as_slice(),
            18,
            None,
        ),
        (b"".as_slice(), b"bar".as_slice(), 3, None),
        (b"".as_slice(), b"".as_slice(), 0, Some(0)),
    ];

    for (haystack, needle, limit, expected) in strnstr_cases {
        assert_eq!(strnstr(haystack, needle, limit), expected);
    }
}

#[test]
fn ports_kstring_c_shaped_api_items() {
    let mut s = KString::new();

    ks_resize(&mut s, 32);
    assert!(s.capacity() >= 32);
    ks_expand(&mut s, 4);

    assert_eq!(kputs("abc", &mut s), 3);
    assert_eq!(kputsn(b"def", &mut s), 3);
    assert_eq!(kputc(b'g', &mut s), b'g');
    assert_eq!(kputc_(b'h', &mut s), 1);
    assert_eq!(kputsn_(b"ij", &mut s), 2);
    assert_eq!(ks_str(&s), b"abcdefghij");
    assert_eq!(ks_c_str(&s), b"abcdefghij");
    assert_eq!(ks_len(&s), 10);

    assert_eq!(kputuw(u32::MAX, &mut s), u32::MAX.to_string().len());
    assert_eq!(kputw(i32::MIN, &mut s), i32::MIN.to_string().len());
    assert_eq!(kputll(i64::MIN, &mut s), i64::MIN.to_string().len());
    assert_eq!(kputl(-123, &mut s), 4);

    kinsert_char(b'X', 0, &mut s).unwrap();
    kinsert_str("YZ", 1, &mut s).unwrap();
    assert!(s.as_bytes().starts_with(b"XYZ"));

    assert_eq!(kmemmem(b"a\0bc", b"\0b"), Some(1));
    assert_eq!(kstrstr("foo", "oo"), Some(1));
    assert_eq!(kstrnstr(b"foo\0bar", b"bar", 7), None);

    let released = ks_release(&mut s);
    assert!(!released.is_empty());
    assert!(s.is_empty());

    kputs("clear", &mut s);
    ks_clear(&mut s);
    assert!(s.is_empty());

    kputs("free", &mut s);
    ks_free(&mut s);
    assert!(s.is_empty());

    kputs("reset", &mut s);
    ks_initialize(&mut s);
    assert!(s.is_empty());

    let mut chunks = [b"abc\r\n".as_slice()].into_iter();
    assert!(kgetline_from_chunks(&mut s, &mut chunks));
    assert_eq!(s.as_bytes(), b"abc");
}
