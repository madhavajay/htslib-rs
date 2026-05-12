use htslib_rs::text::{IntParse, str_to_int, str_to_uint, strprint};

#[test]
fn ports_test_str2int_integer_boundaries() {
    for bits in 1..64 {
        let max = (1u64 << bits) - 1;
        let min_offset = if bits < 5 {
            -((1i128) << (bits - 1))
        } else {
            -16
        };

        for offset in min_offset..=30 {
            let expected_failed = offset > 0;
            let input_value = (max as i128 + offset) as u64;
            let input = format!("{input_value}#");
            let actual = str_to_uint(&input, bits);

            assert_eq!(
                actual,
                IntParse {
                    value: if expected_failed { max } else { input_value },
                    end: input.len() - 1,
                    failed: expected_failed,
                },
                "str_to_uint {bits} bits {input}"
            );
        }

        for offset in min_offset..=30 {
            let expected_failed = offset > 0;
            let input_value = (max as i128 + offset) as u64;
            let input = format!("{input_value}#");
            let actual = str_to_int(&input, bits + 1);

            assert_eq!(
                actual,
                IntParse {
                    value: if expected_failed {
                        max as i64
                    } else {
                        input_value as i64
                    },
                    end: input.len() - 1,
                    failed: expected_failed,
                },
                "positive str_to_int {} bits {input}",
                bits + 1
            );
        }

        for offset in min_offset..=30 {
            let expected_failed = offset > 0;
            let input_magnitude = (max as i128 + offset + 1) as u64;
            let input = format!("-{input_magnitude}#");
            let actual = str_to_int(&input, bits + 1);
            let expected_magnitude = if expected_failed {
                max + 1
            } else {
                input_magnitude
            };

            assert_eq!(actual.end, input.len() - 1, "negative end {bits} bits");
            assert_eq!(
                actual.failed, expected_failed,
                "negative failed {bits} bits"
            );
            assert_eq!(
                actual.value.unsigned_abs(),
                expected_magnitude,
                "negative magnitude {} bits {input}",
                bits + 1
            );
        }
    }
}

#[test]
fn ports_test_str2int_uint64_max_special_case() {
    for offset in 0..=999 {
        let expected_failed = offset > 615;
        let input = format!("18446744073709551{offset:03}#");
        let actual = str_to_uint(&input, 64);

        assert_eq!(
            actual,
            IntParse {
                value: if expected_failed {
                    u64::MAX
                } else {
                    18446744073709551000u64 + offset
                },
                end: input.len() - 1,
                failed: expected_failed,
            },
            "uint64 max offset {offset}"
        );
    }
}

#[test]
fn ports_test_strprint_expectations() {
    assert_eq!(strprint(b"chr10", 9, None), "chr10");
    assert_eq!(strprint(b"chr10", 6, None), "chr10");
    assert_eq!(strprint(b"chr10", 5, None), "c...");
    assert_eq!(strprint(b"chr10", 4, None), "...");
    assert_eq!(strprint(b"tab\twxyz", 10, None), "tab\\twxyz");
    assert_eq!(strprint(b"tab\twxyz", 9, None), "tab\\t...");
    assert_eq!(strprint(b"tab\twxyz", 8, None), "tab\\...");
    assert_eq!(strprint(b"tab\twxyz", 7, None), "tab...");
    assert_eq!(strprint(b"tab\twxyz", 6, None), "ta...");
    assert_eq!(strprint(b"\xab", 5, None), "\\xAB");
    assert_eq!(strprint(b"\xab", 4, None), "...");
    assert_eq!(strprint(b"hello\xff", 40, None), "hello\\xFF");
    assert_eq!(strprint(b"hello\xff", 10, None), "hello\\xFF");
    assert_eq!(strprint(b"hello\xff", 9, None), "hello...");
    assert_eq!(strprint(b"hello\t", 40, None), "hello\\t");
    assert_eq!(strprint(b"hello\t", 8, None), "hello\\t");
    assert_eq!(strprint(b"hello\t", 7, None), "hel...");
    assert_eq!(strprint(b"\t", 40, None), "\\t");
    assert_eq!(strprint(b"", 40, None), "");

    assert_eq!(strprint(b"chr10", 9, Some(b'\'')), "'chr10'");
    assert_eq!(strprint(b"chr10", 8, Some(b'\'')), "'chr10'");
    assert_eq!(strprint(b"chr10", 7, Some(b'\'')), "'c'...");
    assert_eq!(strprint(b"chr10", 6, Some(b'\'')), "''...");
    assert_eq!(strprint(b"quo'wxyz", 12, Some(b'\'')), "'quo\\'wxyz'");
    assert_eq!(strprint(b"quo'wxyz", 11, Some(b'\'')), "'quo\\''...");
    assert_eq!(strprint(b"quo'wxyz", 10, Some(b'\'')), "'quo\\'...");

    assert_eq!(strprint(b"foo", 10, None), "foo");
    assert_eq!(strprint(b"foo\0bar", 10, None), "foo\\0bar");
    assert_eq!(strprint(b"foo\0bar", 9, None), "foo\\0bar");
    assert_eq!(strprint(b"foo\0bar", 8, None), "foo\\...");
}
