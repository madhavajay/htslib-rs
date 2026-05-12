use htslib_rs::kroundup::{
    kroundup_size_t, kroundup32, kroundup64, roundup_i32, roundup_size_t, roundup_u32, roundup_u64,
};

#[test]
fn ports_kroundup_unsigned_cases() {
    for exp in 0..u64::BITS {
        let expected = 1u64 << exp;
        let deltas: &[i128] = if exp > 1 { &[-1, 0, 1] } else { &[0] };

        for &delta in deltas {
            let input = (expected as i128 + delta) as u64;
            let actual = roundup_u64(input);
            let alias_actual = kroundup64(input);

            if delta <= 0 {
                assert_eq!(actual, expected, "roundup_u64({input:#x})");
                assert_eq!(alias_actual, expected, "kroundup64({input:#x})");
            } else if exp < 63 {
                assert_eq!(actual, expected * 2, "roundup_u64({input:#x})");
                assert_eq!(alias_actual, expected * 2, "kroundup64({input:#x})");
            } else {
                assert_eq!(actual, u64::MAX, "roundup_u64({input:#x})");
                assert_eq!(alias_actual, u64::MAX, "kroundup64({input:#x})");
            }
        }
    }

    for exp in 0..u32::BITS {
        let expected = 1u32 << exp;
        let deltas: &[i64] = if exp > 1 { &[-1, 0, 1] } else { &[0] };

        for &delta in deltas {
            let input = (expected as i64 + delta) as u32;
            let actual = roundup_u32(input);
            let alias_actual = kroundup32(input);

            if delta <= 0 {
                assert_eq!(actual, expected, "roundup_u32({input:#x})");
                assert_eq!(alias_actual, expected, "kroundup32({input:#x})");
            } else if exp < 31 {
                assert_eq!(actual, expected * 2, "roundup_u32({input:#x})");
                assert_eq!(alias_actual, expected * 2, "kroundup32({input:#x})");
            } else {
                assert_eq!(actual, u32::MAX, "roundup_u32({input:#x})");
                assert_eq!(alias_actual, u32::MAX, "kroundup32({input:#x})");
            }
        }
    }
}

#[test]
fn ports_kroundup_size_t_and_signed_cases() {
    assert_eq!(roundup_size_t(0), 0);
    assert_eq!(roundup_size_t(3), 4);
    assert_eq!(kroundup_size_t(3), 4);
    assert_eq!(
        roundup_size_t((1usize << (usize::BITS - 1)) + 1),
        usize::MAX
    );
    assert_eq!(
        kroundup_size_t((1usize << (usize::BITS - 1)) + 1),
        usize::MAX
    );

    assert_eq!(roundup_i32(0), 0);
    assert_eq!(roundup_i32(-1), 0);
    assert_eq!(roundup_i32(3), 4);
    assert_eq!(roundup_i32((1 << 30) + 1), i32::MAX as u32);
}
