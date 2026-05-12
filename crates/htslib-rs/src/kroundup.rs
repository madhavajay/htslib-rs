//! Rust replacements for HTSlib kroundup helpers.

/// Rounds an unsigned 64-bit value up to the next power of two.
///
/// Values that would overflow saturate to `u64::MAX`, matching HTSlib's
/// `kroundup64` macro.
pub fn roundup_u64(value: u64) -> u64 {
    if value == 0 {
        0
    } else {
        value.checked_next_power_of_two().unwrap_or(u64::MAX)
    }
}

/// C-shaped alias for HTSlib's `kroundup64` macro.
pub fn kroundup64(value: u64) -> u64 {
    roundup_u64(value)
}

/// Rounds an unsigned 32-bit value up to the next power of two.
///
/// Values that would overflow saturate to `u32::MAX`.
pub fn roundup_u32(value: u32) -> u32 {
    if value == 0 {
        0
    } else {
        value.checked_next_power_of_two().unwrap_or(u32::MAX)
    }
}

/// C-shaped alias for HTSlib's `kroundup32` macro.
pub fn kroundup32(value: u32) -> u32 {
    roundup_u32(value)
}

/// Rounds a `usize` value up to the next power of two.
///
/// Values that would overflow saturate to `usize::MAX`.
pub fn roundup_size_t(value: usize) -> usize {
    if value == 0 {
        0
    } else {
        value.checked_next_power_of_two().unwrap_or(usize::MAX)
    }
}

/// C-shaped alias for HTSlib's `kroundup_size_t` macro.
pub fn kroundup_size_t(value: usize) -> usize {
    roundup_size_t(value)
}

/// Rounds a signed 32-bit value to the next allocation boundary.
///
/// Nonpositive values return zero. Values above `1 << 30` saturate to
/// `i32::MAX`, preserving the historical signed 32-bit behavior used by
/// HTSlib's tests.
pub fn roundup_i32(value: i32) -> u32 {
    if value <= 0 {
        0
    } else if value <= 1 << 30 {
        (value as u32).next_power_of_two()
    } else {
        i32::MAX as u32
    }
}

#[cfg(test)]
mod tests {
    use super::{
        kroundup_size_t, kroundup32, kroundup64, roundup_i32, roundup_size_t, roundup_u32,
        roundup_u64,
    };

    #[test]
    fn test_unsigned_roundup() {
        assert_eq!(roundup_u32(0), 0);
        assert_eq!(kroundup32(0), 0);
        assert_eq!(roundup_u32(1), 1);
        assert_eq!(kroundup32(1), 1);
        assert_eq!(roundup_u32(3), 4);
        assert_eq!(kroundup32(3), 4);
        assert_eq!(roundup_u32((1 << 31) + 1), u32::MAX);
        assert_eq!(kroundup32((1 << 31) + 1), u32::MAX);

        assert_eq!(roundup_u64(0), 0);
        assert_eq!(kroundup64(0), 0);
        assert_eq!(roundup_u64(3), 4);
        assert_eq!(kroundup64(3), 4);
        assert_eq!(roundup_u64((1 << 63) + 1), u64::MAX);
        assert_eq!(kroundup64((1 << 63) + 1), u64::MAX);
        assert_eq!(roundup_size_t(3), 4);
        assert_eq!(kroundup_size_t(3), 4);
    }

    #[test]
    fn test_signed_roundup() {
        assert_eq!(roundup_i32(-1), 0);
        assert_eq!(roundup_i32(3), 4);
        assert_eq!(roundup_i32((1 << 30) + 1), i32::MAX as u32);
    }
}
