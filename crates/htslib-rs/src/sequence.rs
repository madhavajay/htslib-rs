//! HTSlib-compatible sequence encoding helpers.

/// BAM nt16 code to base lookup string.
pub const SEQ_NT16_STR: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Returns the BAM nt16 code at base index `i` from a four-bit-packed sequence.
pub fn bam_seqi(sequence: &[u8], i: usize) -> u8 {
    (sequence[i >> 1] >> (((!i) & 1) << 2)) & 0x0f
}

/// Converts a BAM four-bit-packed sequence to HTSlib nt16 base characters.
pub fn nibble_to_bases(sequence: &[u8], len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| SEQ_NT16_STR[usize::from(bam_seqi(sequence, i))])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{bam_seqi, nibble_to_bases};

    #[test]
    fn test_bam_seqi() {
        let src = [0x12, 0x48];

        assert_eq!(bam_seqi(&src, 0), 1);
        assert_eq!(bam_seqi(&src, 1), 2);
        assert_eq!(bam_seqi(&src, 2), 4);
        assert_eq!(bam_seqi(&src, 3), 8);
    }

    #[test]
    fn test_nibble_to_bases() {
        assert_eq!(nibble_to_bases(&[0x12, 0x48], 4), b"ACGT");
    }
}
