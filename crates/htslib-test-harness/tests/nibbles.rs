use htslib_rs::sequence::{SEQ_NT16_STR, bam_seqi, nibble_to_bases};

fn nibble_to_bases_single(nibble: &[u8], len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| SEQ_NT16_STR[usize::from(bam_seqi(nibble, i))])
        .collect()
}

#[test]
fn ports_test_nibbles_validation_loop() {
    let nibble: Vec<_> = (0..5000).map(|i| (i % 256) as u8).collect();

    for start in 0..80 {
        for len in 0..400 {
            let expected = nibble_to_bases_single(&nibble[start..], len);
            let actual = nibble_to_bases(&nibble[start..], len);

            assert_eq!(actual, expected, "start={start}, len={len}");
        }
    }
}

#[test]
fn decodes_all_nt16_symbols() {
    assert_eq!(nibble_to_bases(&[0x01, 0x23, 0x45, 0x67], 8), b"=ACMGRSV");
    assert_eq!(nibble_to_bases(&[0x89, 0xab, 0xcd, 0xef], 8), b"TWYHKDBN");
}
