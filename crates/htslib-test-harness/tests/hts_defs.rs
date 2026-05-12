use htslib_rs::hts_defs::hts_prefetch;

#[test]
fn ports_hts_defs_prefetch_call_shape() {
    let value = 42_u32;
    hts_prefetch(&value);

    let values = [1_u8, 2, 3, 4];
    hts_prefetch(values.as_slice());
}
