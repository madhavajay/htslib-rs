use htslib_rs::hts_os::{Rand48, hts_drand48, hts_erand48, hts_lrand48, hts_srand48};

#[test]
fn ports_hts_os_rand48_deterministic_sequences() {
    let mut rand = Rand48::seeded(0);

    assert_eq!(rand.erand48(), 0.17082803610628972);
    assert_eq!(rand.lrand48(), 1_610_402_240);

    let mut rand = Rand48::seeded(1);

    assert_eq!(rand.erand48(), 0.041630344771878214);
    assert_eq!(rand.lrand48(), 976_015_093);
}

#[test]
fn ports_hts_os_global_and_erand48_state() {
    hts_srand48(0);
    assert_eq!(hts_drand48(), 0.17082803610628972);
    assert_eq!(hts_lrand48(), 1_610_402_240);

    let mut seed = [0x330e, 0, 0];
    assert_eq!(hts_erand48(&mut seed), 0.17082803610628972);
    assert_eq!(seed, [0x5101, 0x62dc, 0x2bbb]);
}
