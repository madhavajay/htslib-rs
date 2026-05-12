use htslib_rs::variant::{VariantType, bcf_acgt2int, bcf_ij2g, bcf_int2acgt, classify_variant};

#[test]
fn ports_test_bcf_set_variant_type_cases() {
    assert_eq!(classify_variant("A", "T").variant_type, VariantType::SNP);

    assert_eq!(
        classify_variant("A", "AA").variant_type,
        VariantType::INDEL | VariantType::INS
    );
    assert_eq!(
        classify_variant("AA", "A").variant_type,
        VariantType::INDEL | VariantType::DEL
    );

    for (reference, alternate) in [
        ("N", "N]16:33625444]"),
        ("N", "N[16:33625444["),
        ("N", "]16:33625444]N"),
        ("N", "[16:33625444[N"),
        ("T", "]chrB:123]AGTNNNNNCAT"),
        ("C", "CAGTNNNNNCA[chrA:321["),
    ] {
        assert_eq!(
            classify_variant(reference, alternate).variant_type,
            VariantType::BND
        );
    }

    assert_eq!(
        classify_variant("A", "<NON_REF>").variant_type,
        VariantType::REF
    );
    assert_eq!(classify_variant("A", "<*>").variant_type, VariantType::REF);

    assert_eq!(classify_variant("AA", "TT").variant_type, VariantType::MNP);
    assert_eq!(
        classify_variant("A", "*").variant_type,
        VariantType::OVERLAP
    );
    assert_eq!(classify_variant("A", ".").variant_type, VariantType::REF);
}

#[test]
fn ports_vcfutils_small_helpers() {
    assert_eq!(bcf_acgt2int('A'), Some(0));
    assert_eq!(bcf_acgt2int('c'), Some(1));
    assert_eq!(bcf_acgt2int('G'), Some(2));
    assert_eq!(bcf_acgt2int('t'), Some(3));
    assert_eq!(bcf_acgt2int('N'), None);

    assert_eq!(bcf_int2acgt(0), Some('A'));
    assert_eq!(bcf_int2acgt(1), Some('C'));
    assert_eq!(bcf_int2acgt(2), Some('G'));
    assert_eq!(bcf_int2acgt(3), Some('T'));
    assert_eq!(bcf_int2acgt(4), None);

    assert_eq!(
        (0..=3)
            .flat_map(|j| (0..=j).map(move |i| bcf_ij2g(i, j)))
            .collect::<Vec<_>>(),
        (0..10).collect::<Vec<_>>()
    );
}
