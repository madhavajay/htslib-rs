use htslib_rs::variant_io_compat::{
    SyncedPairLogic, SyncedVariantGroup, pair_synced_variant_groups,
};

fn canonical_rows(rows: Vec<Vec<Option<String>>>) -> Vec<Vec<Option<String>>> {
    let mut rows = rows;
    rows.sort();
    rows
}

#[test]
fn ports_bcf_sr_exact_pairing_logic() {
    let groups = [
        SyncedVariantGroup {
            variants: vec!["A>C".into()],
            input_indexes: vec![0],
        },
        SyncedVariantGroup {
            variants: vec!["A>C".into()],
            input_indexes: vec![1],
        },
        SyncedVariantGroup {
            variants: vec!["A>G".into()],
            input_indexes: vec![2],
        },
    ];

    let rows = pair_synced_variant_groups(&groups, SyncedPairLogic::Exact);

    assert_eq!(
        rows,
        vec![
            vec![Some("C".into()), Some("C".into()), None],
            vec![None, None, Some("G".into())],
        ]
    );
}

#[test]
fn ports_bcf_sr_basic_pairing_modes() {
    let groups = [
        SyncedVariantGroup {
            variants: vec!["A>C".into()],
            input_indexes: vec![0],
        },
        SyncedVariantGroup {
            variants: vec!["A>G".into()],
            input_indexes: vec![1],
        },
        SyncedVariantGroup {
            variants: vec!["A>AT".into()],
            input_indexes: vec![2],
        },
        SyncedVariantGroup {
            variants: vec!["A>AC".into()],
            input_indexes: vec![3],
        },
        SyncedVariantGroup {
            variants: vec!["A>.".into()],
            input_indexes: vec![4],
        },
    ];

    assert_eq!(
        canonical_rows(pair_synced_variant_groups(&groups, SyncedPairLogic::Snps)),
        canonical_rows(vec![
            vec![Some("C".into()), Some("G".into()), None, None, None],
            vec![None, None, Some("AT".into()), None, None],
            vec![None, None, None, Some("AC".into()), None],
            vec![None, None, None, None, Some(".".into())],
        ])
    );

    assert_eq!(
        canonical_rows(pair_synced_variant_groups(&groups, SyncedPairLogic::Indels)),
        canonical_rows(vec![
            vec![Some("C".into()), None, None, None, None],
            vec![None, Some("G".into()), None, None, None],
            vec![None, None, Some("AT".into()), Some("AC".into()), None],
            vec![None, None, None, None, Some(".".into())],
        ])
    );

    assert_eq!(
        canonical_rows(pair_synced_variant_groups(&groups, SyncedPairLogic::Both)),
        canonical_rows(vec![
            vec![Some("C".into()), Some("G".into()), None, None, None],
            vec![None, None, Some("AT".into()), Some("AC".into()), None],
            vec![None, None, None, None, Some(".".into())],
        ])
    );
}

#[test]
fn ports_bcf_sr_subset_and_logic_specific_pairing() {
    let groups = [
        SyncedVariantGroup {
            variants: vec!["A>C,A>G".into()],
            input_indexes: vec![0],
        },
        SyncedVariantGroup {
            variants: vec!["A>C".into()],
            input_indexes: vec![1],
        },
        SyncedVariantGroup {
            variants: vec!["A>AT".into()],
            input_indexes: vec![2],
        },
    ];

    let some_rows = pair_synced_variant_groups(&groups, SyncedPairLogic::Some);
    assert_eq!(
        some_rows,
        vec![
            vec![Some("C,G".into()), Some("C".into()), None],
            vec![None, None, Some("AT".into())],
        ]
    );

    let all_rows = pair_synced_variant_groups(&groups, SyncedPairLogic::All);
    assert_eq!(
        all_rows,
        vec![vec![
            Some("C,G".into()),
            Some("C".into()),
            Some("AT".into())
        ]]
    );
}

#[test]
fn ports_bcf_sr_reference_pairing_modes() {
    let groups = [
        SyncedVariantGroup {
            variants: vec!["A>C".into()],
            input_indexes: vec![0],
        },
        SyncedVariantGroup {
            variants: vec!["A>.".into()],
            input_indexes: vec![1],
        },
        SyncedVariantGroup {
            variants: vec!["A>AT".into()],
            input_indexes: vec![2],
        },
    ];

    let snps_ref_rows = pair_synced_variant_groups(&groups, SyncedPairLogic::SnpsAndReference);
    assert_eq!(
        snps_ref_rows,
        vec![
            vec![Some("C".into()), Some(".".into()), None],
            vec![None, None, Some("AT".into())],
        ]
    );

    let indels_ref_rows = pair_synced_variant_groups(&groups, SyncedPairLogic::IndelsAndReference);
    assert_eq!(
        indels_ref_rows,
        vec![
            vec![Some("C".into()), None, None],
            vec![None, Some(".".into()), Some("AT".into())],
        ]
    );

    let both_ref_rows = pair_synced_variant_groups(&groups, SyncedPairLogic::BothAndReference);
    assert_eq!(
        both_ref_rows,
        vec![
            vec![Some("C".into()), Some(".".into()), None],
            vec![None, None, Some("AT".into())],
        ]
    );
}

#[test]
fn ports_bcf_sr_randomized_script_shape_deterministically() {
    let groups = randomized_script_shape_groups();

    let cases = [
        (
            SyncedPairLogic::Snps,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "G", "-", "-", "C", "T,C"]),
                row(["-", "-", ".", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::Indels,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "-", "-", "-", "C", "T,C"]),
                row(["-", "G", "-", "-", "-", "-"]),
                row(["-", "-", ".", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::Both,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "G", "-", "-", "C", "T,C"]),
                row(["-", "-", ".", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::SnpsAndReference,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "G", ".", "-", "C", "T,C"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::IndelsAndReference,
            vec![
                row(["AG", "-", ".", "AG", "AG", "-"]),
                row(["C", "-", "-", "-", "C", "T,C"]),
                row(["-", "G", "-", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::BothAndReference,
            vec![
                row(["AG", "-", ".", "AG", "AG", "-"]),
                row(["C", "G", "-", "-", "C", "T,C"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::Exact,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "-", "-", "-", "C", "-"]),
                row(["-", "G", "-", "-", "-", "-"]),
                row(["-", "-", ".", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
                row(["-", "-", "-", "-", "-", "T,C"]),
            ],
        ),
        (
            SyncedPairLogic::Some,
            vec![
                row(["AG", "-", "-", "AG", "AG", "-"]),
                row(["C", "-", "-", "-", "C", "T,C"]),
                row(["-", "G", "-", "-", "-", "-"]),
                row(["-", "-", ".", "-", "-", "-"]),
                row(["-", "-", "-", "AT", "-", "-"]),
            ],
        ),
        (
            SyncedPairLogic::All,
            vec![
                row(["AG", "G", ".", "AG", "AG", "T,C"]),
                row(["C", "-", "-", "AT", "C", "-"]),
            ],
        ),
    ];

    for (logic, expected) in cases {
        assert_eq!(
            pair_synced_variant_groups(&groups, logic),
            expected,
            "{logic:?}"
        );
    }
}

#[test]
fn ports_bcf_sr_randomized_script_shape_is_stable_across_group_order() {
    let groups = randomized_script_shape_groups();
    let shuffled_groups = [
        groups[3].clone(),
        groups[1].clone(),
        groups[4].clone(),
        groups[0].clone(),
        groups[2].clone(),
    ];

    for logic in [
        SyncedPairLogic::Snps,
        SyncedPairLogic::Indels,
        SyncedPairLogic::Both,
        SyncedPairLogic::SnpsAndReference,
        SyncedPairLogic::IndelsAndReference,
        SyncedPairLogic::BothAndReference,
        SyncedPairLogic::Exact,
        SyncedPairLogic::Some,
        SyncedPairLogic::All,
    ] {
        assert_eq!(
            canonical_rows(pair_synced_variant_groups(&groups, logic)),
            canonical_rows(pair_synced_variant_groups(&shuffled_groups, logic)),
            "{logic:?}"
        );
    }
}

#[test]
fn ports_bcf_sr_randomized_script_shape_is_stable_across_variant_and_input_order() {
    let groups = randomized_script_shape_groups();
    let shuffled_within_groups = [
        SyncedVariantGroup {
            variants: vec!["A>AG".into(), "A>C".into()],
            input_indexes: vec![4, 0],
        },
        groups[1].clone(),
        groups[2].clone(),
        SyncedVariantGroup {
            variants: vec!["A>AT".into(), "A>AG".into()],
            input_indexes: vec![3],
        },
        groups[4].clone(),
    ];

    for logic in [
        SyncedPairLogic::Snps,
        SyncedPairLogic::Indels,
        SyncedPairLogic::Both,
        SyncedPairLogic::SnpsAndReference,
        SyncedPairLogic::IndelsAndReference,
        SyncedPairLogic::BothAndReference,
        SyncedPairLogic::Exact,
        SyncedPairLogic::Some,
        SyncedPairLogic::All,
    ] {
        assert_eq!(
            canonical_rows(pair_synced_variant_groups(&groups, logic)),
            canonical_rows(pair_synced_variant_groups(&shuffled_within_groups, logic)),
            "{logic:?}"
        );
    }
}

#[test]
fn ports_bcf_sr_randomized_script_duplicate_input_order_with_multiple_shapes() {
    let cases = [
        vec![
            SyncedVariantGroup {
                variants: vec!["C>T".into(), "C>CA".into(), "C>.".into()],
                input_indexes: vec![0, 6],
            },
            SyncedVariantGroup {
                variants: vec!["C>G".into(), "C>CT".into()],
                input_indexes: vec![1],
            },
            SyncedVariantGroup {
                variants: vec!["C>T,C>G".into()],
                input_indexes: vec![2, 7],
            },
            SyncedVariantGroup {
                variants: vec!["C>CA".into(), "C>CG".into()],
                input_indexes: vec![3],
            },
            SyncedVariantGroup {
                variants: vec!["C>.".into(), "C>CT".into()],
                input_indexes: vec![4, 5],
            },
        ],
        vec![
            SyncedVariantGroup {
                variants: vec!["G>A".into(), "G>GT".into()],
                input_indexes: vec![2],
            },
            SyncedVariantGroup {
                variants: vec!["G>.".into()],
                input_indexes: vec![0, 4],
            },
            SyncedVariantGroup {
                variants: vec!["G>GA".into(), "G>A,G>C".into()],
                input_indexes: vec![1],
            },
            SyncedVariantGroup {
                variants: vec!["G>C".into(), "G>GT".into(), "G>GAA".into()],
                input_indexes: vec![3, 5],
            },
        ],
    ];

    for groups in cases {
        let reversed_input_groups = groups
            .iter()
            .map(|group| {
                let mut input_indexes = group.input_indexes.clone();
                input_indexes.reverse();

                SyncedVariantGroup {
                    variants: group.variants.clone(),
                    input_indexes,
                }
            })
            .collect::<Vec<_>>();

        for logic in all_pair_logics() {
            let expected = canonical_rows(pair_synced_variant_groups(&groups, logic));

            assert_eq!(
                canonical_rows(pair_synced_variant_groups(&reversed_input_groups, logic)),
                expected,
                "{logic:?}"
            );
        }
    }
}

fn all_pair_logics() -> [SyncedPairLogic; 9] {
    [
        SyncedPairLogic::Snps,
        SyncedPairLogic::Indels,
        SyncedPairLogic::Both,
        SyncedPairLogic::SnpsAndReference,
        SyncedPairLogic::IndelsAndReference,
        SyncedPairLogic::BothAndReference,
        SyncedPairLogic::Exact,
        SyncedPairLogic::Some,
        SyncedPairLogic::All,
    ]
}

fn randomized_script_shape_groups() -> Vec<SyncedVariantGroup> {
    vec![
        SyncedVariantGroup {
            variants: vec!["A>C".into(), "A>AG".into()],
            input_indexes: vec![0, 4],
        },
        SyncedVariantGroup {
            variants: vec!["A>G".into()],
            input_indexes: vec![1],
        },
        SyncedVariantGroup {
            variants: vec!["A>.".into()],
            input_indexes: vec![2],
        },
        SyncedVariantGroup {
            variants: vec!["A>AG".into(), "A>AT".into()],
            input_indexes: vec![3],
        },
        SyncedVariantGroup {
            variants: vec!["A>T,A>C".into()],
            input_indexes: vec![5],
        },
    ]
}

fn row<const N: usize>(values: [&str; N]) -> Vec<Option<String>> {
    values
        .into_iter()
        .map(|value| (value != "-").then(|| value.to_string()))
        .collect()
}
