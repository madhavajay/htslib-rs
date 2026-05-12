use htslib_rs::kbitset::{
    BitSet, BitSetIter, KBS_ELTBITS, kbs_clear, kbs_delete, kbs_elt, kbs_exists, kbs_init,
    kbs_init2, kbs_insert, kbs_insert_all, kbs_last_mask, kbs_mask, kbs_next, kbs_resize,
    kbs_resize2, kbs_start,
};

#[test]
fn ports_kbitset_operations() {
    let mut bitset = BitSet::new(100);

    bitset.insert(5);
    bitset.insert(68);
    bitset.delete(37);

    assert!(bitset.exists(5));
    assert!(bitset.exists(68));
    assert!(!bitset.exists(37));
    assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![5, 68]);

    bitset.insert_all();
    assert_eq!(bitset.iter().count(), 100);

    bitset.clear();
    assert_eq!(bitset.iter().next(), None);

    bitset.resize_with(103, true);
    assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![100, 101, 102]);

    bitset.resize(2);
    assert_eq!(bitset.len(), 2);
    assert_eq!(bitset.iter().next(), None);
}

#[test]
fn ports_kbitset_c_shaped_api_items() {
    assert_eq!(kbs_elt(KBS_ELTBITS + 1), 1);
    assert_eq!(kbs_mask(KBS_ELTBITS + 1), 2);
    assert_eq!(kbs_last_mask(0), usize::MAX);
    assert_eq!(kbs_last_mask(5), 31);

    let mut bitset = kbs_init2(8, true);
    assert_eq!(
        bitset.iter().collect::<Vec<_>>(),
        (0..8).collect::<Vec<_>>()
    );

    kbs_clear(&mut bitset);
    assert_eq!(bitset.iter().next(), None);

    kbs_insert(&mut bitset, 1);
    kbs_insert(&mut bitset, 6);
    kbs_delete(&mut bitset, 1);
    assert!(!kbs_exists(&bitset, 1));
    assert!(kbs_exists(&bitset, 6));

    let mut iter = BitSetIter::default();
    kbs_start(&mut iter);
    assert_eq!(kbs_next(&bitset, &mut iter), Some(6));
    assert_eq!(kbs_next(&bitset, &mut iter), None);

    kbs_resize2(&mut bitset, 10, true);
    assert!(kbs_exists(&bitset, 8));
    assert!(kbs_exists(&bitset, 9));

    kbs_resize(&mut bitset, 3);
    assert_eq!(bitset.len(), 3);

    let mut empty = kbs_init(3);
    kbs_insert_all(&mut empty);
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![0, 1, 2]);
}
