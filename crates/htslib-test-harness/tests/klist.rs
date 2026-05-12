use htslib_rs::klist::{KList, kl_begin, kl_destroy, kl_init, kl_push, kl_pushp, kl_shift};

#[test]
fn ports_klist_fifo_pushp_and_shift_semantics() {
    let mut list = KList::new();

    list.push(10);
    *list.push_default() = 20;
    list.push(30);

    assert_eq!(list.len(), 3);
    assert_eq!(list.iter().copied().collect::<Vec<_>>(), vec![10, 20, 30]);
    assert_eq!(list.shift(), Some(10));
    assert_eq!(list.shift(), Some(20));
    assert_eq!(list.shift(), Some(30));
    assert_eq!(list.shift(), None);
    assert!(list.is_empty());
}

#[test]
fn ports_klist_c_shaped_api_items() {
    let mut list = kl_init();

    kl_push(&mut list, 7);
    *kl_pushp(&mut list) = 11;
    kl_push(&mut list, 13);

    assert_eq!(
        kl_begin(&list).copied().collect::<Vec<_>>(),
        vec![7, 11, 13]
    );
    assert_eq!(kl_shift(&mut list), Some(7));
    assert_eq!(kl_shift(&mut list), Some(11));
    assert_eq!(kl_shift(&mut list), Some(13));
    assert_eq!(kl_shift(&mut list), None);

    kl_destroy(list);
}
