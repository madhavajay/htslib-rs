use htslib_rs::ksort::{
    combsort_by, heapadjust_by, heapmake_by, heapsort_by, introsort_by, ks_combsort_by,
    ks_heapadjust_by, ks_heapmake_by, ks_heapsort_by, ks_introsort_by, ks_ksmall_by, ks_lt_generic,
    ks_lt_str, ks_mergesort_by, ks_shuffle_by, ksmall_by, ksort_swap, mergesort_by, shuffle_by,
};

#[test]
fn ports_ksort_sorting_variants() {
    let expected = [-5, 0, 1, 2, 3, 3, 8, 13];

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    introsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    combsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    heapsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    mergesort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);
}

#[test]
fn ports_ksort_heap_and_ksmall_helpers() {
    let mut values = [7, 1, 4, 9, 2, 6, 3];

    heapmake_by(&mut values, |a, b| a < b);
    assert_eq!(values[0], 9);

    values[0] = 0;
    let heap_len = values.len();
    heapadjust_by(&mut values, 0, heap_len, |a, b| a < b);
    assert_eq!(values[0], 7);

    let kth = ksmall_by(&mut values, 4, |a, b| a < b);
    assert_eq!(kth, Some(4));
    let len = values.len();
    assert_eq!(ksmall_by(&mut values, len, |a, b| a < b), None);
}

#[test]
fn ports_ksort_string_and_shuffle_helpers() {
    let mut values = ["sam", "bcf", "bam", "cram"];
    introsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, ["bam", "bcf", "cram", "sam"]);

    let mut values = [0, 1, 2, 3];
    let mut draws = [0.0, 0.5, 0.0].into_iter();
    shuffle_by(&mut values, || draws.next().unwrap());
    assert_eq!(values, [2, 3, 1, 0]);
}

#[test]
fn ports_ksort_c_shaped_api_items() {
    assert!(ks_lt_generic(&1, &2));
    assert!(ks_lt_str("bam", "sam"));

    let mut values = [1, 2, 3];
    assert!(ksort_swap(&mut values, 0, 2));
    assert_eq!(values, [3, 2, 1]);
    assert!(!ksort_swap(&mut values, 0, 3));

    let expected = [-5, 0, 1, 2, 3, 3, 8, 13];

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    ks_mergesort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    ks_introsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    ks_combsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut values = [3, 1, 8, -5, 3, 13, 0, 2];
    ks_heapsort_by(&mut values, |a, b| a < b);
    assert_eq!(values, expected);

    let mut heap = [7, 1, 4, 9, 2, 6, 3];
    ks_heapmake_by(&mut heap, |a, b| a < b);
    assert_eq!(heap[0], 9);
    heap[0] = 0;
    let heap_len = heap.len();
    ks_heapadjust_by(&mut heap, 0, heap_len, |a, b| a < b);
    assert_eq!(heap[0], 7);

    let mut values = [7, 1, 4, 9, 2, 6, 3];
    assert_eq!(ks_ksmall_by(&mut values, 4, |a, b| a < b), Some(6));

    let mut values = [0, 1, 2, 3];
    let mut draws = [0.0, 0.5, 0.0].into_iter();
    ks_shuffle_by(&mut values, || draws.next().unwrap());
    assert_eq!(values, [2, 3, 1, 0]);
}
