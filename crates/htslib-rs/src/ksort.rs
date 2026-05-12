//! Rust replacements for selected HTSlib ksort behavior.

use std::cmp::Ordering;

fn ordering_from_less<T>(a: &T, b: &T, is_less: &mut impl FnMut(&T, &T) -> bool) -> Ordering {
    if is_less(a, b) {
        Ordering::Less
    } else if is_less(b, a) {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// HTSlib-style generic less-than predicate.
pub fn ks_lt_generic<T: PartialOrd>(a: &T, b: &T) -> bool {
    a < b
}

/// HTSlib-style string less-than predicate.
pub fn ks_lt_str(a: &str, b: &str) -> bool {
    a < b
}

/// HTSlib-style swap helper.
pub fn ksort_swap<T>(data: &mut [T], a: usize, b: usize) -> bool {
    if a >= data.len() || b >= data.len() {
        return false;
    }

    data.swap(a, b);
    true
}

/// Sorts a slice using HTSlib `ks_mergesort`-style stable ordering.
pub fn mergesort_by<T>(data: &mut [T], mut is_less: impl FnMut(&T, &T) -> bool) {
    data.sort_by(|a, b| ordering_from_less(a, b, &mut is_less));
}

/// C-shaped alias for `ks_mergesort`.
pub fn ks_mergesort_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    mergesort_by(data, is_less);
}

/// Sorts a slice using HTSlib `ks_introsort`-style unstable ordering.
pub fn introsort_by<T>(data: &mut [T], mut is_less: impl FnMut(&T, &T) -> bool) {
    data.sort_unstable_by(|a, b| ordering_from_less(a, b, &mut is_less));
}

/// C-shaped alias for `ks_introsort`.
pub fn ks_introsort_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    introsort_by(data, is_less);
}

/// Sorts a slice using HTSlib `ks_combsort`-style unstable ordering.
///
/// The C implementation uses combsort as an implementation detail for
/// introsort fallback. The Rust compatibility layer exposes the same
/// observable contract: the slice is sorted according to the supplied
/// less-than predicate.
pub fn combsort_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    introsort_by(data, is_less);
}

/// C-shaped alias for `ks_combsort`.
pub fn ks_combsort_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    combsort_by(data, is_less);
}

/// Repairs the max-heap rooted at `index` for the prefix `data[..heap_len]`.
pub fn heapadjust_by<T>(
    data: &mut [T],
    index: usize,
    heap_len: usize,
    mut is_less: impl FnMut(&T, &T) -> bool,
) {
    if index >= heap_len || heap_len > data.len() {
        return;
    }

    sift_down(data, index, heap_len, &mut is_less);
}

/// C-shaped alias for `ks_heapadjust`.
pub fn ks_heapadjust_by<T>(
    data: &mut [T],
    index: usize,
    heap_len: usize,
    is_less: impl FnMut(&T, &T) -> bool,
) {
    heapadjust_by(data, index, heap_len, is_less);
}

/// Turns the slice into a max-heap according to the supplied less-than predicate.
pub fn heapmake_by<T>(data: &mut [T], mut is_less: impl FnMut(&T, &T) -> bool) {
    for index in (0..(data.len() / 2)).rev() {
        sift_down(data, index, data.len(), &mut is_less);
    }
}

/// C-shaped alias for `ks_heapmake`.
pub fn ks_heapmake_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    heapmake_by(data, is_less);
}

/// Sorts a slice using heap-sort semantics.
pub fn heapsort_by<T>(data: &mut [T], mut is_less: impl FnMut(&T, &T) -> bool) {
    heapmake_by(data, &mut is_less);

    for end in (1..data.len()).rev() {
        data.swap(0, end);
        sift_down(data, 0, end, &mut is_less);
    }
}

/// C-shaped alias for `ks_heapsort`.
pub fn ks_heapsort_by<T>(data: &mut [T], is_less: impl FnMut(&T, &T) -> bool) {
    heapsort_by(data, is_less);
}

fn sift_down<T>(
    data: &mut [T],
    mut root: usize,
    heap_len: usize,
    is_less: &mut impl FnMut(&T, &T) -> bool,
) {
    loop {
        let left = root * 2 + 1;

        if left >= heap_len {
            break;
        }

        let right = left + 1;
        let mut child = left;

        if right < heap_len && is_less(&data[left], &data[right]) {
            child = right;
        }

        if !is_less(&data[root], &data[child]) {
            break;
        }

        data.swap(root, child);
        root = child;
    }
}

/// Partitions around and returns the `k`th smallest element.
///
/// This mirrors `ks_ksmall`: the input slice may be reordered, and `None` is
/// returned for out-of-range indexes instead of invoking undefined behavior.
pub fn ksmall_by<T: Clone>(
    data: &mut [T],
    k: usize,
    mut is_less: impl FnMut(&T, &T) -> bool,
) -> Option<T> {
    if k >= data.len() {
        return None;
    }

    let (_, nth, _) = data.select_nth_unstable_by(k, |a, b| ordering_from_less(a, b, &mut is_less));

    Some(nth.clone())
}

/// C-shaped alias for `ks_ksmall`.
pub fn ks_ksmall_by<T: Clone>(
    data: &mut [T],
    k: usize,
    is_less: impl FnMut(&T, &T) -> bool,
) -> Option<T> {
    ksmall_by(data, k, is_less)
}

/// Shuffles a slice using HTSlib `ks_shuffle`-style Fisher-Yates swaps.
///
/// `next_unit` must return values in `[0.0, 1.0)`, matching `hts_drand48`.
/// Values outside that range are clamped so callers cannot produce an
/// out-of-bounds swap index.
pub fn shuffle_by<T>(data: &mut [T], mut next_unit: impl FnMut() -> f64) {
    for i in (2..=data.len()).rev() {
        let raw = (next_unit() * i as f64) as isize;
        let j = raw.clamp(0, i as isize - 1) as usize;
        data.swap(j, i - 1);
    }
}

/// C-shaped alias for `ks_shuffle`.
pub fn ks_shuffle_by<T>(data: &mut [T], next_unit: impl FnMut() -> f64) {
    shuffle_by(data, next_unit);
}

#[cfg(test)]
mod tests {
    use super::{
        combsort_by, heapadjust_by, heapmake_by, heapsort_by, introsort_by, ks_combsort_by,
        ks_heapadjust_by, ks_heapmake_by, ks_heapsort_by, ks_introsort_by, ks_ksmall_by,
        ks_lt_generic, ks_lt_str, ks_mergesort_by, ks_shuffle_by, ksmall_by, ksort_swap,
        mergesort_by, shuffle_by,
    };

    #[test]
    fn test_sort_variants() {
        let mut values = [9, -1, 4, 4, 0, 12, 3];
        introsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, [-1, 0, 3, 4, 4, 9, 12]);

        let mut values = [9, -1, 4, 4, 0, 12, 3];
        combsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, [-1, 0, 3, 4, 4, 9, 12]);

        let mut values = [9, -1, 4, 4, 0, 12, 3];
        heapsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, [-1, 0, 3, 4, 4, 9, 12]);

        let mut values = [9, -1, 4, 4, 0, 12, 3];
        mergesort_by(&mut values, |a, b| a < b);
        assert_eq!(values, [-1, 0, 3, 4, 4, 9, 12]);
    }

    #[test]
    fn test_heap_and_selection_helpers() {
        let mut values = [3, 1, 4, 1, 5, 9, 2];
        heapmake_by(&mut values, |a, b| a < b);
        assert_eq!(values[0], 9);

        values[0] = 0;
        let len = values.len();
        heapadjust_by(&mut values, 0, len, |a, b| a < b);
        assert_eq!(values[0], 5);

        let kth = ksmall_by(&mut values, 3, |a, b| a < b);
        assert_eq!(kth, Some(2));
    }

    #[test]
    fn test_shuffle_by() {
        let mut values = [0, 1, 2, 3];
        let mut draws = [0.0, 0.5, 0.0].into_iter();

        shuffle_by(&mut values, || draws.next().unwrap());

        assert_eq!(values, [2, 3, 1, 0]);
    }

    #[test]
    fn test_c_shaped_aliases() {
        assert!(ks_lt_generic(&1, &2));
        assert!(ks_lt_str("bam", "sam"));

        let mut values = [1, 2, 3];
        assert!(ksort_swap(&mut values, 0, 2));
        assert_eq!(values, [3, 2, 1]);
        assert!(!ksort_swap(&mut values, 0, 3));

        let expected = [1, 2, 3, 4];

        let mut values = [3, 1, 4, 2];
        ks_mergesort_by(&mut values, |a, b| a < b);
        assert_eq!(values, expected);

        let mut values = [3, 1, 4, 2];
        ks_introsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, expected);

        let mut values = [3, 1, 4, 2];
        ks_combsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, expected);

        let mut values = [3, 1, 4, 2];
        ks_heapsort_by(&mut values, |a, b| a < b);
        assert_eq!(values, expected);

        let mut values = [3, 1, 4, 2];
        ks_heapmake_by(&mut values, |a, b| a < b);
        values[0] = 0;
        let heap_len = values.len();
        ks_heapadjust_by(&mut values, 0, heap_len, |a, b| a < b);
        assert_eq!(values[0], 3);

        let mut values = [9, 7, 5, 3, 1];
        assert_eq!(ks_ksmall_by(&mut values, 2, |a, b| a < b), Some(5));

        let mut values = [0, 1, 2, 3];
        let mut draws = [0.0, 0.5, 0.0].into_iter();
        ks_shuffle_by(&mut values, || draws.next().unwrap());
        assert_eq!(values, [2, 3, 1, 0]);
    }
}
