use htslib_rs::khash::{
    Str2Int, StrIntMap, kh_clear_str_int, kh_del_str_int, kh_destroy_str_int, kh_exist_str_int,
    kh_foreach_str_int, kh_get_str_int, kh_init_str_int, kh_put_str_int, kh_size_str_int,
    khash_str2int_destroy, khash_str2int_destroy_free, khash_str2int_get, khash_str2int_has_key,
    khash_str2int_inc, khash_str2int_init, khash_str2int_set, khash_str2int_size,
};

fn make_keys(max: usize) -> Vec<String> {
    (0..max).map(|i| format!("test{i}")).collect()
}

fn roundup_power_of_two(mut value: usize) -> usize {
    value = value.saturating_sub(1);
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;

    if usize::BITS > 32 {
        value |= value >> 32;
    }

    value + 1
}

fn delete_victims(max: usize, to_delete: usize) -> Vec<bool> {
    let mut flags = vec![false; max];
    let mask = roundup_power_of_two(max) - 1;
    let mut r = 0x533du32;

    for _ in 0..to_delete {
        let victim = loop {
            r = (r >> 1) ^ ((r & 1).wrapping_mul(0x80000057));
            let raw = ((r as usize) & mask).wrapping_sub(1);

            if raw < max && !flags[raw] {
                break raw;
            }
        };

        flags[victim] = true;
    }

    flags
}

#[test]
fn ports_test_khash_default_str2int_case() {
    let max = 1000;
    let to_delete = max / 4;
    let keys = make_keys(max);
    let mut map = StrIntMap::default();

    for (i, key) in keys.iter().enumerate() {
        assert!(map.insert(key, i as u32), "failed to insert {key}");
    }

    assert_eq!(map.len(), max);

    for (i, key) in keys.iter().enumerate() {
        assert_eq!(map.get(key), Some(i as u32), "missing {key}");
    }

    let deleted = delete_victims(max, to_delete);

    for (i, key) in keys.iter().enumerate() {
        if deleted[i] {
            assert_eq!(map.remove(key), Some(i as u32), "failed to delete {key}");
        }
    }

    assert_eq!(map.len(), max - to_delete);

    for (i, key) in keys.iter().enumerate() {
        let expected = (!deleted[i]).then_some(i as u32);
        assert_eq!(map.get(key), expected, "unexpected state for {key}");
    }

    for (i, key) in keys.iter().enumerate() {
        if deleted[i] {
            assert!(map.insert(key, i as u32), "failed to reinsert {key}");
        }
    }

    assert_eq!(map.len(), max);

    for (i, key) in keys.iter().enumerate() {
        assert_eq!(map.get(key), Some(i as u32), "missing after reinsert {key}");
    }
}

#[test]
fn ports_khash_str2int_wrapper_semantics() {
    let mut map = Str2Int::default();

    assert!(!map.has_key("A"));
    assert_eq!(map.get("A"), None);
    assert_eq!(map.inc("A"), 0);
    assert_eq!(map.inc("B"), 1);
    assert_eq!(map.inc("A"), 0);
    assert!(map.has_key("B"));
    assert_eq!(map.get("B"), Some(1));
    assert_eq!(map.len(), 2);

    assert!(!map.set("B", 9));
    assert_eq!(map.get("B"), Some(9));
    assert!(map.set("C", -1));
    assert_eq!(map.get("C"), Some(-1));
    assert_eq!(map.len(), 3);
}

#[test]
fn ports_khash_c_shaped_str_int_api_items() {
    let mut map = kh_init_str_int();

    assert!(kh_put_str_int(&mut map, "test0", 0));
    assert!(kh_put_str_int(&mut map, "test1", 1));
    assert!(!kh_put_str_int(&mut map, "test0", 10));

    assert_eq!(kh_get_str_int(&map, "test0"), Some(10));
    assert!(kh_exist_str_int(&map, "test1"));
    assert_eq!(kh_size_str_int(&map), 2);

    let mut pairs = kh_foreach_str_int(&map).collect::<Vec<_>>();
    pairs.sort_unstable();
    assert_eq!(pairs, vec![("test0", 10), ("test1", 1)]);

    assert_eq!(kh_del_str_int(&mut map, "test1"), Some(1));
    assert!(!kh_exist_str_int(&map, "test1"));

    kh_clear_str_int(&mut map);
    assert_eq!(kh_size_str_int(&map), 0);

    kh_destroy_str_int(map);
}

#[test]
fn ports_khash_str2int_c_shaped_api_items() {
    let mut map = khash_str2int_init();

    assert!(!khash_str2int_has_key(&map, "A"));
    assert_eq!(khash_str2int_get(&map, "A"), None);
    assert_eq!(khash_str2int_inc(&mut map, "A"), 0);
    assert_eq!(khash_str2int_inc(&mut map, "B"), 1);
    assert_eq!(khash_str2int_inc(&mut map, "A"), 0);
    assert!(khash_str2int_has_key(&map, "B"));
    assert_eq!(khash_str2int_get(&map, "B"), Some(1));
    assert_eq!(khash_str2int_size(&map), 2);

    assert!(!khash_str2int_set(&mut map, "B", 9));
    assert!(khash_str2int_set(&mut map, "C", -1));
    assert_eq!(khash_str2int_get(&map, "B"), Some(9));
    assert_eq!(khash_str2int_get(&map, "C"), Some(-1));

    khash_str2int_destroy(map);
    khash_str2int_destroy_free(khash_str2int_init());
}
