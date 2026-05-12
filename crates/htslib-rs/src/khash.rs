//! Rust replacements for HTSlib khash use cases.

use std::collections::HashMap;

/// A string-to-integer map used by the ported khash tests.
#[derive(Clone, Debug, Default)]
pub struct StrIntMap {
    inner: HashMap<String, u32>,
}

impl StrIntMap {
    /// Inserts a key-value pair.
    ///
    /// Returns `true` when the key was newly inserted and `false` when an
    /// existing value was replaced.
    pub fn insert(&mut self, key: impl Into<String>, value: u32) -> bool {
        self.inner.insert(key.into(), value).is_none()
    }

    /// Gets a value by key.
    pub fn get(&self, key: &str) -> Option<u32> {
        self.inner.get(key).copied()
    }

    /// Deletes a value by key.
    pub fn remove(&mut self, key: &str) -> Option<u32> {
        self.inner.remove(key)
    }

    /// Clears all entries.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Returns whether a key exists.
    pub fn contains_key(&self, key: &str) -> bool {
        self.inner.contains_key(key)
    }

    /// Returns an iterator over key-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&str, u32)> + '_ {
        self.inner.iter().map(|(key, value)| (key.as_str(), *value))
    }

    /// Returns the number of entries in the map.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns whether the map is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

/// HTSlib-style initializer for a string-to-`u32` map.
pub fn kh_init_str_int() -> StrIntMap {
    StrIntMap::default()
}

/// HTSlib-style destroy operation for an owned Rust map.
pub fn kh_destroy_str_int(_map: StrIntMap) {}

/// HTSlib-style clear operation.
pub fn kh_clear_str_int(map: &mut StrIntMap) {
    map.clear();
}

/// HTSlib-style put operation.
///
/// Returns `true` when the key was newly inserted.
pub fn kh_put_str_int(map: &mut StrIntMap, key: impl Into<String>, value: u32) -> bool {
    map.insert(key, value)
}

/// HTSlib-style get operation.
pub fn kh_get_str_int(map: &StrIntMap, key: &str) -> Option<u32> {
    map.get(key)
}

/// HTSlib-style delete operation.
pub fn kh_del_str_int(map: &mut StrIntMap, key: &str) -> Option<u32> {
    map.remove(key)
}

/// HTSlib-style existence check.
pub fn kh_exist_str_int(map: &StrIntMap, key: &str) -> bool {
    map.contains_key(key)
}

/// HTSlib-style size accessor.
pub fn kh_size_str_int(map: &StrIntMap) -> usize {
    map.len()
}

/// HTSlib-style iteration helper over key-value pairs.
pub fn kh_foreach_str_int(map: &StrIntMap) -> impl Iterator<Item = (&str, u32)> + '_ {
    map.iter()
}

/// A Rust replacement for HTSlib's `khash_str2int_*` wrapper API.
#[derive(Clone, Debug, Default)]
pub struct Str2Int {
    inner: HashMap<String, i32>,
}

impl Str2Int {
    /// Returns whether the key exists.
    pub fn has_key(&self, key: &str) -> bool {
        self.inner.contains_key(key)
    }

    /// Gets a value by key.
    ///
    /// This mirrors `khash_str2int_get`: missing keys return `None`.
    pub fn get(&self, key: &str) -> Option<i32> {
        self.inner.get(key).copied()
    }

    /// Adds a key with the next sequential integer if it is missing.
    ///
    /// Existing keys return their current value.
    pub fn inc(&mut self, key: impl Into<String>) -> i32 {
        let len = self.inner.len();
        *self
            .inner
            .entry(key.into())
            .or_insert_with(|| i32::try_from(len).unwrap_or(i32::MAX))
    }

    /// Sets a key to an explicit value.
    ///
    /// Returns `true` when the key was newly inserted and `false` when an
    /// existing value was replaced.
    pub fn set(&mut self, key: impl Into<String>, value: i32) -> bool {
        self.inner.insert(key.into(), value).is_none()
    }

    /// Returns the number of entries in the map.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns whether the map is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

/// HTSlib-style `khash_str2int_init`.
pub fn khash_str2int_init() -> Str2Int {
    Str2Int::default()
}

/// HTSlib-style `khash_str2int_destroy`.
pub fn khash_str2int_destroy(_map: Str2Int) {}

/// HTSlib-style `khash_str2int_destroy_free`.
///
/// Rust owns and drops keys with the map, so this is equivalent to destroy.
pub fn khash_str2int_destroy_free(map: Str2Int) {
    khash_str2int_destroy(map);
}

/// HTSlib-style `khash_str2int_has_key`.
pub fn khash_str2int_has_key(map: &Str2Int, key: &str) -> bool {
    map.has_key(key)
}

/// HTSlib-style `khash_str2int_get`.
pub fn khash_str2int_get(map: &Str2Int, key: &str) -> Option<i32> {
    map.get(key)
}

/// HTSlib-style `khash_str2int_inc`.
pub fn khash_str2int_inc(map: &mut Str2Int, key: impl Into<String>) -> i32 {
    map.inc(key)
}

/// HTSlib-style `khash_str2int_set`.
///
/// Returns `true` when the key was newly inserted.
pub fn khash_str2int_set(map: &mut Str2Int, key: impl Into<String>, value: i32) -> bool {
    map.set(key, value)
}

/// HTSlib-style `khash_str2int_size`.
pub fn khash_str2int_size(map: &Str2Int) -> usize {
    map.len()
}

#[cfg(test)]
mod tests {
    use super::{
        Str2Int, StrIntMap, kh_clear_str_int, kh_del_str_int, kh_destroy_str_int, kh_exist_str_int,
        kh_foreach_str_int, kh_get_str_int, kh_init_str_int, kh_put_str_int, kh_size_str_int,
        khash_str2int_destroy, khash_str2int_destroy_free, khash_str2int_get,
        khash_str2int_has_key, khash_str2int_inc, khash_str2int_init, khash_str2int_set,
        khash_str2int_size,
    };

    #[test]
    fn test_str_int_map() {
        let mut map = StrIntMap::default();

        assert!(map.insert("test0", 0));
        assert!(!map.insert("test0", 1));
        assert_eq!(map.get("test0"), Some(1));
        assert_eq!(map.remove("test0"), Some(1));
        assert_eq!(map.get("test0"), None);
        assert!(map.is_empty());
    }

    #[test]
    fn test_str2int_wrapper_semantics() {
        let mut map = Str2Int::default();

        assert!(!map.has_key("sample0"));
        assert_eq!(map.get("sample0"), None);
        assert_eq!(map.inc("sample0"), 0);
        assert_eq!(map.inc("sample1"), 1);
        assert_eq!(map.inc("sample0"), 0);
        assert!(map.has_key("sample1"));
        assert_eq!(map.get("sample1"), Some(1));
        assert_eq!(map.len(), 2);

        assert!(!map.set("sample1", 42));
        assert_eq!(map.get("sample1"), Some(42));
        assert!(map.set("sample2", -7));
        assert_eq!(map.get("sample2"), Some(-7));
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn test_c_shaped_str_int_aliases() {
        let mut map = kh_init_str_int();

        assert!(kh_put_str_int(&mut map, "a", 1));
        assert!(kh_put_str_int(&mut map, "b", 2));
        assert!(!kh_put_str_int(&mut map, "a", 3));

        assert_eq!(kh_get_str_int(&map, "a"), Some(3));
        assert!(kh_exist_str_int(&map, "b"));
        assert_eq!(kh_size_str_int(&map), 2);

        let mut pairs = kh_foreach_str_int(&map).collect::<Vec<_>>();
        pairs.sort_unstable();
        assert_eq!(pairs, vec![("a", 3), ("b", 2)]);

        assert_eq!(kh_del_str_int(&mut map, "a"), Some(3));
        assert!(!kh_exist_str_int(&map, "a"));

        kh_clear_str_int(&mut map);
        assert_eq!(kh_size_str_int(&map), 0);

        kh_destroy_str_int(map);
    }

    #[test]
    fn test_c_shaped_str2int_aliases() {
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
}
