//! Rust replacement for HTSlib kbitset use cases.

/// Number of bits in the storage word used by HTSlib's `kbitset.h` macros.
pub const KBS_ELTBITS: usize = usize::BITS as usize;

/// A finite bit set over indexes `0..len`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BitSet {
    bits: Vec<bool>,
}

impl BitSet {
    /// Creates an empty bit set with the given capacity.
    pub fn new(len: usize) -> Self {
        Self {
            bits: vec![false; len],
        }
    }

    /// Creates a bit set with all bits set.
    pub fn full(len: usize) -> Self {
        Self {
            bits: vec![true; len],
        }
    }

    /// Returns the number of representable indexes.
    pub fn len(&self) -> usize {
        self.bits.len()
    }

    /// Returns whether the bit set is empty.
    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }

    /// Resizes the bit set, filling newly added bits with `false`.
    pub fn resize(&mut self, new_len: usize) {
        self.resize_with(new_len, false);
    }

    /// Resizes the bit set, filling newly added bits with `fill`.
    pub fn resize_with(&mut self, new_len: usize, fill: bool) {
        self.bits.resize(new_len, fill);
    }

    /// Clears all bits.
    pub fn clear(&mut self) {
        self.bits.fill(false);
    }

    /// Sets all bits.
    pub fn insert_all(&mut self) {
        self.bits.fill(true);
    }

    /// Inserts an index.
    pub fn insert(&mut self, index: usize) -> bool {
        let Some(bit) = self.bits.get_mut(index) else {
            return false;
        };

        let was_set = *bit;
        *bit = true;
        !was_set
    }

    /// Deletes an index.
    pub fn delete(&mut self, index: usize) -> bool {
        let Some(bit) = self.bits.get_mut(index) else {
            return false;
        };

        let was_set = *bit;
        *bit = false;
        was_set
    }

    /// Returns whether an index is present.
    pub fn exists(&self, index: usize) -> bool {
        self.bits.get(index).copied().unwrap_or(false)
    }

    /// Returns indexes present in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.bits
            .iter()
            .enumerate()
            .filter_map(|(i, is_set)| is_set.then_some(i))
    }
}

/// HTSlib-style storage element index for a bit position.
pub fn kbs_elt(index: usize) -> usize {
    index / KBS_ELTBITS
}

/// HTSlib-style storage mask for a bit position.
pub fn kbs_mask(index: usize) -> usize {
    1usize << (index % KBS_ELTBITS)
}

/// HTSlib-style final storage-word mask for a set capacity.
pub fn kbs_last_mask(len: usize) -> usize {
    let mask = kbs_mask(len) - 1;
    if mask == 0 { usize::MAX } else { mask }
}

/// HTSlib-style initializer with optional full-set fill.
pub fn kbs_init2(len: usize, fill: bool) -> BitSet {
    if fill {
        BitSet::full(len)
    } else {
        BitSet::new(len)
    }
}

/// HTSlib-style empty-set initializer.
pub fn kbs_init(len: usize) -> BitSet {
    kbs_init2(len, false)
}

/// HTSlib-style resize with optional fill for newly added indexes.
pub fn kbs_resize2(bitset: &mut BitSet, new_len: usize, fill: bool) {
    bitset.resize_with(new_len, fill);
}

/// HTSlib-style resize that fills newly added indexes as absent.
pub fn kbs_resize(bitset: &mut BitSet, new_len: usize) {
    kbs_resize2(bitset, new_len, false);
}

/// HTSlib-style clear operation.
pub fn kbs_clear(bitset: &mut BitSet) {
    bitset.clear();
}

/// HTSlib-style full-set operation.
pub fn kbs_insert_all(bitset: &mut BitSet) {
    bitset.insert_all();
}

/// HTSlib-style insert operation.
pub fn kbs_insert(bitset: &mut BitSet, index: usize) {
    bitset.insert(index);
}

/// HTSlib-style delete operation.
pub fn kbs_delete(bitset: &mut BitSet, index: usize) {
    bitset.delete(index);
}

/// HTSlib-style membership check.
pub fn kbs_exists(bitset: &BitSet, index: usize) -> bool {
    bitset.exists(index)
}

/// Iterator state matching `kbitset_iter_t` behavior for Rust callers.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BitSetIter {
    next_index: usize,
}

/// HTSlib-style iterator initializer.
pub fn kbs_start(iter: &mut BitSetIter) {
    iter.next_index = 0;
}

/// HTSlib-style next-element operation.
pub fn kbs_next(bitset: &BitSet, iter: &mut BitSetIter) -> Option<usize> {
    while iter.next_index < bitset.len() {
        let index = iter.next_index;
        iter.next_index += 1;

        if bitset.exists(index) {
            return Some(index);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::{
        BitSet, BitSetIter, KBS_ELTBITS, kbs_clear, kbs_delete, kbs_elt, kbs_exists, kbs_init,
        kbs_init2, kbs_insert, kbs_insert_all, kbs_last_mask, kbs_mask, kbs_next, kbs_resize,
        kbs_resize2, kbs_start,
    };

    #[test]
    fn test_bit_set_operations() {
        let mut bitset = BitSet::new(8);

        assert_eq!(bitset.len(), 8);
        assert!(bitset.insert(5));
        assert!(!bitset.insert(5));
        assert!(bitset.insert(2));
        assert!(bitset.exists(5));
        assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![2, 5]);

        assert!(bitset.delete(5));
        assert!(!bitset.exists(5));
        assert!(!bitset.delete(5));

        bitset.resize_with(10, true);
        assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![2, 8, 9]);

        bitset.resize(1);
        assert_eq!(bitset.iter().collect::<Vec<_>>(), Vec::<usize>::new());

        bitset.insert_all();
        assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![0]);
        bitset.clear();
        assert_eq!(bitset.iter().next(), None);
    }

    #[test]
    fn test_c_shaped_aliases() {
        assert_eq!(kbs_elt(KBS_ELTBITS * 2 + 3), 2);
        assert_eq!(kbs_mask(KBS_ELTBITS * 2 + 3), 8);
        assert_eq!(kbs_last_mask(3), 7);
        assert_eq!(kbs_last_mask(KBS_ELTBITS), usize::MAX);

        let mut bitset = kbs_init(10);
        kbs_insert(&mut bitset, 3);
        kbs_insert(&mut bitset, 9);
        kbs_delete(&mut bitset, 4);

        assert!(kbs_exists(&bitset, 3));
        assert!(kbs_exists(&bitset, 9));
        assert!(!kbs_exists(&bitset, 4));

        let mut iter = BitSetIter::default();
        kbs_start(&mut iter);
        assert_eq!(kbs_next(&bitset, &mut iter), Some(3));
        assert_eq!(kbs_next(&bitset, &mut iter), Some(9));
        assert_eq!(kbs_next(&bitset, &mut iter), None);

        kbs_resize2(&mut bitset, 12, true);
        assert!(kbs_exists(&bitset, 10));
        assert!(kbs_exists(&bitset, 11));

        kbs_clear(&mut bitset);
        assert_eq!(bitset.iter().next(), None);

        kbs_resize(&mut bitset, 3);
        kbs_insert_all(&mut bitset);
        assert_eq!(bitset.iter().collect::<Vec<_>>(), vec![0, 1, 2]);

        assert_eq!(kbs_init2(3, true).iter().collect::<Vec<_>>(), vec![0, 1, 2]);
    }
}
