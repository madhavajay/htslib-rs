//! Rust replacement for HTSlib klist use cases.

use std::{collections::VecDeque, iter::FusedIterator};

/// A FIFO list backed by Rust's standard `VecDeque`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct KList<T> {
    inner: VecDeque<T>,
}

impl<T> KList<T> {
    /// Creates an empty list.
    pub fn new() -> Self {
        Self {
            inner: VecDeque::new(),
        }
    }

    /// Returns the number of elements.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Pushes a value to the back of the list.
    pub fn push(&mut self, value: T) {
        self.inner.push_back(value);
    }

    /// Pushes a default value and returns a mutable reference to it.
    ///
    /// This is the Rust analogue of `kl_pushp`, where callers receive a place
    /// to fill the newly appended node.
    pub fn push_default(&mut self) -> &mut T
    where
        T: Default,
    {
        self.inner.push_back(T::default());
        self.inner
            .back_mut()
            .expect("newly pushed list element must exist")
    }

    /// Removes and returns the front element.
    pub fn shift(&mut self) -> Option<T> {
        self.inner.pop_front()
    }

    /// Returns an iterator over list elements.
    pub fn iter(&self) -> Iter<'_, T> {
        Iter {
            inner: self.inner.iter(),
        }
    }
}

impl<T> IntoIterator for KList<T> {
    type Item = T;
    type IntoIter = std::collections::vec_deque::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}

/// An iterator over a `KList`.
pub struct Iter<'a, T> {
    inner: std::collections::vec_deque::Iter<'a, T>,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

impl<T> ExactSizeIterator for Iter<'_, T> {
    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl<T> FusedIterator for Iter<'_, T> {}

/// HTSlib-style list initializer.
pub fn kl_init<T>() -> KList<T> {
    KList::new()
}

/// HTSlib-style destroy operation.
///
/// Rust drops the owned list at the end of this function.
pub fn kl_destroy<T>(_list: KList<T>) {}

/// HTSlib-style push-by-value helper for Rust callers.
pub fn kl_push<T>(list: &mut KList<T>, value: T) {
    list.push(value);
}

/// HTSlib-style push-and-fill helper.
pub fn kl_pushp<T>(list: &mut KList<T>) -> &mut T
where
    T: Default,
{
    list.push_default()
}

/// HTSlib-style shift operation.
pub fn kl_shift<T>(list: &mut KList<T>) -> Option<T> {
    list.shift()
}

/// HTSlib-style begin operation as a Rust iterator.
pub fn kl_begin<T>(list: &KList<T>) -> Iter<'_, T> {
    list.iter()
}

#[cfg(test)]
mod tests {
    use super::{KList, kl_begin, kl_destroy, kl_init, kl_push, kl_pushp, kl_shift};

    #[test]
    fn test_klist_fifo_operations() {
        let mut list = KList::new();

        assert!(list.is_empty());
        list.push(1);
        *list.push_default() = 2;
        list.push(3);

        assert_eq!(list.len(), 3);
        assert_eq!(list.iter().copied().collect::<Vec<_>>(), vec![1, 2, 3]);
        assert_eq!(list.shift(), Some(1));
        assert_eq!(list.shift(), Some(2));
        assert_eq!(list.shift(), Some(3));
        assert_eq!(list.shift(), None);
        assert!(list.is_empty());
    }

    #[test]
    fn test_c_shaped_aliases() {
        let mut list = kl_init();
        kl_push(&mut list, 10);
        *kl_pushp(&mut list) = 20;
        kl_push(&mut list, 30);

        assert_eq!(
            kl_begin(&list).copied().collect::<Vec<_>>(),
            vec![10, 20, 30]
        );
        assert_eq!(kl_shift(&mut list), Some(10));
        assert_eq!(kl_shift(&mut list), Some(20));
        assert_eq!(kl_shift(&mut list), Some(30));
        assert_eq!(kl_shift(&mut list), None);

        kl_destroy(list);
    }
}
