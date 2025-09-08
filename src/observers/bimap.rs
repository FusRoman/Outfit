use std::collections::HashMap;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    forward: HashMap<K, V>,
    reverse: HashMap<V, K>,
}

impl<K, V> Default for BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<K, V> BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    pub fn new() -> Self {
        Self {
            forward: HashMap::new(),
            reverse: HashMap::new(),
        }
    }

    pub fn entry_or_insert_by_key(&mut self, key: K, value: V) -> &mut V {
        self.forward.entry(key.clone()).or_insert_with(|| {
            self.reverse.insert(value.clone(), key);
            value
        })
    }

    pub fn entry_or_insert_by_value(&mut self, value: V, key: K) -> &mut K {
        self.reverse.entry(value.clone()).or_insert_with(|| {
            self.forward.insert(key.clone(), value);
            key
        })
    }

    pub fn insert(&mut self, key: K, value: V) {
        self.forward.insert(key.clone(), value.clone());
        self.reverse.insert(value, key);
    }

    pub fn get_by_key(&self, key: &K) -> Option<&V> {
        self.forward.get(key)
    }

    pub fn get_by_value(&self, value: &V) -> Option<&K> {
        self.reverse.get(value)
    }

    pub fn remove_by_key(&mut self, key: &K) {
        if let Some(val) = self.forward.remove(key) {
            self.reverse.remove(&val);
        }
    }

    pub fn remove_by_value(&mut self, value: &V) {
        if let Some(key) = self.reverse.remove(value) {
            self.forward.remove(&key);
        }
    }

    pub fn len(&self) -> usize {
        self.forward.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    // ---------------------------
    // Iteration helpers (immutable)
    // ---------------------------

    /// Iterate over (&K, &V) using the forward map.
    ///
    /// Return
    /// ----------
    /// * An iterator yielding `(&K, &V)` pairs. Order is not guaranteed.
    ///
    /// See also
    /// ------------
    /// * [`BiMap::iter_rev`]
    /// * [`BiMap::keys`], [`BiMap::values`]
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.forward.iter()
    }

    /// Iterate over (&V, &K) using the reverse map.
    ///
    /// Return
    /// ----------
    /// * An iterator yielding `(&V, &K)` pairs. Order is not guaranteed.
    ///
    /// See also
    /// ------------
    /// * [`BiMap::iter`]
    pub fn iter_rev(&self) -> impl Iterator<Item = (&V, &K)> {
        self.reverse.iter()
    }

    /// Iterate over keys (&K).
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.forward.keys()
    }

    /// Iterate over values (&V).
    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.forward.values()
    }
}

// ---------------------------
// IntoIterator implementations
// ---------------------------

impl<'a, K, V> IntoIterator for &'a BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    type Item = (&'a K, &'a V);
    type IntoIter = std::collections::hash_map::Iter<'a, K, V>;

    /// Consume `&BiMap` into an iterator over `(&K, &V)` on the forward map.
    fn into_iter(self) -> Self::IntoIter {
        self.forward.iter()
    }
}

impl<'a, K, V> IntoIterator for &'a mut BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    type Item = (&'a K, &'a mut V);
    type IntoIter = std::collections::hash_map::IterMut<'a, K, V>;

    /// Consume `&mut BiMap` into an iterator over `(&K, &mut V)` on the forward map.
    ///
    /// Warning
    /// -------
    /// Mutating values can desynchronize the reverse map if you change logical identity.
    /// Use with care; prefer removing/re-inserting pairs instead.
    fn into_iter(self) -> Self::IntoIter {
        self.forward.iter_mut()
    }
}

impl<K, V> IntoIterator for BiMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    type Item = (K, V);
    type IntoIter = std::collections::hash_map::IntoIter<K, V>;

    /// Consume the bimap and iterate over owned `(K, V)` pairs using the forward map.
    fn into_iter(self) -> Self::IntoIter {
        self.forward.into_iter()
    }
}

#[cfg(test)]
mod bimap_iter_tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn iter_and_iter_rev_cover_same_pairs() {
        let mut m = BiMap::new();
        m.insert("a", 1);
        m.insert("b", 2);
        m.insert("c", 3);

        let fwd: HashSet<_> = m.iter().map(|(k, v)| ((*k).to_string(), *v)).collect();
        let rev: HashSet<_> = m.iter_rev().map(|(v, k)| ((*k).to_string(), *v)).collect();
        assert_eq!(fwd, rev);
    }

    #[test]
    fn into_iterator_by_ref_and_by_value() {
        let mut m = BiMap::new();
        m.insert("x", 10);
        m.insert("y", 20);

        // &BiMap
        let pairs_ref: HashSet<_> = (&m)
            .into_iter()
            .map(|(k, v)| ((*k).to_string(), *v))
            .collect();
        assert!(pairs_ref.contains(&("x".to_string(), 10)));
        assert!(pairs_ref.contains(&("y".to_string(), 20)));

        // BiMap (by value)
        let pairs_val: HashSet<_> = m.into_iter().collect();
        assert!(pairs_val.contains(&("x", 10)));
        assert!(pairs_val.contains(&("y", 20)));
        // m is moved here; no further use
    }

    #[test]
    fn keys_and_values_match_len() {
        let mut m = BiMap::new();
        m.insert(1, "one");
        m.insert(2, "two");
        m.insert(3, "three");

        assert_eq!(m.len(), m.keys().count());
        assert_eq!(m.len(), m.values().count());
    }
}
