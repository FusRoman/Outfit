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
}
