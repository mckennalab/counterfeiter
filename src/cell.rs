use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use crate::genome::{Genome, GenomeEventKey};

// from https://stackoverflow.com/questions/32935808/generate-sequential-ids-for-each-instance-of-a-struct
// and https://doc.rust-lang.org/std/sync/atomic/
// we want to have a unique ID for each cell, so we can relate them
static OBJECT_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A cell is a collection of events at target sites. This struct is kept as minimal
/// as possible to reduce memory usage, and events and sites are referred to by index.
/// The cells are passed the relevant look-up object when they 'split' to create new cells.
#[derive(Debug)]
pub struct Cell {
    pub id: usize,
    pub events: HashSet<GenomeEventKey>,
    pub visible_in_output: bool,
}

impl Cell {
    pub fn new() -> Cell {
        let id = OBJECT_COUNTER.fetch_add(1, Ordering::SeqCst);
        Cell{id, events: HashSet::new(), visible_in_output: true }
    }

    /// we want to be able to fully clone a cell when we need to but have the canonical clone
    /// protect the ID by making it unique.
    pub fn pure_clone(&self) -> Cell {
        Cell{id: self.id.clone(), events: self.events.clone(), visible_in_output: true}
    }

    pub fn increment_id_clone(&self) -> Self {
        let id = OBJECT_COUNTER.fetch_add(1, Ordering::SeqCst);
        Cell{id, events: self.events.clone(), visible_in_output: true}
    }
}

#[cfg(test)]
mod tests {
    use crate::lineagemodels::model::SimpleDivision;
    use super::*;

    #[test]
    fn fetch_add_test() {
        let cell = Cell::new();
        let cell2 = Cell::new();
        assert_eq!(cell.id + 1, cell2.id);
    }
}