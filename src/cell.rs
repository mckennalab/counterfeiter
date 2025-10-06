use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use crate::genome::{Genome, GenomeDescription, GenomeEventKey};

// from https://stackoverflow.com/questions/32935808/generate-sequential-ids-for-each-instance-of-a-struct
// and https://doc.rust-lang.org/std/sync/atomic/
// we want to have a unique ID for each cell, so we can relate them
static OBJECT_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A cell is a collection of events at target sites. This struct is kept as minimal
/// as possible to reduce memory usage, and events and sites are referred to by index.
/// The cells are passed the relevant look-up object when they 'split' to create new cells.
#[derive(Debug, Clone)]
pub struct Cell {
    pub id: usize,
    pub parent_id: usize,
    pub events: HashSet<GenomeEventKey>,
    pub optional_genomes: Vec<usize>,
    pub visible_in_output: bool,
}

impl Cell {
    /// Creates a new cell with a unique identifier.
    ///
    /// Generates a globally unique ID using an atomic counter and initializes
    /// the cell with no parent (parent_id = usize::MAX), no editing events,
    /// and default visibility in output.
    ///
    /// # Returns
    ///
    /// A new `Cell` instance with:
    /// - Unique ID from global atomic counter
    /// - No parent (parent_id = usize::MAX)
    /// - Empty event set
    /// - Visible in output by default
    ///
    /// # Examples
    ///
    /// ```rust
    /// let cell = Cell::new();
    /// assert!(cell.events.is_empty());
    /// assert_eq!(cell.parent_id, usize::MAX);
    /// assert!(cell.visible_in_output);
    /// ```
    pub fn new() -> Cell {
        let id = OBJECT_COUNTER.fetch_add(1, Ordering::SeqCst);
        Cell{id, parent_id: usize::MAX, events: HashSet::new(), optional_genomes: vec![], visible_in_output: true }
    }

    /// Creates an exact copy of the cell preserving the same ID.
    ///
    /// This method performs a deep clone of all cell data while maintaining
    /// the original cell's unique identifier. Used when you need multiple
    /// references to the same logical cell or for checkpoint/snapshot purposes.
    ///
    /// # Returns
    ///
    /// A `Cell` with identical data to the original, including:
    /// - Same unique ID
    /// - Same parent ID
    /// - Clone of the events set
    /// - Reset visibility to true
    ///
    /// # Examples
    ///
    /// ```rust
    /// let original = Cell::new();
    /// let clone = original.pure_clone();
    /// assert_eq!(original.id, clone.id);
    /// assert_eq!(original.parent_id, clone.parent_id);
    /// ```
    pub fn pure_clone(&self) -> Cell {
        Cell{id: self.id.clone(), 
            parent_id: self.parent_id, 
            events: self.events.clone(), 
            optional_genomes: vec![], 
            visible_in_output: true}
    }

    /// Creates a copy of the cell with a new unique ID and specified parent.
    ///
    /// This method is used during cell division to create offspring cells.
    /// Each offspring gets a new unique identifier while inheriting the
    /// editing events from the parent cell.
    ///
    /// # Arguments
    ///
    /// * `parent_id` - The ID of the parent cell (typically self.id)
    ///
    /// # Returns
    ///
    /// A new `Cell` with:
    /// - New unique ID from global atomic counter
    /// - Specified parent ID
    /// - Clone of parent's events
    /// - Default visibility in output
    ///
    /// # Examples
    ///
    /// ```rust
    /// let parent = Cell::new();
    /// let child = parent.increment_id_clone(parent.id);
    /// assert_ne!(parent.id, child.id);
    /// assert_eq!(child.parent_id, parent.id);
    /// ```
    pub fn increment_id_clone(&self, parent_id: usize) -> Self {
        let id = OBJECT_COUNTER.fetch_add(1, Ordering::SeqCst);
        Cell{id, parent_id, events: self.events.clone(), optional_genomes: vec![], visible_in_output: true}
    }
    
    
    pub fn filter_to_genome_events(&self, genome_id: &u16) -> Vec<GenomeEventKey> {
        let mut ret = vec![];
        self.events.iter().for_each(|x| {
            if x.genome_index == *genome_id {
                ret.push(x.clone());
            } 
        });
        ret
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