use crate::cell::Cell;

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventPosition = u16;

// our reserved mutational 'space'
pub type EventOutcomeIndex = u16;

pub trait CellFactory {

    fn estimated_event_space(&self) -> usize;

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell>;

    fn to_mix_array(&mut self, drop_rate: &f64, input_cell: &mut Cell) -> Option<Vec<u8>>;

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String;
}
pub struct SimpleDivision {
    pub offspring_count: usize,
}

/// A simple cell division model where each cell divides into a fixed number of offspring cells.
///
/// The `SimpleDivision` struct represents a basic lineage model where each input cell divides
/// into a specified number of offspring cells. This model is deterministic and does not consider
/// stochastic events or additional cell states.
impl CellFactory for SimpleDivision {
    /// Estimates the event space for the lineage model.
    ///
    /// Returns `0` because the `SimpleDivision` model does not involve any stochastic events
    /// or varying outcomes; each cell division is deterministic with a fixed number of offspring.
    ///
    /// # Returns
    ///
    /// * `usize` - The estimated size of the event space (always `0` for this model).
    fn estimated_event_space(&self) -> usize {
        0
    }

    /// Divides an input cell into multiple offspring cells.
    ///
    /// Creates `offspring_count` new `Cell` instances, each cloned from the input cell with
    /// an incremented identifier. This simulates the division process where each cell gives
    /// rise to multiple identical offspring.
    ///
    /// # Arguments
    ///
    /// * `input_cell` - A reference to the parent `Cell` that will be divided.
    ///
    /// # Returns
    ///
    /// * `Vec<Cell>` - A vector containing the newly created offspring cells.
    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = Vec::new();
        for _i in 0..self.offspring_count {
            new_cells.push(input_cell.increment_id_clone());
        }
        new_cells
    }

    fn to_mix_array(&mut self, drop_rate: &f64, input_cell: &mut Cell) -> Option<Vec<u8>> {
        None
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        format!("SimpleDivision:{}", index)
    }
}

