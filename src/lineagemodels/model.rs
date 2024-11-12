use crate::cell::Cell;
use crate::genome::GenomeEventCollection;
use rand_distr::{Poisson, Distribution};
use rand_distr::num_traits::ToPrimitive;

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventPosition = u16;

// our reserved mutational 'space'
pub type EventOutcomeIndex = u16;

pub trait CellFactory {
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell>;

    fn to_mix_array(&mut self, genome: &GenomeEventCollection, drop_rate: &f64, input_cell: &mut Cell) -> Option<Vec<u8>>;

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String;

    fn from_file(file: &String) -> Self;
}

pub trait DivisionModel {
    fn divide(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell>;
}

pub struct SimpleDivision {
    pub offspring_count: usize,
}


/// A simple cell division model where each cell divides into a fixed number of offspring cells.
///
/// The `SimpleDivision` struct represents a basic lineage model where each input cell divides
/// into a specified number of offspring cells. This model is deterministic and does not consider
/// stochastic events or additional cell states.
impl DivisionModel for SimpleDivision {

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
    fn divide(&self, input_cell: &mut Cell, _genome: &mut GenomeEventCollection) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = Vec::new();
        for _i in 0..self.offspring_count {
            new_cells.push(input_cell.increment_id_clone(input_cell.parent_id));
        }
        new_cells
    }

}


pub struct StochasticDivision {
    pub sampler: Poisson<f64>,
}

impl StochasticDivision {
    pub fn new(lambda: &f64) -> StochasticDivision {
        let poi = Poisson::new(*lambda).unwrap();
        StochasticDivision{
            sampler: poi,
        }
    }
}
/// A simple cell division model where each cell divides into a fixed number of offspring cells.
///
/// The `SimpleDivision` struct represents a basic lineage model where each input cell divides
/// into a specified number of offspring cells. This model is deterministic and does not consider
/// stochastic events or additional cell states.
impl DivisionModel for StochasticDivision {

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
    fn divide(&self, input_cell: &mut Cell, _genome: &mut GenomeEventCollection) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = Vec::new();
        for _i in 0..self.sampler.sample(&mut rand::thread_rng()).round().to_i64().unwrap() {
            new_cells.push(input_cell.increment_id_clone(input_cell.parent_id));
        }
        new_cells
    }

}

