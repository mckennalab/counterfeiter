use crate::cell::Cell;
use crate::genome::GenomeEventCollection;
use rand_distr::{Poisson, Distribution};
use rand_distr::num_traits::ToPrimitive;

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventPosition = u16;

// our reserved mutational 'space'
pub type EventOutcomeIndex = u16;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DroppedAllele {
    Sampled,
    Dropped
}

/// Trait for implementing cellular editing systems and mutation models.
///
/// This trait defines the interface for different editing technologies
/// (e.g., CRISPR-Cas systems, base editors) that can modify cells during
/// simulation. Implementors handle the stochastic application of editing
/// events and conversion to analysis formats.
pub trait CellFactory {
    /// Applies editing events to a cell and returns the modified cell(s).
    ///
    /// This is the core function that simulates the editing process on a
    /// cell. It may introduce new mutations, modifications, or other editing
    /// events based on the specific editing system's characteristics. Implementations
    /// should generally worry about modifying their event outcomes, or providing
    /// mechanisms for dividing the cells (say biforcating), but not both.
    ///
    /// # Arguments
    ///
    /// * `input_cell` - The cell to apply editing events to
    /// * `genome` - Global collection of genome events for deduplication
    ///
    /// # Returns
    ///
    /// Vector of cells (typically just the modified input cell) after
    /// applying editing events.
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell>;

    /// Converts a cell's editing state to MIX format array representation.
    ///
    /// Generates a binary representation of the cell's editing state suitable
    /// for phylogenetic analysis. Includes handling of barcode dropout based
    /// on the specified drop rate.
    ///
    /// # Arguments
    ///
    /// * `genome` - Global genome event collection for lookups
    /// * `drop_rate` - Probability of barcode dropout (0.0 to 1.0)
    /// * `input_cell` - Cell to convert to MIX format
    ///
    /// # Returns
    ///
    /// * `Some(Vec<u8>)` - Binary array representing editing state
    /// * `None` - If cell should be excluded (e.g., all barcodes dropped)
    fn to_mix_array(&self, genome: &GenomeEventCollection, input_cell: &mut Cell, force_retention: &bool) -> (DroppedAllele,Option<Vec<u8>>);

    /// Maps an outcome index to its string representation.
    ///
    /// Converts internal numeric outcome indices to human-readable strings
    /// for output and debugging purposes.
    ///
    /// # Arguments
    ///
    /// * `index` - The outcome index to convert
    ///
    /// # Returns
    ///
    /// String representation of the outcome (e.g., "0", "1", "deletion").
    fn get_mapping(&self, index: &EventOutcomeIndex) -> String;

}

/// Trait for implementing cell division models.
///
/// This trait defines how cells divide during the simulation, allowing
/// for different division behaviors such as deterministic binary division,
/// stochastic division with variable offspring counts, or more complex
/// division rules based on cell state or environmental conditions.
pub trait DivisionModel {
    /// Divides a cell into offspring cells.
    ///
    /// Takes a parent cell and generates one or more offspring cells
    /// according to the specific division model. Each offspring inherits
    /// the parent's editing events and receives a new unique identifier.
    ///
    /// # Arguments
    ///
    /// * `input_cell` - The parent cell to divide
    /// * `genome` - Global genome event collection (may be used by some models)
    ///
    /// # Returns
    ///
    /// Vector of offspring cells generated from the parent. The number and
    /// characteristics of offspring depend on the specific division model.
    ///
    /// # Examples
    ///
    /// For binary division (2 offspring):
    /// ```rust
    /// let offspring = division_model.divide(&mut parent_cell, &mut genome);
    /// assert_eq!(offspring.len(), 2);
    /// ```
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


/// A stochastic cell division model with Poisson-distributed offspring count.
///
/// This division model generates a variable number of offspring based on
/// a Poisson distribution. This can model more realistic biological scenarios
/// where cell division rates vary stochastically around a mean value.
pub struct StochasticDivision {
    pub sampler: Poisson<f64>,
}

impl StochasticDivision {
    /// Creates a new stochastic division model with specified mean offspring count.
    ///
    /// Initializes the Poisson distribution sampler with the given lambda parameter,
    /// which represents the mean number of offspring per division event.
    ///
    /// # Arguments
    ///
    /// * `lambda` - Mean number of offspring (must be positive)
    ///
    /// # Returns
    ///
    /// A new `StochasticDivision` instance with the configured Poisson sampler.
    ///
    /// # Panics
    ///
    /// Panics if lambda is not positive or finite.
    ///
    /// # Examples
    ///
    /// ```rust
    /// // Mean of 2.0 offspring per division
    /// let division_model = StochasticDivision::new(&2.0);
    /// ```
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

