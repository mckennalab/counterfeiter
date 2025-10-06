use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io;
use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, IntervalTree};

use std::sync::atomic::{AtomicU16, Ordering};
use rand_distr::num_traits::ToPrimitive;
use regex::internal::Input;
use crate::cell::Cell;
use crate::lineagemodels::model::CellFactory;
use io::Write;

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub enum Genome {
    CRISPRBits,
    ABECas12a,
    ECDNA,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub enum DecayType {
    Permanent,
    StochasticDrift(u64), // x / 100,000 == prop of retention; capped at 1,000,000,000 (1B) by asserts. Using u64 fo hash without headache
}

// A global atomic counter
static GLOBAL_ID: AtomicU16 = AtomicU16::new(1);

pub fn next_genome_id() -> u16 {
    GLOBAL_ID.fetch_add(1, Ordering::Relaxed)
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Probability(u64);

impl Probability {
    pub const MIN: u64 = 0;
    pub const MAX: u64 = 999999;

    pub fn new(value: u64) -> Option<Self> {
        if value >= Self::MIN && value <= Self::MAX {
            Some(Probability(value))
        } else {
            None
        }
    }

    pub fn new_from_f64(value: &f64) -> Option<Self> {
        if value >= &0.0 && value <= &1.0 {
            Some(Probability((*value * Probability::MAX.to_f64().unwrap()).round() as u64))
        } else {
            None
        }
    }

    pub fn get(self) -> f64 {
        self.0.to_f64().unwrap() / Probability::MAX.to_f64().unwrap()
    }
}


#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomeID(u16);

impl GenomeID {
    pub fn new() -> u16 { next_genome_id() }
    pub fn get(self) -> u16 { self.0 }
}


#[derive(Debug, PartialEq, Eq, Hash,  Clone)]
pub struct GenomeDescription {
    pub genome: Genome,
    pub name: String,
    pub allows_overlap: bool,
    pub decay_type: DecayType,
    pub drop_rate: Probability, 
    pub id: u16,
}

impl PartialOrd for GenomeDescription {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.id.partial_cmp(&other.id)
    }
}
impl Ord for GenomeDescription {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.cmp(&other.id)
    }
}

/// Compactly store the outcome of mutating the genome. This stores indexes into
/// the global registry of events; the goal is simply to store an offset into a much richer
/// database. Ordered to avoid splitting the u16s and getting weird alignment in the worst case
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct GenomeEventKey {
    pub position_index: u32,
    pub genome_index: u16,
    pub outcome_index: u16,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub enum Modification {
    Substitution,
    Insertion,
    Deletion,
    Inversion
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct EditingOutcome {
    pub start: u32,
    pub stop: u32,
    pub change: Modification,
    pub nucleotides: Vec<u8>,
    pub internal_outcome_id: u16,
}


/// Generates MIX format output file from cell population.
///
/// Converts the editing states of all cells into MIX format suitable for
/// phylogenetic analysis. Applies barcode dropout and filters cells based
/// on the cell_ids_to_keep mapping. Creates a tab-separated file with
/// binary character states.
///
/// # Arguments
///
/// * `genome` - Global genome event collection for event lookups
/// * `cells` - Vector of cells to convert to MIX format
/// * `drop_rate` - Probability of barcode dropout (0.0 to 1.0)
/// * `cell_ids_to_keep` - Mutable map of cell IDs to include in output
/// * `output` - Path to output MIX file
///
/// # Output Format
/// - Header line: number of cells and total characters
/// - Data lines: cell name followed by binary string (0=unedited, 1=edited, ?=missing)
///
/// # Side Effects
/// - Creates/overwrites the output file
/// - Modifies cell_ids_to_keep by removing cells with all barcodes dropped
/// - May set cell.visible_in_output to false for excluded cells
pub fn create_mix_file(
    genomes: &GenomeEventCollection,
    ordered_editors: &Vec<Box<dyn CellFactory>>,
    cells: &mut Vec<Cell>,
    cell_ids_to_keep: &mut HashMap<usize, bool>,
    output_file: &String,
) {
    let mut cell_to_output: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut target_count = 0; 
    
    cells.iter_mut().enumerate().for_each(|(index,cell)| {
        if cell_ids_to_keep.contains_key(&cell.id) {
            ordered_editors.iter().for_each(|(genome)| {
                match genome.to_mix_array(genomes, cell) {
                    Some(x) => {
                        println!("x {:?} {}",x, cell_to_output.contains_key(&cell.id));
                        cell_to_output.entry(cell.id).and_modify(|v| v.extend(x)).or_insert(vec![]);
                        println!("TTTT {}",cell_to_output.get(&cell.id).unwrap().len());
                        if target_count < cell_to_output.get(&cell.id).unwrap().len() {
                            target_count = cell_to_output.get(&cell.id).unwrap().len();
                        }
                    }
                    None => {
                        error!("Unable to process cell {}",cell.id);
                        cell_ids_to_keep.remove(&cell.id);
                    }
                }
            });
        }
    });

    let mut out = File::create(output_file).unwrap();
    write!(
        out,
        "\t{}\t{}\n",
        cell_to_output.len(),
        target_count
    )
        .unwrap();
    cell_to_output.iter().for_each(|(k, v)| {
        println!("CCCCCCCCCCCCCCC len {}",v.len());
        write!(
            out,
            "{:<10}\t{}\n",
            format!("n{}", k),
            String::from_utf8(v.clone()).unwrap()
        );
    });
}


#[derive(Clone)]
pub struct GenomeEventCollection {
    pub genomes: BTreeMap<GenomeDescription,HashMap<u32,EditingOutcome>>,
    genomes_offsets: HashMap<GenomeDescription,u16>,
    key_to_outcome: HashMap<GenomeEventKey,EditingOutcome>,
    outcome_to_key: HashMap<EditingOutcome,GenomeEventKey>,
}

impl GenomeEventCollection {
    /// Creates a new empty genome event collection.
    ///
    /// Initializes all internal data structures for storing genome events,
    /// their outcomes, and the mappings between them. This serves as the
    /// global registry for all editing events that occur during simulation.
    ///
    /// # Returns
    ///
    /// A new `GenomeEventCollection` with empty collections for:
    /// - Genome descriptions and their associated events
    /// - Genome type offsets for efficient indexing
    /// - Bidirectional mappings between event keys and outcomes
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut genome_events = GenomeEventCollection::new();
    /// assert_eq!(genome_events.current_genome_offset, 0);
    /// ```
    pub fn new() -> GenomeEventCollection {
        GenomeEventCollection{
            genomes: BTreeMap::default(),
            genomes_offsets: HashMap::new(),
            key_to_outcome: HashMap::new(),
            outcome_to_key: HashMap::new() }
    }

    /// Creates a compact event key for a given genome and editing outcome.
    ///
    /// Generates a `GenomeEventKey` that serves as a compact reference to
    /// a specific editing event. The key combines position, genome type,
    /// and outcome information into a space-efficient identifier.
    ///
    /// # Arguments
    ///
    /// * `genome` - The genome description for the editing system
    /// * `outcome` - The editing outcome to create a key for
    ///
    /// # Returns
    ///
    /// A `GenomeEventKey` containing:
    /// - Position index from the outcome's start position
    /// - Genome index from the genome type offset
    /// - Outcome index from the outcome's internal ID
    ///
    /// # Panics
    ///
    /// Panics if the genome has not been registered in `genomes_offsets`.
    fn create_event(&mut self, genome: &GenomeDescription, outcome: &EditingOutcome) -> GenomeEventKey {

        GenomeEventKey{
            position_index: outcome.start.clone(),
            genome_index: *self.genomes_offsets.get(genome).unwrap(),
            outcome_index: outcome.internal_outcome_id.clone(),
        }
    }

    /// Filters events for a specific genome and returns their outcomes.
    ///
    /// Takes a set of event keys and returns only the editing outcomes
    /// that belong to the specified genome type. This is used to extract
    /// genome-specific events from a cell's complete event set.
    ///
    /// # Arguments
    ///
    /// * `genome` - The genome description to filter by
    /// * `set` - Set of event keys to filter
    ///
    /// # Returns
    ///
    /// Vector of `EditingOutcome`s that match the specified genome type.
    /// Returns empty vector if the genome type is not registered.
    ///
    /// # Panics
    ///
    /// Panics if an event key in the set is not found in the collection.
    pub fn filter_for_genome_and_matching_events(&self, genome: &GenomeDescription, set: &HashSet<GenomeEventKey>) -> Vec<EditingOutcome> {
        if !self.genomes_offsets.contains_key(genome) {
            warn!("Unable to find genome {} {:?}",genome.id, self.genomes_offsets.len());
            Vec::new()
        } else {
            let genome_offset = self.genomes_offsets.get(genome).unwrap();
            set.iter().map(|it| {
                match self.key_to_outcome.get(it) {
                    None => {
                        panic!("We dont know about event {:?}", it);
                    }
                    Some(x) => {
                        match *genome_offset == it.genome_index {
                            true => {
                                Some(x.clone())
                            }
                            false => {
                                None
                            }
                        }
                    }
                }
            }).flatten().collect::<Vec<EditingOutcome>>()
        }
    }


    /// Adds a new editing event to the collection.
    ///
    /// Registers a new editing outcome for a specific genome type. If this is
    /// the first event for the genome type, initializes the genome's data
    /// structures. Prevents duplicate events at the same position.
    ///
    /// # Arguments
    ///
    /// * `genome` - The genome description for the editing system
    /// * `outcome` - The editing outcome to add
    ///
    /// # Returns
    ///
    /// * `Some(GenomeEventKey)` - Key for the added or existing event
    /// * `None` - If the event could not be added (should not occur in current implementation)
    ///
    /// # Behavior
    ///
    /// - For new genome types: Creates new data structures and assigns offset
    /// - For existing genomes: Checks for position conflicts
    /// - Updates bidirectional mappings between keys and outcomes
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut collection = GenomeEventCollection::new();
    /// let genome_desc = GenomeDescription { /* ... */ };
    /// let outcome = EditingOutcome { /* ... */ };
    /// let key = collection.add_event(&genome_desc, outcome);
    /// ```
    pub fn add_event(&mut self, genome: &GenomeDescription, outcome: EditingOutcome) -> Option<GenomeEventKey> {

        match self.genomes.contains_key(genome) {
            false => {
                // we dont have a record of this genome -- set it up and add the new outcome
                let mut lt = HashMap::new();
                lt.insert(outcome.start, outcome.clone());
                self.genomes.insert(genome.clone(), lt);
                println!("adding to find key {:?} {:?}",genome,outcome);

                self.genomes_offsets.insert(genome.clone(), genome.id);
                let id = self.create_event(genome,&outcome);

                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());
                Some(id)
            }
            true => {
                // we have a record -- check (1) for overlap and (2) if we allow that.
                let id = self.create_event(genome,&outcome);
                println!("adding to prev key {:?} {:?} {:?}",id,genome,outcome);
                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());

                match self.genomes.get(genome).unwrap().contains_key(&outcome.start) {
                    false => {
                        self.genomes.get_mut(genome).unwrap().insert(outcome.start,outcome.clone());
                        let id = self.create_event(genome,&outcome);

                        self.key_to_outcome.insert(id.clone(), outcome.clone());
                        self.outcome_to_key.insert(outcome.clone(), id.clone());
                        Some(id)
                    },
                    true => {
                        Some(self.outcome_to_key.get(&outcome).unwrap().clone())
                    }
                }
            }
        }
    }
    
}


