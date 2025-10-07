/// Cas9 wild-type (WT) lineage tracing system implementation.
///
/// This module implements a Cas9-based genome editing system that simulates
/// wild-type Cas9 behavior for lineage tracing applications. Cas9 creates
/// double-strand breaks that result in deletions and insertions through
/// non-homologous end joining (NHEJ) repair.

use crate::rand_distr::Distribution;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::atomic::{AtomicU16, Ordering};
use rand::{rng, Rng};
use rand_distr::Poisson;
use crate::cell::Cell;
use crate::genome::{next_genome_id, DecayType, EditingOutcome, Genome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification, Probability};
use crate::lineagemodels::cas12a_abe::Cas12aABE;
use crate::lineagemodels::model::{CellFactory, DroppedAllele, EventOutcomeIndex};


/// Cas9 wild-type editing system for lineage tracing.
///
/// Models Cas9-induced double-strand breaks and NHEJ repair outcomes including
/// deletions and insertions. Each Cas9WT instance represents a single integrated
/// barcode with multiple target sites.
///
/// # Fields
///
/// * `edit_rate` - Per-target editing probability for each target site
/// * `insertion_rate` - Probability of insertion vs deletion when editing occurs
/// * `edit_width` - Poisson mean for deletion size distribution at each target
/// * `cut_positions` - Genomic positions where Cas9 creates double-strand breaks
/// * `size` - Total barcode length in base pairs
/// * `targets_per_barcode` - Number of target sites per barcode
/// * `description` - Human-readable identifier for this editing system
/// * `genome` - Genome description containing decay and dropout properties
/// * `interdependent_rate` - Rate of correlated editing between sites (currently unused)
#[derive(Clone, Debug)]
pub struct Cas9WT {
    edit_rate: Vec<f64>,
    insertion_rate: f64,
    edit_width: Vec<f64>,
    cut_positions: Vec<u32>,
    size: u32,
    pub targets_per_barcode: u32,
    description: String,
    genome: GenomeDescription,
    interdependent_rate: f64,
}


/// Global atomic counter for unique editing event IDs.
///
/// Used to assign unique identifiers to each editing outcome across all
/// Cas9 instances and simulation runs.
static GLOBAL_ID: AtomicU16 = AtomicU16::new(1);

/// Generates the next unique editing event ID.
///
/// Thread-safe counter using atomic operations to ensure unique IDs
/// even in concurrent contexts.
///
/// # Returns
///
/// A unique 16-bit identifier for an editing event
pub fn next_editing_id() -> u16 {
    GLOBAL_ID.fetch_add(1, Ordering::Relaxed)
}


impl Cas9WT {
    /// Creates multiple Cas9WT instances with uniform editing rates.
    ///
    /// Constructs a vector of Cas9 editing systems, one for each integrated barcode.
    /// Each barcode has identical editing properties and multiple target sites.
    ///
    /// # Barcode Structure
    ///
    /// Each barcode consists of:
    /// - 20bp flanking region (start)
    /// - N target sites, each 24bp apart
    /// - 20bp flanking region (end)
    /// - Total length: 40 + (target_count × 24) bp
    /// - Cut position per target: flank + (target_index × 24) + 16
    ///
    /// # Arguments
    ///
    /// * `rate` - Base editing probability per target site per generation
    /// * `target_count` - Number of target sites per barcode
    /// * `integration_count` - Number of independent barcodes to create
    /// * `insertion_rate` - Probability of insertion vs deletion upon editing
    /// * `description` - Name/identifier for this editing system
    /// * `drop_rate` - Probability of complete barcode dropout per generation
    /// * `interdependent_rate` - Correlation between target sites (currently unused)
    ///
    /// # Returns
    ///
    /// Vector of `integration_count` independent Cas9WT instances, each with
    /// unique genome IDs but identical editing parameters
    pub fn from_editing_rate(
        rate: &f64,
        target_count: &u32,
        integration_count: &u32,
        insertion_rate: f64,
        description: String,
        drop_rate: &f64,
        interdependent_rate: &f64,
    ) -> Vec<Cas9WT> {
        // our length is going to be the target: (count * 24) + 20 + 20
        let size = 40 + (*target_count * 24);
        let cut_positions = (0..*target_count).map(|x| 20 + (x * 24) + 16).collect::<Vec<u32>>();
        let edit_rate = (0..*target_count).map(|x| rate.clone()).collect::<Vec<f64>>();
        let edit_width = (0..*target_count).map(|x| 8.0f64).collect::<Vec<f64>>();

        let mut return_vec = Vec::new();
        for i in (0..*integration_count) {
            return_vec.push(
                Cas9WT {
                    edit_rate: edit_rate.clone(),
                    insertion_rate: insertion_rate.clone(),
                    edit_width: edit_width.clone(),
                    cut_positions: cut_positions.clone(),
                    size: size,
                    targets_per_barcode: *target_count,
                    description: description.clone(),
                    genome: GenomeDescription {
                        genome: Genome::ABECas12a,
                        name: description.clone(),
                        allows_overlap: true,
                        decay_type: DecayType::Permanent,
                        drop_rate: Probability::new_from_f64(drop_rate).unwrap(),
                        id: next_genome_id(),
                    },
                    interdependent_rate: *interdependent_rate,
                })
        }
        return_vec
    }

    /// Determines which target sites remain active after previous editing events.
    ///
    /// Target sites become inactive if prior editing events (deletions/insertions)
    /// overlap with the target's effective editing window. This models the biological
    /// reality that Cas9 cannot bind to disrupted target sequences.
    ///
    /// # Inactivation Window
    ///
    /// A target at position `cutsite` is inactivated if any existing event overlaps
    /// the region `[cutsite - 16, cutsite + 6]`, representing the ~22bp target sequence.
    ///
    /// # Arguments
    ///
    /// * `event_outcome_index` - Cell's current editing events to check
    /// * `genome` - Global genome collection for looking up event details
    ///
    /// # Returns
    ///
    /// Boolean vector of length `self.edit_rate.len()` where `true` indicates
    /// the target site is still active for editing
    pub fn events_to_active_target_sites(&self, event_outcome_index: &Vec<GenomeEventKey>, genome: &mut GenomeEventCollection) -> Vec<bool> {
        let mut ret = vec![true; self.edit_rate.len()];
        event_outcome_index.iter().for_each(|x| {
            let mut hash: HashSet<GenomeEventKey> = HashSet::default();
            hash.insert(x.clone());
            genome.filter_for_genome_and_matching_events(&self.genome, &hash).iter().for_each(|x| {
                self.cut_positions.iter().enumerate().for_each(|(index, cutsite)| {
                    if (x.start > *cutsite - 16 && x.start < *cutsite + 6) ||
                        (x.stop > *cutsite - 16 && x.stop < *cutsite + 6) {
                        ret[index] = false;
                    }
                });
            });
        });
        ret
    }
}

impl CellFactory for Cas9WT {
    /// Applies Cas9 editing to a cell, generating deletions and/or insertions.
    ///
    /// Simulates one generation of Cas9-mediated editing by:
    /// 1. Checking which target sites are still active
    /// 2. Drawing editing events stochastically for active targets
    /// 3. Determining insertion vs deletion outcomes
    /// 4. Creating merged deletion spans when multiple targets are edited
    /// 5. Adding valid insertions outside deletion regions
    ///
    /// # Editing Logic
    ///
    /// - For each active target, draw editing with probability `edit_rate[i]`
    /// - If editing occurs, draw insertion vs deletion with probability `insertion_rate`
    /// - Deletions: Poisson-distributed size centered at cut position
    /// - Insertions: Random nucleotide sequence, Poisson(mean=3) length
    /// - Multiple deletions are merged into a single deletion spanning all cuts
    /// - Insertions overlapping deletions are discarded
    ///
    /// # Arguments
    ///
    /// * `input_cell` - Cell to edit (modified in place with new events)
    /// * `genome` - Global genome collection for tracking all editing outcomes
    ///
    /// # Returns
    ///
    /// Single-element vector containing the edited cell (1-to-1 conversion)
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell> {
        let current_events = input_cell.filter_to_genome_events(&self.genome.id);
        
        let mut cut_start = self.size;
        let mut cut_end = 0;
        let mut event_drawn = false;
        let mut insertion_drawn = false;
        let mut insertions: Vec<(u32, String)> = Vec::new();

        let active_targets = self.events_to_active_target_sites(&current_events, genome);
        self.edit_rate.iter().zip(active_targets).enumerate().for_each(|(barcode_index, (barcode_edit_rate, active))| {
            if active {
                // draw an event with the proportion provided
                // Draw a float between 0.0 and 1.0 (continuous uniform)
                let rnd: f64 = rand::random::<f64>();

                let cut_position = self.cut_positions.get(barcode_index).unwrap();

                if rnd < *barcode_edit_rate {
                    let rnd2: f64 = rand::random::<f64>();
                    if rnd2 < self.insertion_rate {
                        let inserted_bases = random_nucleotides(Poisson::new(3.0f64).unwrap().sample(&mut rand::thread_rng()).round() as usize);
                        insertions.push((*cut_position, inserted_bases));
                        insertion_drawn = true;
                    } else {
                        let edit_start_offset = cut_position - Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                        let edit_stop_offset = cut_position + Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                        cut_start = cut_start.min(edit_start_offset);
                        cut_end = cut_end.max(edit_stop_offset);
                        event_drawn = true;
                    }
                }
            }
        });

        if event_drawn || insertion_drawn {
            println!("mutate {} {} genome {}",event_drawn,insertion_drawn,self.genome.id);
            match (event_drawn, insertion_drawn) {
                (false, false) => { panic!("we had at least one event") }
                (true, false) => {
                    let outcome = EditingOutcome {
                        start: cut_start,
                        stop: cut_end + 1, // [x,y) intervals
                        change: Modification::Deletion,
                        nucleotides: vec![],
                        internal_outcome_id: next_editing_id(),
                    };
                    let ky = genome.add_event(&self.genome, outcome).unwrap();
                    input_cell.events.insert(ky);
                }
                (false, true) => {
                    // just an insertion(s)
                    for insertion in insertions {
                        let outcome = EditingOutcome {
                            start: insertion.0,
                            stop: insertion.0 + 1, // [x,y) intervals
                            change: Modification::Insertion,
                            nucleotides: insertion.1.into_bytes(),
                            internal_outcome_id: next_editing_id(),
                        };
                        let ky = genome.add_event(&self.genome, outcome).unwrap();
                        input_cell.events.insert(ky);
                    }
                }
                (true, true) => {
                    // we have to do the deletion
                    let outcome = EditingOutcome {
                        start: cut_start,
                        stop: cut_end + 1, // [x,y) intervals
                        change: Modification::Deletion,
                        nucleotides: vec![],
                        internal_outcome_id: next_editing_id(),
                    };
                    let ky = genome.add_event(&self.genome, outcome).unwrap();
                    input_cell.events.insert(ky);

                    // now for each insertion, check that if's outside any deletions
                    for insertion in insertions {
                        if insertion.0 < cut_start || insertion.0 > cut_end {
                            let outcome = EditingOutcome {
                                start: insertion.0,
                                stop: insertion.0 + 1, // [x,y) intervals
                                change: Modification::Insertion,
                                nucleotides: insertion.1.into_bytes(),
                                internal_outcome_id: next_editing_id(),
                            };
                            let ky = genome.add_event(&self.genome, outcome).unwrap();
                            input_cell.events.insert(ky);
                        }
                    }
                }
            }
        }
        vec![input_cell.pure_clone()]
    }

    fn to_mix_array(&self, genome_lookup_object: &GenomeEventCollection, input_cell: &mut Cell, force_retention: &bool) -> (DroppedAllele,Option<Vec<u8>>) {
        let existing_events: &HashMap<u32,EditingOutcome> = match genome_lookup_object.genomes.get(&self.genome) {
            None => {
                error!("genome does not exist {:?}, likely because we haven't edited that genome",&self.genome);
                &HashMap::default()
            }
            Some(x) => {x}
        };
        println!("drop rate {}",self.genome.drop_rate.get());
        let drawed = rand::rng().random::<f64>();
        if drawed < self.genome.drop_rate.get() && !force_retention {
            println!("dropping {}",drawed);
            (DroppedAllele::Dropped,Some(existing_events.iter().map(|(x, y)| b'?').collect()))
        } else {
            println!("existing len {}", existing_events.len());
            let mut has_it = 0;
            let ret = existing_events.iter().map(|(x, y)| {
                let mut has_event = false;
                input_cell.events.iter().for_each(|cell_key| {
                    if cell_key.outcome_index == y.internal_outcome_id {
                        has_event = true;
                    }
                });
                if has_event {
                    has_it += 1;
                    b'1'
                } else {
                    b'0'
                }
            }).collect::<Vec<u8>>();
            println!("ret {:?} has it {}", ret,has_it);
            (DroppedAllele::Sampled,Some(ret))
        }
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        todo!()
    }
}


fn random_nucleotides(length: usize) -> String {
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut rng = rand::thread_rng();

    (0..length)
        .map(|_| nucleotides[rng.gen_range(0..4)])
        .collect()
}

#[cfg(test)]
mod tests {
    use std::io::BufWriter;
    use super::*;

    #[test]
    fn test_iterative_example() {}

    #[test]
    fn test_load_from_pileup_file() {}
}
