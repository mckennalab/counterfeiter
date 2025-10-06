use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::atomic::{AtomicU16, Ordering};
use rand::distributions::Distribution;
use rand::Rng;
use rand_distr::Poisson;
use crate::cell::Cell;
use crate::genome::{next_genome_id, DecayType, EditingOutcome, Genome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification, Probability};
use crate::lineagemodels::cas12a_abe::Cas12aABE;
use crate::lineagemodels::model::{CellFactory, EventOutcomeIndex};

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


// A global atomic counter
static GLOBAL_ID: AtomicU16 = AtomicU16::new(1);

pub fn next_editing_id() -> u16 {
    GLOBAL_ID.fetch_add(1, Ordering::Relaxed)
}


impl Cas9WT {
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
                    interdependent_rate: 0.0,
                })
        }
        return_vec
    }

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
    // we make the new cell! Just a 1-to-1 conversion
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

    fn to_mix_array(&self, genome_lookup_object: &GenomeEventCollection, input_cell: &mut Cell) -> Option<Vec<u8>> {
        let existing_events = genome_lookup_object.genomes.get(&self.genome).unwrap();
        
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
        Some(ret)
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
