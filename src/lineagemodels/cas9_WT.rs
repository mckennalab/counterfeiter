use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::atomic::{AtomicU16, Ordering};
use rand::distributions::Distribution;
use rand::Rng;
use rand_distr::Poisson;
use crate::cell::Cell;
use crate::genome::{EditingOutcome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification};
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
    pub fn events_to_active_target_sites(&self, event_outcome_index: &Vec<GenomeEventKey>, genome: &mut GenomeEventCollection) -> Vec<bool> {
        let mut ret = vec![true; self.edit_rate.len()];
        event_outcome_index.iter().for_each(|x| {
            let mut hash: HashSet<GenomeEventKey> = HashSet::default();
            hash.insert(x.clone());
            genome.filter_events_and_get_outcomes(&self.genome,&hash).iter().for_each(|x| {
                self.cut_positions.iter().enumerate().for_each(|(index,cutsite)| {
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
        let mut insertions : Vec<(u32, String)> = Vec::new();
        
        let active_targets = self.events_to_active_target_sites(&current_events,genome);
        self.edit_rate.iter().zip(active_targets).enumerate().for_each(|(barcode_index, (barcode_edit_rate,active))| {
            if active {
                // draw an event with the proportion provided
                // Draw a float between 0.0 and 1.0 (continuous uniform)
                let rnd: f64 = rand::random::<f64>();

                let cut_position = self.cut_positions.get(barcode_index).unwrap();
                
                if rnd < *barcode_edit_rate {
                    let rnd2: f64 = rand::random::<f64>();
                    if rnd2 < self.insertion_rate { 
                        let inserted_bases = random_nucleotides(Poisson::new(3.0f64).unwrap().sample(&mut rand::thread_rng()).round() as usize);
                        insertions.push((*cut_position,inserted_bases));
                        event_drawn = true;
                    } else {
                        
                        let edit_start_offset = cut_position - Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                        let edit_stop_offset = cut_position + Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                        cut_start = cut_start.min(edit_start_offset);
                        cut_end = cut_end.min(edit_stop_offset);
                        event_drawn = true;
                    }
                }
            }
        });
        
        if event_drawn {
            match (cut_end  == 0, insertions.len() == 0) {
                (true, true) => {panic!("we had at least one event")},
                (false, true) => {
                    let outcome = EditingOutcome {
                        start: cut_start,
                        stop: cut_end + 1, // [x,y) intervals
                        change: Modification::Deletion,
                        nucleotides: vec![],
                        internal_outcome_id: next_editing_id(),
                    };
                    let ky = genome.add_event(&self.genome, outcome).unwrap();
                    input_cell.events.insert(ky);
                    
                },
                (true,false) => {
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
                },
                (false, false) => {
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

    fn to_mix_array(&self, genome: &GenomeEventCollection, input_cell: &mut Cell) -> Option<Vec<u8>> {
        let mut existing_events = genome
            .filter_events_and_get_outcomes(&self.genome, &input_cell.events)
            .iter()
            .map(|x| (x.start, x.clone()))
            .collect::<Vec<(u32, EditingOutcome)>>();
        existing_events.sort_by(|a, b| b.0.cmp(&a.0));
        
        let ret = existing_events.iter().map(|(x,y)| {
            let mut has_event = false;
            input_cell.events.iter().for_each(|cell_key| {
                if cell_key.outcome_index == y.internal_outcome_id {
                    has_event = true;
                }
            });
            if has_event {
                b'1'
            } else {
                b'0'
            }
        }).collect::<Vec<u8>>();
        
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
    fn test_iterative_example() {
        
    }

    #[test]
    fn test_load_from_pileup_file() {
        
    }
}
