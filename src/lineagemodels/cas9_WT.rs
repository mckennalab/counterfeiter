use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::atomic::{AtomicU16, Ordering};
use rand::distributions::Distribution;
use rand_distr::Poisson;
use crate::cell::Cell;
use crate::genome::{EditingOutcome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification};
use crate::lineagemodels::cas12a_abe::Cas12aABE;
use crate::lineagemodels::model::{CellFactory, EventOutcomeIndex};

#[derive(Clone, Debug)]
pub struct Cas9WT {
    edit_rate: Vec<f64>,
    edit_width: Vec<f64>,
    positions: Vec<(u32,u32)>,
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
                self.positions.iter().enumerate().for_each(|(index,(tstart,tstop))| {
                    if (x.start > *tstart && x.start < *tstop) ||
                        (x.stop > *tstart && x.stop < *tstop) {
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
        let active_targets = self.events_to_active_target_sites(&current_events,genome);
        self.edit_rate.iter().zip(active_targets).enumerate().for_each(|(barcode_index, (barcode_edit_rate,active))| {
            if active {
                // draw an event with the proportion provided
                // Draw a float between 0.0 and 1.0 (continuous uniform)
                let rnd: f64 = rand::random::<f64>();

                if rnd < *barcode_edit_rate {
                    event_drawn = true;
                    let edit_start_offset = Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                    let edit_stop_offset = Poisson::new(self.edit_width.get(barcode_index).unwrap().clone()).unwrap().sample(&mut rand::thread_rng()).round() as u32;
                    cut_start = cut_start.min(edit_start_offset);
                    cut_end = cut_end.min(edit_stop_offset);
                }
            }
        });
        
        if event_drawn {
            let outcome = EditingOutcome {
                start: cut_start,
                stop: cut_end + 1, // [x,y) intervals
                change: Modification::Deletion,
                nucleotides: vec![],
                internal_outcome_id: 1,
            };
            let evt = genome.add_event(&self.genome, outcome);
            match evt {
                None => {
                    panic!("unable to add event");
                }
                Some(x) => {
                    input_cell.events.insert(x);
                }
            }
        }
        vec![input_cell.pure_clone()]
    }

    fn to_mix_array(&self, genome: &GenomeEventCollection, input_cell: &mut Cell) -> Option<Vec<u8>> {
        todo!()
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        todo!()
    }

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
