use std::collections::HashMap;
use rand::Rng;
use crate::cell::Cell;
use crate::genome::{Genome, GenomeEventLookup};
use weighted_rand::table::WalkerTable;
use weighted_rand::builder::WalkerTableBuilder;
use weighted_rand::builder::NewBuilder;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, LineageModel};

#[derive(Clone, Debug)]
pub struct CRISPRBitRate {
    none: f64,
    left: f64,
    right: f64,
    both: f64,
    bind_modified_target: f64, // how often a modified target is bound by CRISPR and re-edited
    weighted_draw: WalkerTable,
}

impl CRISPRBitRate {
    pub fn new(none: f32, left: f32, right: f32, both: f32, bind_modified_target: f64) -> CRISPRBitRate {
        let mut weighted_draw = WalkerTableBuilder::new(&vec![none, left, right, both]).build();
        CRISPRBitRate{none: (none as f64), left: (left as f64), right: (right as f64), both: (both as f64), bind_modified_target, weighted_draw}
    }

}

pub struct CRISPRBits {
    integration_count: usize,
    targets_per_integration: usize,
    rates: Vec<CRISPRBitRate>,
}

impl CRISPRBits {
    const NEITHER: EventOutcomeIndex = 0;
    const LEFT: EventOutcomeIndex = 1;
    const RIGHT: EventOutcomeIndex = 2;
    const BOTH: EventOutcomeIndex = 3;

    pub fn new(int_count: &usize, targets_per_integration: &usize, rates: Vec<CRISPRBitRate>) -> CRISPRBits {
        CRISPRBits{
            integration_count: *int_count,
            targets_per_integration: *targets_per_integration,
            rates,
        }
    }

    pub fn index_to_outcome(index: &usize) -> EventOutcomeIndex {
        match index {
            0 => CRISPRBits::NEITHER,
            1 => CRISPRBits::LEFT,
            2 => CRISPRBits::RIGHT,
            3 => CRISPRBits::BOTH,
            _ => panic!("Invalid index for CRISPRBits"),
        }
    }

    pub fn outcome_to_index(outcome: &EventOutcomeIndex) -> usize {
        match outcome {
            &CRISPRBits::NEITHER => 0,
            &CRISPRBits::LEFT => 1,
            &CRISPRBits::RIGHT => 2,
            &CRISPRBits::BOTH => 3,
            _ => panic!("Invalid outcome for CRISPRBits"),
        }
    }

    fn draw_new_event(&self, current_state: &EventOutcomeIndex, position: &usize) -> EventOutcomeIndex {
        match current_state {
            &CRISPRBits::NEITHER => {
                CRISPRBits::index_to_outcome(&self.rates[*position].weighted_draw.next())
            }
            &CRISPRBits::BOTH => {
                CRISPRBits::BOTH
            }
            &CRISPRBits::LEFT => {
                let nxt = CRISPRBits::index_to_outcome(&self.rates[*position].weighted_draw.next());
                let mut rng: f64 = rand::thread_rng().gen();
                match (nxt, rng) {
                    (CRISPRBits::BOTH,x) if x < self.rates[*position].bind_modified_target => CRISPRBits::BOTH,
                    _ => CRISPRBits::LEFT,
                }
            }
            &CRISPRBits::RIGHT => {
                let nxt = CRISPRBits::index_to_outcome(&self.rates[*position].weighted_draw.next());
                let mut rng: f64 = rand::thread_rng().gen();
                match (nxt, rng) {
                    (CRISPRBits::BOTH,x) if x < self.rates[*position].bind_modified_target => CRISPRBits::BOTH,
                    _ => CRISPRBits::RIGHT,
                }
            }
            _ => {
                panic!("Invalid CRISPRBits state");
            }
        }
    }
}

impl LineageModel for CRISPRBits {

    fn estimated_event_space(&self) -> usize {
        self.targets_per_integration * self.integration_count
    }

    fn transform(&self, input_cell: &mut Cell) -> Vec<Cell> {
        assert_eq!(self.targets_per_integration,self.rates.len());
        (0..self.integration_count).for_each(|i| {
            let mut existing_events = input_cell.events.entry(
                Genome::CRISPRBits(i)).or_insert(GenomeEventLookup::new());

            (0..self.targets_per_integration).for_each(|t| {
                existing_events.events.
                    insert(t as EventPosition,
                           self.draw_new_event(
                               &existing_events.events.get(&(t as EventPosition)).
                                   map_or(CRISPRBits::NEITHER, |x| *x), &t));
            });
        });
        vec![input_cell.clone()]
    }


    fn get_mapping(&self, x: &EventOutcomeIndex) -> String {
        String::new()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crispr_bits_transform_cell() {
        
        let equal_rates = CRISPRBitRate::new(0.98, 0.0099, 0.0099, 0.0002, 0.01);

        let cBits = CRISPRBits::new(&1, &1, vec![equal_rates.clone()]);

        let mut counts = HashMap::new();

        for i in 0..100 {
            let mut cell = Cell::new(&10);
            for j in 0..100 {
                let cell = cBits.transform(&mut cell);
            };
            let outcome = cell.events.get(&Genome::CRISPRBits(0)).unwrap().events.get(&0).unwrap().clone();
            *counts.entry(outcome).or_insert(0) += 1;

        };
        println!("{:?}",counts);
    }
}