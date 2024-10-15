use std::collections::HashMap;
use std::fs::File;
use rand::Rng;
use crate::cell::Cell;
use crate::genome::{Genome, GenomeEventLookup};
use weighted_rand::table::WalkerTable;
use weighted_rand::builder::WalkerTableBuilder;
use weighted_rand::builder::NewBuilder;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, LineageModel};
use std::io::Write;

#[derive(Clone, Debug)]
pub struct CRISPRBitRate {
    none: f32,
    left: f32,
    right: f32,
    both: f32,
    bind_modified_target: f32,
    // how often a modified target is bound by CRISPR and re-edited
    weighted_draw: WalkerTable,
}

impl CRISPRBitRate {
    pub fn new(none: f32, left: f32, right: f32, both: f32, bind_modified_target: f32) -> CRISPRBitRate {
        let weighted_draw = WalkerTableBuilder::new(&vec![none, left, right, both]).build();
        CRISPRBitRate { none, left, right, both, bind_modified_target, weighted_draw }
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
        CRISPRBits {
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
                let rng: f32 = rand::thread_rng().gen();
                match (nxt, rng) {
                    (CRISPRBits::BOTH, x) if x < self.rates[*position].bind_modified_target => CRISPRBits::BOTH,
                    _ => CRISPRBits::LEFT,
                }
            }
            &CRISPRBits::RIGHT => {
                let nxt = CRISPRBits::index_to_outcome(&self.rates[*position].weighted_draw.next());
                let rng: f32 = rand::thread_rng().gen();
                match (nxt, rng) {
                    (CRISPRBits::BOTH, x) if x < self.rates[*position].bind_modified_target => CRISPRBits::BOTH,
                    _ => CRISPRBits::RIGHT,
                }
            }
            _ => {
                panic!("Invalid CRISPRBits state");
            }
        }
    }


    pub fn to_mix_input(cells: &Vec<Cell>, cBits: &CRISPRBits, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        for cell in cells {
            cell_to_output.insert(cell.id, cBits.to_mix_array(&cell));
            if event_len.is_none() {
                event_len = Some(cell_to_output.get(&cell.id).unwrap().len());
            } else {
                assert_eq!(cell_to_output.get(&cell.id).unwrap().len(), event_len.unwrap());
            }
        }

        // filter out columns that are all one value
        let mut kept_columns = Vec::new();
        for i in 0..event_len.unwrap() {
            let mut counts = HashMap::new();
            let mut values = cell_to_output.iter().for_each(|(k, v)| *counts.entry(v[i]).or_insert(0) += 1);
            if counts.len() > 1 {
                kept_columns.push(i);
            }
        }
        let mut out = File::create(output).unwrap();
        write!(out, "\t{}\t{}\n", cell_to_output.len(), kept_columns.len()).unwrap();
        cell_to_output.iter().for_each(|(k, v)| {
            write!(out, "{:<10}\t{}\n", format!("n{}", k), kept_columns.iter().map(|x| if v[*x] == 0 { "0".to_string() } else { "1".to_string() }).collect::<Vec<String>>().join("")).unwrap();
        });
    }

    pub fn to_mix_input_full(generations: HashMap<usize, Vec<Cell>>, cBits: &CRISPRBits, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        for (generation, cells) in generations {
            for cell in cells {
                cell_to_output.insert(cell.id, cBits.to_mix_array(&cell));
                if event_len.is_none() {
                    event_len = Some(cell_to_output.get(&cell.id).unwrap().len());
                } else {
                    assert_eq!(cell_to_output.get(&cell.id).unwrap().len(), event_len.unwrap());
                }
            }
        }
        // filter out columns that are all one value
        let mut kept_columns = Vec::new();
        for i in 0..event_len.unwrap() {
            let mut counts = HashMap::new();
            let mut values = cell_to_output.iter().for_each(|(k, v)| *counts.entry(v[i]).or_insert(0) += 1);
            if counts.len() > 1 {
                kept_columns.push(i);
            }
        }
        let mut out = File::create(output).unwrap();
        write!(out, "\t{}\t{}\n", cell_to_output.len(), kept_columns.len()).unwrap();
        cell_to_output.iter().for_each(|(k, v)| {
            write!(out, "{:<10}\t{}\n", format!("n{}", k), kept_columns.iter().map(|x| if v[*x] == 0 { "0".to_string() } else { "1".to_string() }).collect::<Vec<String>>().join("")).unwrap();
        });
    }

    pub fn to_newick_tree(parent_child_map: &HashMap<usize, Vec<usize>>, output: &String) {
        let mut out = File::create(output).unwrap();
        write!(out, "{};\n", CRISPRBits::recursive_tree_builder(parent_child_map, &0)).expect("Unable to write file");
    }

    pub fn recursive_tree_builder(parent_child_map: &HashMap<usize, Vec<usize>>, current_index: &usize) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                //format!("(({})n{})", children.iter().map(|x| recursive_tree_builder(parent_child_map, x)).collect::<Vec<String>>().join(","), current_index)
                format!("({})", children.iter().map(|x| CRISPRBits::recursive_tree_builder(parent_child_map, x)).collect::<Vec<String>>().join(","))
            }
            false => {
                //println!("Done at {}",current_index);
                format!("n{}", current_index)
            }
        }
    }

}

impl LineageModel for CRISPRBits {
    fn estimated_event_space(&self) -> usize {
        self.targets_per_integration * self.integration_count
    }

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell> {
        let mut ic = input_cell.pure_clone();
        assert_eq!(self.targets_per_integration, self.rates.len());
        (0..self.integration_count).for_each(|i| {
            let existing_events = ic.events.entry(
                Genome::CRISPRBits(i)).or_insert(GenomeEventLookup::new());

            (0..self.targets_per_integration).for_each(|t| {
                existing_events.events.
                    insert(t as EventPosition,
                           self.draw_new_event(
                               &existing_events.events.get(&(t as EventPosition)).
                                   map_or(CRISPRBits::NEITHER, |x| *x), &t));
            });
        });
        vec![ic]
    }

    fn to_mix_array(&self, input_cell: &Cell) -> Vec<usize> {
        let mut ret = Vec::new();
        for i in 0..self.integration_count {
            let existing_events = input_cell.events.get(&Genome::CRISPRBits(i));
            match existing_events {
                None => {
                    println!("No CRISPRBits events found for index {:?}", input_cell);
                    panic!("No CRISPRBits events found for index {}", i);
                }
                Some(x) => {
                    for j in 0..(self.targets_per_integration as u16) {
                        for k in (0 as u16)..(4 as u16) {
                            if x.events.get(&j).unwrap() == &k {
                                ret.push(1 as usize);
                            } else {
                                ret.push(0 as usize);
                            }
                        }
                    }
                }
            }
        }
        ret
    }


    fn get_mapping(&self, x: &EventOutcomeIndex) -> String {
        String::new()
    }

}


#[cfg(test)]
mod tests {
    use std::fs;
    use std::fs::File;
    use crate::lineagemodels::model::SimpleDivision;
    use super::*;

    #[test]
    fn crispr_bits_transform_cell() {
        let equal_rates = CRISPRBitRate::new(0.98, 0.0099, 0.0099, 0.0002, 0.01);

        let cBits = CRISPRBits::new(&1, &1, vec![equal_rates.clone()]);

        let mut counts = HashMap::new();

        for _i in 0..100 {
            let mut cell = Cell::new();
            for _j in 0..100 {
                cell = cBits.divide(&mut cell).iter().next().unwrap().increment_id_clone();
            };
            let outcome = cell.events.get(&Genome::CRISPRBits(0)).unwrap().events.get(&0).unwrap().clone();
            *counts.entry(outcome).or_insert(0) += 1;
        };
        println!("{:?}", counts);
    }



}