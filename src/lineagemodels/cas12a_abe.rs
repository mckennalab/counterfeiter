use crate::cell::Cell;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, CellFactory};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;
use rand::prelude::StdRng;
use rand::prelude::*;
use crate::genome::{Genome, GenomeEventLookup};
use std::io::Write;

#[derive(Clone, Debug)]
pub struct Cas12aABE {
    edit_rate: Vec<Vec<f64>>,
    positions: Vec<u32>,
    pub targets_per_barcode: usize,
    random: StdRng,
    description: String,
}


impl Cas12aABE {


    pub fn from_editing_rate(rate: &f64,
                             target_count: &usize,
                             integration_count: &usize,
                             description: String,
    ) -> Cas12aABE {
        let edit_rates = (0..(*integration_count)).map(|e|
        (0..*target_count).map(|x| *rate).collect::<Vec<f64>>()).collect::<Vec<Vec<f64>>>();

        Cas12aABE {
            edit_rate: edit_rates,
            positions: (0..(integration_count * target_count)).map(|x| x as u32).collect::<Vec<u32>>(),
            targets_per_barcode: *target_count,
            random: StdRng::from_entropy(),
            description,
        }
    }
    fn draw_new_event(&mut self, current_state: EventOutcomeIndex, position: &usize) -> EventOutcomeIndex {
        let target = position % self.targets_per_barcode;
        let barcode = position / self.targets_per_barcode;
        match current_state {
            0 => {
                let proportion = self.edit_rate.get(barcode).unwrap();
                let proportion = proportion.get(target).unwrap();

                let rando = self.random.gen::<f64>();
                if rando <= *proportion {
                    1 as EventOutcomeIndex
                } else {
                    0 as EventOutcomeIndex
                }
            }
            _ => {
                //println!("current state {} ",current_state);
                current_state
            }
        }
    }

    pub fn to_mix_input(&mut self, cells: &Vec<Cell>, drop_rate: &f64, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        for cell in cells {
            cell_to_output.insert(cell.id, self.to_mix_array(drop_rate, &cell));
            if event_len.is_none() {
                event_len = Some(cell_to_output.get(&cell.id).unwrap().len());
            } else {
                assert_eq!(cell_to_output.get(&cell.id).unwrap().len(), event_len.unwrap());
            }
        }

        let mut out = File::create(output).unwrap();
        write!(out, "\t{}\t{}\n", cell_to_output.len(), self.targets_per_barcode * self.edit_rate.len()).unwrap();
        cell_to_output.iter().for_each(|(k, v)| {
            write!(out, "{:<10}\t{}\n", format!("n{}", k), String::from_utf8(v.clone()).unwrap());
        });
    }


    pub fn to_newick_tree(parent_child_map: &HashMap<usize, Vec<usize>>, output: &String) {
        let mut out = File::create(output).unwrap();
        write!(out, "{};\n", Cas12aABE::recursive_tree_builder(parent_child_map, &0)).expect("Unable to write file");
    }

    pub fn recursive_tree_builder(parent_child_map: &HashMap<usize, Vec<usize>>, current_index: &usize) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                format!("({})", children.iter().map(|x| Cas12aABE::recursive_tree_builder(parent_child_map, x)).collect::<Vec<String>>().join(","))
            }
            false => {
                format!("n{}", current_index)
            }
        }
    }
}


impl CellFactory for Cas12aABE {
    fn estimated_event_space(&self) -> usize {
        self.positions.len()
    }

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell> {
        let mut ic = input_cell.pure_clone();

        let existing_events = ic.events.entry(
            Genome::ABECas12a(self.description.clone())).or_insert(GenomeEventLookup::new());

        let edit_rates = self.edit_rate.clone();
        edit_rates.iter().enumerate().for_each(|(barcode_index,barcode)| {
            barcode.iter().enumerate().for_each(| (position, edit_rate)| {
                let pos = (position + (barcode_index * self.targets_per_barcode)) as EventPosition;
                let new_event = self.draw_new_event(existing_events.events.get(&pos).map_or(0 as EventOutcomeIndex, |x| *x), &(pos as usize));
                //println!("pos {} new event {}",pos,new_event);
                existing_events.events.
                    insert(pos, new_event);
            });
        });
        vec![ic]
    }

    fn to_mix_array(&mut self, drop_rate: &f64, input_cell: &Cell) -> Vec<u8> {
        let mut ret = Vec::new();
        let existing_events = input_cell.events.get(&Genome::ABECas12a(self.description.clone())).unwrap();
        //println!("Event size space {}", existing_events.events.len());
        for integration in 0..self.edit_rate.len() {
            let draw = self.random.gen::<f64>();
            //println!("draw {} < {} threshold ",draw,drop_rate);
            if draw > *drop_rate {
                for i in 0..self.edit_rate.get(integration).unwrap().len() {
                    //println!("size {}",existing_events.events.len());

                    let position_outcome = existing_events.events.get(&((i + (integration * self.targets_per_barcode)) as EventPosition));
                    match position_outcome {
                        None => {
                            println!("No Cas12aABE events found for index {:?}", input_cell);
                            panic!("No Cas12aABE events found for index {}", i);
                        }
                        Some(x) => {
                            match x {
                                0 => ret.push(b'0'),
                                1 => ret.push(b'1'),
                                _ => panic!("unknown symbol")
                            }
                        }
                    }
                }
            } else {
                for i in 0..self.edit_rate.get(integration).unwrap().len() {
                    ret.push(b'?');
                }
            }
        }
        //println!("Ret {:?}",ret);
        ret
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        match index {
            0 => "0".to_string(),
            _ => "1".to_string()
        }
    }
}