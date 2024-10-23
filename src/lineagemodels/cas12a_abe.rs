use crate::cell::Cell;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, CellFactory};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;
use crate::genome::{EditingOutcome, Genome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification};
use std::io::Write;
use rand::prelude::*;

#[derive(Clone, Debug)]
pub struct Cas12aABE {
    edit_rate: Vec<Vec<f64>>,
    positions: Vec<u32>,
    pub targets_per_barcode: u32,
    description: String,
    genome: GenomeDescription,
}


impl Cas12aABE {
    pub fn from_editing_rate(rate: &f64,
                             target_count: &u32,
                             integration_count: &u32,
                             description: String,
    ) -> Cas12aABE {
        let edit_rates = (0..(*integration_count)).map(|e|
        (0..*target_count).map(|x| *rate).collect::<Vec<f64>>()).collect::<Vec<Vec<f64>>>();

        Cas12aABE {
            edit_rate: edit_rates,
            positions: (0..(integration_count * target_count)).map(|x| x as u32).collect::<Vec<u32>>(),
            targets_per_barcode: *target_count,
            description: description.clone(),
            genome: GenomeDescription {
                genome: Genome::ABECas12a,
                name: description,
                allows_overlap: false,
            },
        }
    }
    fn draw_new_event(&self, position: &u32, genome: &mut GenomeEventCollection) -> Option<GenomeEventKey> {
        let target = position % self.targets_per_barcode;
        let barcode = position / self.targets_per_barcode;
        let proportion = self.edit_rate.get(barcode as usize).unwrap();
        let proportion = proportion.get(target as usize).unwrap();

        //let mut rng = rand::rng();
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();

        let rando = rng.gen::<f64>();
        if rando <= *proportion {
            let outcome = EditingOutcome {
                start: *position,
                stop: *position + 1, // [x,y) intervals
                change: Modification::Substitution,
                nucleotides: vec![b'A'],
                internal_outcome_id: 1,
            };
            genome.add_event(&self.genome, outcome)

        } else {
            None
        }
    }

    pub fn to_mix_input(&mut self, genome: &GenomeEventCollection, cells: &mut Vec<Cell>, drop_rate: &f64, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        cells.iter_mut().for_each(|cell| {
            match self.to_mix_array(genome, drop_rate, cell) {
                Some(x) => {
                    cell_to_output.insert(cell.id, x);
                    if event_len.is_none() {
                        event_len = Some(cell_to_output.get(&cell.id).unwrap().len());
                    } else {
                        assert_eq!(cell_to_output.get(&cell.id).unwrap().len(), event_len.unwrap());
                    }
                }
                None => {}
            }
        });


        let mut out = File::create(output).unwrap();
        write!(out, "\t{}\t{}\n", cell_to_output.len(), self.targets_per_barcode as usize * self.edit_rate.len()).unwrap();
        cell_to_output.iter().for_each(|(k, v)| {
            write!(out, "{:<10}\t{}\n", format!("n{}", k), String::from_utf8(v.clone()).unwrap());
        });
    }


    pub fn to_newick_tree(cells: &Vec<Cell>, parent_child_map: &HashMap<usize, Vec<usize>>, output: &String) {
        let mut out = File::create(output).unwrap();
        write!(out, "{};\n", Cas12aABE::recursive_tree_builder(cells, parent_child_map, &0)).expect("Unable to write file");
    }

    pub fn recursive_tree_builder(cells: &Vec<Cell>, parent_child_map: &HashMap<usize, Vec<usize>>, current_index: &usize) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                let child_map = children.iter().map(|x| Cas12aABE::recursive_tree_builder(cells, parent_child_map, x))
                    .filter(|x| x != &"".to_string())
                    .collect::<Vec<String>>().join(",");
                if child_map.len() > 0 { format!("({})", child_map) } else { format!("") }
            }
            false => {
                let mut is_silent = true;
                for cell in cells {
                    if cell.id == *current_index && cell.visible_in_output {
                        is_silent = false;
                    }
                }
                if is_silent {
                    format!("")
                } else {
                    format!("n{}", current_index)
                }
            }
        }
    }
}


impl CellFactory for Cas12aABE {
    fn estimated_event_space(&self) -> usize {
        self.positions.len()
    }

    fn divide(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell> {
        &self.edit_rate.iter().enumerate().for_each(|(barcode_index, barcode)| {
            barcode.iter().enumerate().for_each(|(position, edit_rate)| {
                let pos = (position + (barcode_index * self.targets_per_barcode as usize)) as EventPosition;
                match self.draw_new_event(&(pos as u32), genome) {
                    None => {}
                    Some(x) => {
                        input_cell.events.insert(x);
                    }
                }
            });
        });
        vec![input_cell.pure_clone()]
    }

    fn to_mix_array(&mut self, genome: &GenomeEventCollection, drop_rate: &f64, input_cell: &mut Cell) -> Option<Vec<u8>> {
        let mut ret = Vec::new();

        let existing_events = genome.filter_events_and_get_Outcomes(&self.genome, &input_cell.events).iter().map(|x| {
            (x.start,x.clone())
        }).collect::<HashMap<u32, EditingOutcome>>();


        let mut all_empty = true;
        for integration in 0..self.edit_rate.len() {
            let mut rng: rand::rngs::ThreadRng = rand::thread_rng();

            let draw = rng.gen::<f64>();
            //println!("draw {} < {} threshold ",draw,drop_rate);
            if draw > *drop_rate {
                for i in 0..self.edit_rate.get(integration).unwrap().len() {
                    let i = i as u32;
                    //println!("size {}",existing_events.events.len());
                    all_empty = false;
                    let position_outcome = existing_events.get(&(i + (integration as u32 * self.targets_per_barcode)));
                    match position_outcome {
                        None => {
                            ret.push(b'0');
                        }
                        Some(x) => {
                            match x.internal_outcome_id {
                                0 => {panic!("We shouldn't store 0 outcomes")}
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
        if all_empty {
            input_cell.visible_in_output = false;
            None
        } else {
            Some(ret)
        }
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        match index {
            0 => "0".to_string(),
            _ => "1".to_string()
        }
    }
}