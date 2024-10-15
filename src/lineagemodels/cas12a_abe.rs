use crate::cell::Cell;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, LineageModel};
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
    edit_rate: Vec<f64>,
    positions: Vec<u32>,
    pub original_size: usize,
    random: StdRng,
    description: String,
}


impl Cas12aABE {

    pub fn total_edit_sites(&self) -> usize {
        self.edit_rate.len()
    }

    // Function to process the pileup bases string and count variant bases
    fn process_bases(bases: &str, ref_base: u8) -> HashMap<u8, usize> {
        let mut counts: HashMap<u8, usize> = HashMap::new();
        let mut i = 0;
        let chars: Vec<char> = bases.chars().collect();
        let mut skip_next = 0;

        while i < chars.len() {
            if skip_next > 0 {
                skip_next -= 1;
                i += 1;
                continue;
            }

            match chars[i] {
                '^' => {
                    // Start of a read segment; skip the next character (mapping quality)
                    i += 2;
                }
                '$' => {
                    // End of a read segment; move to the next character
                    i += 1;
                }
                '+' | '-' => {
                    // Insertion or deletion
                    let sign = chars[i];
                    i += 1;
                    // Extract the number of bases inserted/deleted
                    let mut num_str = String::new();
                    while i < chars.len() && chars[i].is_digit(10) {
                        num_str.push(chars[i]);
                        i += 1;
                    }
                    let num_bases: usize = num_str.parse().unwrap_or(0);
                    // Skip the inserted or deleted bases
                    skip_next = num_bases;
                }
                '.' | ',' => {
                    // Match to the reference base
                    let base = ref_base.to_ascii_uppercase();
                    *counts.entry(base).or_insert(0) += 1;
                    i += 1;
                }
                'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n' => {
                    // Observed base; convert to uppercase
                    let base = chars[i].to_ascii_uppercase();
                    *counts.entry(base as u8).or_insert(0) += 1;
                    i += 1;
                }
                _ => {
                    // Other symbols; skip
                    i += 1;
                }
            }
        }

        counts
    }
    pub fn from_editing_rate(rate: &f64,
        target_count: &usize,
                             description: String,
    ) -> Cas12aABE {
        Cas12aABE{
            edit_rate: vec![*rate; *target_count],
            positions: (0..*target_count).map(|x| x as u32).collect::<Vec<u32>>(),
            original_size: *target_count,
            random: StdRng::from_entropy(),
            description,
        }
    }
    pub fn from_mpileup_file(filename: &String,
                             minimum_coverage: &usize,
                             minimum_mutation_rate: &f64,
                             allowed_mutations: HashMap<u8, u8>,
                             description: String,
                             generations: &usize,
                            duplicate_barcodes: &usize,
                             output_file: &mut File,
    ) -> Cas12aABE {

        let mut mutation_counts: Vec<f64> = Vec::new();
        let mut mutation_positions: Vec<u32> = Vec::new();

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);

        // Process each line of the Mpileup file
        for (index,line) in reader.lines().enumerate() {
            let line = line.unwrap();
            if line.trim().is_empty() {
                //eprintln!("Warning: Skipping empty line: {}", line);
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                // eprintln!("Warning: Skipping malformed line: {}", line);
                continue;
            }

            let pos = u32::from_str(fields[1]).unwrap();
            let ref_base = fields[2].as_bytes()[0];
            if allowed_mutations.contains_key(&ref_base) {
                let read_count: usize = fields[3].parse().unwrap_or(0);
                let bases = fields[4];

                let counts = Cas12aABE::process_bases(bases, ref_base);

                // Calculate frequencies
                let total_bases: usize = counts.values().sum();
                let mutated_prop = *counts.get(allowed_mutations.get(&ref_base).unwrap()).unwrap_or(&0) as f64 / total_bases as f64;
                //let mutated_prop = mutated_prop / (*generations as f64);
                if total_bases >= *minimum_coverage && mutated_prop >= *minimum_mutation_rate {
                    mutation_counts.push(mutated_prop);
                    mutation_positions.push(pos);
                }
            } else {
                //eprintln!("Warning: Skipping non-targeted base line: {}", line);
            }
        }
        // now duplicate out to the number of barcodes we'll use

        let mut final_edit_rates : Vec<f64> = Vec::new();
        let mut final_edit_rates_orig : Vec<f64> = Vec::new();
        let mut final_edit_positions : Vec<u32> = Vec::new();

        println!("edits {:?}",mutation_counts);

        for i in (0..*duplicate_barcodes) {
            final_edit_rates.append(&mut mutation_counts.iter().enumerate().map(|(index,x)| {
                let new_rate = 1.0 - f64::exp(f64::ln((1.0-x))/7.0);
                new_rate
            }).collect());
            final_edit_rates_orig.append(&mut mutation_counts.clone());
            final_edit_positions.append(&mut mutation_positions.iter().map(|x|*x + (i as u32 * mutation_positions.len() as u32)).collect());
        }

        // header
        output_file.write_all(format!("run").as_bytes());
        for i in 0..final_edit_rates_orig.len() {
            output_file.write_all(format!("\tindex{}", i).as_bytes());
        }
        output_file.write_all(format!("\n").as_bytes());

        // summary info
        let output_str = final_edit_rates_orig.iter().map(|x| {
            format!("{:.3}", *x)
        }).collect::<Vec<String>>().join("\t");

        //println!("output string {}",output_str);
        output_file.write_all(format!("0\t{}\n", output_str).as_bytes()).unwrap();

            //println!("Editing rate size {}",mutation_positions.len());
        Cas12aABE { edit_rate: final_edit_rates, positions: final_edit_positions, original_size: mutation_counts.len(), random: StdRng::from_entropy(), description: description.clone() }
    }

    fn draw_new_event(&mut self, current_state: EventOutcomeIndex, position: &usize) -> EventOutcomeIndex {
        match current_state {
            0 => {
                let proportion = self.edit_rate.get(*position).unwrap();
                let rando = self.random.gen::<f64>();
                //println!("proportion {} rate {} state {}",rando, *proportion, rando  <= *proportion);
                if rando <= *proportion {
                    //println!("1 proportion {} rate {} state {}",rando, *proportion, rando  <= *proportion);
                    1 as EventOutcomeIndex
                } else {
                    //println!("0 proportion {} rate {} state {}",rando, *proportion, rando  <= *proportion);
                    0 as EventOutcomeIndex
                }
            }
            _ => {
                //println!("current state {} ",current_state);
                current_state
            }
        }
    }

    pub fn to_mix_input(&self, cells: &Vec<Cell>, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        for cell in cells {
            cell_to_output.insert(cell.id, self.to_mix_array(&cell));
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
            if counts.len() >= 1 {
                kept_columns.push(i);
            } else {
                println!("drop");
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
        write!(out, "{};\n", Cas12aABE::recursive_tree_builder(parent_child_map, &0)).expect("Unable to write file");
    }

    pub fn recursive_tree_builder(parent_child_map: &HashMap<usize, Vec<usize>>, current_index: &usize) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                //format!("(({})n{})", children.iter().map(|x| recursive_tree_builder(parent_child_map, x)).collect::<Vec<String>>().join(","), current_index)
                format!("({})", children.iter().map(|x| Cas12aABE::recursive_tree_builder(parent_child_map, x)).collect::<Vec<String>>().join(","))
            }
            false => {
                //println!("Done at {}",current_index);
                format!("n{}", current_index)
            }
        }
    }
}


impl LineageModel for Cas12aABE {
    fn estimated_event_space(&self) -> usize {
        self.positions.len()
    }

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell> {
        let mut ic = input_cell.pure_clone();
        let existing_events = ic.events.entry(
            Genome::ABECas12a(self.description.clone())).or_insert(GenomeEventLookup::new());

        (0..self.edit_rate.len()).for_each(|t| {
            let pos = t as EventPosition;
            let new_event = self.draw_new_event(existing_events.events.get(&(t as EventPosition)).map_or(0 as EventOutcomeIndex, |x| *x), &t);
            //println!("pos {} new event {}",pos,new_event);
            existing_events.events.
                insert(pos, new_event);
        });
        vec![ic]
    }

    fn to_mix_array(&self, input_cell: &Cell) -> Vec<usize> {
        let mut ret = Vec::new();
        let existing_events = input_cell.events.get(&Genome::ABECas12a(self.description.clone())).unwrap();
        //println!("Event size space {}", existing_events.events.len());
        for i in 0..self.edit_rate.len() {
            let position_outcome = existing_events.events.get(&(i as EventPosition));
            match position_outcome {
                None => {
                    println!("No Cas12aABE events found for index {:?}", input_cell);
                    panic!("No Cas12aABE events found for index {}", i);
                }
                Some(x) => {
                    ret.push(*x as usize);
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