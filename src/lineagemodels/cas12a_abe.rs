use crate::cell::Cell;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition, CellFactory};
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;
use crate::genome::{EditingOutcome, Genome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification};
use std::io::Write;
use rustc_hash::FxHashMap;
use rand::prelude::*;
use regex::{Match, Regex};


#[derive(Clone, Debug)]
pub struct Cas12aABE {
    edit_rate: Vec<Vec<f64>>,
    positions: Vec<u32>,
    pub targets_per_barcode: u32,
    description: String,
    genome: GenomeDescription,
    interdependent_rate: f64,
}


impl Cas12aABE {
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

    pub fn from_mpileup_file(filename: &String,
                             minimum_coverage: &usize,
                             minimum_mutation_rate: &f64,
                             allowed_mutations: &HashMap<u8, u8>,
                             description: &String,
                             generations: &usize,
                             duplicate_barcodes: &usize,
                             output_file: &mut File,
                             interdependent_rate: &f64,

    ) -> Cas12aABE {
        let mut mutation_counts: Vec<f64> = Vec::new();
        let mut mutation_positions: Vec<u32> = Vec::new();

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);

        // Process each line of the Mpileup file
        for (index, line) in reader.lines().enumerate() {
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

        let mut final_edit_rates: Vec<Vec<f64>> = Vec::new();
        let mut final_edit_rates_orig: Vec<f64> = Vec::new();
        let mut final_edit_positions: Vec<u32> = Vec::new();

        println!("edits {:?}", mutation_counts);
        println!("edits len {}", mutation_counts.len());

        for i in (0..*duplicate_barcodes) {
            let final_edit_rate_barcode: Vec<f64> = mutation_counts.iter().enumerate().map(|(index, x)| {
                let new_rate = 1.0 - f64::exp(f64::ln((1.0 - x)) / (*generations as f64));
                new_rate
            }).collect();
            final_edit_rates.push(final_edit_rate_barcode);
            final_edit_rates_orig.append(&mut mutation_counts.clone());
            final_edit_positions.append(&mut mutation_positions.iter().map(|x| *x + (i as u32 * mutation_positions.len() as u32)).collect());
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
        Cas12aABE {
            edit_rate: final_edit_rates,
            positions: final_edit_positions,
            targets_per_barcode: mutation_counts.len() as u32,
            description: description.clone(),
            genome: GenomeDescription {
                genome: Genome::ABECas12a,
                name: "Cas12aFrommPileup".to_string(),
                allows_overlap: false,
            },
            interdependent_rate: *interdependent_rate,
        }
    }

    pub fn from_editing_rate(rate: &f64,
                             target_count: &u32,
                             integration_count: &u32,
                             description: String,
                             interdependent_rate: &f64,
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
                allows_overlap: true,
            },
            interdependent_rate: *interdependent_rate,
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
                nucleotides: vec![b'G'],
                internal_outcome_id: 1,
            };
            genome.add_event(&self.genome, outcome)
        } else {
            None
        }
    }

    pub fn to_mix_input(&mut self, genome: &GenomeEventCollection, cells: &mut Vec<Cell>, drop_rate: &f64, cell_ids_to_keep: &mut HashMap<usize, bool>, output: &String) {
        let mut cell_to_output = HashMap::new();
        let mut event_len: Option<usize> = None;
        cells.iter_mut().for_each(|cell| {
            if cell_ids_to_keep.contains_key(&cell.id) {
                match self.to_mix_array(genome, drop_rate, cell) {
                    Some(x) => {
                        cell_to_output.insert(cell.id, x);
                        if event_len.is_none() {
                            event_len = Some(cell_to_output.get(&cell.id).unwrap().len());
                        } else {
                            assert_eq!(cell_to_output.get(&cell.id).unwrap().len(), event_len.unwrap());
                        }
                    }
                    None => {
                        cell_ids_to_keep.remove(&cell.id);
                    }
                }
            }
        });


        let mut out = File::create(output).unwrap();
        write!(out, "\t{}\t{}\n", cell_to_output.len(), self.targets_per_barcode as usize * self.edit_rate.len()).unwrap();
        cell_to_output.iter().for_each(|(k, v)| {
            write!(out, "{:<10}\t{}\n", format!("n{}", k), String::from_utf8(v.clone()).unwrap());
        });
    }


    pub fn to_newick_tree(cells: &Vec<Cell>, parent_child_map: &HashMap<usize, Vec<usize>>, cell_ids_to_keep: &HashMap<usize, bool>, output: &String) {
        let mut out = File::create(output).unwrap();
        write!(out, "{};\n", Cas12aABE::recursive_tree_builder(parent_child_map, cell_ids_to_keep, &0)).expect("Unable to write file");
    }

    pub fn recursive_tree_builder(parent_child_map: &HashMap<usize, Vec<usize>>,
                                  cell_ids_to_keep: &HashMap<usize, bool>,
                                  current_index: &usize) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                let child_map = children.iter().map(|x| Cas12aABE::recursive_tree_builder(parent_child_map, cell_ids_to_keep, x))
                    .filter(|x| x != &"".to_string())
                    .collect::<Vec<String>>();
                if child_map.len() == 1 && child_map.get(0).unwrap() == "" {
                    "".to_string()
                } else if child_map.len() == 1 {
                    let re = Regex::new(r"^(.*):(\d+)$").unwrap();
                    if let Some(caps) = re.captures(child_map.get(0).unwrap()) {
                        // Extract the node and value parts as strings
                        let node = caps.get(1).unwrap().as_str();
                        let value: i32 = caps.get(2).unwrap().as_str().parse().unwrap();
                        // Format the result back into the desired string
                        format!("{}:{}", node, value + 1)
                    } else {
                        // Return None if the input string doesn't match the pattern
                        panic!("no match for {}", child_map.get(0).unwrap());
                    }
                } else if child_map.len() > 0 {
                    format!("({})n{}:1", child_map.join(","), current_index)
                } else { "".to_string() }
            }
            false => {
                if cell_ids_to_keep.contains_key(current_index) && cell_ids_to_keep.get(current_index).unwrap() == &true {
                    format!("n{}:1", current_index)
                } else {
                    "".to_string()
                }
            }
        }
    }
}

impl CellFactory for Cas12aABE {
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell> {
        &self.edit_rate.iter().enumerate().for_each(|(barcode_index, barcode)| {
            barcode.iter().enumerate().for_each(|(position, edit_rate)| {
                let pos = (position + (barcode_index * self.targets_per_barcode as usize)) as EventPosition;
                let previous_pos_edited = (position + (barcode_index * self.targets_per_barcode as usize)) as EventPosition;
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

        let existing_events = genome.filter_events_and_get_outcomes(&self.genome, &input_cell.events).iter().map(|x| {
            (x.start, x.clone())
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
                                0 => { panic!("We shouldn't store 0 outcomes") }
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

    fn from_file(file_string: &String) -> Self {
        let tokens: Vec<&str> = file_string.split(",").collect();
        assert!(tokens.len() == 5);
        let generations = tokens.get(0).unwrap().parse::<usize>().unwrap();
        let minimum_coverage = tokens.get(1).unwrap().parse::<usize>().unwrap();
        let minimum_mutation_rate = tokens.get(2).unwrap().parse::<f64>().unwrap();
        let interdependent_rate = tokens.get(3).unwrap().parse::<f64>().unwrap();
        let mut output_file = File::open(tokens.get(4).unwrap()).unwrap();
        let mut allowed_mutations = HashMap::new();
        allowed_mutations.insert(b'A', b'G');
        allowed_mutations.insert(b'T', b'C');
        Cas12aABE::from_mpileup_file(&tokens.get(3).unwrap().to_string(),
                                     &minimum_coverage,
                                     &minimum_mutation_rate,
                                     &allowed_mutations,
                                     &"FROMFILE".to_string(),
                                     &generations,
                                     &0,
                                     &mut output_file,
                                     &interdependent_rate)
    }
}


#[cfg(test)]
mod tests {
    use crate::lineagemodels::model::SimpleDivision;
    use super::*;

    #[test]
    fn test_iterative_example() {
        /*

    pub fn iterative_tree_builder(
        cells: &Vec<Cell>,
        parent_child_map: &HashMap<usize, Vec<usize>>,
        cell_ids_to_keep: &HashMap<usize,bool>,
        root_index: &usize
         */
        let mut parent_child_map: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut cell_ids_to_keep: HashMap<usize, bool> = HashMap::new();
        let root_index: usize = 0;

        parent_child_map.insert(0, vec!(1, 2));
        parent_child_map.insert(1, vec!(3, 4));
        parent_child_map.insert(2, vec!(5, 6));
        cell_ids_to_keep.insert(3, true);
        cell_ids_to_keep.insert(5, true);
        cell_ids_to_keep.insert(4, true);
        cell_ids_to_keep.insert(6, true);
        println!("tree {} ", Cas12aABE::recursive_tree_builder(&parent_child_map, &cell_ids_to_keep, &root_index));
    }
}