use crate::cell::Cell;
use crate::genome::{next_genome_id, DecayType, EditingOutcome, Genome, GenomeDescription, GenomeEventCollection, GenomeEventKey, Modification, Probability};
use crate::lineagemodels::model::{CellFactory, EventOutcomeIndex, EventPosition};
use rand::prelude::*;
use regex::{Match, Regex};
use rustc_hash::FxHashMap;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

#[derive(Clone, Debug)]
pub struct Cas12aABE {
    edit_rate: Vec<Vec<f64>>,
    proportions: Option<Vec<f64>>,
    positions: Vec<u32>,
    pub targets_per_barcode: u32,
    description: String,
    genome: GenomeDescription,
    interdependent_rate: f64,
}

impl Cas12aABE {
    /// Processes mpileup bases string and counts variant bases.
    ///
    /// Parses the complex mpileup bases format to count occurrences of each
    /// nucleotide at a given position. Handles special mpileup symbols like
    /// insertions (+), deletions (-), read starts (^), read ends ($), and
    /// reference matches (. and ,).
    ///
    /// # Arguments
    ///
    /// * `bases` - The mpileup bases string to process
    /// * `ref_base` - The reference base at this position
    ///
    /// # Returns
    ///
    /// HashMap mapping nucleotides (A, C, G, T, N) to their observed counts.
    ///
    /// # Mpileup Format
    /// - `.` or `,`: Match to reference (forward/reverse strand)
    /// - `A`, `C`, `G`, `T`: Observed bases (case indicates strand)
    /// - `^X`: Start of read (X is mapping quality)
    /// - `$`: End of read
    /// - `+N[bases]`: Insertion of N bases
    /// - `-N[bases]`: Deletion of N bases
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

    pub fn from_mpileup_file(
        filename: &String,
        minimum_coverage: &usize,
        minimum_mutation_rate: &f64,
        allowed_mutations: &HashMap<u8, u8>,
        description: &String,
        generations: &usize,
        barcode_count: &usize,
        drop_rate: &f64,
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
                let mutated_prop = *counts
                    .get(allowed_mutations.get(&ref_base).unwrap())
                    .unwrap_or(&0) as f64
                    / total_bases as f64;
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

        for i in (0..*barcode_count) {
            let final_edit_rate_barcode: Vec<f64> = mutation_counts
                .iter()
                .enumerate()
                .map(|(index, x)| {
                    let new_rate = 1.0 - f64::exp(f64::ln((1.0 - x)) / (*generations as f64));
                    new_rate
                })
                .collect();
            final_edit_rates.push(final_edit_rate_barcode);
            final_edit_rates_orig.append(&mut mutation_counts.clone());
            final_edit_positions.append(
                &mut mutation_positions
                    .iter()
                    .map(|x| *x + (i as u32 * mutation_positions.len() as u32))
                    .collect(),
            );
        }


        //println!("Editing rate size {}",mutation_positions.len());
        Cas12aABE {
            edit_rate: final_edit_rates,
            proportions: Some(mutation_counts.clone()),
            positions: final_edit_positions,
            targets_per_barcode: mutation_counts.len() as u32,
            description: description.clone(),
            genome: GenomeDescription {
                genome: Genome::ABECas12a,
                name: "Cas12aFrommPileup".to_string(),
                allows_overlap: false,
                decay_type: DecayType::Permanent,
                drop_rate: Probability::new_from_f64(drop_rate).unwrap(),
                id: next_genome_id(),
            },
            interdependent_rate: *interdependent_rate,
        }
    }

    /// Creates a Cas12a ABE system with uniform editing rates.
    ///
    /// Constructs a Cas12a adenine base editor with identical editing rates
    /// across all target sites and integrations. This is useful for controlled
    /// simulations where you want to test the effect of specific editing rates.
    ///
    /// # Arguments
    ///
    /// * `rate` - Per-generation editing rate (0.0 to 1.0) for each site
    /// * `target_count` - Number of target sites per barcode integration
    /// * `integration_count` - Number of barcode integrations in the system
    /// * `description` - Description string for this editing system
    /// * `interdependent_rate` - Rate of interdependent editing events (0.0 to 1.0)
    ///
    /// # Returns
    ///
    /// A new `Cas12aABE` instance with uniform editing rates across all sites.
    ///
    /// # Examples
    ///
    /// ```rust
    /// // Create system with 1% editing rate, 20 sites per barcode, 5 barcodes
    /// let cas12a = Cas12aABE::from_editing_rate(
    ///     &0.01, &20, &5, "uniform_1pct".to_string(), &0.0
    /// );
    /// ```
    pub fn from_editing_rate(
        rate: &f64,
        target_count: &u32,
        integration_count: &u32,
        description: String,
        drop_rate: &f64,
        interdependent_rate: &f64,
    ) -> Cas12aABE {
        let edit_rates = (0..(*integration_count))
            .map(|e| (0..*target_count).map(|x| *rate).collect::<Vec<f64>>())
            .collect::<Vec<Vec<f64>>>();

        Cas12aABE {
            edit_rate: edit_rates,
            proportions: None,
            positions: (0..(integration_count * target_count))
                .map(|x| x as u32)
                .collect::<Vec<u32>>(),
            targets_per_barcode: *target_count,
            description: description.clone(),
            genome: GenomeDescription {
                genome: Genome::ABECas12a,
                name: description,
                allows_overlap: true,
                decay_type: DecayType::Permanent,
                drop_rate: Probability::new_from_f64(drop_rate).unwrap(),
                id: next_genome_id(),
            },
            interdependent_rate: *interdependent_rate,
        }
    }
    /// Stochastically generates a new editing event at a given position.
    ///
    /// Uses the configured editing rate for the specified position to determine
    /// whether an A→G editing event occurs. Creates an appropriate EditingOutcome
    /// and registers it with the genome event collection if editing occurs.
    ///
    /// # Arguments
    ///
    /// * `position` - Global position index (across all barcodes and sites)
    /// * `genome` - Global genome event collection for registration
    ///
    /// # Returns
    ///
    /// * `Some(GenomeEventKey)` - If editing occurs, returns the event key
    /// * `None` - If no editing occurs at this position
    ///
    /// # Editing Logic
    /// 1. Calculates barcode and site indices from global position
    /// 2. Retrieves position-specific editing rate
    /// 3. Performs random sampling to determine if editing occurs
    /// 4. Creates A→G substitution outcome if editing occurs
    fn draw_new_event(
        &self,
        position: &u32,
        genome: &mut GenomeEventCollection,
    ) -> Option<GenomeEventKey> {
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

    /// Generates a Newick format phylogenetic tree from parent-child relationships.
    ///
    /// Constructs a phylogenetic tree in Newick format representing the lineage
    /// relationships between cells. Only includes cells specified in cell_ids_to_keep
    /// and uses recursive tree building to handle complex branching patterns.
    ///
    /// # Arguments
    ///
    /// * `cells` - Vector of all cells (used for validation/reference)
    /// * `parent_child_map` - Map from parent IDs to vectors of child IDs
    /// * `cell_ids_to_keep` - Map specifying which cells to include in the tree
    /// * `output` - Path to output Newick tree file
    ///
    /// # Output Format
    /// Newick format: `((child1,child2)parent1,(child3,child4)parent2)root;`
    /// - Leaf nodes: `nCELL_ID:1`
    /// - Internal nodes: `(children)nNODE_ID:1`
    /// - Tree terminated with semicolon
    ///
    /// # Side Effects
    /// Creates/overwrites the output file with the Newick tree representation.
    pub fn to_newick_tree(
        cells: &Vec<Cell>,
        parent_child_map: &HashMap<usize, Vec<usize>>,
        cell_ids_to_keep: &HashMap<usize, bool>,
        output: &String,
    ) {
        let mut out = File::create(output).unwrap();
        write!(
            out,
            "{};\n",
            Cas12aABE::recursive_tree_builder(parent_child_map, cell_ids_to_keep, &0)
        )
        .expect("Unable to write file");
    }

    /// Recursively builds Newick tree representation from parent-child relationships.
    ///
    /// Core recursive function for constructing Newick format trees. Handles
    /// various tree topology scenarios including single-child nodes, multi-child
    /// nodes, and leaf nodes. Automatically collapses single-child internal nodes
    /// and adjusts branch lengths accordingly.
    ///
    /// # Arguments
    ///
    /// * `parent_child_map` - Map from parent IDs to vectors of child IDs
    /// * `cell_ids_to_keep` - Map specifying which cells to include in output
    /// * `current_index` - Current node ID being processed
    ///
    /// # Returns
    ///
    /// String representation of the subtree rooted at current_index in Newick format.
    /// Returns empty string if the node should not be included in the final tree.
    ///
    /// # Tree Building Logic
    /// - Leaf nodes: Return `nID:1` if in cell_ids_to_keep
    /// - Single child: Collapse node and increment branch length
    /// - Multiple children: Return `(child1,child2,...)nID:1`
    /// - Filtered nodes: Return empty string to exclude from tree
    pub fn recursive_tree_builder(
        parent_child_map: &HashMap<usize, Vec<usize>>,
        cell_ids_to_keep: &HashMap<usize, bool>,
        current_index: &usize,
    ) -> String {
        match parent_child_map.contains_key(current_index) {
            true => {
                let children = parent_child_map.get(current_index).unwrap();
                let child_map = children
                    .iter()
                    .map(|x| {
                        Cas12aABE::recursive_tree_builder(parent_child_map, cell_ids_to_keep, x)
                    })
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
                } else {
                    "".to_string()
                }
            }
            false => {
                if cell_ids_to_keep.contains_key(current_index)
                    && cell_ids_to_keep.get(current_index).unwrap() == &true
                {
                    format!("n{}:1", current_index)
                } else {
                    "".to_string()
                }
            }
        }
    }

    fn from_file(file_string: &String, drop_rate: &f64) -> Self {
        let tokens: Vec<&str> = file_string.split(",").collect();
        assert!(tokens.len() == 5);
        let generations = tokens.get(0).unwrap().parse::<usize>().unwrap();
        let minimum_coverage = tokens.get(1).unwrap().parse::<usize>().unwrap();
        let minimum_mutation_rate = tokens.get(2).unwrap().parse::<f64>().unwrap();
        let interdependent_rate = tokens.get(3).unwrap().parse::<f64>().unwrap();
        let mut allowed_mutations = HashMap::new();
        allowed_mutations.insert(b'A', b'G');
        allowed_mutations.insert(b'T', b'C');
        Cas12aABE::from_mpileup_file(
            &tokens.get(3).unwrap().to_string(),
            &minimum_coverage,
            &minimum_mutation_rate,
            &allowed_mutations,
            &"FROMFILE".to_string(),
            &generations,
            &0,
            drop_rate,
            &interdependent_rate,
        )
    }
}

impl CellFactory for Cas12aABE {
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell> {
        &self
            .edit_rate
            .iter()
            .enumerate()
            .for_each(|(barcode_index, barcode)| {
                barcode
                    .iter()
                    .enumerate()
                    .for_each(|(position, edit_rate)| {
                        let pos = (position + (barcode_index * self.targets_per_barcode as usize))
                            as EventPosition;
                        let previous_pos_edited = (position
                            + (barcode_index * self.targets_per_barcode as usize))
                            as EventPosition;
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

    fn to_mix_array(
        &self,
        genome: &GenomeEventCollection,
        input_cell: &mut Cell,
    ) -> Option<Vec<u8>> {
        let mut ret = Vec::new();

        let existing_events = genome
            .filter_events_and_get_outcomes(&self.genome, &input_cell.events)
            .iter()
            .map(|x| (x.start, x.clone()))
            .collect::<HashMap<u32, EditingOutcome>>();

        let mut all_empty = true;
        for integration in 0..self.edit_rate.len() {
            let mut rng: rand::rngs::ThreadRng = rand::thread_rng();

            let draw = rng.gen::<f64>();
            //println!("draw {} < {} threshold ",draw,drop_rate);
            if draw > self.genome.drop_rate.get() {
                for i in 0..self.edit_rate.get(integration).unwrap().len() {
                    let i = i as u32;
                    //println!("size {}",existing_events.events.len());
                    all_empty = false;
                    let position_outcome =
                        existing_events.get(&(i + (integration as u32 * self.targets_per_barcode)));
                    match position_outcome {
                        None => {
                            ret.push(b'0');
                        }
                        Some(x) => match x.internal_outcome_id {
                            0 => {
                                panic!("We shouldn't store 0 outcomes")
                            }
                            1 => ret.push(b'1'),
                            _ => panic!("unknown symbol"),
                        },
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
            _ => "1".to_string(),
        }
    }

}

#[cfg(test)]
mod tests {
    use std::io::BufWriter;
    use super::*;
    use crate::lineagemodels::model::SimpleDivision;

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

        parent_child_map.insert(0, vec![1, 2]);
        parent_child_map.insert(1, vec![3, 4]);
        parent_child_map.insert(2, vec![5, 6]);
        cell_ids_to_keep.insert(3, true);
        cell_ids_to_keep.insert(5, true);
        cell_ids_to_keep.insert(4, true);
        cell_ids_to_keep.insert(6, true);
        println!(
            "tree {} ",
            Cas12aABE::recursive_tree_builder(&parent_child_map, &cell_ids_to_keep, &root_index)
        );
    }

    #[test]
    fn test_load_from_pileup_file() {
        let mut allowed_mutations = HashMap::new();
        allowed_mutations.insert(b'A', b'G');
        allowed_mutations.insert(b'T', b'C');
        let from_pileup = Cas12aABE::from_mpileup_file(
            &"test_data/bc20_mpileup.txt".to_string(),
            &20,
            &0.01,
            &allowed_mutations,
            &"testing".to_string(),
            &30,
            &1,
            &0.0,
            &0.0
        );

        let mut output_file = File::create("test_data/mileup_conversion_to_rates.txt").unwrap();
        let mut writer = BufWriter::new(output_file);
        writeln!(writer, "site\tprop\trate").unwrap();
        println!("len {}",from_pileup.edit_rate.len());
        println!("len {}",from_pileup.positions.len());
        from_pileup.edit_rate.get(0).unwrap().iter().enumerate().for_each(|(index,x)| {
            let mutation_rate = from_pileup.proportions.as_ref().unwrap().get(index).unwrap();
            writeln!(writer, "{}\t{}\t{}",index,mutation_rate,x).unwrap();    
        });
    }
}
