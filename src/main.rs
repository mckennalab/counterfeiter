mod cell;
mod genome;

mod lineagemodels {
    pub mod model;
    // pub mod crispr_bit;
    pub mod cas12a_abe;
    pub mod ec_dna;
    
    pub mod cas9_WT;
}

extern crate rand;
extern crate array_tool;
extern crate bio;
extern crate rustc_hash;
extern crate rand_distr;
#[macro_use]
extern crate log;

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufRead, Write};
use std::error::Error;
use rustc_hash::FxHashMap;
use rand::seq::SliceRandom;
use rand::thread_rng;


use std::collections::HashMap;
use clap::{Parser, Subcommand};
use pretty_trace::PrettyTrace;
use crate::cell::Cell;
use crate::genome::{create_mix_file, GenomeEventCollection};
use crate::lineagemodels::cas12a_abe::Cas12aABE;
use crate::lineagemodels::cas9_WT::Cas9WT;
use crate::lineagemodels::model::{DivisionModel, SimpleDivision};
use crate::lineagemodels::model::CellFactory;

#[derive(Subcommand, Debug)]
enum Cmd {
    Pileup {
        #[clap(long)]
        mpileup: String,

        #[clap(long)]
        output_mix: String,

        #[clap(long)]
        output_tree: String,

        #[clap(long)]
        generations: usize,

        #[clap(long)]
        mpileupgenerations: usize,

        #[clap(long)]
        integrated_barcodes: u32,

        #[clap(long)]
        barcode_drop_rate: f64,

        #[clap(long)]
        trials: usize,

        #[clap(long)]
        subsampling_number: usize,

        #[clap(long)]
        interdependent_rate: f64,

    },
    Rate {
        #[clap(long)]
        output_mix: String,

        #[clap(long)]
        output_tree: String,

        #[clap(long)]
        generations: usize,

        #[clap(long)]
        sites_per_barcode: u32,

        #[clap(long)]
        integrated_barcodes: u32,

        #[clap(long)]
        barcode_drop_rate: f64,

        #[clap(long)]
        editrate: f64,

        #[clap(long)]
        subsampling_number: usize,

        #[clap(long)]
        interdependent_rate: f64,
    },
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    cmd: Cmd,
}

/// Main entry point for the counterfeiter simulation application.
///
/// Initializes logging, parses command-line arguments, and orchestrates the entire
/// simulation workflow including cell division, editing events, and output generation.
///
/// The function supports two main simulation modes:
/// - Rate: Uses fixed editing rates across all sites
/// - Pileup: Derives editing rates from mpileup file data
///
/// # Workflow
/// 1. Setup logging and error handling
/// 2. Parse command-line arguments
/// 3. Initialize editing system (Cas12a ABE)
/// 4. Run simulation for specified generations
/// 5. Generate MIX and Newick tree outputs
fn main() {
    PrettyTrace::new().ctrlc().on();

    if let Err(_) = std::env::var("RUST_LOG") {
        std::env::set_var("RUST_LOG", "info");
    }

    pretty_env_logger::init_timed();

    let parameters = Args::parse();
    trace!("{:?}", &parameters.cmd);


    //     pub fn from_mpileup_file(filename: &String, minimum_coverage: &usize, minimum_mutation_rate: &f64, allowed_mutations: HashMap<u8, u8>) -> Cas12aABE {
    let mut mp: HashMap<u8, u8> = HashMap::new();
    mp.insert(b'A', b'G');


    let mut allowed_mutations = HashMap::new();
    allowed_mutations.insert(b'A', b'G');
    allowed_mutations.insert(b'T', b'C');
    let mut fl = File::create("pileup_output_txt").unwrap();

    match &parameters.cmd {
        Cmd::Rate {
            output_mix,
            output_tree,
            generations,
            sites_per_barcode,
            integrated_barcodes,
            barcode_drop_rate,
            editrate,
            subsampling_number,
            interdependent_rate
        } => {
            /*let cas12a = Box::new(Cas12aABE::from_editing_rate(
                editrate,
                sites_per_barcode,
                integrated_barcodes,
                "50mer".to_string(),
                barcode_drop_rate,
                interdependent_rate,
            ));*/
            let cas9 = Box::new(Cas9WT::from_editing_rate(
                editrate,
                sites_per_barcode,
                integrated_barcodes,
                0.01,
                "Cas9WT".to_string(),
                barcode_drop_rate,
                &0.0,
            ));
            println!("Cas9 length {}",cas9.len());
            let simple_division = Box::new(SimpleDivision { offspring_count: 2 });

            let mutators: Vec<Box<dyn CellFactory>> = cas9.into_iter().map(|x| {let xx : Box<dyn CellFactory> = Box::new(x); xx}).collect();
            let dividers: Vec<Box<dyn DivisionModel>>= vec![simple_division];
            
            let mut genome = GenomeEventCollection::new();

            let mut simulation_results = run_simulation(&mut genome, &1, generations, &mutators, &dividers);
            
            // subsample cells that we're interested in the final generation
            let cell_ids: Vec<usize> = simulation_results.cells.iter().map(|c| c.id.clone()).collect();
            let mut cell_ids_to_keep: HashMap<usize, bool> = subsample(&cell_ids, *subsampling_number).iter().map(|x| (*x, true)).collect();

            simulation_results.cell_generations.insert(*generations, simulation_results.cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

            println!("generating mix input file");
            create_mix_file(&genome, &mutators, &mut simulation_results.cells, &mut cell_ids_to_keep, output_mix);

            println!("generating tree file");
            Cas12aABE::to_newick_tree(&simulation_results.cells, &simulation_results.parent_child_map, &cell_ids_to_keep, &(*output_tree).to_string());
            println!("generating summary");
        }
        Cmd::Pileup {
            mpileup,
            output_mix,
            output_tree,
            generations,
            mpileupgenerations,
            integrated_barcodes,
            barcode_drop_rate,
            trials,
            subsampling_number,
            interdependent_rate
        } => {
            /*let mut cas12a = Cas12aABE::from_mpileup_file(mpileup,
                                                          &100,
                                                          &0.001,
                                                          &allowed_mutations,
                                                          &"MPILEUP".to_string(),
                                                          mpileupgenerations,
                                                          &(*integrated_barcodes as usize),
                                                          interdependent_rate);


            
            let mut genome = GenomeEventCollection::new();

            let mut current_cells = vec![Cell::new()].into_iter().map(|mut x| cas12a.mutate(&mut x, &mut genome)).next().unwrap();
            let mut cell_generations: FxHashMap<usize, Vec<Cell>> = FxHashMap::default();

            for i in 0..*generations {
                let mut next_cells = Vec::new();

                cell_generations.insert(i, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

                current_cells.into_iter().for_each(|mut cell| {
                    let parent_id = cell.id;

                    let generated_children = cas12a.mutate(&mut cell, &mut genome).into_iter().
                        flat_map(|mut x| simple_division.divide(&mut x, &mut genome)).collect::<Vec<Cell>>();

                    for child in generated_children {
                        let child_id = child.id;
                        parent_child_map.entry(parent_id).or_insert(Vec::new()).push(child_id);
                        next_cells.push(child);
                    }
                });

                println!("{}: {}, {}", i, next_cells.len(), next_cells.iter().next().unwrap().id);
                current_cells = next_cells;
            }

            // subsample cells that we're interested in the final generation
            let cell_ids: Vec<usize> = current_cells.iter().map(|c| c.id.clone()).collect();
            let mut cell_ids_to_keep: HashMap<usize, bool> = subsample(&cell_ids, *subsampling_number).iter().map(|x| (*x, true)).collect();

            cell_generations.insert(*generations, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

            println!("generating mix input file");
            cas12a.to_mix_input(&genome, &mut current_cells, barcode_drop_rate, &mut cell_ids_to_keep, output_mix);

            println!("generating tree file");
            Cas12aABE::to_newick_tree(&current_cells, &parent_child_map, &cell_ids_to_keep, &output_tree.to_string());
            println!("generating summary");
*/
        }
    }
}

pub struct SimulationResult {
    pub cells: Vec<Cell>,
    pub cell_generations: FxHashMap<usize, Vec<Cell>>,
    pub parent_child_map: HashMap<usize, Vec<usize>>,
}

/// Runs a cell division simulation over multiple generations with modular cell factories.
///
/// This function simulates cell division and mutation events over a specified number of
/// generations using a composable set of cell factories (dividers and mutators). Each
/// factory in the pipeline processes cells sequentially, allowing for complex behaviors
/// like editing, division, and barcode dropout.
///
/// # Arguments
///
/// * `initial_cell_count` - Number of cells to start the simulation with
/// * `generations` - Number of division cycles to simulate
/// * `cell_dividers` - Vector of trait objects implementing `CellFactory` that define
///                     the mutation and division behaviors. These are applied sequentially
///                     to each cell in the order provided.
///
/// # Returns
///
/// A vector of all cells present in the final generation after completing all
/// division cycles.
///
/// # Cell Factory Pipeline
///
/// The function applies each `CellFactory` in sequence to every cell:
/// 1. Start with parent cell
/// 2. Apply first factory (e.g., editing events)
/// 3. Apply second factory (e.g., cell division) to all outputs from step 2
/// 4. Continue through all factories
/// 5. Collect final daughter cells for next generation
///
/// # Internal Tracking
///
/// The function maintains:
/// - `genome`: Global event collection tracking all mutation events
/// - `cell_generations`: Historical snapshot of cells at each generation
/// - `parent_child_map`: Lineage relationships between parent and daughter cells
///
/// # Example
///
/// ```rust
/// let factories: Vec<Box<dyn CellFactory>> = vec![
///     Box::new(editing_system),
///     Box::new(simple_division),
/// ];
/// let final_cells = run_simulation(&1, &10, factories);
/// println!("Final population: {} cells", final_cells.len());
/// ```
///
/// # Performance Note
///
/// Prints progress to stdout showing generation number, cell count, and first cell ID
/// for each generation.
pub fn run_simulation(genome: &mut GenomeEventCollection,
                      initial_cell_count: &usize,
                      generations: &usize,
                      cell_mutators: &Vec<Box<dyn CellFactory>>,
                      cell_dividers: &Vec<Box<dyn DivisionModel>>,
) -> SimulationResult {

    let mut current_cells = vec![Cell::new(); *initial_cell_count];

    let mut cell_generations: FxHashMap<usize, Vec<Cell>> = FxHashMap::default();

    // store the relationship between children and parents
    let mut parent_child_map: HashMap<usize, Vec<usize>> = HashMap::new();


    for i in 0..*generations {
        let mut next_cells : Vec<Cell> = Vec::new();

        cell_generations.insert(i, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

        current_cells.into_iter().for_each(|mut cell| {

            let parent_id = cell.id;


            let mut daughters = vec![cell];

            cell_mutators.iter().for_each(|divider| {
                let mut new_daughters = vec![];
                daughters.iter().for_each( |child_id| {
                    let mut child = child_id.clone();
                    new_daughters.extend(divider.as_ref().mutate(&mut child, genome).into_iter());
                });
                daughters = new_daughters;
            });

            cell_dividers.iter().for_each(|divider| {
                let mut new_daughters = vec![];
                daughters.iter().for_each( |child_id| {
                    let mut child = child_id.clone();
                    new_daughters.extend(divider.as_ref().divide(&mut child, genome).into_iter());
                });
                daughters = new_daughters;
            });

            daughters.into_iter().for_each(|daughter| {
                let child_id = daughter.id;
                parent_child_map.entry(parent_id).or_insert(Vec::new()).push(child_id);
                next_cells.push(daughter);
            });
        });

        println!("{}: {}, {}", i, next_cells.len(), next_cells.iter().next().unwrap().id);
        current_cells = next_cells;
    }
    
    SimulationResult {
        cells: current_cells,
        cell_generations,
        parent_child_map,
    }
    
}



/// Randomly subsamples a vector to a specified size.
///
/// Creates a random subset of the input vector by shuffling all elements
/// and taking the first `sample_size` items. This is used to select a
/// representative subset of cells from the final generation for analysis.
///
/// # Arguments
///
/// * `vec` - The input vector to subsample from
/// * `sample_size` - The desired number of elements in the output
///
/// # Returns
///
/// A vector containing `sample_size` randomly selected elements from the input.
/// If `sample_size` is larger than the input vector, returns all elements.
///
/// # Type Parameters
///
/// * `T` - The type of elements in the vector, must implement `Clone`
fn subsample<T: Clone>(vec: &Vec<T>, sample_size: usize) -> Vec<T> {
    let mut rng = thread_rng();

    // Shuffle the original vector and take the first `sample_size` elements
    let mut vec_clone = vec.clone();
    vec_clone.shuffle(&mut rng);

    vec_clone.into_iter().take(sample_size).collect()
}

/// Reads a MIX format file and calculates the proportions of 1s in each column.
///
/// The MIX format is a phylogenetic analysis format where each row represents
/// a taxon (cell) and each column represents a character (editing site).
/// This function analyzes the editing patterns across all cells.
///
/// # MIX File Format
/// - Header line: two integers (number of rows, number of columns)
/// - Data lines: taxon name followed by binary character string
/// - Characters: '0' (unedited), '1' (edited), '?' (missing data)
///
/// # Arguments
///
/// * `filename` - Path to the MIX format file to analyze
///
/// # Returns
///
/// * `Ok(Vec<f64>)` - Vector of proportions (0.0 to 1.0) for each column
/// * `Err(Box<dyn Error>)` - File I/O or parsing error
///
/// # Errors
///
/// This function will return an error if:
/// - The file cannot be opened
/// - The header format is invalid
/// - The number of rows doesn't match the header
/// - Invalid characters are encountered (other than 0, 1, ?)
///
/// # Example
///
/// ```rust
/// let proportions = calculate_column_proportions(&"data.mix".to_string())?;
/// println!("Edit proportions: {:?}", proportions);
/// ```
fn calculate_column_proportions(filename: &String) -> Result<Vec<f64>, Box<dyn Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut lines = reader.lines();

    // Read the header line
    let header_line = match lines.next() {
        Some(line) => line?,
        None => return Err(From::from("Empty file")),
    };

    // Parse the header line to get nrows and ncols
    let header_parts: Vec<&str> = header_line.trim().split_whitespace().collect();
    if header_parts.len() != 2 {
        return Err(From::from("Header line must contain exactly two numbers"));
    }
    
    let nrows: usize = header_parts[0].parse()?;
    let ncols: usize = header_parts[1].parse()?;

    // Initialize 'counts' vector
    let mut counts = vec![0usize; ncols];
    let mut actual_rows = 0usize;

    // Process each line
    for line_result in lines {
        let line = line_result?;
        // Split the line into words
        let words: Vec<&str> = line.trim().split_whitespace().collect();

        // The name is the first word
        let _name = words[0];

        for (i, val_str) in words[1].as_bytes().iter().enumerate() {
            match *val_str {
                b'1' => counts[i] += 1,
                b'0' => {}
                _ => {} //return Err(From::from(format!("Invalid data value '{}'", val_str))),
            }
        }

        actual_rows += 1;
    }

    if actual_rows != nrows {
        return Err(From::from(format!("Expected {} rows, found {}", nrows, actual_rows)));
    }

    // Calculate proportions
    let proportions: Vec<f64> = counts.iter()
        .map(|&count| count as f64 / f64::max(nrows as f64, 1.0))
        .collect();

    Ok(proportions)
}
