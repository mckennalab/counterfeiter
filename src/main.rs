mod cell;
mod genome;

mod lineagemodels {
    pub mod model;
    // pub mod crispr_bit;
    pub mod cas12a_abe;
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
use crate::genome::GenomeEventCollection;
use crate::lineagemodels::cas12a_abe::Cas12aABE;
//use crate::lineagemodels::crispr_bit::{CRISPRBitRate, CRISPRBits};
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

fn main() {
    PrettyTrace::new().ctrlc().on();

    if let Err(_) = std::env::var("RUST_LOG") {
        std::env::set_var("RUST_LOG", "info");
    }

    pretty_env_logger::init_timed();

    let parameters = Args::parse();
    trace!("{:?}", &parameters.cmd);


    // store the relationship between children and parents
    let mut parent_child_map: HashMap<usize, Vec<usize>> = HashMap::new();

    //     pub fn from_mpileup_file(filename: &String, minimum_coverage: &usize, minimum_mutation_rate: &f64, allowed_mutations: HashMap<u8, u8>) -> Cas12aABE {
    let mut mp: HashMap<u8, u8> = HashMap::new();
    mp.insert(b'A', b'G');

    let mut simple_division = SimpleDivision { offspring_count: 2 };

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
            let mut cas12a = Cas12aABE::from_editing_rate(
                editrate,
                sites_per_barcode,
                integrated_barcodes,
                "50mer".to_string(),
                interdependent_rate,
            );

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
            Cas12aABE::to_newick_tree(&current_cells, &parent_child_map, &cell_ids_to_keep, &(*output_tree).to_string());
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
            let mut cas12a = Cas12aABE::from_mpileup_file(mpileup,
                                                          &100,
                                                          &0.001,
                                                          &allowed_mutations,
                                                          &"MPILEUP".to_string(),
                                                          mpileupgenerations,
                                                          &(*integrated_barcodes as usize),
                                                          &mut fl,
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

        }
    }
}

fn subsample<T: Clone>(vec: &Vec<T>, sample_size: usize) -> Vec<T> {
    let mut rng = thread_rng();

    // Shuffle the original vector and take the first `sample_size` elements
    let mut vec_clone = vec.clone();
    vec_clone.shuffle(&mut rng);

    vec_clone.into_iter().take(sample_size).collect()
}

/// Reads a MIX format file and calculates the proportions of 1s in each column.
///
/// The MIX format file should have:
/// - A header line containing two integers: the number of rows and the number of columns.
/// - Subsequent lines containing a whitespace-padded name followed by a series of 0s and 1s.
///
/// # Arguments
///
/// * `filename` - The path to the MIX format file.
///
/// # Returns
///
/// A `Result` containing a vector of proportions (as `f64`) if successful, or an error.
///
/// # Example
///
/// ```rust
/// let proportions = calculate_column_proportions("data.mix")?;
/// println!("{:?}", proportions);
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

    // Initialize counts vector
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
