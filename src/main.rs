mod cell;
mod genome;

mod lineagemodels {
    pub mod model;
    // pub mod crispr_bit;
    pub mod cas12a_abe;
}

extern crate lazy_static;
extern crate rand;
extern crate array_tool;

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufRead, Write};
use std::error::Error;

use std::collections::HashMap;
use clap::Parser;
use crate::cell::Cell;
use crate::lineagemodels::cas12a_abe::Cas12aABE;
//use crate::lineagemodels::crispr_bit::{CRISPRBitRate, CRISPRBits};
use crate::lineagemodels::model::SimpleDivision;
use crate::lineagemodels::model::LineageModel;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(long)]
    mpileup: String,

    #[clap(long)]
    output_mix: String,

    #[clap(long)]
    output_tree: String,

    #[clap(long)]
    generations: usize,

    #[clap(long)]
    integrations: usize,

    #[clap(long)]
    mpileupgenerations: usize,

    #[clap(long)]
    integrated_barcodes: usize,

    #[clap(long)]
    trials: usize,
}

fn main() {
    let parameters = Args::parse();

    // store the relationship between children and parents
    let mut parent_child_map: HashMap<usize, Vec<usize>> = HashMap::new();

    //     pub fn from_mpileup_file(filename: &String, minimum_coverage: &usize, minimum_mutation_rate: &f64, allowed_mutations: HashMap<u8, u8>) -> Cas12aABE {
    let mut mp: HashMap<u8, u8> = HashMap::new();
    mp.insert(b'A', b'G');


    let mut simple_division = SimpleDivision { offspring_count: 2 };

    // Open the file in append mode, creating it if it doesn't exist
    let mut file = File::create("simulation_summary.txt").unwrap();
    let mut cBits = Cas12aABE::from_mpileup_file(&parameters.mpileup,
                                                 &20,
                                                 &0.005,
                                                 mp,
                                                 "50mer".to_string(),
                                                 &parameters.mpileupgenerations,
                                                 &parameters.integrated_barcodes,
                                                 &mut file);


    file.flush();


    for trial in 0..parameters.trials {
        let mut current_cells = vec![Cell::new()].iter().map(|x| cBits.divide(x)).next().unwrap();
        let mut generations: std::collections::HashMap<usize, Vec<Cell>> = HashMap::new();

        for i in 0..parameters.generations {
            let mut next_cells = Vec::new();

            generations.insert(i, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

            for mut cell in current_cells {
                let parent_id = cell.id;

                let generated_children = cBits.divide(&mut cell).iter().
                    flat_map(|x| simple_division.divide(x)).collect::<Vec<Cell>>();

                for child in generated_children {
                    let child_id = child.id;
                    parent_child_map.entry(parent_id).or_insert(Vec::new()).push(child_id);
                    next_cells.push(child);
                }
            }

            println!("{}: {}, {}", i, next_cells.len(), next_cells.iter().next().unwrap().id);
            current_cells = next_cells;
        }

        generations.insert(parameters.generations, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());
        cBits.to_mix_input(&current_cells, &parameters.output_mix);

        let proportions = calculate_column_proportions(&parameters.output_mix).unwrap();

        // Open the file in append mode, creating it if it doesn't exist
        let mut file = OpenOptions::new()
            .append(true) // Open in append mode
            .create(true) // Create the file if it doesn't exist
            .open("simulation_summary.txt").unwrap();

        // Write the line to the file
        let output_str = proportions.iter().map(|x| {
            format!("{:.3}", *x)
        }).collect::<Vec<String>>().join("\t");
        //println!("output string {}",output_str);
        file.write_all(format!("{}\t", trial+1).as_bytes()).unwrap();

        file.write_all(output_str.as_bytes()).unwrap();
        file.write_all("\n".as_bytes()).unwrap();

        // Optionally, flush to ensure all data is written
        file.flush().unwrap();
    }
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
                _ => return Err(From::from(format!("Invalid data value '{}'", val_str))),
            }
        }

        actual_rows += 1;
    }

    if actual_rows != nrows {
        return Err(From::from(format!("Expected {} rows, found {}", nrows, actual_rows)));
    }

    // Calculate proportions
    let proportions: Vec<f64> = counts.iter()
        .map(|&count| count as f64 / nrows as f64)
        .collect();

    Ok(proportions)
}
