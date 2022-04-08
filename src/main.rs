
mod allreads;
mod cell;

extern crate lazy_static;
extern crate rand;
extern crate array_tool;

use std::collections::HashSet;

use clap::Parser;

use allreads::read_all_read_counts;
use cell::Cell;
use std::fs::File;
use std::str::FromStr;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(long)]
    targets: usize,

    #[clap(long)]
    allreads: String,

    #[clap(long)]
    output: String,

    #[clap(long)]
    generations: usize,

    #[clap(long)]
    enumerate: String,

    #[clap(long)]
    rate: f64,
}

fn main() {
    let parameters = Args::parse();

    let known_events = read_all_read_counts(&parameters.allreads).unwrap();

    let _output = File::create(parameters.output).unwrap();

    let enumerated = parameters.enumerate.split(",").collect::<Vec<&str>>().iter().map(|x| usize::from_str(x).unwrap()).collect();
    let known_events = known_events.emulate_sites(enumerated);

    let mut cells : Vec<Cell> = Vec::new();
    cells.push(Cell::new(parameters.targets));

    for _i in 0..parameters.generations {
        let mut new_cells : Vec<Cell> = Vec::new();
        for cell in &cells {
            cell.split(2,parameters.rate, &known_events).iter().for_each(|ncell| new_cells.push(ncell.to_owned()));
        }
        let mut books = HashSet::new();
        let mut avg_editing_rate = 0.0;
        for cl in new_cells.clone().iter() {
            books.insert(cl.to_comp_string());
            avg_editing_rate += cl.edited_rate();
        };
        println!("Number of cells: {}, number of unique alleles: {}, editing rate: {}",new_cells.len(),books.len(), avg_editing_rate/(new_cells.len() as f64));
        println!("Book 1: {}",books.iter().next().unwrap());
        cells = new_cells;
    }
}
