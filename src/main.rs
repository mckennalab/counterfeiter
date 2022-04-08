#![feature(iter_advance_by)]

mod allreads;
mod cell;
//mod cell;

#[macro_use]
extern crate lazy_static;
extern crate rand;
extern crate array_tool;

use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::str;
use std::sync::mpsc;
use std::time::{Duration, SystemTime};
use std::io::Write;
use std::sync::mpsc::sync_channel;
use std::thread;
use std::fs;
use std::collections::HashSet;

use std::io;
use clap::Parser;
use std::sync::{Arc, Mutex};

use allreads::read_all_read_counts;
use cell::Cell;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(long)]
    targets: usize,

    #[clap(long)]
    allreads: String,

    #[clap(long)]
    generations: usize,

    #[clap(long)]
    rate: f64,
}

fn main() {
    let parameters = Args::parse();

    let known_events = read_all_read_counts(&parameters.allreads, parameters.targets).unwrap();

    //let output = Arc::new(Mutex::new(File::create(parameters.output).unwrap()));
    let mut cells : Vec<Cell> = Vec::new();
    cells.push(Cell::new(parameters.targets));
    for i in 0..parameters.generations {
        let mut new_cells : Vec<Cell> = Vec::new();
        for cell in &cells {
            cell.split(2,parameters.rate, &known_events).iter().for_each(|ncell| new_cells.push(ncell.to_owned()));
        }
        let mut books = HashSet::new();
        println!("NCell size {}",new_cells.len());
        for cl in new_cells.clone().iter() {books.insert(cl.to_comp_string());};
        println!("Book size {}",books.len());
        for b in books {
            println!("b size {}",b);
        }

        cells = new_cells;

        println!("Cell size {}",cells.len());
    }
}
