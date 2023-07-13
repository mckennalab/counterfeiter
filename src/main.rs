
mod cell;
mod genome;
mod lineagemodels {
    pub mod model;
    pub mod crispr_bit;
}

extern crate lazy_static;
extern crate rand;
extern crate array_tool;

use std::collections::HashMap;
use clap::Parser;
use crate::cell::Cell;
use crate::lineagemodels::crispr_bit::{CRISPRBitRate, CRISPRBits};
use crate::lineagemodels::model::SimpleDivision;
use crate::lineagemodels::model::LineageModel;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {

    #[clap(long)]
    left: f32,

    #[clap(long)]
    right: f32,

    #[clap(long)]
    both: f32,

    #[clap(long)]
    reedit: f32,

    #[clap(long)]
    output: String,

    #[clap(long)]
    generations: usize,

    #[clap(long)]
    intergrations: usize,
}

fn main() {
    let parameters = Args::parse();

// store the relationship between children and parents
    let mut parent_child_map: HashMap<usize, Vec<usize>> = HashMap::new();

    // our starting balance rate
    let total = parameters.left + parameters.right + parameters.both;
    assert!(total < 1.0);
    let neither = 1.0 - total;

    let rates = CRISPRBitRate::new(neither, parameters.left, parameters.right, parameters.both, parameters.reedit);
    let cBits = CRISPRBits::new(&parameters.intergrations, &1, vec![rates.clone()]);
    let simple_division = SimpleDivision { offspring_count: 2 };

    let mut current_cells = vec![Cell::new()].iter().map(|x| cBits.transform(x)).next().unwrap();
    let mut generations: std::collections::HashMap<usize, Vec<Cell>> = HashMap::new();
    for i in 0..parameters.generations {
        let mut next_cells = Vec::new();

        generations.insert(i, current_cells.iter().map(|x| x.pure_clone()).collect::<Vec<Cell>>());

        for mut cell in current_cells {
            let parent_id = cell.id;

            let generated_children = cBits.transform(&mut cell).iter().
                flat_map(|x| simple_division.transform(x)).collect::<Vec<Cell>>();

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

    CRISPRBits::to_newick_tree(&parent_child_map, &"test_tree.newick".to_string());
    CRISPRBits::to_mix_input(generations.get(&parameters.generations).unwrap(), &cBits, &"test_mix_input.txt".to_string());
}
