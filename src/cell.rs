use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use nohash_hasher::BuildNoHashHasher;
use std::sync::Arc;
use nohash_hasher::NoHashHasher;
use crate::genome::{Genome, GenomeEventLookup};


// from https://stackoverflow.com/questions/32935808/generate-sequential-ids-for-each-instance-of-a-struct
// we want to have a unique ID for each cell, so we can relate them
use std::{
    sync::atomic::{AtomicUsize, Ordering},
    thread,
};

static OBJECT_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A cell is a collection of events at target sites. This struct is kept as minimal
/// as possible to reduce memory usage, and events and sites are referred to by index.
/// The cells are passed the relevent look-up object when they 'split' to create new cells.
#[derive(Clone, Debug)]
pub struct Cell {
    id: usize,
    pub events: HashMap<Genome,GenomeEventLookup>,//, BuildHasherDefault<NoHashHasher<Genome>>>,
}

impl Cell {
    pub fn new(guessed_genome_count: &usize) -> Cell {
        //let events: HashMap::<Genome, GenomeEventLookup, BuildHasherDefault<NoHashHasher<Genome>>> =
        //    HashMap::with_capacity_and_hasher(*guessed_genome_count, BuildNoHashHasher::default());
        let id = OBJECT_COUNTER.fetch_add(1, Ordering::SeqCst);
        Cell{id, events: HashMap::new() }
    }
}

