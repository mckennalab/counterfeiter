use std::collections::HashMap;
use crate::lineagemodels::model::{EventOutcomeIndex, EventPosition};

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
#[repr(u8)]
pub enum Genome {
    MITOCHONDRIAL(usize),
    NUCLEAR(usize),
    CRISPRBits(usize),
    THREE_D_POSITION(usize),
    FOUR_D_POSITION(usize),
}

#[derive(Clone, Debug)]
pub struct GenomeEventLookup {
    pub events: HashMap<EventPosition, EventOutcomeIndex>,
}

impl GenomeEventLookup {
    pub fn new() -> GenomeEventLookup {
        GenomeEventLookup{events: HashMap::new()}
    }
}

