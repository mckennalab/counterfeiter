use std::collections::{HashMap, HashSet};
use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, IntervalTree};


#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub enum Genome {
    CRISPRBits,
    ABECas12a,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct GenomeDescription {
    pub genome: Genome,
    pub name: String,
    pub allows_overlap: bool,
}

/// Compactly store the outcome of mutating the genome. This stores indexes into
/// the global registry of events; the goal is simply to store an offset into a much richer
/// database. Ordered to avoid splitting the u16s and getting weird alignment in the worst case
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct GenomeEventKey {
    pub position_index: u32,
    pub genome_index: u16,
    pub outcome_index: u16,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub enum Modification {
    Substitution,
    Insertion,
    Deletion,
    Inversion
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct EditingOutcome {
    pub start: u32,
    pub stop: u32,
    pub change: Modification,
    pub nucleotides: Vec<u8>,
    pub internal_outcome_id: u16,
}

/**
This struct can be somewhat heavyweight
*
#[derive(Clone)]
pub struct GenomeEventCollection {
    genomes: HashMap<GenomeDescription,ArrayBackedIntervalTree<u32,EditingOutcome>>,
    genomes_offsets: HashMap<GenomeDescription,u16>,
    current_genome_offset: u16,
    key_to_outcome: HashMap<GenomeEventKey,EditingOutcome>,
    outcome_to_key: HashMap<EditingOutcome,GenomeEventKey>,
}

impl GenomeEventCollection {
    pub fn new() -> GenomeEventCollection {
        GenomeEventCollection{
            genomes: HashMap::new(),
            genomes_offsets: HashMap::new(),
            current_genome_offset: 0,
            key_to_outcome: HashMap::new(),
            outcome_to_key: HashMap::new() }
    }

    fn create_event(&mut self, genome: &GenomeDescription, outcome: &EditingOutcome) -> GenomeEventKey {

        GenomeEventKey{
            position_index: outcome.start.clone(),
            genome_index: *self.genomes_offsets.get(genome).unwrap(),
            outcome_index: outcome.internal_outcome_id.clone(),
        }
    }

    pub fn filter_events_and_get_Outcomes(&self, genome: &GenomeDescription, set: &HashSet<GenomeEventKey>) -> Vec<EditingOutcome> {
        let genome_offset = self.genomes_offsets.get(genome).unwrap();
        set.iter().map(|it| {
            match self.key_to_outcome.get(it) {
                None => {
                    panic!("We dont know about event {:?}",it);
                }
                Some(x) => {
                    match *genome_offset == it.genome_index {
                        true => {
                            Some(x.clone())
                        }
                        false => {
                            None
                        }
                    }
                }
            }
        }).flatten().collect::<Vec<EditingOutcome>>()
    }

    pub fn add_event(&mut self, genome: &GenomeDescription, outcome: EditingOutcome) -> Option<GenomeEventKey> {

        match self.genomes.get_mut(genome) {
            None => {
                // we dont have a record of this genome -- set it up and add the new outcome
                let mut lt = ArrayBackedIntervalTree::new();
                lt.insert(outcome.start..outcome.stop, outcome.clone());
                lt.index();
                self.genomes.insert(genome.clone(), lt);
                self.genomes_offsets.insert(genome.clone(), self.current_genome_offset.clone());
                self.current_genome_offset += 1;
                let id = self.create_event(genome,&outcome);

                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());
                Some(id)
            }
            Some(x) => {
                // we have a record -- check (1) for overlap and (2) if we allow that.
                let range = outcome.start..outcome.stop;
                let id = self.create_event(genome,&outcome);
                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());

                match genome.allows_overlap {
                    true => {
                        // just add it
                        let tr = self.genomes.get_mut(genome).unwrap();
                        tr.insert(range, outcome);
                        tr.index();
                        Some(id)
                    }
                    false => {
                        // check if an overlapping event exists
                        match self.genomes.get_mut(genome).unwrap().find(range.clone()).len() {
                            0 => {
                                let tr = self.genomes.get_mut(genome).unwrap();
                                tr.insert(range, outcome);
                                tr.index();
                                Some(id)
                            }
                            _ => {
                                None
                            }
                        }
                    }
                }
            }
        }
    }
}
*/
#[derive(Clone)]
pub struct GenomeEventCollection {
    genomes: HashMap<GenomeDescription,HashMap<u32,EditingOutcome>>,
    genomes_offsets: HashMap<GenomeDescription,u16>,
    current_genome_offset: u16,
    key_to_outcome: HashMap<GenomeEventKey,EditingOutcome>,
    outcome_to_key: HashMap<EditingOutcome,GenomeEventKey>,
}

impl GenomeEventCollection {
    pub fn new() -> GenomeEventCollection {
        GenomeEventCollection{
            genomes: HashMap::new(),
            genomes_offsets: HashMap::new(),
            current_genome_offset: 0,
            key_to_outcome: HashMap::new(),
            outcome_to_key: HashMap::new() }
    }

    fn create_event(&mut self, genome: &GenomeDescription, outcome: &EditingOutcome) -> GenomeEventKey {

        GenomeEventKey{
            position_index: outcome.start.clone(),
            genome_index: *self.genomes_offsets.get(genome).unwrap(),
            outcome_index: outcome.internal_outcome_id.clone(),
        }
    }

    pub fn filter_events_and_get_outcomes(&self, genome: &GenomeDescription, set: &HashSet<GenomeEventKey>) -> Vec<EditingOutcome> {
        if !self.genomes.contains_key(genome) {
            Vec::new()
        } else {
            let genome_offset = self.genomes_offsets.get(genome).unwrap();
            set.iter().map(|it| {
                match self.key_to_outcome.get(it) {
                    None => {
                        panic!("We dont know about event {:?}", it);
                    }
                    Some(x) => {
                        match *genome_offset == it.genome_index {
                            true => {
                                Some(x.clone())
                            }
                            false => {
                                None
                            }
                        }
                    }
                }
            }).flatten().collect::<Vec<EditingOutcome>>()
        }
    }

    pub fn add_event(&mut self, genome: &GenomeDescription, outcome: EditingOutcome) -> Option<GenomeEventKey> {

        match self.genomes.contains_key(genome) {
            false => {
                // we dont have a record of this genome -- set it up and add the new outcome
                let mut lt = HashMap::new();
                lt.insert(outcome.start, outcome.clone());
                self.genomes.insert(genome.clone(), lt);
                println!("adding to find key {:?}",genome);

                self.genomes_offsets.insert(genome.clone(), self.current_genome_offset.clone());
                self.current_genome_offset += 1;
                let id = self.create_event(genome,&outcome);

                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());
                Some(id)
            }
            true => {
                // we have a record -- check (1) for overlap and (2) if we allow that.
                let id = self.create_event(genome,&outcome);
                self.key_to_outcome.insert(id.clone(), outcome.clone());
                self.outcome_to_key.insert(outcome.clone(), id.clone());

                match self.genomes.get(genome).unwrap().contains_key(&outcome.start) {
                    false => {
                        self.genomes.get_mut(genome).unwrap().insert(outcome.start,outcome.clone());
                        let id = self.create_event(genome,&outcome);

                        self.key_to_outcome.insert(id.clone(), outcome.clone());
                        self.outcome_to_key.insert(outcome.clone(), id.clone());
                        Some(id)
                    },
                    true => {
                        Some(self.outcome_to_key.get(&outcome).unwrap().clone())
                    }
                }
            }
        }
    }
}


