use crate::cell::Cell;
use crate::genome::{GenomeDescription, GenomeEventCollection};
use crate::lineagemodels::model::{CellFactory, DroppedAllele, EventOutcomeIndex};

#[derive(Clone, Debug)]
pub struct EcDNA {
    edit_rate: Vec<Vec<f64>>,
    positions: Vec<u32>,
    pub targets_per_barcode: u32,
    description: String,
    genome: GenomeDescription,
    interdependent_rate: f64,
}

impl CellFactory for EcDNA {
    fn mutate(&self, input_cell: &mut Cell, genome: &mut GenomeEventCollection) -> Vec<Cell> {
        todo!()
    }

    fn to_mix_array(&self, genome: &GenomeEventCollection, input_cell: &mut Cell, force_retention: &bool) -> (DroppedAllele,Option<Vec<u8>>) {
        todo!()
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        todo!()
    }

}