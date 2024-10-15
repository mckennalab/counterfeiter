use crate::cell::Cell;

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventPosition = u16;

// our reserved mutational 'space'
pub type EventOutcomeIndex = u16;

pub trait LineageModel {

    fn estimated_event_space(&self) -> usize;

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell>;

    fn to_mix_array(&self, input_cell: &Cell) -> Vec<usize>;

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String;
}

pub struct SimpleDivision {
    pub offspring_count: usize,
}

impl LineageModel for SimpleDivision {
    fn estimated_event_space(&self) -> usize {
        0
    }

    fn divide(&mut self, input_cell: &Cell) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = Vec::new();
        for _i in 0..self.offspring_count {
            new_cells.push(input_cell.increment_id_clone());
        }
        new_cells
    }

    fn to_mix_array(&self, input_cell: &Cell) -> Vec<usize> {
        Vec::new()
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        format!("SimpleDivision:{}", index)
    }
}

