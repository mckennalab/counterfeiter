use std::rc::Rc;
use std::sync::Arc;
use crate::cell::Cell;

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventPosition = u16;
// our reserved mutational 'space'
pub type EventOutcomeIndex = u16;

pub(crate) trait LineageModel {

    fn estimated_event_space(&self) -> usize;

    fn transform(&self, input_cell: &mut Cell) -> Vec<Cell>;

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String;
}

pub struct SimpleDivision {
    pub offspring_count: usize,
}

impl LineageModel for SimpleDivision {
    fn estimated_event_space(&self) -> usize {
        0
    }

    fn transform(&self, input_cell: &mut Cell) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = Vec::new();
        for _i in 0..self.offspring_count {
            new_cells.push(input_cell.clone());
        }
        new_cells
    }

    fn get_mapping(&self, index: &EventOutcomeIndex) -> String {
        format!("SimpleDivision:{}", index)
    }
}

