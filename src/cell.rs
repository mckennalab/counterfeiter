use crate::allreads::{EditEvent, AllEvents};
use rand::Rng;

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone)]
pub struct Cell {
    pub events: Vec<EditEvent>,
}


impl Cell {
    pub fn new(sites: usize) -> Cell {
        let mut events: Vec<EditEvent> = Vec::new();
        for _i in 0..sites {
            events.push(EditEvent::new_WT());
        }
        assert_eq!(sites, events.len());
        Cell{events}
    }

    pub fn split(&self, cell_count: usize, mutation_rate: f64, events: &AllEvents) -> Vec<Cell> {
        let mut ret_cell: Vec<Cell> = Vec::new();
        let mut rng = rand::thread_rng();

        // for each cell we've been asked to make
        for i in 0..cell_count {
            let mut new_cell_events= self.events.clone();
            let mut cnt = 0;

            // for each event location
            self.events.iter().enumerate().for_each(|(index,event)| {
                cnt += 1;
                let rando: f64 = rng.gen();
                let occupied_sites = Cell::occupied_sites(&new_cell_events);
                // if it's been mutated, leave it be. Otherwise try to mutate it
                if rando <= mutation_rate && !event.is_mutated() {
                    let new_event = events.targets_to_events[&i].draw_compatible_event(index, &occupied_sites, 10);
                    match new_event {
                        Some(x) => {
                            let low_values = vec![index as i32 + x.offset,0];
                            let high_values = vec![(index as i32 + x.offset) + x.length as i32,self.events.len() as i32];

                            for i in (*low_values.iter().max().unwrap() as usize)..(*high_values.iter().min().unwrap() as usize) {
                                new_cell_events[i as usize] = x.clone();
                            }
                        },
                        None => {
                            // do nothing for now
                        },
                    }
                } else {
                    new_cell_events[index] = event.clone();
                }
            });
            assert_eq!(self.events.len(), cnt);
            assert_eq!(self.events.len(), new_cell_events.len());
            ret_cell.push(Cell{events: new_cell_events})
        }
        ret_cell
    }

    pub fn to_comp_string(&self) -> String {
        self.events.iter().map(|k| {k.event_string.to_owned()}).collect::<Vec<String>>().join("_")
    }

    pub fn edited_rate(&self) -> f64 {
        self.events.iter().map(|k| {k.is_mutated() as u8 as f64}).sum::<f64>() / (self.events.len() as f64)
    }

    pub fn occupied_sites(current_edits: &Vec<EditEvent>) -> Vec<usize> {
        current_edits.iter().enumerate().filter(|(index,edit)| edit.is_mutated()).map(|(index,edit)| index).collect()
    }
}




#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashSet, HashMap};
    use crate::allreads::TargetToEvents;

    fn create_stub_all_reads() -> AllEvents {
        let sites = vec![0];
        let ed = EditEvent{event_string: "Test".to_string(), occupied_sites: sites};
        let event_to_count = vec![(ed,1)];
        let tg = TargetToEvents { event_to_count};
        let mut targets_to_events: HashMap<usize, TargetToEvents> = HashMap::new();
        targets_to_events.insert(0usize, tg);
        AllEvents{targets_to_events}
    }

    #[test]
    fn test_basic_load() {
        let fake_all_reads = create_stub_all_reads();
        let cell = Cell::new(1);
        let new_cells = cell.split(1, 1.0, &fake_all_reads);
        assert_eq!(1, new_cells.len());
        assert_eq!(new_cells[0].events[0].event_string,"Test".to_string());

    }
}