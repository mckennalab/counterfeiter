use crate::allreads::{EditEvent, AllEvents, UNEDITED};
use std::borrow::Borrow;
use rand::Rng;

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone)]
pub struct Cell {
    pub events: Vec<EditEvent>,
}


impl Cell {
    pub fn new(sites: usize) -> Cell {
        let mut events: Vec<EditEvent> = Vec::new();
        for i in 0..sites {
            events.push(EditEvent::new_none());
        }
        assert_eq!(sites, events.len());
        Cell{events}
    }

    pub fn split(&self, cell_count: usize, mutation_rate: f64, events: &AllEvents) -> Vec<Cell> {
        let mut ret_cell: Vec<Cell> = Vec::new();
        let mut rng = rand::thread_rng();

        // for each cell we've been asked to make
        for i in 0..cell_count {
            let mut new_cell_events: Vec<EditEvent> = Vec::new();

            // for each event location
            self.events.iter().for_each(|event| {
                let rando: f64 = rng.gen();

                // if it's been mutated, leave it be. Otherwise try to mutate it
                if event.is_mutated() {
                    new_cell_events.push(event.clone());
                }
                else if rando <= mutation_rate {
                    let new_event = events.targets_to_events[&i].draw_compatible_event(event.occupied_sites.borrow(), 10);
                    match new_event {
                        Some(x) => {
                            println!("{}",x.event_string);
                            new_cell_events.push(x)
                        },
                        None => {
                            println!("NO DRAW ");
                            new_cell_events.push(EditEvent::new(i))
                        },
                    }
                }
            });
            ret_cell.push(Cell{events: new_cell_events})
        }
        ret_cell
    }

    pub fn to_comp_string(&self) -> String {
        let ret = self.events.iter().map(|k| {k.event_string.to_owned()}).collect::<Vec<String>>().join("_");
        //println!("ret {}",ret);
        ret

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