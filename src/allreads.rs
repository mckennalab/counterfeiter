use std::collections::{HashMap};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::str;
use rand::prelude::SliceRandom;
use std::hash::{Hasher, Hash};

// the current simulation caps this at 65K unique events -- we can raise this in the future
pub type EventIndex = u16;

// two reserved indices -- the wildtype sequence and an unknown event
pub const WT_INDEX: EventIndex = 0;
pub const UNKNOWN_INDEX: EventIndex = 1;

pub type EventPosition = u16;

#[derive(PartialOrd, Ord, Clone)]
pub struct EditEvent {
    pub event_string: String,

    /// these are both zero for NONE (WT edits), or other indicators of reversible editing
    pub offset: i32,
    pub length: usize,
}

impl Hash for EditEvent {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.event_string.hash(state)
    }
}

impl PartialEq for EditEvent {
    fn eq(&self, other: &Self) -> bool {
        self.event_string == other.event_string
    }
}
impl Eq for EditEvent {}

pub(crate) static UNEDITED: &str = "NONE";

impl EditEvent {
    pub fn new_none_vec(size: usize) -> Vec<EditEvent> {
        let mut newvec: Vec<EditEvent> = Vec::new();
        for _i in 0..size {
            newvec.push(EditEvent::new_WT());
        }
        newvec
    }

    pub fn new_WT() -> EditEvent {
        EditEvent{event_string: UNEDITED.to_string(), offset: 0, length: 0 }
    }

    pub fn is_mutated(&self) -> bool {
        self.event_string != UNEDITED.to_string()
    }
}

// ************************************************************
//
// A collection of editing events observed at a target site, with
// their counts
//
#[derive(Clone)]
pub struct TargetToEvents {
    pub event_to_count: Vec<(EditEvent, usize)>,
    pub event_to_eventID: HashMap<EditEvent, EventIndex>,
}

impl TargetToEvents {
    pub fn draw_weighted_event(&self) -> EditEvent {
        let mut rng = rand::thread_rng();
        self.event_to_count.choose_weighted(&mut rng, |item| item.1).unwrap().to_owned().0
    }

    pub fn draw_compatible_event(&self, site: usize, occupied_sites: &Vec<usize>, max_tries: usize) -> Option<EditEvent> {
        for _i in 0..max_tries {
            let evt = self.draw_weighted_event();

            let min_site = (site as i32 + evt.offset) as usize;
            let max_mite = min_site + evt.length;
            if !occupied_sites.contains(&min_site) && !occupied_sites.contains(&max_mite) {
                return Some(evt);
            }
        }
        None
    }
}

// ************************************************************
//
// over all targets, what's the collection of editing outcomes
//
pub struct AllEvents {
    pub targets_to_events: HashMap<usize, TargetToEvents>,
}

impl AllEvents {
    /// Create a new AllEvents with the old_sites vector enumerating each new site in the
    /// produced AllEvents.
    pub fn emulate_sites(&self, old_sites: Vec<usize>) -> AllEvents {
        let mut targets_to_events : HashMap<usize, TargetToEvents> = HashMap::new();
        old_sites.iter().enumerate().for_each(|(index,old_site)| {
            targets_to_events.insert(index,self.targets_to_events[old_site].clone());
        });
        AllEvents{targets_to_events}
    }
}
// ************************************************************
//
// handle reading an allEventCounts file into an AllEvents object
//
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_all_read_counts(filename: &String) -> Option<AllEvents> {

    let mut targets_to_events: HashMap<usize, HashMap<EditEvent, usize>> = HashMap::new();

    if let Ok(mut lines) = read_lines(filename) {
        lines.next();
        for line in lines {
            if let Ok(ip) = line {
                // split on the columns, then split the edits by their separator "_"
                let split = ip.split("\t");
                let columns: Vec<&str> = split.collect();
                let edits: Vec<&str> = columns[0].split("_").collect();

                // if our target_to_events hasn't been setup, use the first line to create the mappings
                if targets_to_events.len() == 0 {
                    for i in 0..edits.len() {
                        let event_to_count: HashMap<EditEvent, usize> = HashMap::new();
                        let target_to_events = event_to_count;
                        targets_to_events.insert(i,target_to_events);
                    }
                }

                let count: usize = columns[2].parse::<usize>().unwrap();

                let mut event_to_targets: HashMap<String, Vec<usize>> = HashMap::new();
                edits.iter().enumerate().for_each(|(index, edit)| {
                    let default_vec: Vec<usize> = Vec::new();
                    event_to_targets.entry(edit.to_string()).or_insert(default_vec).push(index);
                });

                event_to_targets.iter().for_each(|(event_string, sites)| {
                    sites.iter().enumerate().for_each(|(index,site)| {
                        let evt = EditEvent { event_string: event_string.to_string(), offset: -1 * index as i32, length: sites.len()};
                        let contains_key = targets_to_events[site].contains_key(&evt);
                        let new_value = count + if contains_key { targets_to_events[site][&evt]} else { 0 };
                        targets_to_events.get_mut(site).unwrap().insert(evt, new_value);
                    });
                });
            }
        }
    }

    // now make this into our structured data
    let mut targets_to_events_final: HashMap<usize, TargetToEvents> = HashMap::new();
    for i in 0..targets_to_events.len() {
        let tte = TargetToEvents{event_to_count: targets_to_events.get_mut(&i).unwrap().into_iter().map(|(k, v)| (k.to_owned(),v.clone())).collect()};
        targets_to_events_final.insert(i, tte);
    }
    Some(AllEvents{targets_to_events: targets_to_events_final})
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_basic_load() {
        let events = read_all_read_counts(&"data/Lba1.allReadCounts".to_string(),15);
        let events_clone = events.unwrap();
        assert_eq!(292, events_clone.targets_to_events[&14].event_to_count.len());
    }

    #[test]
    fn test_random_draw() {
        let events = read_all_read_counts(&"data/Lba1.allReadCounts".to_string(),15).unwrap();
        let mut books: HashMap<EditEvent,usize> = HashMap::new();
        for i in 0..1000 {
            let evt = events.targets_to_events[&14].draw_weighted_event();
            if books.contains_key(&evt.to_owned()) {
                books.insert(evt.to_owned(), books[&evt.to_owned()] + 1usize);
            } else {
                books.insert(evt.to_owned(),1usize);
            }
        }
        for b in books {
            println!("b {}->{}",b.0.event_string, b.1);
        }

    }
}