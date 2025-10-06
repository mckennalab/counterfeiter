use std::collections::{HashMap};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::str;
use rand::prelude::SliceRandom;
use std::hash::{Hasher, Hash};

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
    /// Creates a vector of wild-type (unedited) events.
    ///
    /// Generates a vector containing the specified number of unedited
    /// EditEvent instances. Used for initializing target sites before
    /// applying editing events.
    ///
    /// # Arguments
    ///
    /// * `size` - Number of unedited events to create
    ///
    /// # Returns
    ///
    /// Vector of `EditEvent`s, all representing wild-type (unedited) state.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let events = EditEvent::new_none_vec(5);
    /// assert_eq!(events.len(), 5);
    /// assert!(events.iter().all(|e| !e.is_mutated()));
    /// ```
    pub fn new_none_vec(size: usize) -> Vec<EditEvent> {
        let mut newvec: Vec<EditEvent> = Vec::new();
        for _i in 0..size {
            newvec.push(EditEvent::new_WT());
        }
        newvec
    }

    /// Creates a new wild-type (unedited) event.
    ///
    /// Generates an EditEvent representing the unedited state of a target site.
    /// This is the default state before any editing events occur.
    ///
    /// # Returns
    ///
    /// An `EditEvent` with:
    /// - event_string: "NONE" (indicating no editing)
    /// - offset: 0
    /// - length: 0
    ///
    /// # Examples
    ///
    /// ```rust
    /// let wt_event = EditEvent::new_WT();
    /// assert!(!wt_event.is_mutated());
    /// assert_eq!(wt_event.event_string, "NONE");
    /// ```
    pub fn new_WT() -> EditEvent {
        EditEvent{event_string: UNEDITED.to_string(), offset: 0, length: 0 }
    }

    /// Checks if this event represents a mutation.
    ///
    /// Determines whether the event represents an edited state by comparing
    /// the event string to the unedited constant. Returns false for wild-type
    /// events and true for any editing outcomes.
    ///
    /// # Returns
    ///
    /// * `true` - If the event represents an editing outcome
    /// * `false` - If the event is wild-type (unedited)
    ///
    /// # Examples
    ///
    /// ```rust
    /// let wt_event = EditEvent::new_WT();
    /// assert!(!wt_event.is_mutated());
    ///
    /// let edited_event = EditEvent {
    ///     event_string: "DELETION_5BP".to_string(),
    ///     offset: -2,
    ///     length: 5
    /// };
    /// assert!(edited_event.is_mutated());
    /// ```
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
    /// Randomly selects an editing event weighted by observed frequency.
    ///
    /// Uses weighted random sampling to select from the available editing
    /// events at this target site. Events with higher observed counts
    /// are more likely to be selected, reflecting experimental frequencies.
    ///
    /// # Returns
    ///
    /// An `EditEvent` selected according to the weighted probability
    /// distribution based on observed event counts.
    ///
    /// # Panics
    ///
    /// Panics if no events are available or if all weights are zero.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let target_events = TargetToEvents { /* ... */ };
    /// let event = target_events.draw_weighted_event();
    /// // Event selected based on experimental frequencies
    /// ```
    pub fn draw_weighted_event(&self) -> EditEvent {
        let mut rng = rand::thread_rng();
        self.event_to_count.choose_weighted(&mut rng, |item| item.1).unwrap().to_owned().0
    }

    /// Attempts to draw an editing event that doesn't conflict with occupied sites.
    ///
    /// Repeatedly samples weighted events until finding one that doesn't
    /// overlap with already occupied genomic positions. This prevents
    /// overlapping edits that might be impossible in reality.
    ///
    /// # Arguments
    ///
    /// * `site` - The target site position to edit
    /// * `occupied_sites` - Vector of already occupied positions
    /// * `max_tries` - Maximum number of sampling attempts
    ///
    /// # Returns
    ///
    /// * `Some(EditEvent)` - Compatible event that doesn't overlap
    /// * `None` - If no compatible event found within max_tries
    ///
    /// # Examples
    ///
    /// ```rust
    /// let target_events = TargetToEvents { /* ... */ };
    /// let occupied = vec![10, 11, 15];
    /// let event = target_events.draw_compatible_event(5, &occupied, 100);
    /// if let Some(evt) = event {
    ///     // Event doesn't conflict with positions 10, 11, or 15
    /// }
    /// ```
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
    /// Creates a new AllEvents collection by remapping sites from an existing collection.
    ///
    /// This function allows creating a subset or rearrangement of target sites
    /// by specifying which sites from the original collection should be included
    /// and in what order. The new collection will have sites numbered sequentially
    /// starting from 0, regardless of the original site indices.
    ///
    /// # Arguments
    ///
    /// * `old_sites` - Vector of site indices from the current collection to include
    ///                 in the new collection
    ///
    /// # Returns
    ///
    /// A new `AllEvents` collection where:
    /// - Site 0 contains events from `self.targets_to_events[old_sites[0]]`
    /// - Site 1 contains events from `self.targets_to_events[old_sites[1]]`
    /// - And so on...
    ///
    /// # Panics
    ///
    /// Panics if any index in `old_sites` is not present in the current collection.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let all_events = AllEvents { /* ... */ };
    /// // Create subset with only sites 5, 10, and 15 from original
    /// let subset = all_events.emulate_sites(vec![5, 10, 15]);
    /// // New collection has 3 sites (0, 1, 2) with events from original sites 5, 10, 15
    /// ```
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
/// Opens a file and returns an iterator over its lines.
///
/// A utility function that opens a file and creates a buffered reader
/// with line-by-line iteration capability. This is used for reading
/// large data files efficiently without loading the entire file into memory.
///
/// # Arguments
///
/// * `filename` - Path to the file to read (anything that implements `AsRef<Path>`)
///
/// # Returns
///
/// * `Ok(Lines<BufReader<File>>)` - Iterator over file lines
/// * `Err(io::Error)` - If the file cannot be opened
///
/// # Examples
///
/// ```rust
/// if let Ok(lines) = read_lines("data.txt") {
///     for line in lines {
///         if let Ok(content) = line {
///             println!("{}", content);
///         }
///     }
/// }
/// ```
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// Reads an allEventCounts file and constructs an AllEvents collection.
///
/// Parses a tab-separated file containing editing event data and builds
/// a comprehensive collection of events organized by target site. The file
/// format expects one header line followed by data lines with event patterns,
/// frequencies, and counts.
///
/// # File Format
/// - Header line: ignored (typically column names)
/// - Data lines: `event1_event2_...` `\t` `frequency` `\t` `count`
/// - Events separated by underscores represent multi-site editing patterns
/// - Each event applies to consecutive target sites
///
/// # Arguments
///
/// * `filename` - Path to the allEventCounts file to read
///
/// # Returns
///
/// * `Some(AllEvents)` - Successfully parsed event collection
/// * `None` - If file cannot be read or parsing fails
///
/// # File Processing
/// 1. Skips header line
/// 2. For each data line, splits events by underscore
/// 3. Creates EditEvent objects with position offsets and lengths
/// 4. Accumulates counts for identical events at each site
/// 5. Converts to final TargetToEvents structure
///
/// # Examples
///
/// ```rust
/// if let Some(events) = read_all_read_counts(&"data.allEventCounts".to_string()) {
///     // Use events for simulation
/// }
/// ```
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