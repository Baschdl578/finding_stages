use std::env;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::sync::Mutex;
use std::thread;
use std::collections::{HashMap, HashSet};

extern crate time;
#[macro_use]
extern crate lazy_static;
extern crate num_cpus;
extern crate pbr;


use pbr::ProgressBar;

const L1_D_CACHE_SIZE: u64 = 32 * 1024;
const L1_I_CACHE_SIZE: u64 = 32 * 1024;
const L2_CACHE_SIZE: u64 = 256 * 1024;
const LINE_SIZE: u8 = 64;
const L1_ASSOCIATIVITY: u8 = 8;
const L2_ASSOCIATIVITY: u8 = 8;

const SLICE_FACTOR: usize = 500; //How many points to pick from heuristics

lazy_static! {
    static ref SWITCH_TIMES : Mutex<Vec<u64>> = Mutex::new(Vec::new()); //Already determined switch times
    static ref NEW_SWITCHES: Mutex<Vec<u64>> = Mutex::new(Vec::new()); //New switch times for simulation
    static ref SWITCH_COSTS : Mutex<u64> = Mutex::new(0); //Current lowest switch costs
    static ref PROGRESSBAR : Mutex<ProgressBar<std::io::Stdout>> = Mutex::new({
      let mut bar = ProgressBar::new(SLICE_FACTOR as u64);
      bar.format("╢▌▌░╟");
      bar.tick_format("▀▐▄▌");
      bar
    });
    static ref ACCESSES: Vec<(u64, u64, u8, u8)> = {  //List of all accesses, each element a tuple of (time, address, size, operation)
      let args: Vec<String> = env::args().collect();
      read_file(&args[1])
    };
    static ref REFERENCE_MAP: HashMap<u64, Vec<u64>> = build_reference_hashmap(&ACCESSES); //HashMap with each block number as key and a list of access times as value
}



fn main() {
    let start = current_time_millis();
    let args: Vec<String> = env::args().collect();
    let threads = num_cpus::get(); //For multithreaded computation
    let cpus = num_cpus::get_physical(); //For number of stages
    println!("Detected {} CPUs and {} threads", cpus, threads);

    let mut max_switches: usize = 0;

    let accesscount;
    let unoptimized_costs: u64;
    let unoptimized_misses: u64;

    accesscount = ACCESSES.len();

    println!("Parsed {} accesses", accesscount);

    let sim_start = current_time_millis();

    unoptimized_misses = simulate_simple(&ACCESSES);

    let simulation_length = current_time_millis() - sim_start;
    println!(
        "Simple simulation took {} ms, got {} L2 misses",
        simulation_length,
        unoptimized_misses
    );
    unoptimized_costs = unoptimized_misses * 40;

    let mut i = 1;
    while i < args.len() {
        max_switches = match args[i].parse::<usize>() {
            Ok(v) => v,
            Err(_) => {
                i += 1;
                continue;
            }
        };
        break;
    }
    if max_switches > 0 {
        max_switches -= 1;
    } else {
        max_switches = cpus - 1;
    }
    if max_switches == 0 {
        max_switches = 1;
    }

    let len = REFERENCE_MAP.len(); //Build the reference HashMap now
    println!("Got {} unique blocks", len);

    let mut misses: Vec<(u64, u64, u64)> = Vec::new(); //Contains for each number of switches the number of total, shared and dirty L3 accesses
    for point in 0..max_switches {
        //Outer loop: each run adds one switch point
        println!("Finding point no. {} of {}", point + 1, max_switches);

        {
            //Initialize progressbar and reset list of new switch points and cost
            let mut pb = PROGRESSBAR.lock().unwrap();
            *pb = ProgressBar::new(SLICE_FACTOR as u64);
            pb.tick_format("▀▐▄▌");
            pb.format("╢▌▌░╟");
            *SWITCH_COSTS.lock().unwrap() = u64::max_value();
            NEW_SWITCHES.lock().unwrap().clear();
        }
        let mut guards = Vec::with_capacity(threads); //List of threads
        if args.contains(&format!("brute-force")) {
            //Simulate for EVERY possible switch time, no heuristics
            let mut switches = NEW_SWITCHES.lock().unwrap();
            for &(time, _, _, _) in &*ACCESSES {
                switches.push(time.clone());
            }
        } else {
            for t in 0..threads {
                let mut thread_slices: Vec<u64> = Vec::new();
                let mut i: u64 = 0;

                let factor: u64 = SLICE_FACTOR as u64; //Assign all slices to the available threads
                while i < factor {
                    if i as usize % threads == t {
                        thread_slices.push(i);
                    }
                    i += 1;
                }

                let mut thread_switches: Vec<u64>; //Copy the already determined switches to a thread-local list to avoid constant (re)locking
                {
                    let switches = SWITCH_TIMES.lock().unwrap();
                    thread_switches = Vec::with_capacity(switches.len());
                    let mut i = 0;
                    while i < switches.len() {
                        thread_switches.push(switches.get(i).unwrap().clone());
                        i += 1;
                    }
                }

                guards.push(thread::spawn(
                    move || for slice in thread_slices {
                        let mut index: usize = slice as usize * accesscount / SLICE_FACTOR; //Index for inserted switch; starts as first index in the slice
                        let end: usize;
                        if slice == SLICE_FACTOR as u64 - 1 {
                            end = accesscount - 1; //Last slice, end is highest access index
                        } else {
                            end = (slice as usize + 1) * accesscount / SLICE_FACTOR; //End is index of beginning of next slice
                        }

                        let mut slice_overlaps: Vec<(u64, u64)> =
                            Vec::with_capacity(accesscount / SLICE_FACTOR + 1); //The list of all overlaps in this slice
                        let mut segments: Vec<HashMap<u64, u64>> = Vec::new(); //List of HashMaps with block numbers as key and referece count as value
                        let mut overlaps = 0;
                        let mut last_i = usize::max_value();
                        let mut last_index = usize::max_value() - 1;

                        while index < end {
                            let (time, _, _, _) = ACCESSES[index];
                            let mut temp_switches: Vec<u64> = Vec::new(); //Build list of found switch times, add the current one and sort
                            for time2 in &thread_switches {
                                temp_switches.push(time2.clone());
                            }
                            if temp_switches.contains(&time) { //Skip if two switches at the same time
                                index += 1;
                                continue;
                            }
                            temp_switches.push(time);
                            temp_switches.sort_unstable();

                            let mut i: usize = 0;
                            while i < temp_switches.len() { //Find current switch time in local list of switches
                                if temp_switches[i] == time {
                                    break;
                                }
                                i += 1;
                            }
                            if last_i != i || last_index + 1 != index { //Rebuild segments and overlap from scratch on first run or if we skipped
                                segments = block_sets(&temp_switches);
                                overlaps = segment_overlap(&segments);
                                last_i = i;
                            } else { //Calculate next overlap incrementally
                                let (_, addr, size, _) = ACCESSES[index];
                                let block1 = block_nr(addr);
                                let block2 = block_nr(addr + size as u64 - 1); //block2 equals block1 if we only referenced one block

                                //From here on, Segment A (A) refers to the segment before the moved switch.
                                //Segment B (B) is the segment after the switch
                                //<Segment A>|<SWITCH>|<Segment B>

                                //The following section inserts block1 and block2 into Segment A (if neccessary)
                                let mut insert1 = true;
                                let mut insert2 = true;
                                if segments[i].contains_key(&block1) { //Segment A already has block1, increase reference count
                                    insert1 = false;
                                    if let Some(x) = segments[i].get_mut(&block1) {
                                        *x += 1;
                                    }
                                } else {
                                    segments[i].insert(block1, 1); //Insert block1
                                }
                                if block1 != block2 { //Same procedure for second block (if neccessary)
                                    if segments[i].contains_key(&block2) {
                                        insert2 = false;
                                        if let Some(x) = segments[i].get_mut(&block2) {
                                            *x += 1;
                                        }
                                    } else {
                                        segments[i].insert(block2, 1);
                                    }
                                } else {
                                    insert2 = false;
                                }

                                //The following section checks for both blocks in Segment B, decreases their reference count or deletes them
                                let mut duplicate_in_next1 = false;
                                let mut duplicate_in_next2 = false;
                                if segments[i + 1].contains_key(&block1) {
                                    let appearances = segments[i + 1].get(&block1).unwrap().clone();
                                    if appearances > 1 {
                                        if let Some(x) = segments[i + 1].get_mut(&block1) {
                                            *x -= 1;
                                        }
                                        duplicate_in_next1 = true;
                                    } else {
                                        segments[i + 1].remove(&block1);
                                    }
                                }
                                if block1 != block2 {
                                    if segments[i + 1].contains_key(&block2) {
                                        let appearances =
                                            segments[i + 1].get(&block2).unwrap().clone();
                                        if appearances > 1 {
                                            if let Some(x) = segments[i + 1].get_mut(&block2) {
                                                *x -= 1;
                                            }
                                            duplicate_in_next2 = true;
                                        } else {
                                            segments[i + 1].remove(&block2);
                                        }
                                    }
                                }


                                if (insert1 && duplicate_in_next1) ||
                                    (!insert1 && !duplicate_in_next1) ||
                                    (block1 != block2 &&
                                         ((insert2 && duplicate_in_next2) ||
                                              (!insert2 && !duplicate_in_next2)))
                                { //Here, we will change the overlap count
                                    //This is only necessary if we inserted or removed a block in Segments A or B


                                    //The following section counts the other segments blocks 1 or 2 appear in
                                    let mut block1_other_appearances = 0;
                                    let mut block2_other_appearances = 0;
                                    let mut s = 0;
                                    while s < segments.len() {
                                        if s == i || s == i + 1 {
                                            s += 1;
                                            continue;
                                        }

                                        if segments[s].contains_key(&block1) {
                                            block1_other_appearances +=
                                                *segments[s].get(&block1).unwrap();
                                        }
                                        if block1 != block2 {
                                            if segments[s].contains_key(&block2) {
                                                block2_other_appearances +=
                                                    *segments[s].get(&block2).unwrap();
                                            }
                                        }
                                        s += 1;
                                    }


                                    //This changes the overlap count from the previous one
                                    if insert1 && duplicate_in_next1 {
                                        overlaps += block1_other_appearances + 1;
                                    } else {
                                        if !insert1 && !duplicate_in_next1 &&
                                            overlaps > block1_other_appearances
                                        {
                                            overlaps -= block1_other_appearances + 1;
                                        } else if overlaps <= block1_other_appearances {
                                            overlaps = 0;
                                        }
                                    }

                                    if block1 != block2 {
                                        if insert2 && duplicate_in_next2 {
                                            overlaps += block2_other_appearances + 1;
                                        } else {
                                            if !insert2 && !duplicate_in_next2 &&
                                                overlaps > block2_other_appearances
                                            {
                                                overlaps -= block2_other_appearances + 1;
                                            } else if overlaps <= block2_other_appearances {
                                                overlaps = 0;
                                            }
                                        }
                                    }
                                }
                            }
                            last_index = index;
                            slice_overlaps.push((overlaps, time)); //Insert new overlap into list
                            index += 1;

                            if index % 10000 == 0 {
                                PROGRESSBAR.lock().unwrap().tick();
                            }
                        }
                        //All the overlaps are calculated now

                        let args: Vec<String> = env::args().collect();
                        let short_run = !args.contains(&format!("long"));
                        //If parameter "long" is specified, we simulate for all the smallest
                        //overlaps of the section (if multiple points have that same overlap)
                        //Note: This will take A LONG TIME
                        //If that parameter is not specified, we just simulate the first occurance of that overlap

                        //The following section looks for points with the lowest overlap and puts
                        //them into a list for simulation
                        let mut smallest_overlap: u64 = u64::max_value();
                        let mut best_time = 0;
                        let mut best_times: Vec<u64> = Vec::new();
                        for &(overlap, time2) in &slice_overlaps {
                            if smallest_overlap > overlap {
                                smallest_overlap = overlap;
                                if short_run {
                                    best_time = time2;
                                } else {
                                    best_times.clear();
                                }
                            }
                            if !short_run {
                                if smallest_overlap == overlap {
                                    best_times.push(time2);
                                }
                            }
                        }
                        slice_overlaps.clear();
                        if short_run {
                            NEW_SWITCHES.lock().unwrap().push(best_time);
                        } else {
                            let mut new_switches = NEW_SWITCHES.lock().unwrap();
                            for time2 in best_times {
                                new_switches.push(time2);
                            }
                        }
                        PROGRESSBAR.lock().unwrap().inc();
                    },
                ));
            }
        }

        for g in guards {
            let _ = g.join();
        }
        println!("");

        println!("Simulating each point:");
        let mut guards = Vec::with_capacity(threads);
        {
            //Reset progress bar
            let mut pb = PROGRESSBAR.lock().unwrap();
            *pb = ProgressBar::new(NEW_SWITCHES.lock().unwrap().len() as u64);
            pb.format("╢▌▌░╟");
            pb.tick_format("▀▐▄▌");
        }

        for t in 0..threads {
            let mut thread_times = Vec::new();
            {
                //Assign new switch times to threads
                let times = NEW_SWITCHES.lock().unwrap();
                let mut i = 0;
                while i < times.len() {
                    if i % threads == t {
                        thread_times.push(times[i].clone());
                    }
                    i += 1;
                }
            }

            let mut thread_switches: Vec<u64> = Vec::new();
            {
                //Copy determined switches for each tread
                let timevec = SWITCH_TIMES.lock().unwrap();
                let mut i = 0;
                while i < timevec.len() {
                    thread_switches.push(timevec[i].clone());
                    i += 1;
                }
            }
            guards.push(thread::spawn(move || {
                let mut thread_costs = u64::max_value();
                let mut switch_time = 0;
                for time in thread_times {
                    let mut temp_switches: Vec<u64> = Vec::new(); //Build list of switches with the new switch inserted and sort
                    for time2 in &thread_switches {
                        temp_switches.push(time2.clone());
                    }
                    temp_switches.push(time);
                    temp_switches.sort_unstable();

                    let (sim_misses, shared, transfers) = simulate_with_switches(temp_switches);

                    let costs = calc_costs((sim_misses, shared, transfers));

                    if costs < thread_costs { //Check if new simulation produces less cycles for L3 accesses and save costs
                        switch_time = time;
                        thread_costs = costs;
                    }
                    PROGRESSBAR.lock().unwrap().inc();

                }

                {

                    let mut overall_costs = SWITCH_COSTS.lock().unwrap(); //Put the calculated cheapest point into the global list
                    if *overall_costs > thread_costs {
                        let mut times = SWITCH_TIMES.lock().unwrap();
                        if times.len() <= point {
                            times.push(switch_time);
                        } else {
                            (*times)[point] = switch_time;
                        }
                        *overall_costs = thread_costs;
                    }
                }
            }));
        }

        NEW_SWITCHES.lock().unwrap().clear();
        for g in guards {
            let _ = g.join();
        }
        println!("");

        //This will do one final simulation to determine the costs and print detailed statistics
        let mut times: Vec<u64> = Vec::new();
        {
            let timevec = SWITCH_TIMES.lock().unwrap();
            let mut i: usize = 0;
            while i < timevec.len() {
                times.push(timevec[i].clone());
                i += 1;
            }
        }
        times.sort_unstable();
        let optimized_misses: (u64, u64, u64) = simulate_with_switches(times);
        let (miss, _, transf) = optimized_misses;
        let costs = calc_costs(optimized_misses);
        println!(
            "New switch saved {} cycles ({}%), {} misses and transferred {} dirty cache lines",
            unoptimized_costs as i64 - costs as i64,
            (100.0 - (costs as f64 * 100.0 / unoptimized_costs as f64)) as f32,
            unoptimized_misses as i64 - miss as i64,
            transf
        );
        misses.push(optimized_misses);
    }

    let mut optimized_costs: u64 = u64::max_value();
    let mut optimized_misses: u64 = u64::max_value();
    let mut i = 0;
    let mut switchcount = 0;
    let mut transfers = 0;
    for (miss, shared, transf) in misses {
        //Find out the best number of switches
        i += 1;

        let costs = calc_costs((miss, shared, transf));
        if costs < optimized_costs {
            switchcount = i;
            optimized_costs = costs;
            optimized_misses = miss;
            transfers = transf;
        }
    }
    let mut final_switches: Vec<u64> = Vec::new(); //List of final switches
    i = 0;
    {
        let timevec = SWITCH_TIMES.lock().unwrap();
        while i < switchcount && i < timevec.len() {
            final_switches.push(timevec.get(i).unwrap().clone());
            i += 1;
        }
    }
    final_switches.sort_unstable();

    //This is just the output
    println!("Switching {} times is optimal", final_switches.len());
    if !final_switches.is_empty() {
        let mut i: usize = 0;
        println!("Switch times: ");
        while i < final_switches.len() {
            let time = final_switches.get(i).unwrap().clone();
            let &(time2, _, _, _) = &ACCESSES[accesscount - 1];
            println!(
                "At {}ns ({}% of the total runtime)",
                time,
                (time * 100) as f32 / time2 as f32
            );

            let mut j: usize = 0;
            while j < accesscount {
                let &(time2, _, _, _) = &ACCESSES[j];
                if time2 == time {
                    let (file, function, line) = code_pos(&args[1], time);
                    println!(
                        "\tSwitch after position number {} of {}\n\tFile: {}\n\tFunction: {}\n\tLine: {}",
                        j,
                        accesscount,
                        file,
                        function,
                        line
                    );
                    break;
                }
                j += 1;
            }
            i += 1;
        }
        println!(
            "Savings: {}%, {} less misses, {} less cycles",
            (100.0 - (optimized_costs as f64 * 100.0 / unoptimized_costs as f64)) as f32,
            unoptimized_misses as i64 - optimized_misses as i64,
            unoptimized_costs as i64 - optimized_costs as i64
        );
        println!("This will transfer {} dirty cache lines", transfers);
    }

    let end = current_time_millis();
    println!("Took: {}mins", ((end - start) / 1000 / 60));
}

/// Builds a HashMap with the blocks as key and thre reference count as value for each segment
/// resulting from the given switch configuration
///
/// WARNING: The given switch point configuration must be sorted!
fn block_sets(switches: &Vec<u64>) -> Vec<HashMap<u64, u64>> {
    let mut blocksets: Vec<HashMap<u64, u64>> = Vec::with_capacity(switches.len());

    let mut i: usize = 0;
    while i <= switches.len() {
        blocksets.push(HashMap::new());
        i += 1;
    }
    let mut i: usize = 0;
    for &(time, addr, size, _) in &*ACCESSES {
        let block1 = block_nr(addr);
        let block2 = block_nr(addr + size as u64 - 1);

        if !(blocksets[i].contains_key(&block1)) {
            (blocksets[i]).insert(block1, 1); //Insert new block
        } else {
            if let Some(x) = blocksets[i].get_mut(&block1) {
                *x += 1; //Increase reference counter
            }
        }
        if block1 != block2 {
            if !(blocksets[i].contains_key(&block2)) {
                (blocksets[i]).insert(block2, 1);
            } else {
                if let Some(x) = blocksets[i].get_mut(&block2) {
                    *x += 1;
                }
            }
        }
        if i < switches.len() && switches[i] == time {
            i += 1; //Go to next segment
        }
    }
    blocksets
}

/// Calculates the working set overlap of the given list of segments
fn segment_overlap(segments: &Vec<HashMap<u64, u64>>) -> u64 {
    let mut overlap: u64 = 0;
    let mut i: usize = 0;
    while i < segments.len() {
        let segment = &segments[i];
        for block in segment.keys() {
            let mut j: usize = i + 1;
            while j < segments.len() {
                if segments[j].contains_key(&block) {
                    overlap += 1;
                }
                j += 1;
            }
        }
        i += 1;
    }
    overlap
}

/// Gets the file, function and line number of a reference identified by its timestamp
fn code_pos(path: &str, time: u64) -> (String, String, u64) {

    let file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);

    let mut line: String;
    loop {
        line = String::new();
        match reader.read_line(&mut line) {
            Err(_) => {
                return (String::new(), String::new(), 0);
            }
            Ok(_) => {
                if line.starts_with("Operation") {
                    continue;
                }
                if line == "" {
                    return (String::new(), String::new(), 0);
                }
                let mut string = String::new();
                string.push_str(&line);
                let mut strings: Vec<String> = Vec::new();
                for s in string.split(';') {
                    strings.push(format!("{}", s));
                }
                if strings.len() == 8 {
                    let parse_time = match strings[7].trim().parse::<u64>() {
                        Ok(v) => v,
                        Err(_) => {
                            println!("Could not parse to u64: {}", strings[7]);
                            panic!("Parse error");
                        }
                    };
                    if parse_time != time {
                        continue;
                    }

                    let parse_line = match strings[4].trim().parse::<u64>() {
                        Ok(v) => v,
                        Err(_) => {
                            println!("Could not parse to u64: {}", strings[4]);
                            panic!("Parse error");
                        }
                    };

                    let my_file = format!("{}", strings[2].trim());
                    let my_function = format!("{}", strings[3].trim());

                    return (my_file, my_function, parse_line);

                } else {
                    continue;
                }
            }
        }

    }
}

/// This will parse the input file that contains all references into a list where each element
/// is a tuple containing (Timestamp, Address, Size, Operation)
fn read_file(path: &str) -> Vec<(u64, u64, u8, u8)> {
    println!("Parsing file: {}", path);

    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines = reader.lines().count();
    let mut pb = ProgressBar::new(lines as u64);
    pb.format("╢▌▌░╟");
    pb.tick_format("▀▐▄▌");

    let mut filter: HashSet<(String, String)> = HashSet::new();
    let args: Vec<String> = env::args().collect();
    let mut i = 0;
    while i < args.len() {
        let mut arg = format!("{}", &args[i]);
        if arg.starts_with("filter=") {
            //This will filter out all the references caused by
            // funtions in a seperate filter file that has the same format as the main input file
            // For that purpose we build a HashSet that contains a tuple of (File, Function)
            let offset = arg.find('=').unwrap_or(arg.len()) + 1;
            let _ = arg.drain(..offset);
            let file = File::open(arg).unwrap();
            let mut reader = BufReader::new(file);
            let mut line: String;
            loop {
                line = String::new();
                match reader.read_line(&mut line) {
                    Err(_) => {
                        break;
                    }
                    Ok(_) => {
                        if line.starts_with("Operation") {
                            continue;
                        }
                        if line == "" {
                            break;
                        }
                        let mut string = String::new();
                        string.push_str(&line);
                        let mut strings: Vec<String> = Vec::new();
                        for s in string.split(';') {
                            strings.push(format!("{}", s));
                        }
                        if strings.len() == 8 {
                            filter.insert((
                                format!("{}", strings[2].trim()),
                                format!("{}", strings[3].trim()),
                            ));
                        } else {
                            continue;
                        }
                    }
                }
            }
        }
        i += 1;
    }


    let mut accesses: Vec<(u64, u64, u8, u8)> = Vec::with_capacity(lines);
    let file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);
    let mut line: String;
    let mut i: u64 = 0;
    loop {
        i += 1;
        line = String::new();
        match reader.read_line(&mut line) {
            Err(_) => {
                break;
            }
            Ok(_) => {
                if line.starts_with("Operation") {
                    continue;
                }
                if line == "" {
                    break;
                }
                let mut string = String::new();
                string.push_str(&line);
                let mut strings: Vec<String> = Vec::new();
                for s in string.split(';') {
                    strings.push(format!("{}", s));
                }
                if strings.len() == 8 {
                    let parse_time = match strings[7].trim().parse::<u64>() {
                        Ok(v) => v,
                        Err(_) => {
                            println!("Could not parse to u64: {}", strings[7]);
                            panic!("Parse error");
                        }
                    };
                    let parse_size = match strings[6].trim().parse::<u8>() {
                        Ok(v) => v,
                        Err(_) => {
                            println!("Could not parse to u8: {}", strings[6]);
                            panic!("Parse error");
                        }
                    };
                    let parse_addr = match strings[5].trim().parse::<u64>() {
                        Ok(v) => v,
                        Err(_) => {
                            println!("Could not parse to u64: {}", strings[5]);
                            panic!("Parse error");
                        }
                    };

                    let op_in = &strings[0].trim();

                    let parse_op: u8 = match op_in {
                        &"Ir" => 0,
                        &"Dr" => 1,
                        &"Dw" => 2,
                        _ => 3,
                    };
                    let file = format!("{}", strings[2]);
                    let function = format!("{}", strings[3]);
                    if !filter.contains(&(file, function)) {
                        let access = (parse_time, parse_addr, parse_size, parse_op);
                        accesses.push(access);
                    }
                } else {
                    continue;
                }
            }
        }
        if i % 50000 == 0 {
            pb.tick();
        }
        if i % 500000 == 0 {
            pb.set(i);
        }
    }
    pb.set(lines as u64);
    println!("");
    accesses
}


/// Simple simulation that calculates the number of L2 misses without switching processors
fn simulate_simple(accesses: &Vec<(u64, u64, u8, u8)>) -> u64 {
    let mut current_state = Box::new(State::new());
    let mut misses: u64 = 0;

    for &(time, addr, size, op) in accesses {
        current_state.time = time;
        if current_state.reference(addr, size, op) {
            misses += 1;
        }
    }

    let mut previous_misses: u64 = 0;
    while previous_misses != misses {
        previous_misses = misses;
        misses = 0;
        for &(time, addr, size, op) in accesses {
            current_state.time = time;
            if current_state.reference(addr, size, op) {
                misses += 1;
            }

        }
    }
    misses
}


/// Builds a HashMap of all referenced blocks and their access times
fn build_reference_hashmap(accesses: &Vec<(u64, u64, u8, u8)>) -> HashMap<u64, Vec<u64>> {
    let mut out: HashMap<u64, Vec<u64>> = HashMap::new();
    for &(time, addr, size, _) in accesses {
        let block1 = block_nr(addr);
        let block2 = block_nr(addr + size as u64 - 1);

        if out.contains_key(&block1) {
            out.get_mut(&block1).unwrap().push(time);
        } else {
            let mut times: Vec<u64> = Vec::new();
            times.push(time);
            times.sort_unstable();
            out.insert(block1, times);
        }
        if block1 != block2 {
            if out.contains_key(&block2) {
                out.get_mut(&block2).unwrap().push(time);
            } else {
                let mut times: Vec<u64> = Vec::new();
                times.push(time);
                times.sort();
                out.insert(block2, times);
            }
        }
    }
    return out;
}

/// Simulate with switching processors
/// Returns the following tuple: (Total Misses, Shared Accesses, Dirty Accesses)
/// Important: For the shared accessses, we assume a system running multiple queries, so every block
/// that is refenced in multiple segments will count as a shared access (worst case)
fn simulate_with_switches(switch_times: Vec<u64>) -> (u64, u64, u64) {
    let mut switch_states: Vec<Box<State>> = Vec::with_capacity(switch_times.len() + 1);
    for _ in 0..(switch_times.len() + 1) {
        switch_states.push(Box::new(State::new()));
    }
    let mut current_state = Box::new(State::new());

    let mut count: usize;
    let mut switch_count: usize;
    let mut shared_accesses: u64 = 0;
    let mut access: &(u64, u64, u8, u8);
    let accesscount = ACCESSES.len();

    let mut misses: u64 = 0;
    let mut previous_misses: u64 = 1;
    while previous_misses != misses {
        //Run this until we reach a steady state (at lease 3 times)
        previous_misses = misses;
        misses = 0;
        count = 0;
        switch_count = 0;
        shared_accesses = 0;
        while count < accesscount {
            access = ACCESSES.get(count).unwrap();
            let &(time, addr, size, op) = access;
            let miss = current_state.reference(addr, size, op);
            if miss {

                misses += 1;
                let block1 = block_nr(addr);
                let block2 = block_nr(addr + size as u64 - 1);

                //Check if shared access
                let mut found = false;
                for ref_time in REFERENCE_MAP.get(&block1).unwrap() {
                    if switch_count > 0 && switch_count < switch_times.len() {
                        if ref_time <= switch_times.get(switch_count - 1).unwrap() ||
                            ref_time > switch_times.get(switch_count).unwrap()
                        {
                            found = true;
                            shared_accesses += 1;
                            break;
                        }
                    } else {
                        if switch_count > 0 {
                            if ref_time <= switch_times.get(switch_count - 1).unwrap() {
                                found = true;
                                shared_accesses += 1;
                                break;
                            }
                        } else {
                            if switch_count < switch_times.len() &&
                                ref_time > switch_times.get(switch_count).unwrap()
                            {
                                found = true;
                                shared_accesses += 1;
                                break;
                            }
                        }
                    }
                }
                if !found && block1 != block2 {
                    for ref_time in REFERENCE_MAP.get(&block2).unwrap() {
                        if switch_count > 0 && switch_count < switch_times.len() {
                            if ref_time < switch_times.get(switch_count - 1).unwrap() ||
                                ref_time > switch_times.get(switch_count).unwrap()
                            {
                                shared_accesses += 1;
                                break;
                            }
                        } else {
                            if switch_count > 0 {
                                if ref_time < switch_times.get(switch_count - 1).unwrap() {
                                    shared_accesses += 1;
                                    break;
                                }
                            } else {
                                if switch_count < switch_times.len() &&
                                    ref_time > switch_times.get(switch_count).unwrap()
                                {
                                    shared_accesses += 1;
                                    break;
                                }
                            }
                        }
                    }
                }

            }
            current_state.time = time;
            if switch_times.contains(&time) && switch_count + 1 < switch_states.len() {
                switch_states[switch_count] = current_state.clone();
                current_state = switch_states[switch_count + 1].clone();;

                switch_count += 1;
            }
            count += 1;
            if count == accesscount {
                switch_states[switch_times.len()] = current_state.clone();
            }
            if count % (accesscount / 10) == 0 {
                PROGRESSBAR.lock().unwrap().tick();
            }
        }
    }
    let mut dirty_count: u64 = 0;
    // CHeck for dirty accesses
    let mut i: usize = 0;
    while i < switch_states.len() - 1 {
        let state: Box<State> = switch_states.get(i).unwrap().clone();
        let mut j: usize = 0;
        while j < state.dirty_blocks.len() {
            let block = state.dirty_blocks.get(j).unwrap();

            for ref_time in REFERENCE_MAP.get(&block).unwrap() {
                if ref_time > &state.time {
                    dirty_count += 1;
                    break;
                }
            }
            j += 1;
        }
        i += 1;
    }

    (misses, shared_accesses, dirty_count)
}

/// Calculate the block nr of a given address
fn block_nr(address: u64) -> u64 {
    return address / (LINE_SIZE as u64);
}

/// Calculate the L2 set number of a given block
fn get_l2_set(block: u64) -> u16 {
    return (block % (L2_ASSOCIATIVITY as u64)) as u16;
}

/// Calculate the L1 set number of a given block
fn get_l1_set(block: u64) -> u16 {
    return (block % (L1_ASSOCIATIVITY as u64)) as u16;
}

/// Calculate the current time in ms
fn current_time_millis() -> u64 {
    let time = time::get_time();
    let now: u64 = (time.sec as u64 * 1000) + (time.nsec as u64 / 1000 / 1000);
    now
}

/// Calculate the cycles spent on L3 accesses
fn calc_costs(values: (u64, u64, u64)) -> u64 {
    let (sim_misses, shared, transfers) = values;
    let temp_misses;
    let temp_shared;
    if sim_misses > shared {
        temp_misses = sim_misses - shared;
    } else {
        temp_misses = 0; //This should never happen
    }
    if shared > transfers {
        temp_shared = shared - transfers;
    } else {
        temp_shared = 0; //This should never happen
    }
    let costs = temp_misses * 40 + temp_shared * 65 + transfers * 75;
    costs
}

/// Cache state
struct State {
    pub time: u64, //Timestamp
    pub l1_i_cache:
        [L1Set; ((L1_I_CACHE_SIZE / LINE_SIZE as u64) / L1_ASSOCIATIVITY as u64) as usize],
    pub l1_d_cache:
        [L1Set; ((L1_D_CACHE_SIZE / LINE_SIZE as u64) / L1_ASSOCIATIVITY as u64) as usize],
    pub l2_cache: [L2Set; ((L2_CACHE_SIZE / LINE_SIZE as u64) / L2_ASSOCIATIVITY as u64) as usize],
    pub dirty_blocks: Vec<u64>, //List of dirty blocks
}

struct L1Set {
    pub addresses: [u64; L1_ASSOCIATIVITY as usize],
}

struct L2Set {
    pub addresses: [u64; L2_ASSOCIATIVITY as usize],
}

impl State {
    fn new() -> State {
        State {
            time: 0,
            l1_i_cache: [L1Set::new();
                ((L1_I_CACHE_SIZE / LINE_SIZE as u64) / L1_ASSOCIATIVITY as u64) as usize],
            l1_d_cache: [L1Set::new();
                ((L1_D_CACHE_SIZE / LINE_SIZE as u64) / L1_ASSOCIATIVITY as u64) as usize],
            l2_cache: [L2Set::new();
                ((L2_CACHE_SIZE / LINE_SIZE as u64) / L2_ASSOCIATIVITY as u64) as usize],
            dirty_blocks: Vec::new(),
        }
    }

    /// This updates the cache state when a ressource is referenced
    fn reference(&mut self, addr: u64, size: u8, op: u8) -> bool {
        let block1 = block_nr(addr);
        let block2 = block_nr(addr + size as u64 - 1);
        let mut out: bool;

        if op == 2 {
            //Data write, add blocks to dirty list
            if !self.dirty_blocks.contains(&(block1 as u64)) {
                self.dirty_blocks.push(block1 as u64);
            }
            if block2 != block1 {
                if !self.dirty_blocks.contains(&(block2 as u64)) {
                    self.dirty_blocks.push(block2 as u64);
                }
            }
        }

        if op == 0 {
            //Instruction read, reference L1 Instruction cache
            let (miss, evict) = self.l1_i_cache[get_l1_set(block1) as usize].reference(block1);
            out = miss;
            if miss && evict != 0 {
                if self.dirty_blocks.contains(&evict) {
                    let mut new_dirty_blocks: Vec<u64> =
                        Vec::with_capacity(self.dirty_blocks.len() - 1);
                    let mut i: usize = 0;
                    while i < self.dirty_blocks.len() {
                        let block = self.dirty_blocks.get(i).unwrap().clone();
                        if block != evict {
                            new_dirty_blocks.push(block);
                        }
                        i += 1;
                    }
                    self.dirty_blocks = new_dirty_blocks;
                }
            }
            if block1 + 1 == block2 {
                if out {
                    self.l1_i_cache[get_l1_set(block2) as usize].reference(block2);
                } else {
                    let (miss, evict) =
                        self.l1_i_cache[get_l1_set(block2) as usize].reference(block2);
                    out = miss;
                    if miss && evict != 0 {
                        if self.dirty_blocks.contains(&evict) {
                            let mut new_dirty_blocks: Vec<u64> =
                                Vec::with_capacity(self.dirty_blocks.len() - 1);
                            let mut i: usize = 0;
                            while i < self.dirty_blocks.len() {
                                let block = self.dirty_blocks.get(i).unwrap().clone();
                                if block != evict {
                                    new_dirty_blocks.push(block);
                                }
                                i += 1;
                            }
                            self.dirty_blocks = new_dirty_blocks;
                        }
                    }
                }
            } else {
                if block1 != block2 {
                    panic!("Cache Access too large")
                }
            }
        } else {
            //Data referece, reference L1 Data cache first
            let (miss, evict) = self.l1_d_cache[get_l1_set(block1) as usize].reference(block1);
            out = miss;
            if miss && evict != 0 {
                if self.dirty_blocks.contains(&evict) {
                    let mut new_dirty_blocks: Vec<u64> =
                        Vec::with_capacity(self.dirty_blocks.len() - 1);
                    let mut i: usize = 0;
                    while i < self.dirty_blocks.len() {
                        let block = self.dirty_blocks.get(i).unwrap().clone();
                        if block != evict {
                            new_dirty_blocks.push(block);
                        }
                        i += 1;
                    }
                    self.dirty_blocks = new_dirty_blocks;
                }
            }
            if block1 + 1 == block2 {
                if out {
                    self.l1_d_cache[get_l1_set(block2) as usize].reference(block2);
                } else {
                    let (miss, evict) =
                        self.l1_d_cache[get_l1_set(block2) as usize].reference(block2);
                    out = miss;
                    if miss && evict != 0 {
                        if self.dirty_blocks.contains(&evict) {
                            let mut new_dirty_blocks: Vec<u64> =
                                Vec::with_capacity(self.dirty_blocks.len() - 1);
                            let mut i: usize = 0;
                            while i < self.dirty_blocks.len() {
                                let block = self.dirty_blocks.get(i).unwrap().clone();
                                if block != evict {
                                    new_dirty_blocks.push(block);
                                }
                                i += 1;
                            }
                            self.dirty_blocks = new_dirty_blocks;
                        }
                    }
                }
            } else {
                if block1 != block2 {
                    panic!("Cache Access too large")
                }
            }
        }

        if out {
            //Miss in L1 cache, reference L2 cache now
            let (miss, evict) = self.l2_cache[get_l2_set(block1) as usize].reference(block1);
            out = miss;
            if miss && evict != 0 {
                if self.dirty_blocks.contains(&evict) {
                    let mut new_dirty_blocks: Vec<u64> =
                        Vec::with_capacity(self.dirty_blocks.len() - 1);
                    let mut i: usize = 0;
                    while i < self.dirty_blocks.len() {
                        let block = self.dirty_blocks.get(i).unwrap().clone();
                        if block != evict {
                            new_dirty_blocks.push(block);
                        }
                        i += 1;
                    }
                    self.dirty_blocks = new_dirty_blocks;
                }
            }
            if block1 + 1 == block2 {
                if out {
                    self.l2_cache[get_l2_set(block2) as usize].reference(block2);
                } else {
                    let (miss, evict) =
                        self.l2_cache[get_l2_set(block2) as usize].reference(block2);
                    out = miss;
                    if miss && evict != 0 {
                        if self.dirty_blocks.contains(&evict) {
                            let mut new_dirty_blocks: Vec<u64> =
                                Vec::with_capacity(self.dirty_blocks.len() - 1);
                            let mut i: usize = 0;
                            while i < self.dirty_blocks.len() {
                                let block = self.dirty_blocks.get(i).unwrap().clone();
                                if block != evict {
                                    new_dirty_blocks.push(block);
                                }
                                i += 1;
                            }
                            self.dirty_blocks = new_dirty_blocks;
                        }
                    }
                }
            } else {
                if block1 != block2 {
                    panic!("Cache Access too large")
                }
            }
        }
        return out;
    }
}

impl L1Set {
    fn new() -> L1Set {
        L1Set { addresses: [0 as u64; L1_ASSOCIATIVITY as usize] }
    }

    /// This will simulate referencing a block
    /// The set consists of a list of blocks.
    /// If the referenced block is already present in the set, we move it to the first spot in the list.
    /// All blocks that were previously before the referenced one are moved a spot to the back.
    ///
    /// If the block is not referneced yet, all blocks are moved to the back by one and the new one is inserted in the front.
    /// If the Set is full, the last block is evicted.
    /// Because we always move referenced blocks to the front, this implements an LRU eviction policy
    fn reference(&mut self, address: u64) -> (bool, u64) {
        let mut pos: usize = (L1_ASSOCIATIVITY - 1) as usize;
        let mut i: usize = 0;
        let mut found = false;
        while i < L1_ASSOCIATIVITY as usize {
            if self.addresses[i] == address {
                pos = i;
                found = true;
                break;
            }
            if self.addresses[i] == 0 {
                pos = i;
                break;
            }
            i += 1;
        }
        let mut evict: u64 = 0;
        if !found && pos == (L1_ASSOCIATIVITY - 1) as usize {
            evict = self.addresses[pos].clone();
        }

        while pos > 0 {
            self.addresses[pos] = self.addresses[pos - 1];
            pos -= 1;
        }
        self.addresses[0] = address;
        return (!found, evict);
    }
}

impl Clone for L1Set {
    fn clone(&self) -> L1Set {
        let mut out = L1Set::new();
        let mut i: usize = 0;
        while i < L1_ASSOCIATIVITY as usize {
            out.addresses[i] = self.addresses[i].clone();
            i += 1;
        }
        out
    }
}

impl Copy for L1Set {}

impl L2Set {
    fn new() -> L2Set {
        L2Set { addresses: [0 as u64; L2_ASSOCIATIVITY as usize] }
    }

    /// This will simulate referencing a block
    /// The set consists of a list of blocks.
    /// If the referenced block is already present in the set, we move it to the first spot in the list.
    /// All blocks that were previously before the referenced one are moved a spot to the back.
    ///
    /// If the block is not referneced yet, all blocks are moved to the back by one and the new one is inserted in the front.
    /// If the Set is full, the last block is evicted.
    /// Because we always move referenced blocks to the front, this implements an LRU eviction policy
    fn reference(&mut self, address: u64) -> (bool, u64) {
        let mut pos: usize = (L2_ASSOCIATIVITY - 1) as usize;
        let mut i: usize = 0;
        let mut found = false;
        while i < L2_ASSOCIATIVITY as usize {
            if self.addresses[i] == address {
                pos = i;
                found = true;
                break;
            }
            if self.addresses[i] == 0 {
                pos = i;
                break;
            }
            i += 1;
        }
        let mut evict: u64 = 0;
        if !found && pos == (L1_ASSOCIATIVITY - 1) as usize {
            evict = self.addresses[pos].clone();
        }

        while pos > 0 {
            self.addresses[pos] = self.addresses[pos - 1];
            pos -= 1;
        }
        self.addresses[0] = address;
        return (!found, evict);
    }
}

impl Clone for L2Set {
    fn clone(&self) -> L2Set {
        let mut out = L2Set::new();
        let mut i: usize = 0;
        while i < L2_ASSOCIATIVITY as usize {
            out.addresses[i] = self.addresses[i].clone();
            i += 1;
        }
        out
    }
}

impl Copy for L2Set {}

impl Clone for State {
    fn clone(&self) -> State {
        let mut out = State::new();
        out.time = self.time;
        let mut i: usize = 0;
        while i < self.l1_i_cache.len() {
            out.l1_i_cache[i] = self.l1_i_cache[i].clone();
            i += 1;
        }
        let mut i: usize = 0;
        while i < self.l1_d_cache.len() {
            out.l1_d_cache[i] = self.l1_d_cache[i].clone();
            i += 1;
        }
        let mut i: usize = 0;
        while i < self.l2_cache.len() {
            out.l2_cache[i] = self.l2_cache[i].clone();
            i += 1;
        }
        let mut i: usize = 0;
        while i < self.dirty_blocks.len() {
            out.dirty_blocks.push(
                self.dirty_blocks.get(i).unwrap().clone(),
            );
            i += 1;
        }
        out
    }
}
