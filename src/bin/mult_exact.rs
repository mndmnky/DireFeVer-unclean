use std::time::Instant;
use clap::{Arg, App};
use std::error;
use std::path::PathBuf;
use std::fs::File;
use std::thread;
use std::thread::JoinHandle;
use std::io::{BufReader, Write};
use dfvs_solver::{digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule, cust_errors::{ImportError, GraphError}, statistics::RuleStats};
use regex::Regex;

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let digits = Regex::new(r"\d+").unwrap();
    let m = App::new("statistics")
        .arg(Arg::new("files")
             .takes_value(true)
             .multiple_values(true)
             .short('f'))
        .arg(Arg::new("csv")
             .required(true)
             .takes_value(true)
             .short('c'))
        .arg(Arg::new("time_out")
             .takes_value(true)
             .validator_regex(&digits, "only numbers are allowed")
             .short('t'))
        .get_matches();
    // Get, as input all public instances. 
    let files: Vec<PathBuf> = m.values_of("files").unwrap().map(|p| PathBuf::from(p)).collect();
    let csv: &str = m.value_of("csv").unwrap();
    let time_out: Option<u128> = m.value_of("time_out").map(|val| val.parse::<u128>().expect("input checked with regex"));
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }
    let mut threads: Vec<JoinHandle<Result<_,GraphError>>> = Vec::new();
    for (graph, name) in graphs.clone() {
        threads.push(thread::spawn(move || {
            let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);
            let (n, m, _, _) = dfvsi.graph.graph.n_m_max_min_stats();
            // Set timer if `time_out` was set 
            if let Some(time) = time_out {
                dfvsi.set_time_interrupter(time);
            }
            // initial reduction
            let priority = vec![Rule::SimpleRules, Rule::Dome, Rule::Clique, Rule::Core, Rule::SCC, Rule::LinkNode];
            let start = Instant::now();
            if !dfvsi.exhaustive_reductions(&priority){
                // 0: timeout at reductions
                return Ok((name, (n,m), (-1, -1), -1, -1, -1, -1, -1, -1))
            }
            let exRedDur = start.elapsed().as_millis();

            let (kn, km, _, _) = dfvsi.graph.graph.n_m_max_min_stats();

            if !dfvsi.compute_and_set_upper_lower(true){
                return Ok((name, (n,m), (kn as isize, km as isize), exRedDur as i128, -1, -1, -1, -1, -1))
            }
            let upper = dfvsi.upper_bound.expect("was set");
            let lower = dfvsi.lower_bound.expect("was set");
            let subgraphs: isize = dfvsi.graph.graph.split_into_connected_components_alt().len() as isize;
            // Solve 
            match dfvsi.advanced_branching_stat(&priority, true){
                Ok((time,resu)) => return Ok((name, (n, m), (kn as isize, km as isize), exRedDur as i128, upper as isize, lower as isize, subgraphs, resu.len() as isize, time as isize)),
                Err(_) => return Ok((name, (n, m), (kn as isize, km as isize), exRedDur as i128, upper as isize, lower as isize, subgraphs, -1, -1)),
            };
        }));
    }
    let mut out_file = File::create(csv)?;
    writeln!(&mut out_file, "instance, n(g), m(g), n(k), m(k), t(init_rules),\
    number_scc, lower_bound, upper_bound, exact, time")?;

    for join_handle in threads {
        match join_handle.join() {
            Ok(Ok((name, (n,m), (kn, km), exRedDur, upper, lower, sccs, exact, time))) => {
                eprintln!("joined {:?}", name);
                writeln!(&mut out_file, "{:?}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}", name, n, m, kn, km, exRedDur, sccs, lower, upper, exact, time)?;
            },
            Ok(Err(_)) => eprintln!("Some thread paniced"),
            Err(_) => eprintln!("Some thread paniced"),
        }
    }
    Ok(())
}
