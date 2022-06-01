use std::collections::HashMap;
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
            let mut results = HashMap::new();
            let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);
            let (n, m, _, _) = dfvsi.graph.graph.n_m_max_min_stats();
            // Set timer if `time_out` was set 
            if let Some(time) = time_out {
                dfvsi.set_time_interrupter(time);
            }
            // initial reduction
            let priority = vec![Rule::SimpleRules, Rule::Dome, Rule::Clique, Rule::Core, Rule::SCC];
            let start = Instant::now();
            if !dfvsi.exhaustive_reductions(&priority){
                // 0: timeout at reductions
                results.insert("\"Time-out\"".to_owned(), (0,0));
                return Ok((results, name, (n,m), (0, 0), 0, 0))
            }
            let exRedDur = start.elapsed().as_millis();

            let (kn, km, _, _) = dfvsi.graph.graph.n_m_max_min_stats();

            let start = Instant::now();
            let chr = dfvsi.clique_heuristic(&priority, true);
            if chr.is_none() {
                // print results so far if any
                // 1: timeout after reductions
                results.insert("\"Time-out\"".to_owned(), (1,1));
                return Ok((results, name, (n,m), (kn, km), exRedDur, 0))
            }
            let (lower, upper, res) = chr.expect("is not none");
            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                let duration = start.elapsed();
                results.insert("\"CliqueHeuristic/wloc\"".to_owned(), (better_res.len(), duration.as_millis()));
            } else {
                let duration = start.elapsed();
                results.insert("\"CliqueHeuristic\"".to_owned(), (res.len(), duration.as_millis()));
            }
            if dfvsi.graph.num_nodes()<50 {
                // Solve 
                // TODO: for the actual casewe need a branching that returns an intermediate solution in
                // case of a sigint.
                match dfvsi.advanced_branching_stat(&priority, true){
                    Ok((time,resu)) => results.insert("Exact".to_owned(), (resu.len(), time)),
                    Err(e) => results.insert("Error on exact".to_owned(), (1,1)),
                };
            } else if dfvsi.graph.num_nodes() < 10000 {
                // Apply lots of upper bounds 
                // * If kernel contains cliques biggern than 15 perform pseudo branch.
                if dfvsi.graph.greedy_max_clique().len() > 5 {
                    let start = Instant::now();
                    let res = dfvsi.advanced_clique_heuristic(&DFVSInstance::bottom_up_weight_heuristic, &Digraph::cai_weight, (0.3,0f64), &priority)?;
                    if res.is_none() {
                        return Ok((results, name, (n,m), (kn, km), exRedDur, lower))
                    }
                    let res = res.expect("is not none");
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        let duration = start.elapsed();
                        results.insert("\"BranchHeuristic/cai0.3/wloc\"".to_owned(), (better_res.len(), duration.as_millis()));
                    } else {
                        let duration = start.elapsed();
                        results.insert("\"BranchHeuristic/cai0.3\"".to_owned(), (res.len(), duration.as_millis()));
                    }
                }
                // * Perfrom bottom-up/top-down/switch weight heuristics with different weight functions
                // and different weights
                // Cai:
                for i in vec![0.2,0.3,0.4] {
                    let start = Instant::now();
                    if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::cai_weight, (i, 0f64), &priority, true) {
                        if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                            let duration = start.elapsed();
                            results.insert(format!("\"BottomUp/cai{}/wloc\"", i), (better_res.len(),duration.as_millis()));
                        } else {
                            let duration = start.elapsed();
                            results.insert(format!("\"BottomUp/cai{}\"", i), (res.len(),duration.as_millis()));
                        }
                    } else {
                        break;
                    }
                    let start = Instant::now();
                    if let Some(res) = dfvsi.top_down_weight_heuristic(&Digraph::cai_weight, (i, 0f64), &priority, true) {
                        if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                            let duration = start.elapsed();
                            results.insert(format!("\"TopDown/cai{}/wloc\"", i), (better_res.len(), duration.as_millis()));
                        } else {
                            let duration = start.elapsed();
                            results.insert(format!("\"TopDown/cai{}\"", i), (res.len(), duration.as_millis()));
                        }
                    } else {
                        break;
                    }
                    let start = Instant::now();
                    if let Some(res) = dfvsi.top_bottom_switch_weight_heuristic(&Digraph::cai_weight, (i, 0f64), &priority, true) {
                        if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                            let duration = start.elapsed();
                            results.insert(format!("\"Switch/cai{}/wloc\"", i), (better_res.len(), duration.as_millis()));
                        } else {
                            let duration = start.elapsed();
                            results.insert(format!("\"Switch/cai{}\"", i), (res.len(), duration.as_millis()));
                        }
                    } else {
                        break;
                    }
                }
                // Lin:
                'outer: for i in vec![0.25,0.5,0.75] {
                    for j in vec![0.25,0.5,0.75] {
                        if i == j && i != 0.25 {
                            continue 'outer
                        }
                        let start = Instant::now();
                        if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                                let duration = start.elapsed();
                                results.insert(format!("\"BottomUp/lin({},{})/wloc\"", i, j), (better_res.len(),duration.as_millis()));
                            } else {
                                let duration = start.elapsed();
                                results.insert(format!("\"BottomUp/lin({},{})\"", i, j), (res.len(), duration.as_millis()));
                            }
                        } else {
                            break;
                        }
                        let start = Instant::now();
                        if let Some(res) = dfvsi.top_down_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                                let duration = start.elapsed();
                                results.insert(format!("\"TopDown/lin({},{})/wloc\"", i, j), (better_res.len(), duration.as_millis()));
                            } else {
                                let duration = start.elapsed();
                                results.insert(format!("\"TopDown/lin({},{})\"", i, j), (res.len(), duration.as_millis()));
                            }
                        } else {
                            break;
                        }
                        let start = Instant::now();
                        if let Some(res) = dfvsi.top_bottom_switch_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                                let duration = start.elapsed();
                                results.insert(format!("\"Switch/lin({},{})/wloc\"", i, j), (better_res.len(), duration.as_millis()));
                            } else {
                                let duration = start.elapsed();
                                results.insert(format!("\"Switch/lin({},{})\"", i, j), (res.len(), duration.as_millis()));
                            }
                        } else {
                            break;
                        }
                    }
                }
            } else {
                if dfvsi.graph.greedy_max_clique().len() > 10 {
                    let start = Instant::now();
                    let res = dfvsi.advanced_clique_heuristic(&DFVSInstance::bottom_up_weight_heuristic, &Digraph::cai_weight, (0.3,0f64), &vec![Rule::SimpleRules])?;
                    if res.is_none() {
                        // 1: timeout after reductions
                        results.insert("\"Time-out\"".to_owned(), (1,1));
                        return Ok((results, name, (n,m), (kn, km), exRedDur, lower))
                    }
                    let res = res.expect("is not none");
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        let duration = start.elapsed();
                        results.insert("\"BranchHeuristic/cai0.3/wloc\"".to_owned(), (better_res.len(),duration.as_millis()));
                    } else {
                        let duration = start.elapsed();
                        results.insert("\"BranchHeuristic/cai0.3\"".to_owned(), (res.len(), duration.as_millis()));
                    }
                }
                // Cai:
                let start = Instant::now();
                if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::cai_weight, (0.3, 0f64), &vec![Rule::SimpleRules], true) {
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        let duration = start.elapsed();
                        results.insert("\"BottomUp/cai0.3/wloc\"".to_owned(), (better_res.len(),duration.as_millis()));
                    } else {
                        let duration = start.elapsed();
                        results.insert("\"BottomUp/cai0.3\"".to_owned(), (res.len(), duration.as_millis()));
                    }
                }
            }
            Ok((results, name, (n,m), (kn, km), exRedDur, lower))
        }));
    }
    let mut out_file = File::create(csv)?;
    writeln!(&mut out_file, "instance, n(g), m(g), n(k), m(k), t(init_rules),\
    lower_bound, heuristic, result, time")?;

    // TODO join
    for join_handle in threads {
        match join_handle.join() {
            Ok(Ok((results, name, (n,m), (kn, km), exRedDur, lower))) => {
                eprintln!("joined {:?}", name);
                for (heuristic, (result, time)) in results {
                    writeln!(&mut out_file, "{:?}, {}, {}, {}, {}, {}, {}, {}, {}, {}", name, n, m, kn, km, exRedDur, lower, heuristic, result, time)?;
                }
            },
            Ok(Err(_)) => eprintln!("Some thread paniced"),
            Err(_) => eprintln!("Some thread paniced"),
        }
    }

    // print results: 
    // TODO: Check quality of both heuristic results with the lower bound

    Ok(())
}
