use clap::{Arg, App};
use std::error;
use std::path::PathBuf;
use std::fs::File;
use std::thread;
use std::io::{BufReader, Write};
use dfvs_solver::{digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule, cust_errors::ImportError, statistics::RuleStats};
use regex::Regex;

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let digits = Regex::new(r"\d+").unwrap();
    let m = App::new("statistics")
        .arg(Arg::new("files")
             .takes_value(true)
             .multiple_values(true)
             .short('f'))
        .arg(Arg::new("seperate_rules")
             .takes_value(true)
             .short('s'))
        .arg(Arg::new("core_compare")
             .takes_value(true)
             .short('c'))
        .arg(Arg::new("all")
             .takes_value(true)
             .short('a'))
        .arg(Arg::new("time_out")
             .takes_value(true)
             .validator_regex(&digits, "only numbers are allowed")
             .short('t'))
        .get_matches();
    // Get, as input all public instances. 
    let files: Vec<PathBuf> = m.values_of("files").unwrap().map(|p| PathBuf::from(p)).collect();
    let sep_rules: bool = m.is_present("seperate_rules");
    let core_compare: bool = m.is_present("core_compare");
    let all: bool = m.is_present("all");
    let time_out: Option<u128> = m.value_of("time_out").map(|val| val.parse::<u128>().expect("input checked with regex"));
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }
    if all {
        let csv: &str = m.value_of("all").unwrap();
        // Only record rules and kernel size:
        let mut threads = Vec::new();
        for (graph, name) in graphs.clone() {
            threads.push(thread::spawn(move || {
                let prio = vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
                let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph.clone()), None, None);
                if let Some(time) = time_out {
                    dfvsi.set_time_interrupter(time);
                }
                let (prio_stats, _) = dfvsi.priority_rule_stats(prio);
                (prio_stats, dfvsi, graph, name)
            }));
        }
        let mut out_file = File::create(csv)?;
        writeln!(&mut out_file, "instance, org_n, org_m, reduction, prio_nr, nodes_reduced, edges_reduced, time_took")?;
        // TODO: join first finished
        for join_handle in threads {
            match join_handle.join() {
                Ok((stats, df, g, name)) => {
                    eprintln!("joined {:?}", name);
                    for (i, stat) in stats.iter().enumerate() {
                        writeln!(&mut out_file, "{:?}, {}, {}, {:?}, {}, {}, {}, {}", name, g.num_nodes(), g.num_edges(), stat.rule, i, stat.reduced_nodes, stat.reduced_edges, stat.time_took);
                    }
                },
                Err(_) => eprintln!("Some thread paniced"),
            }
        }
    }
    if core_compare {
        let csv: &str = m.value_of("core_compare").unwrap();
        // Only record rules and kernel size:
        let mut threads = Vec::new();
        for (graph, name) in graphs.clone() {
            threads.push(thread::spawn(move || {
                // TODO print for Rule
                // TODO k flower in
                let prio1 = vec![Rule::SimpleRules, Rule::Dome, Rule::Clique, Rule::SCC, Rule::KFlower];
                let prio2 = vec![Rule::SimpleRules, Rule::Dome, Rule::Clique, Rule::Core, Rule::SCC, Rule::KFlower];
                let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph.clone()), None, None);
                let mut dfvsi2 = dfvsi.clone();
                if let Some(time) = time_out {
                    dfvsi.set_time_interrupter(time);
                }
                let prio_stats_1 = dfvsi.priority_rule_stats(prio1);
                if let Some(time) = time_out {
                    dfvsi2.set_time_interrupter(time);
                }
                let prio_stats_2 = dfvsi2.priority_rule_stats(prio2);
                (prio_stats_1, prio_stats_2, dfvsi, dfvsi2, graph, name)
            }));
        }
        let mut out_file = File::create(csv)?;
        writeln!(&mut out_file, "instance, r1_sr_n, r1_sr_m, r1_sr_t, r1_d_n, r1_d_m, r1_d_t,\
         r1_cl_n, r1_cl_m, r1_cl_t, r1_scc_n, r1_scc_m, r1_scc_t, r1_kf_n, r1_kf_m, r1_kf_t,\
         r1_finished,\
         r2_sr_n, r2_sr_m, r2_sr_t, r2_d_n, r2_d_m, r2_d_t, r2_cl_n, r2_cl_m, r2_cl_t,\
         r2_co_n, r2_co_m, r2_co_t, r2_scc_n, r2_scc_m, r2_scc_t, r2_kf_n, r2_kf_m, r2_kf_t,\
         r2_finished, n, m, k1_n, k1_m, k2_n, k2_m")?;
        // TODO: join first finished
        for join_handle in threads {
            match join_handle.join() {
                Ok((rs1, rs2, df1, df2, graph, name)) => {
                    eprintln!("joined {:?}", name);
                    let red1 = rs1.0.iter().flat_map(|red| format!("{}, {}, {}, ", red.reduced_nodes, red.reduced_edges, red.time_took).chars().collect::<Vec<_>>()).collect::<String>();
                    let red2 = rs2.0.iter().flat_map(|red| format!("{}, {}, {}, ", red.reduced_nodes, red.reduced_edges, red.time_took).chars().collect::<Vec<_>>()).collect::<String>();
                    writeln!(&mut out_file, "{:?}, {}{}, {}{}, {}, {}, {}, {}, {}, {}", name, red1, rs1.1, red2, rs2.1, graph.num_nodes(), graph.num_edges(), df1.graph.graph.num_nodes(), df1.graph.graph.num_edges(), df2.graph.graph.num_nodes(), df2.graph.graph.num_edges())?;
                },
                Err(_) => eprintln!("Some thread paniced"),
            }
        }
    }
    if sep_rules {
        let csv: &str = m.value_of("seperate_rules").unwrap();
        // Only record rules and kernel size:
        let mut threads = Vec::new();
        for (graph, name) in graphs {
            threads.push(thread::spawn(move || {
                let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph.clone()), None, None);
                dfvsi.apply_simple_rules();
                let dome_prio = vec![Rule::SimpleRules, Rule::Dome];
                let scc_prio = vec![Rule::SimpleRules, Rule::SCC];
                let core_prio = vec![Rule::SimpleRules, Rule::Clique, Rule::Core];
                let kflower_prio = vec![Rule::SimpleRules, Rule::KFlower];
                let mut dome_ins = dfvsi.clone();
                if let Some(time) = time_out {
                    dome_ins.set_time_interrupter(time);
                }
                let dome_stats = dome_ins.overall_rule_stats(dome_prio);
                let mut scc_ins = dfvsi.clone();
                if let Some(time) = time_out {
                    scc_ins.set_time_interrupter(time);
                }
                let scc_stats = scc_ins.overall_rule_stats(scc_prio);
                let mut core_ins = dfvsi.clone();
                if let Some(time) = time_out {
                    core_ins.set_time_interrupter(time);
                }
                let core_stats = core_ins.priority_rule_stats(core_prio);
                let mut kflower_ins = dfvsi.clone();
                if let Some(time) = time_out {
                    kflower_ins.set_time_interrupter(time);
                }
                let kflower_stats = kflower_ins.overall_rule_stats(kflower_prio);
                (dome_stats, scc_stats, core_stats, kflower_stats, dfvsi, dome_ins, scc_ins, core_ins, kflower_ins, graph, name)
            }));
        }
        let mut out_file = File::create(csv)?;
        writeln!(&mut out_file, "instance, \
        dome_n, dome_m, dome_t, dome_finished, \
        scc_n, scc_m, scc_t, scc_finished, \
        core_n, clique_of_n, core_of_n, core_m, clique_of_m, core_of_m, core_t, clique_of_t, core_of_t, core_finished, \
        kflower_n, kflower_m, kflower_t, kflower_finished, \
         n, m, sr_n, sr_m")?;
        // TODO: join first finished
        for join_handle in threads {
            match join_handle.join() {
                Ok((
                    (dome_t, dome_n, dome_m, dome_finished),
                    (scc_t, scc_n, scc_m, scc_finished),
                    core_rule_set,
                    (kflower_t, kflower_n, kflower_m, kflower_finished),
                    org_sr_ins,
                    dome_ins,
                    scc_ins,
                    core_ins,
                    kflower_ins,
                    graph,
                    name))
                    => {
                    eprintln!("joined {:?}", name);
                    let core_n: usize = core_rule_set.0
                        .iter()
                        .map(|red| red.reduced_nodes)
                        .sum();
                    let core_m: usize = core_rule_set.0
                        .iter()
                        .map(|red| red.reduced_edges)
                        .sum();
                    let core_t: u128 = core_rule_set.0
                        .iter()
                        .map(|red| red.time_took)
                        .sum();
                    let core_of_n = core_rule_set.0[2].reduced_nodes;
                    let core_of_m = core_rule_set.0[2].reduced_edges;
                    let core_of_t = core_rule_set.0[2].time_took;
                    let clique_of_n = core_rule_set.0[1].reduced_nodes;
                    let clique_of_m = core_rule_set.0[1].reduced_edges;
                    let clique_of_t = core_rule_set.0[1].time_took;
                    writeln!(&mut out_file, "{:?}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}", name, 
                             dome_n, dome_m, dome_t, dome_finished,
                             scc_n, scc_m, scc_t, scc_finished,
                             core_n, clique_of_n, core_of_n, core_m, clique_of_m, core_of_m, core_t, clique_of_t, core_of_t, core_rule_set.1,
                             kflower_n, kflower_m, kflower_t, kflower_finished,
                             graph.num_nodes(), graph.num_edges(),
                             org_sr_ins.graph.graph.num_nodes(), org_sr_ins.graph.graph.num_edges())?;
                },
                Err(_) => eprintln!("Some thread paniced"),
            }
        }
    }
    Ok(())
}


