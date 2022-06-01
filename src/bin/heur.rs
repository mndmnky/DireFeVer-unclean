use std::error;
use std::io;
use std::collections::HashMap;
use dfvs_solver::{digraph::{Digraph,RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Input graph:\t{:>4} Nodes, {:>5} Edges, max degree: {:>3}, min degree: {:>2}", n, m, max, min);

    // Set timer for 30 minutes, which is three times the actual time.
    dfvsi.set_time_interrupter(1800000);

    // initial reduction
    let priority = vec![Rule::SimpleRules, Rule::Dome, Rule::Clique, Rule::Core, Rule::SCC];
    if !dfvsi.exhaustive_reductions(&priority){
        eprintln!("Time didn't suffice to finish the execution of the rules");
        return Ok(())
    }

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Kernel graph:\t{:>4} Nodes, {:>5} Edges, max degree: {:>3}, min degree: {:>2}", n, m, max, min);

    let mut results = HashMap::new();

    if n<300 {
        eprintln!{"Instance can be solved efficiently."}
        // Solve 
        // TODO: for the actual case we need a branching that returns an intermediate solution in
        // case of a sigint.
        let resu = dfvsi.advanced_branching(&priority, true)?;
        eprintln!{"Exact result contains {} nodes.", resu.len()}
        eprintln!("Write solution to stdout.");
        DFVSInstance::write_solution(&resu, &mut stdout)?;
    } else if n < 10000 {
        // Apply lots of upper bounds 
        // * If kernel contains cliques biggern than 5 perform pseudo branch.
        if dfvsi.graph.greedy_max_clique().len() > 5 {
            let res = dfvsi.advanced_clique_heuristic(&DFVSInstance::bottom_up_weight_heuristic, &Digraph::cai_weight, (0f64,0.3), &priority)?;
            if res.is_none() {
                eprintln!("Time interrupt (or sigint)");
                // print results so far if any
                return Ok(())
            }
            eprintln!("Act finished bheuristic");
            let res = res.expect("is not none");
            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                results.insert("BranchHeuristic/cai0.3/wloc".to_owned(), better_res.len());
            } else {
                results.insert("BranchHeuristic/cai0.3".to_owned(), res.len());
            }
        }
        // * Perfrom bottom-up/top-down/switch weight heuristics with different weight functions
        // and different weights
        // Cai:
        for i in vec![0.2,0.3,0.4] {
            if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::cai_weight, (0f64, i), &priority, true) {
                if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                    results.insert(format!("BottomUp/cai{}/wloc", i), better_res.len());
                } else {
                    results.insert(format!("BottomUp/cai{}", i), res.len());
                }
            } else {
                break;
            }
            if let Some(res) = dfvsi.top_down_weight_heuristic(&Digraph::cai_weight, (0f64, i), &priority, true) {
                if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                    results.insert(format!("TopDown/cai{}/wloc", i), better_res.len());
                } else {
                    results.insert(format!("TopDown/cai{}", i), res.len());
                }
            } else {
                break;
            }
            if let Some(res) = dfvsi.top_bottom_switch_weight_heuristic(&Digraph::cai_weight, (0f64, i), &priority, true) {
                if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                    results.insert(format!("Switch/cai{}/wloc", i), better_res.len());
                } else {
                    results.insert(format!("Switch/cai{}", i), res.len());
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
                if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        results.insert(format!("BottomUp/lin({},{})/wloc", i, j), better_res.len());
                    } else {
                        results.insert(format!("BottomUp/lin({},{})", i, j), res.len());
                    }
                } else {
                    break;
                }
                if let Some(res) = dfvsi.top_down_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        results.insert(format!("TopDown/lin({},{})/wloc", i, j), better_res.len());
                    } else {
                        results.insert(format!("TopDown/lin({},{})", i, j), res.len());
                    }
                } else {
                    break;
                }
                if let Some(res) = dfvsi.top_bottom_switch_weight_heuristic(&Digraph::lin_weight, (i, j), &priority, true) {
                    if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                        results.insert(format!("Switch/lin({},{})/wloc", i, j), better_res.len());
                    } else {
                        results.insert(format!("Switch/lin({},{})", i, j), res.len());
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        if dfvsi.graph.greedy_max_clique().len() > 10 {
            let res = dfvsi.advanced_clique_heuristic(&DFVSInstance::bottom_up_weight_heuristic, &Digraph::cai_weight, (0f64,0.3), &vec![Rule::SimpleRules])?;
            if res.is_none() {
                eprintln!("Time interrupt (or sigint)");
                // print results so far if any
                return Ok(())
            }
            let res = res.expect("is not none");
            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                results.insert("BranchHeuristic/cai0.3/wloc".to_owned(), better_res.len());
            } else {
                results.insert("BranchHeuristic/cai0.3".to_owned(), res.len());
            }
        }
        // Cai:
        if let Some(res) = dfvsi.bottom_up_weight_heuristic(&Digraph::cai_weight, (0f64, 0.3), &vec![Rule::SimpleRules], true) {
            if let Some(better_res) = dfvsi.exhaustive_local_search(&res) {
                results.insert("BottomUp/cai0.3/wloc".to_owned(), better_res.len());
            } else {
                results.insert("BottomUp/cai0.3".to_owned(), res.len());
            }
        }
    }
    // print results: 
    eprintln!("{:?}", results);
    // TODO: Check quality of both heuristic results with the lower bound

    Ok(())
}
