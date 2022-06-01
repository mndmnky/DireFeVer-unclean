use std::error;
use std::io;

use dfvs_solver::{digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance,
                  reduction_rules::Rule};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Input graph:\t{} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);

    let (mut simple_time, mut simple_nkc, mut simple_mkc) = (0, 0, 0);
    let (mut scc_time, mut scc_nkc, mut scc_mkc) = (0, 0, 0);
    let (mut adv_scc_time, mut adv_scc_nkc, mut adv_scc_mkc) = (0, 0, 0);
    let (mut dome_time, mut dome_nkc, mut dome_mkc) = (0, 0, 0);
    let (mut clique_time, mut clique_nkc, mut clique_mkc) = (0, 0, 0);
    let (mut adv_clique_time, mut adv_clique_nkc, mut adv_clique_mkc) = (0, 0, 0);
    let (mut exh_clique_time, mut exh_clique_nkc, mut exh_clique_mkc) = (0, 0, 0);
    let (mut core_clique_time, mut core_clique_nkc, mut core_clique_mkc) = (0, 0, 0);
    let (mut core_daisy_time, mut core_daisy_nkc, mut core_daisy_mkc) = (0, 0, 0);
    let (mut core_min_deg_time, mut core_min_deg_nkc, mut core_min_deg_mkc) = (0, 0, 0);
    let (mut local_k_flower_time, mut local_k_flower_nkc, mut local_k_flower_mkc) = (0, 0, 0);
    loop {
        let (time, nkc, mkc) = dfvsi.simple_reduction_stats();
        simple_time += time;
        simple_nkc += nkc;
        simple_mkc += mkc;
        //let (time, nkc, mkc) = dfvsi.scc_reduction_stats();
        //scc_time += time;
        //scc_nkc += nkc;
        //scc_mkc += mkc;
        //if nkc>0 || mkc>0 {
        //    continue;
        //}
        let (time, nkc, mkc) = dfvsi.dome_reduction_stats();
        dome_time += time;
        dome_nkc += nkc;
        dome_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.clique_reduction_stats();
        clique_time += time;
        clique_nkc += nkc;
        clique_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.advanced_clique_reduction_stats();
        adv_clique_time += time;
        adv_clique_nkc += nkc;
        adv_clique_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.exhaustive_clique_reduction_stats();
        exh_clique_time += time;
        exh_clique_nkc += nkc;
        exh_clique_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.core_clique_reduction_stats();
        core_clique_time += time;
        core_clique_nkc += nkc;
        core_clique_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.core_daisy_reduction_stats();
        core_daisy_time += time;
        core_daisy_nkc += nkc;
        core_daisy_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.core_min_deg_reduction_stats();
        core_min_deg_time += time;
        core_min_deg_nkc += nkc;
        core_min_deg_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        let (time, nkc, mkc) = dfvsi.adv_scc_reduction_stats();
        adv_scc_time += time;
        adv_scc_nkc += nkc;
        adv_scc_mkc += mkc;
        if nkc>0 || mkc>0 {
            continue;
        }
        //let (_, upper1) = dfvsi.highest_degree_heuristic_stats(); 
        //dfvsi.update_upper_bound(upper1);
        //let (time, nkc, mkc) = dfvsi.local_k_flower_reduction_stats();
        //local_k_flower_time += time;
        //local_k_flower_nkc += nkc;
        //local_k_flower_mkc += mkc;
        //changed = nkc>0 || mkc>0 || changed;
        break;
    }    

    eprintln!("Simple rules reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", simple_nkc, simple_mkc, simple_time);
    eprintln!("SCC rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", scc_nkc, scc_mkc, scc_time);
    eprintln!("Adv SCC rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", adv_scc_nkc, adv_scc_mkc, adv_scc_time);
    eprintln!("Dome rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", dome_nkc, dome_mkc, dome_time);
    eprintln!("Cluster rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", clique_nkc, clique_mkc, clique_time);
    eprintln!("Adv. cliq. rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", adv_clique_nkc, adv_clique_mkc, adv_clique_time);
    eprintln!("Exh. cliq. rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", exh_clique_nkc, exh_clique_mkc, exh_clique_time);
    eprintln!("Core cliq. rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", core_clique_nkc, core_clique_mkc, core_clique_time);
    eprintln!("Core daisy rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", core_daisy_nkc, core_daisy_mkc, core_daisy_time);
    eprintln!("Core min. rule reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", core_min_deg_nkc, core_min_deg_mkc, core_min_deg_time);
    eprintln!("Local k-flower reduced\t({:>4} Nodes, {:>5} Edges, and took {:?}ms)", local_k_flower_nkc, local_k_flower_mkc, local_k_flower_time);

    let (time, upper1) = dfvsi.highest_degree_heuristic_stats().unwrap(); 
    eprintln!("Highest degree heuristic found a solution of size {} and took {}ms", upper1, time);
    let (time, upper3) = dfvsi.highest_strong_degree_heuristic_stats().unwrap(); 
    eprintln!("Highest strong degree heuristic found a solution of size {} and took {}ms", upper3, time);
    let (time, lower, upper2) = dfvsi.clique_heuristic_stats(&vec![Rule::SimpleRules]).unwrap(); 
    eprintln!("Clique heuristic found a lower bound of {}, a upper bound of {} and took {}ms", lower, upper2, time);

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Kernel graph: {} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);
    let min_direct = dfvsi.graph.graph.min_direct_degree_stat();
    eprintln!("Kernel graph has a min. direct degree of {}", min_direct);
    let daisy_len = dfvsi.graph.graph.max_daisy_stat();
    eprintln!("Kernel graph contains a {} daisy", daisy_len);
    let (time, greedy_clique_len) = dfvsi.graph.graph.greedy_max_clique_stat();
    eprintln!("Kernel graph contains a {} clique, which took {} ms to find", greedy_clique_len, time);
    eprintln!("Current solution size: {}", dfvsi.solution.len());
    let (n, m, solution_length) = dfvsi.instance_stats();
    eprintln!("Branching instance graph: {} nodes, {} edges, current solution size: {}", n, m, solution_length);
    // Panics on interrupt
    let (time, solution) = dfvsi.advanced_branching_stat(&vec![Rule::SimpleRules], true).unwrap();
    eprintln!("The branching algorithm found a solution of size {}, which took {} ms to find.", solution.len(), time);
    eprintln!("Write current solution.");
    DFVSInstance::write_solution(&solution, &mut stdout)?;
    Ok(())
}
