use std::error;
use std::io;

use dfvs_solver::{digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);
    let dfvsi_org = dfvsi.clone();

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Input graph:\t{} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);

    let priority = vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
    dfvsi.exhaustive_reductions(&priority);
    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Kernel graph:\t{} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);
    eprintln!("Current solution size: {}", dfvsi.solution.len());

    dfvsi.set_time_interrupter(3600000);
    //let mut dfvsi_branch = dfvsi.clone();

    eprintln!("Current rebuild size: {}", dfvsi.merge_nodes.len());
    let mut resu = dfvsi.advanced_branching(&priority, true)?;
    eprintln!("found solution of size {}",resu.len());
    //if let Some(improved_solution) = dfvsi_org.exhaustive_local_search(&resu){
    //    eprintln!("local search could improve the solution to size: {}", improved_solution.len());
    //    //let dif = resu.difference(&improved_solution);
    //    //eprintln!("not actually in the solution were: {:?}", dif);
    //    resu = improved_solution;
    //} 

    eprintln!("The branching algorithm found a solution of size {}.", resu.len());
    // validate 
    if dfvsi_org.validate(&resu) {
        eprintln!("Solution could be correct.");
    }
    eprintln!("Write current solution.");
    //DFVSInstance::write_solution(&resu, &mut stdout)?;
    Ok(())
}
