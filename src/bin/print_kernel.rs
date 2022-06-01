use std::error;
use std::io;
use std::path::PathBuf;
use std::fs::File;
use clap::{Arg, App};

use dfvs_solver::{digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule};

pub fn main() -> Result<(), Box<dyn error::Error>> {

    let m = App::new("print_kernel")
        .arg(Arg::new("scc_folder")
             .takes_value(true)
             .required(true)
             .short('s'))
        .get_matches();
    // Get, as input all public instances. 
    let out_dir: PathBuf = PathBuf::from(m.value_of("scc_folder").unwrap());

    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);

    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Input graph:\t{} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);

    let priority = vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
    let (stats, _) = dfvsi.priority_rule_stats(priority);
    for stat in stats {
        eprintln!("{:?}:\t{:>5} nodes,\t{:>6} edges,\t{:>10} ms", stat.rule, stat.reduced_nodes, stat.reduced_edges, stat.time_took);
    }
    let (n, m, max, min) = dfvsi.graph.graph.n_m_max_min_stats();
    eprintln!("Kernel graph:\t{} nodes, {} edges, max degree {}, min degree {}", n, m, max, min);
    eprintln!("Min out {} max out {}",dfvsi.graph.graph.min_out_degree().unwrap(), dfvsi.graph.graph.max_out_degree().unwrap());

    if dfvsi.graph.graph.scc_graph_is_disconnected() {
        // create subgraph for each (keep node ids)
        let subgraphs: Vec<Digraph> = dfvsi.graph.graph.split_into_connected_components_alt();
        for (i,instance) in subgraphs.into_iter().map(|graph| DFVSInstance::new(RebuildDigraph::new(graph), None, None)).enumerate(){
            let mut file = out_dir.clone();
            file.push(format!("scc{}",i));
            instance.graph.graph.write_graph(File::create(file)?)?;
        }
    }

    dfvsi.graph.graph.write_graph(stdout)?;

    Ok(())
}
