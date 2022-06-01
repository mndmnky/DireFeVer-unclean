use std::error;
use std::io::BufReader;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::process::Command;

use dfvs_solver::{vertex_cover::read_vc_solution, digraph::{Digraph, RebuildDigraph}, dfvs_instance::DFVSInstance, reduction_rules::Rule, cust_errors::InterruptError};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(RebuildDigraph::new(graph), None, None);
    let priority = vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC];
    dfvsi.exhaustive_reductions(&priority);
    let mut resu;
    let mut run_vc = Command::new("sh vc_solver ./into_vc.gr > vc_sol");
    if dfvsi.graph.graph.scc_graph_is_disconnected() {
        let mut solution = dfvsi.solution.clone();
        // create subgraph for each (keep node ids)
        let subgraphs: Vec<Digraph> = dfvsi.graph.graph.split_into_connected_components_alt();
        for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(RebuildDigraph::new(graph), None, None)){
            // check size and other important attributes to see if it makes sense to switch to vc.
            eprintln!("st0");
            instance.graph.graph.write_strong_graph_to_gr(BufWriter::new(File::create("./into_vc.gr")?))?;
            eprintln!("st1");
            run_vc.current_dir("./");
            run_vc.status().expect("failed");
            eprintln!("st2");
            let sol = read_vc_solution(BufReader::new(File::open("vc_sol")?))?;

            if instance.validate_left_over(&sol) {
                solution.extend(&sol);
                continue
            } else {
                let lower = sol.len();
                // also get remaining upperbound.
                let mut clone = instance.clone();
                clone.add_all_to_solution(sol);
                let opt_upper = clone.get_good_upper(false);
                instance.update_lower_bound(lower);
                instance.set_current_best(&opt_upper.ok_or(InterruptError::Unclear)?);
            }
            // Do advanced_branching and add solution 
            let sol = instance.advanced_branching(&priority, true)?;
            solution.extend(&sol);
        }
        // If we use the link node rule, we have to finallize the found solution.
        resu = dfvsi.finallize_given_solution_temp(&solution);
    } else {
        eprintln!("st0");
        dfvsi.graph.graph.write_strong_graph_to_gr(BufWriter::new(File::create("./into_vc.gr")?))?;
        eprintln!("st1");
        run_vc.current_dir("./");
        run_vc.status().expect("failed");
        eprintln!("st2");
        let sol = read_vc_solution(BufReader::new(File::open("vc_sol")?))?;

        if dfvsi.validate_left_over(&sol) {
            resu = dfvsi.solution.clone();
            resu.extend(&sol);
            resu = dfvsi.finallize_given_solution_temp(&resu);
        } else {
            let lower = sol.len();
            // also get remaining upperbound.
            let mut clone = dfvsi.clone();
            clone.add_all_to_solution(sol);
            let opt_upper = clone.get_good_upper(false);
            dfvsi.update_lower_bound(lower);
            dfvsi.set_current_best(&opt_upper.ok_or(InterruptError::Unclear)?);
            resu = dfvsi.advanced_branching(&priority, true)?;
        }
    }
    DFVSInstance::write_solution(&resu, &mut stdout)?;
    Ok(())
}
