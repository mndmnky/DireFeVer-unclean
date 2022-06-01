use crate::dfvs_instance::DFVSInstance;
use crate::digraph::{Digraph};
use crate::cust_errors::GraphError;
use crate::reduction_rules::Rule;
use std::collections::HashSet;
use fxhash::{FxHashSet};
use std::time::Instant;
use std::error;

pub struct Statistics {
    pub num_nodes: usize,
    pub num_edges: usize,
    pub best_heuristic: String, 
    pub best_heuristic_size: usize,
    pub rules: Vec<RuleStats>,
    pub kernel_num_nodes: usize,
    pub kernel_num_edges: usize,
}

pub struct RuleStats {
    pub rule: Rule,
    pub reduced_nodes: usize,
    pub reduced_edges: usize,
    pub time_took: u128,
}

impl RuleStats {

    pub fn new(rule: Rule) -> Self {
        RuleStats {
            rule,
            reduced_nodes: 0,
            reduced_edges: 0,
            time_took: 0,
        }
    }

    pub fn add(&mut self, n: usize, m: usize, time: u128) {
        self.reduced_nodes += n;
        self.reduced_edges += m;
        self.time_took += time;
    }

}


impl Digraph {
    /// Returns the number of nodes, the number of edges, the maximal degree and the minimal degree
    /// of `self`.
    pub fn n_m_max_min_stats(&self) -> (usize, usize, usize, usize) {
        let n = self.num_nodes();
        if n == 0{
            (0, 0, 0, 0)
        } else {
        (self.num_nodes(), self.num_edges(), self.max_degree().expect("n > 0"), self.min_degree().expect("n > 0"))
        }
    }

    /// Returns the size of the largest daisy or `0` if no daisy exists.
    pub fn max_daisy_stat(&self) -> usize {
        if let Some((_, neighs)) = self.get_max_strong_degree_node_and_neighbors() {
            neighs.len()
        } else {
            0
        }
    }

    /// Returns the execution time of the greedy max clique algorithm and the size of the largest clique the algorithm finds.
    pub fn greedy_max_clique_stat(&self) -> (u128, usize) {
        let start_time = Instant::now();
        let size = self.greedy_max_clique().len();
        let duration = start_time.elapsed();
        (duration.as_millis(), size)
    }

    /// Returns the minimum direct degree of `self`.
    pub fn min_direct_degree_stat(&self) -> usize {
        self.nodes().map(|node| self.min_direct_degree(node).expect("`node` is in `.nodes()`")).min().unwrap_or(0)
    }

    /// Returns the maximum degree of `self`.
    pub fn max_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.degree(node).expect("`node` does not exist")).max()
    }

    /// Returns the maximum incoming degree of `self`.
    pub fn max_in_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.in_degree(node).expect("`node` does not exist")).max()
    }

    /// Returns the maximum outgoing degree of `self`.
    pub fn max_out_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.out_degree(node).expect("`node` does not exist")).max()
    }

    /// Returns the minimum degree of `self`.
    pub fn min_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.degree(node).expect("`node` does not exist")).min()
    }

    /// Returns the minimum incoming degree of `self`.
    pub fn min_in_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.in_degree(node).expect("`node` does not exist")).min()
    }

    /// Returns the minimum outgoing degree of `self`.
    pub fn min_out_degree(&self) -> Option<usize> {
        self.nodes().map(|node| self.out_degree(node).expect("`node` does not exist")).min()
    }

}

impl DFVSInstance {

    /// Returns execution time, node- and edge kill count of exhaustive application of the simple reduction rules. 
    /// If kill count = 0, no changes have been made.
    pub fn simple_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_simple_rules();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time, node- and edge kill count of the application of the scc reduction rules. 
    /// If kill count = 0, no changes have been made.
    pub fn scc_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_scc_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time, node- and edge kill count of the application of the advanced scc reduction rules. 
    /// If kill count = 0, no changes have been made.
    pub fn adv_scc_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_advanced_scc_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time, node- and edge kill count of the application of the dome reduction rules. 
    /// If kill count = 0, no changes have been made.
    pub fn dome_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_dome_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the local k flower reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn local_k_flower_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_local_k_flower();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the clique reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn clique_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_clique_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the advanced clique reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn advanced_clique_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_advanced_clique_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the exhaustive application of the clique reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn exhaustive_clique_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_exhaustive_clique_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the core clique reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn core_clique_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_core_clique_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the core daisy reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn core_daisy_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_exhaustive_daisy_core_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and kill count of the application of the core min. directed degree reduction rule. 
    /// If kill count = 0, no changes have been made.
    pub fn core_min_deg_reduction_stats(&mut self) -> (u128, usize, usize) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        self.apply_exhaustive_min_direct_core_rule();
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count)
    }

    /// Returns execution time and the solution size of the highest degree heuristic.
    pub fn highest_degree_heuristic_stats(&self) -> Option<(u128, usize)> {
        let start = Instant::now();
        if let Some(heur) = self.highest_degree_heuristic() {
            let duration = start.elapsed();
            return Some((duration.as_millis(), heur.len()))
        }
        None
    }

    /// Returns execution time and the solution size of the highest strong degree heuristic, or
    /// `None` if execution was interrupted by `self.interrupter`.
    pub fn highest_strong_degree_heuristic_stats(&self) -> Option<(u128, usize)> {
        let start = Instant::now();
        if let Some(heur) = self.highest_strong_degree_heuristic() {
            let duration = start.elapsed();
            return Some((duration.as_millis(), heur.len()))
        }
        None
    }

    /// Returns execution time and both bounds of the clique heuristic, or `None` if execution was
    /// interrupted by `self.interrupter`.
    /// Skips initial reductions.
    pub fn clique_heuristic_stats(&self, rule_priority: &Vec<Rule>) -> Option<(u128, usize, usize)> {
        let start = Instant::now();
        if let Some((lower_bound, upper_bound, _)) = self.clique_heuristic(rule_priority, true) {
        let duration = start.elapsed();
        return Some((duration.as_millis(), lower_bound, upper_bound))
        } 
        None
    }

    /// Applies different weighted heuristics with different weights to the graph, picks the best,
    /// uses exhaustive local search and prints the size of the solution.
    pub fn weight_heuristic_stats(&self, rule_priority: &Vec<Rule>, skip_initial_rules: bool) {
        let (mut alpha, mut beta) = (0f64,0f64);
        let mut fun = "";
        let rem_nodes: FxHashSet<usize> = self.graph.nodes().collect();
        let mut best_solution: FxHashSet<usize> = rem_nodes.union(&self.solution).copied().collect();
        // Bottom up: 
        // cai-weight:
        for f in vec![0.2,0.3,0.4,0.5,0.6,0.7] {
            if let Some(sol) = self.bottom_up_weight_heuristic(&Digraph::cai_weight,(f,0f64),rule_priority, skip_initial_rules){
                if sol.len() < best_solution.len(){
                    alpha = f;
                    beta = 0f64;
                    fun = "bottom_up_cai";
                    best_solution = sol;
                }
            } else {
                // TODO: if we add the other rules we have to do more here.
                break
            }
        }
        //// lin-weight:
        //for a in vec![0.0,0.3,0.6,1.0] {
        //    for b in vec![0.0,0.3,0.6,1.0] {
        //        if a == 0.0 && b == 0.0 {
        //            continue;
        //        }
        //        let sol = self.bottom_up_weight_heuristic(&Digraph::lin_weight,(a,b));
        //        if sol.len() < best_solution.len() {
        //            alpha = a;
        //            beta = b;
        //            fun = "bottom_up_lin";
        //            best_solution = sol;
        //        }
        //    }
        //}
        //// Top down:
        //// cai-weight:
        //for f in vec![0.2,0.3,0.4,0.5,0.6,0.7] {
        //    let sol = self.top_down_weight_heuristic(&Digraph::cai_weight,(f,0f64));
        //    if sol.len() < best_solution.len() {
        //        alpha = f;
        //        beta = 0f64;
        //        fun = "top_down_cai";
        //        best_solution = sol;
        //    }
        //}
        //// lin-weight
        //for a in vec![0.0,0.3,0.6,1.0] {
        //    for b in vec![0.0,0.3,0.6,1.0] {
        //        if a == 0.0 && b == 0.0 {
        //            continue;
        //        }
        //        let sol = self.top_down_weight_heuristic(&Digraph::lin_weight,(a,b));
        //        if sol.len() < best_solution.len() {
        //            alpha = a;
        //            beta = b;
        //            fun = "top_down_lin";
        //            best_solution = sol;
        //        }
        //    }
        //}
        //// Top-bottom-switch:
        //// cai-weight:
        //for f in vec![0.2,0.3,0.4,0.5,0.6,0.7] {
        //    let sol = self.top_down_weight_heuristic(&Digraph::cai_weight,(f,0f64));
        //    if sol.len() < best_solution.len() {
        //        alpha = f;
        //        beta = 0f64;
        //        fun = "switch_cai";
        //        best_solution = sol;
        //    }
        //}
        //// lin-weight:
        //for a in vec![0.0,0.3,0.6,1.0] {
        //    for b in vec![0.0,0.3,0.6,1.0] {
        //        if a == 0.0 && b == 0.0 {
        //            continue;
        //        }
        //        let sol = self.top_down_weight_heuristic(&Digraph::lin_weight,(a,b));
        //        if sol.len() < best_solution.len() {
        //            alpha = a;
        //            beta = b;
        //            fun = "switch_lin";
        //            best_solution = sol;
        //        }
        //    }
        //}
        if fun == "" {
            eprintln!("No heuristic could finish before interrupt");
        } else {
            eprintln!("Best heuristic: {}, with alpha {} and beta {}, yielded a solution of size: {}", fun, alpha, beta, best_solution.len());
            if let Some(improved) = self.exhaustive_local_search(&best_solution) {
                eprintln!("Local search improved to: {}", improved.len());
            }
        }
    }

    /// Executes the clique heuristic to find an upper bound and uses exhaustive local search to
    /// fix the solution as far as possible.
    /// Skips initial reductions.
    pub fn clique_and_rebuild_heuristic_stats(&self, rule_priority: &Vec<Rule>) {
        if let Some((_,_,sol)) = self.clique_heuristic(rule_priority, true) {
            eprintln!("Clique heuristic yielded a solution of size: {}", sol.len());
            if let Some(improved) = self.exhaustive_local_search(&sol) {
                eprintln!("Local search improved to: {}", improved.len());
            }
        } else {
            eprintln!("The clique heuristic ran into an interrupt");
        }
    }

    /// Returns execution time and the size of the solution of the advanced branching algorithm, or
    /// `None` if the execution was interrupted.
    pub fn advanced_branching_stat(&mut self, rule_priority: &Vec<Rule>, skip_initial_rules: bool) -> Result<(u128, FxHashSet<usize>), Box<dyn error::Error>> {
        let start = Instant::now();
        // TODO: handle possible errors:
        let resu = self.advanced_branching(rule_priority, skip_initial_rules)?;
        let duration = start.elapsed();
        Ok((duration.as_millis(), resu))
    }
 
    /// Returns the number of nodes, the number of edges and the size of the current solution size
    /// of the branch instance.
    pub fn instance_stats(&self) -> (usize, usize, usize) {
        (self.graph.num_nodes(), self.graph.num_edges(), self.solution.len())
    }

    /// 
    pub fn advanced_clique_heuristic_stats(&self, priority_rules: &Vec<Rule>) -> Result<(),GraphError> {
        if let Some(sol) = self.advanced_clique_heuristic(&DFVSInstance::bottom_up_weight_heuristic, &Digraph::cai_weight, (0f64,0.3), priority_rules)? {
            eprintln!("Advanced clique heuristic yielded a solution of size: {}", sol.len());
            let mut delta_solution = sol.difference(&self.solution).copied().collect();
            loop {
                let new_sol = self.local_search(&delta_solution);
                if let Some(sol) = new_sol {
                    delta_solution = sol;
                } else {
                    break
                }
            }
            eprintln!("Local search improved to: {}", delta_solution.len() + self.solution.len());
            return Ok(())
        }
        eprintln!("Execution was interrupted");
        Ok(())
    }

    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// The priority order is roughly chosen by the time consumption of the respective rules. 
    /// Records the indipendent effect of each rule and returns those stats, together with an
    /// indicator whether this function finished or was interrupted.
    /// TODO: Only record time if something changed?
    pub fn priority_rule_stats(&mut self, priority_list: Vec<Rule>) -> (Vec<RuleStats>, bool) {
        let mut rule_stats: Vec<RuleStats> = priority_list.iter().map(|rule| RuleStats::new(*rule)).collect();
        'outer: while !self.check_interrupt() {
            for rule_stat in &mut rule_stats {
                match rule_stat.rule {
                    Rule::SimpleRules => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_simple_rules();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::SCC => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        // TODO: simple scc can be removed
                        self.apply_advanced_scc_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::Dome => {
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_dome_rule();
                        let duration = start_time.elapsed();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(0, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::Clique => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        // TODO: one rule that conquers all?
                        self.apply_exhaustive_clique_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::Core => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        // TODO: one rule that conquers all?
                        self.apply_exhaustive_core_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::KFlower => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.compute_and_set_simple_upper_lower();
                        self.apply_local_k_daisy() ;
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::LinkNode => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_link_node_rules();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::Crown => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_crown_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::TwinNodes => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_twin_nodes_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                    Rule::Dominion => {
                        let nodes_before = self.graph.num_nodes();
                        let edges_before = self.graph.num_edges();
                        let start_time = Instant::now();
                        self.apply_dominion_rule();
                        let duration = start_time.elapsed();
                        let nodes_kill_count = nodes_before - self.graph.num_nodes();
                        let edges_kill_count = edges_before - self.graph.num_edges();
                        rule_stat.add(nodes_kill_count, edges_kill_count, duration.as_millis());
                        if edges_kill_count > 0 {
                            continue 'outer
                        }
                    },
                }
            }
            return (rule_stats, true)
        }
        return (rule_stats, false)
    }

    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// The priority order is roughly chosen by the time consumption of the respective rules. 
    /// Does not record indipendent changes, but overall changes.
    pub fn overall_rule_stats(&mut self, priority_list: Vec<Rule>) -> (u128, usize, usize, bool) {
        let nodes_before = self.graph.num_nodes();
        let edges_before = self.graph.num_edges();
        let start_time = Instant::now();
        let finished = self.exhaustive_reductions(&priority_list);
        let duration = start_time.elapsed();
        let nodes_kill_count = nodes_before - self.graph.num_nodes();
        let edges_kill_count = edges_before - self.graph.num_edges();
        (duration.as_millis(), nodes_kill_count, edges_kill_count, finished)
    }

    // TODO:
    // `time_limit` in seconds
    pub fn overview_stats(&self, time_limit: u64) -> Result<Statistics, GraphError> {
        let now = Instant::now();
        todo!();

    }

}
