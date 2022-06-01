use std::collections::HashSet;
use fxhash::{FxHashSet};
use std::error;
use crate::cust_errors::{GraphError, InterruptError};
use crate::dfvs_instance::DFVSInstance;
use crate::reduction_rules::Rule;
use crate::digraph::{Digraph, RebuildDigraph};

impl DFVSInstance {

    /// Advanced branch and reduce algorithm that does the following: 
    /// 1. Exhaustive initial reduction and computation of upper and lower bounds.
    /// 1.1. The upper bound also computes an initial solution.
    /// 1.2. If the upper bound is equal to the lower bound, the initial solution is returned.
    /// 3. Priority branching:
    /// 3.1. If a cluster with more than 2 nodes is found, branch on the remaining nodes (the other
    ///   nodes are removed and added to the solution).
    /// 3.2. If a daisy with at least one petal exists, branch on the core and the petals of the
    ///   maximal daisy. 
    /// 3.3. Branch on a node with the hightest cai weight (0.3).
    /// See .advanced_branching_recursion() to see what happens in the invidual branches. 
    /// `self.changes` is reset at the start of the execution of this function.
    /// Returns the optimal solution if one was found or `None` if `self.interrupter` was set.
    ///
    /// 1.3.0 Note: With the link node rule branching on double paths has become redundent.
    ///
    /// # Arguments
    /// * `priority_rules` - A list of `Rule`s specifying which rule is applied during the
    /// execution and in which order they are applied.
    ///
    /// TODO make sure that you dont redo changes too far!!!
    pub fn advanced_branching(&mut self, priority_rules: &Vec<Rule>, skip_initial_rules: bool) -> Result<FxHashSet<usize>, Box<dyn error::Error>> {
        self._reset_changes();
        let mut best_solution;
        if !skip_initial_rules {
            // Exhaustively perform initial reductions w.o. upper bounds.
            if !self.exhaustive_reductions(priority_rules) {
                return Err(Box::new(InterruptError::Unclear))
            }
        }
        
        if self.graph.graph.scc_graph_is_disconnected() {
            let mut solution = self.solution.clone();
            // create subgraph for each (keep node ids)
            let subgraphs: Vec<Digraph> = self.graph.graph.split_into_connected_components_alt();
            for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(RebuildDigraph::new(graph), None, None)){
                instance.copy_interrupter(&self);
                // start advanced_branching on each component
                solution.extend(instance.advanced_branching(priority_rules, true)?);
            }
            // If we use the link node rule, we have to finallize the found solution.
            solution = self.finallize_given_solution_temp(&solution);
            return Ok(solution)
        }

        // If `self.upper_bound` is `None` while `self.lower_bound` is some, we recomputed the
        // lower bound never the less.
        if self.current_best.is_none() {
            // compute (simple?) upper bound (and lower bound)
            if !self.compute_and_set_upper_lower(true) {
                return Err(Box::new(InterruptError::Unclear))
            }
        }
        if self.lower_bound.is_none() {
            if !self.compute_and_set_lower(true) {
                return Err(Box::new(InterruptError::Unclear))
            }
        }
        best_solution = self.current_best.clone().expect("is some");
        // Return solution if best upper == best lower. 
        if self.upper_bound == self.lower_bound {
            // If we use the link node rule, we have to finallize the found solution.
            best_solution = self.finallize_given_solution_temp(&best_solution);
            return Ok(best_solution)
        }
        // Start recursive branching (one recursion of .advanced_branching_recursion)
        let greedy_cluster = self.graph.greedy_max_clique();
        if greedy_cluster.len() > 2 {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            // Start with the whole cluster:
            self.add_all_to_solution(greedy_cluster.clone());
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register and sets a new (second) one.
            self.redo_changes()?;
            self.start_new_changes();
            // Add all but one, contract the one:
            for node in &greedy_cluster {
                // Remove all but node. 
                self.add_all_to_solution(greedy_cluster.clone().iter().copied().filter(|clu| clu !=node).collect());
                self.graph.contract_node(*node);
                // Branch 
                if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                    // Return solution if best upper == best lower. 
                    // TODO: did we even change upper bound to reach this?
                    if self.upper_bound == self.lower_bound {
                        // Recovers up to very first register (pre reductions) and returns
                        self.redo_changes()?;
                        self.redo_changes()?;
                        return Ok(solution)
                    }
                    // Any found solution is better, than the old, due to the upper bound updates.
                    best_solution = solution;
                }
                // Recovers 2nd register and sets a new (second) one.
                self.redo_changes()?;
                self.start_new_changes();
            }
            // Recovers up to very first register (pre reductions) and returns
            self.redo_changes()?;
            self.redo_changes()?;
            return Ok(best_solution)
        // Link node branch, for link nodes n which the rule couldn't opperate. 
        } else if let Some((link_node, neighs)) = self.graph.graph.get_any_link_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_all_to_solution(neighs.iter().copied().collect());
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            self.redo_changes()?;
            self.add_to_solution(link_node);
            self.graph.contract_node(neighs[0]);
            self.graph.contract_node(neighs[1]);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers up to very first register (pre reductions) and returns
            self.redo_changes()?;
            return Ok(best_solution)
        } else if let Some((max_core, max_daisy)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_to_solution(max_core);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register.
            self.redo_changes()?;
            self.add_all_to_solution(max_daisy);
            self.graph.contract_node(max_core);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers up to very first register (pre reductions) and returns
            self.redo_changes()?;
            return Ok(best_solution)
        } else {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            let branch_node = self.graph.get_max_weight_node(&Digraph::cai_weight, (0.3, 0f64)).expect("`self.graph` is not empty");
            self.add_to_solution(branch_node);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register, contracts `node` and set a new 2nd register.
            self.redo_changes()?;
            self.graph.contract_node(branch_node);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recover to very first register.
            self.redo_changes()?;
        }
        Ok(best_solution)
    }

    // TODO: what reductions when?
    /// Branching subroutine of the advanced branch and reduce algorithm that does the following: 
    /// 1. Exhaustively apply simple reductions
    /// 2. Check if the graph is empty
    /// 2.1. If so, check if the current solution is below the upper bound.
    /// 2.1.1. If so, return the current solution and update the upper bound. 
    /// 2.1.2. If not, return `None`.
    /// 3. Priority branching:
    /// 3.1. If a cluster with more than 2 nodes is found, branch on the remaining nodes (the other
    ///   nodes are removed and added to the solution).
    /// 3.2. If a daisy with at least one petal exists, branch on the core and the petals of the
    ///   maximal daisy. 
    /// 3.3. Branch on a node with the hightest cai weight (0.3).
    ///
    /// 1.3.0 Note: With the link node rule branching on double paths has become redundent.
    pub fn advanced_branching_recursion(&mut self, priority_rules: &Vec<Rule>) -> Result<Option<FxHashSet<usize>>, Box<dyn error::Error>> {
        self.send_interrupt()?;
        let mut best_solution = None;
        // Perform some reductions.
        if !self.exhaustive_reductions(priority_rules) {
            return Err(Box::new(InterruptError::Unclear))
        }
        if self.graph.num_nodes() == 0 {
            // Check solution size
            if self.solution.len() < self.upper_bound.expect("upper bound was set") {
                // If we use the link node rule, we have to finallize the found solution.
                let sol = self.finallize_solution_temp();
                self.set_current_best(&sol);
                return Ok(Some(sol));
            } else {
                return Ok(None);
            }
        } else if self.solution.len() >= self.upper_bound.expect("upper bound was set") {
            return Ok(None);
        }
        // TODO: we could do a lower bound on all sccs and check the advanced lower bound criteria
        // we do in the next step.
        if self.graph.graph.scc_graph_is_disconnected() {
            let mut solution = self.solution.clone();
            // create subgraph for each (keep node ids)
            let subgraphs: Vec<Digraph> = self.graph.graph.split_into_connected_components_alt();
            for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(RebuildDigraph::new(graph), None, None)){
                instance.copy_interrupter(&self);
                // start advanced_branching on each component
                solution.extend(instance.advanced_branching(priority_rules, true)?)
            }
            if solution.len() < self.upper_bound.expect("a upper bound was set") {
                // If we use the link node rule, we have to finallize the found solution.
                solution = self.finallize_given_solution_temp(&solution);
                self.set_current_best(&solution);
                return Ok(Some(solution));
            } else {
                return Ok(None);
            }
        }
        // If the current lower bound is equal or greater than our current best solution, we won't find a better
        // solution in this branch.
        if let Some(lower) = self.lower_bound_clique_heuristic(priority_rules, true) {
            if lower >= self.upper_bound.expect("`upper_bound` was set") {
                return Ok(None);
            }
         } else {
             return Err(Box::new(InterruptError::Unclear));
         }
        // TODO (Consider): We could recompute upper bounds and check for matches with the lower
        // bound...
        let greedy_cluster = self.graph.greedy_max_clique();
        if greedy_cluster.len() > 2 {
            // In branch first register after reductions.
            self.start_new_changes();
            // Start with the whole cluster:
            self.add_all_to_solution(greedy_cluster.clone());
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recover first in branch register.
                    self.redo_changes()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recover first in branch register and set a new one.
            self.redo_changes()?;
            self.start_new_changes();
            // Add all but one, contract the one:
            for node in &greedy_cluster {
                // Remove all but node. 
                self.add_all_to_solution(greedy_cluster.clone().iter().copied().filter(|clu| clu != node).collect());
                self.graph.contract_node(*node);
                // Branch 
                if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                        // Return solution if best upper == best lower. 
                        if self.upper_bound == self.lower_bound {
                            // Recover first in branch register.
                            self.redo_changes()?;
                            return Ok(Some(solution))
                        }
                    // Any found solution is better, than the old, due to the upper bound updates.
                    best_solution = Some(solution);
                }
                // Recover first in branch register and set new one.
                self.redo_changes()?;
                self.start_new_changes();
            }
            // Recover first in branch register.
            self.redo_changes()?;
            return Ok(best_solution)
        // Link node branch, for link nodes n which the rule couldn't opperate. 
        } else if let Some((link_node, neighs)) = self.graph.graph.get_any_link_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_all_to_solution(neighs.iter().copied().collect());
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            self.redo_changes()?;
            self.add_to_solution(link_node);
            self.graph.contract_node(neighs[0]);
            self.graph.contract_node(neighs[1]);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recovers up to very first register (pre reductions) and returns
            return Ok(best_solution)
        } else if let Some((max_core, max_daisy)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            self.start_new_changes();
            self.add_to_solution(max_core);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // End of for loop to pop last register
                    self.redo_changes()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            self.redo_changes()?;
            self.add_all_to_solution(max_daisy);
            self.graph.contract_node(max_core);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            return Ok(best_solution)
        } else {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            let branch_node = self.graph.get_max_weight_node(&Digraph::cai_weight, (0.3, 0f64)).expect("`self.graph` is not empty");
            self.add_to_solution(branch_node);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.redo_changes()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recovers 2nd register, contracts `node` and set a new 2nd register.
            self.redo_changes()?;
            self.graph.contract_node(branch_node);
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recover to very first register.
        }
        Ok(best_solution)
    }

    /// Cant be trusted:
    #[deprecated(since="1.0.0", note="please use `advanced_branching()` instead")]
    pub fn branch_on_small_cycle(&mut self) -> Result<Option<FxHashSet<usize>>, GraphError> {
        // TODO: Perform pre reductions
        self.compute_and_set_upper_lower(false);
        self.branch_on_small_cycle_recursion()
    }

    #[deprecated(since="1.0.0", note="please use `advanced_branching()` instead")]
    fn branch_on_small_cycle_recursion(&mut self) -> Result<Option<FxHashSet<usize>>, GraphError> {
        // TODO Perform only local reductions:
        // Clique reduction takes to long to be applied at every step
        let mut changed = true;
        while changed {
            changed = self.apply_simple_rules();
            changed = self.apply_scc_rule() || changed;
            //changed = self.apply_k_flower() || changed;
            // k flower should not be used all the time
            if self.solution.len() > self.upper_bound.expect("upper bound was set") {
                return Ok(None)
            }
        }
        // Branch on all unmarked nodes of the cycle:
        //  1. Take first node into the solution and recurse
        //  2. Mark first node and check if first node has a cycle of marked nodes (and return None
        //  if thats the case)
        //  3. Repeat 2. for all nodes of the cycle.
        // TODO: find smallest cycle with nodes with largest degrees
        if let Some(cycle) = self.graph.find_smallest_cycle() {
            if self.solution.len() >= self.upper_bound.expect("upper bound was set") {
                return Ok(None)
            }
            let mut solutions = Vec::new();
            
            for node in cycle {
                self.start_new_changes();
                self.add_to_solution(node);
                if let Some(solution) = self.branch_on_small_cycle_recursion()? {
                    solutions.push(solution);
                }
                self.redo_changes()?;
            }
            // return best solution 
            if let Some(min) = solutions.iter().min_by_key(|list| list.len()) {
                Ok(Some(min.clone()))
            } else {
                Ok(None)
            }
        } else if self.solution.len() <= self.upper_bound.expect("upper bound was set") {
                self.update_upper_bound(self.solution.len());
                Ok(Some(self.solution.clone()))
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::digraph::{Digraph, RebuildDigraph};
    use crate::dfvs_instance::DFVSInstance;

    #[test]
    fn cycle_branch_test() {
        let gr = Cursor::new("5 13 0\n3 5\n3 4\n1 2 4 5\n1 3 5\n2 3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let mut branch_instance = DFVSInstance::new(graph, Some(5), None);
        let solution = branch_instance.branch_on_small_cycle();
        assert!(solution.is_ok());
        assert_eq!(solution.unwrap(), Some(vec![2,4].into_iter().collect::<FxHashSet<usize>>()));
    }

    #[test]
    fn advanced_branch_test() {
        let gr = Cursor::new("48 138 0\n2 3 4 5 12\n1 3 4 5 6 7 8\n1 2 4 5 14 15 16\n\
                             1 2 3 5 9 10 11\n1 2 3 4 9 10 11\n2 7\n2 8\n2 6 13\n\
                             4 5\n4 5\n4 5\n1 13\n6 12\n3 15\n3 16\n3 14\n18 19 20 21\n\
                             17 19\n17 20\n17 21\n17 18\n23 24 25\n22 24 26 28\n22 23 29 30 31\n\
                             22 27\n23 27\n23 28\n25 26\n24 30\n24 31\n24 29\n33 34 35 36\n\
                             32 37 38 39\n32 40 41 42\n32 43 44 45\n32 46 47 48\n\
                             33 38\n33 39\n33 37\n34 41\n34 42\n34 40\n35 44\n35 45\n\
                             35 43\n36 47\n36 48\n36 46\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let mut branch_instance = DFVSInstance::new(graph, None, None);
        let solution = branch_instance.advanced_branching(&vec![Rule::SimpleRules],false);
        assert!(solution.is_ok());
        assert_eq!(solution.unwrap().len(), 22); 
    }

}
