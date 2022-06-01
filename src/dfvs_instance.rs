use crate::digraph::{RebuildDigraph, Digraph};
use crate::reduction_rules::Rule;
use fxhash::FxHashSet;
use crate::cust_errors::{ImportError, GraphError, InterruptError};
use crate::interrupter::Interrupter;
use std::io::prelude::*;
use std::io;
use std::cmp::max;

#[derive(Debug, Default, Clone, Eq, PartialEq)]
pub struct DFVSInstance {
    pub graph: RebuildDigraph,
    pub solution: FxHashSet<usize>,
    pub upper_bound: Option<usize>,
    pub lower_bound: Option<usize>,
    pub current_best: Option<FxHashSet<usize>>,
    pub interrupter: Option<Interrupter>,
    pub merge_nodes: Vec<(usize, (usize, usize))>,
    register: Vec<usize>,
    changes: Vec<usize>,
}

impl DFVSInstance {

    /// Returns a new `DFVSInstance`
    ///
    /// # Arguments
    ///
    /// * `upper_bound` - If not specified (= None) the upper bound is computed by
    /// `DirectedFeedbackVertexSetInstance::high_degree_heuristic`.
    pub fn new(graph: RebuildDigraph, upper_bound: Option<usize>, lower_bound: Option<usize>) -> Self {
        DFVSInstance {
            graph,
            solution: FxHashSet::default(),
            upper_bound,
            lower_bound,
            current_best: None,
            interrupter: None,
            merge_nodes: Vec::new(),
            register: vec![0],
            changes: Vec::new(),
        }
    }

    /// Returns the effective upper bound of the current instance: 
    /// `self.upper_bound - self.solution.len()`
    pub fn effective_upper_bound(&self) -> Option<usize> {
        self.upper_bound.map(|upper| upper - self.solution.len())
    }

    /// Returns the effective lower bound of the current instance: 
    /// `self.lower_bound - self.solution.len()`
    pub fn effective_lower_bound(&self) -> Option<usize> {
        self.lower_bound.map(|lower| lower - self.solution.len())
    }
    
    /// Updates the old upper bound with `new_upper`. Or adds `new_upper` as upper bound if none
    /// exists.
    pub fn update_upper_bound(&mut self, new_upper: usize) {
        if let Some(upper) = self.upper_bound {
            if upper > new_upper {
                self.upper_bound = Some(new_upper);
            }
        } else {
            self.upper_bound = Some(new_upper);
        }
    }

    /// Updates the old lower bound with `new_lower`. Or adds `new_lowet` as lower bound if none
    /// exists.
    pub fn update_lower_bound(&mut self, new_lower: usize) -> bool {
        if let Some(lower) = self.lower_bound {
            if lower < new_lower {
                self.lower_bound = Some(new_lower);
                return true
            }
        } else {
            self.lower_bound = Some(new_lower);
            return true
        }
        false
    }

    /// Sets `set` to `self.current_best` and `self.upper_bound` to `set.len()` if `self.upper_bound` is 
    /// greater or equal to `set.len()` or `self.upper_bound` is not yet set. 
    /// Does nothing otherwise. 
    ///
    /// Returns `true` if `self.current_best` was updated and `false` otherwise.
    pub fn set_current_best(&mut self, set: &FxHashSet<usize>) -> bool {
        if let Some(upper) = self.upper_bound {
            if upper >= set.len() {
                self.current_best = Some(set.clone());
                self.upper_bound = Some(set.len());
                return true
            }
        } else {
            self.current_best = Some(set.clone());
            self.upper_bound = Some(set.len());
            return true
        }
        false
    }

    /// Computes and sets the best currently available upper- and lower bounds and the best current
    /// solution.
    /// Returns `true` if bounds were computed and set, returns `false` if an interrupt signal was send
    /// before any bounds could have been found.
    ///
    /// To access the best current solution use `self.current_best`.
    pub fn compute_and_set_upper_lower(&mut self, skip_initial_rules: bool) -> bool {
        if let Some((lower, _, set)) = self.get_best_bounds(skip_initial_rules) {
            self.set_current_best(&set);
            self.update_lower_bound(lower);
            return true
        }
        false
    }

    /// Computes and sets simple upper- and lower bounds.
    /// Returns the solution for the best upper bound, or `None` if an interrupt signal was send
    /// before any bounds could have been found.
    #[deprecated(since="1.1.1", note="There is currently no need to use this. Use `compute_and_set_upper_lower()` instead.")]
    pub fn compute_and_set_simple_upper_lower(&mut self) -> Option<FxHashSet<usize>> {
        if let Some((lower, upper, set)) = self.get_simple_bounds(){
            self.upper_bound = Some(upper);
            self.lower_bound = Some(lower);
            self.current_best = Some(set.clone());
            // TODO get rid of set return
            return Some(set)
        }
        None
    }

    /// Computes and sets a lower bound, if the computed lower bound is better then the old one. 
    /// Returns `false` if an interrupt signal was send. TODO: this could be changed. But entails
    /// some things at least in branch and bound
    pub fn compute_and_set_lower(&mut self, skip_initial_rules: bool) -> bool {
        let lower = self.get_some_lower(skip_initial_rules);
        if lower.is_some() {
            self.update_lower_bound(lower.expect("is some")); 
            return true
        } 
        false
    }

    pub fn _reset_changes(&mut self) {
        self.register = vec![0];
        self.changes = Vec::new();
        self.graph._reset_reductions();
    }

    /// This fuction is used to fix certain segments that can be recovered indipendently.
    pub fn start_new_changes(&mut self) {
        self.register.push(self.changes.len());
        self.graph.start_new_reduction();
    }

    /// Reverses the changes since the last time `start_new_changes()` was called.
    pub fn redo_changes(&mut self) -> Result<(), GraphError> {
        self.graph.rebuild_section()?;
        if self.register.is_empty() {
            return Err(GraphError::NothingToRebuildError)
        }
        let up_to = self.register.pop().expect("We checked register.");
        while self.changes.len() > up_to {
            let next_change = &self.changes.pop().expect("This can not be empty at this point");
            if next_change >= &self.graph.num_reserved_nodes() {
                if self.merge_nodes.len() != next_change - self.graph.num_reserved_nodes() + 1 {
                    return Err(GraphError::PopPaceholderOutOfOrder)
                }
                self.merge_nodes.pop();
            }
            self.solution.remove(next_change);
        }
        if self.register.is_empty() {
            self.register.push(0);
        }
        Ok(())
    }

    /// Reverses the changes completely.
    pub fn redo_changes_complete(&mut self) -> Result<(), GraphError> {
        self.graph.rebuild_complete()?;
        if self.register.is_empty() {
            return Err(GraphError::NothingToRebuildError)
        }
        while !self.register.is_empty() {
            let up_to = self.register.pop().expect("We checked register.");
            while self.changes.len() > up_to {
                let next_change = &self.changes.pop().expect("This can not be empty at this point");
                if next_change >= &self.graph.num_reserved_nodes() {
                    if self.merge_nodes.len() != next_change - self.graph.num_reserved_nodes() + 1 {
                        return Err(GraphError::PopPaceholderOutOfOrder)
                    }
                    self.merge_nodes.pop();
                }
                self.solution.remove(next_change);
            }
        }
        self.register.push(0);
        Ok(())
    }

    /// Adds `node` to the solution and removes it from the graph.
    pub fn add_to_solution(&mut self, node: usize) {
        self.solution.insert(node);
        self.changes.push(node);
        self.graph.add_node_to_solution(node);
    }

    /// Adds all nodes in `nodes` to the solution and removes it from the graph.
    pub fn add_all_to_solution(&mut self, nodes: FxHashSet<usize>) {
        for node in nodes {
            self.add_to_solution(node);
        }
    }

    /// Contracts `link` and merges `neighbors[0]` into `neighbors[1]`. Push
    /// `self.graph.num_reserved_nodes() + self.merge_nodes.len()` into the solution as a placeholder, and
    /// add information to `self.merge_nodes` to figure how the placeholder will be converted.
    pub fn contract_link_node(&mut self, link: usize, neighbors: &[usize; 2]) {
        self.graph.remove_node(link);
        self.graph.complete_merge(neighbors[0], neighbors[1]);
        let id = self.merge_nodes.len() + self.graph.num_reserved_nodes();
        self.solution.insert(id);
        self.changes.push(id);
        self.merge_nodes.push((neighbors[1], (neighbors[0], link)));
    }

    /// Replaces placeholder from `self.contract_link_node()` with the actual intended nodes.
    /// Can't currently be reversed.
    /// TODO: Make it not fuck up the reverse.
    pub fn finallize_solution(&mut self) {
        while !self.merge_nodes.is_empty() {
            let id = self.graph.num_reserved_nodes() + self.merge_nodes.len() - 1;
            self.solution.remove(&id);
            let rule = self.merge_nodes.pop().expect("`merge_copy` is not empty");
            if self.solution.contains(&rule.0) {
                self.solution.insert(rule.1.0);
            } else {
                self.solution.insert(rule.1.1);
            }
        }
    }

    /// Finalizes the solution as in `self.finallize_solution()` without actually changing the
    /// solution, so the actual instance can still be redone.
    pub fn finallize_solution_temp(&mut self) -> FxHashSet<usize> {
        let mut merge_copy = self.merge_nodes.clone();
        let mut solution_copy = self.solution.clone();
        while !merge_copy.is_empty() {
            let id = self.graph.num_reserved_nodes() + merge_copy.len() - 1;
            solution_copy.remove(&id);
            let rule = merge_copy.pop().expect("`merge_copy` is not empty");
            if solution_copy.contains(&rule.0) {
                solution_copy.insert(rule.1.0);
            } else {
                solution_copy.insert(rule.1.1);
            }
        }
        solution_copy
    }

    /// Finalizes a given solution `sol` as in `self.finallize_solution()` without actually changing
    /// `self.solution`, so the actual instance can still be redone.
    pub fn finallize_given_solution_temp(&self, sol: &FxHashSet<usize>) -> FxHashSet<usize> {
        let mut merge_copy = self.merge_nodes.clone();
        let mut solution_copy = sol.clone();
        while !merge_copy.is_empty() {
            let id = self.graph.num_reserved_nodes() + merge_copy.len() - 1;
            solution_copy.remove(&id);
            let rule = merge_copy.pop().expect("`merge_copy` is not empty");
            if solution_copy.contains(&rule.0) {
                solution_copy.insert(rule.1.0);
            } else {
                solution_copy.insert(rule.1.1);
            }
        }
        solution_copy
    }

    /// Finds the best upper- and lower bound under the currently implemented heuristics (or
    /// approximations).
    ///
    /// Returns the lower bound, the upper bound and the solution of the best upper bound.
    ///
    /// TODO: make more customizable if we have more time. 
    pub fn get_best_bounds(&self, skip_initial_rules: bool) -> Option<(usize, usize, FxHashSet<usize>)> {
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        let rule_priority = &vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
        let bounds = self.clique_heuristic(&rule_priority, skip_initial_rules);
        upper.push(bounds.clone().map(|(_,_,up)| up));
        lower.push(bounds.clone().map(|(lower,_,_)| lower));
        upper.push(self.top_down_weight_heuristic(&Digraph::cai_weight, (0.2,0f64), &rule_priority, skip_initial_rules));
        // TODO: this might be a waste of time. Find better lower bounds:
        //lower.push(self.disjunct_cycle_heuristic(&rule_priority, skip_initial_rules));
        let best_upper = upper.iter().filter(|ub| ub.is_some()).map(|ub| ub.clone().unwrap()).max_by_key(|ub| ub.len());
        if let Some(mut best_upper) = best_upper {
            if let Some(better) = self.exhaustive_local_search(&best_upper) {
                best_upper = better;
            }
            let best_lower = lower.iter().filter(|lb| lb.is_some()).map(|lb| lb.unwrap()).min();
            if let Some(best_lower) = best_lower {
                return Some((best_lower, best_upper.len(), best_upper)) 
            }
        }
        None
    }

    /// Returns a good upper bound for the given instance, or `None` if an interrupt was send.
    pub fn get_good_upper(&self, skip_initial_rules: bool) -> Option<FxHashSet<usize>> {
        let rule_priority = &vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
        if let Some(mut upper) = self.top_down_weight_heuristic(&Digraph::cai_weight, (0.2,0f64), &rule_priority, skip_initial_rules) {
            if let Some(better) = self.exhaustive_local_search(&upper) {
                upper = better;
            }
            return Some(upper)
        } else {
            return None
        }
    }

    /// Finds a simple/fast upper- and lower bound.
    /// Applied heuristics only use simple rules.
    ///
    /// Returns the best lower bound, upper bound and the solution of the best upper bound. If interrupted, either the best current values are
    /// returned, or `None` if none were found in time.
    #[deprecated(since="1.1.1", note="There is currently no need to use this. Use `get_best_bounds()` instead.")]
    pub fn get_simple_bounds(&self) -> Option<(usize, usize, FxHashSet<usize>)> {
        let simple_rules = vec![Rule::SimpleRules];
        if let Some((lower1, _, upper1_sol)) = self.clique_heuristic(&simple_rules, false){
            if let Some(mut upper2_sol) = self.bottom_up_weight_heuristic(&Digraph::cai_weight, (0f64,0.3), &simple_rules, false){
                if let Some(better) = self.exhaustive_local_search(&upper2_sol) {
                    upper2_sol = better;
                }
                let upper_sol = if upper2_sol.len()<upper1_sol.len() {upper2_sol} else {upper2_sol};
                return Some((lower1, upper_sol.len(), upper_sol))
            }
            return Some((lower1, upper1_sol.len(), upper1_sol))
        }
        None
    }

    /// Finds the best lower bound under the currently implemented heuristics (or
    /// approximations).
    /// Applied heuristics only use simple rules.
    /// Does not skip initial reductions.
    ///
    /// Returns the best lower bound found. If interrupted, either the best current lower bound is
    /// returned, or `None` if none was found in time.
    pub fn get_some_lower(&self, skip_initial_rules: bool) -> Option<usize> {
        if let Some((lower1, _, _)) = self.clique_heuristic(&vec![Rule::SimpleRules], false) {
            if let Some(lower2) = self.disjunct_cycle_heuristic(&vec![Rule::SimpleRules], skip_initial_rules) {
                return Some(max(lower2, lower1))
            }
            return Some(lower1)
        }
        None 
    }

    /// Sets a simple (sigterm) interrupter. TODO: not yet implemented.
    pub fn set_simple_interrupter(&mut self) {
        self.interrupter = Some(Interrupter::new(None));
    }

    /// Sets a time interrupt with `duration` milliseconds.
    pub fn set_time_interrupter(&mut self, duration: u128) {
        self.interrupter = Some(Interrupter::new(Some(duration)));
    }

    pub fn copy_interrupter(&mut self, other: &Self) {
        self.interrupter = other.interrupter.clone();
    }

    /// Check if interrupter is set. 
    /// In case there is no interrupter set, this should always return false.
    ///
    /// TODO: needs testing.
    pub fn check_interrupt(&self) -> bool {
        if self.interrupter.is_some() {
            return self.interrupter.as_ref().expect("is some").check_interrupt()
        }
        return false
        //self.interrupter.as_ref().unwrap().or_default().check_interrupt()
    }

    /// Checks if interrupter is set and throws an error in that case.
    pub fn send_interrupt(&self) -> Result<(),InterruptError> {
        if self.interrupter.is_some() {
            self.interrupter.as_ref().expect("is some").send_interrupt()?;
        }
        Ok(())
    }

}

impl DFVSInstance {

    /// Reads the solution of a DFVS instance from a `BufRead` type.
    /// Returns a HashSet of nodes in the solution.
    pub fn read_solution<R: BufRead>(sol: R) -> Result<FxHashSet<usize>, ImportError> {
        sol.lines()
            .map(|line| {
                line.unwrap().parse::<usize>().or(Err(ImportError::InputMalformedError))
            }).collect()
    }

    /// Writes a solution to a `Write` type.
    pub fn write_solution<W: Write>(solution: &FxHashSet<usize>, mut out: W) -> Result<(), io::Error> { 
        for elem in solution {
            writeln!(out, "{}",elem + 1)?;
        }
        Ok(())
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digraph::{Digraph, RebuildDigraph};
    use std::io::Cursor;
    use fxhash::{FxHashSet};
        
    // TODO: Test whether reduced and unreduced lead to the same result
    #[test]
    fn read_write_sol_test() {
        let gr = Cursor::new("5 13 0\n3 5\n3 4\n1 2 4 5\n1 3 5\n2 3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        let sol = "2\n4";
        let sol_curs = Cursor::new(sol);
        let read_sol = DFVSInstance::read_solution(sol_curs);
        assert!(read_sol.is_ok());
        assert_eq!(read_sol.unwrap(), vec![2,4].into_iter().collect::<FxHashSet<usize>>());
        let stdout = io::stdout();
        let stdout = stdout.lock();
        assert!(DFVSInstance::write_solution(&instance.solution, stdout).is_ok());
    }

}
