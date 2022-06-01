use crate::dfvs_instance::DFVSInstance;
use crate::digraph::RebuildDigraph;
use crate::statistics::RuleStats;
use fxhash::{FxHashSet};
use std::collections::HashSet;
use std::cmp::max;
use crate::interrupter::Interrupter;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rule {
    SimpleRules,
    Dome,
    Clique,
    Core,
    SCC,
    KFlower,
    LinkNode,
    Crown,
    TwinNodes,
    Dominion,
}

impl DFVSInstance {

    /// Applies some simple reduction rules. The rules are applied for each node and in no specific
    /// order until no more rules can be applied.
    /// Returns true if at least one reduction has been applied.
    ///
    /// The rules are: 
    /// Rule 0 (implicit): Remove merge multi edges to one. 
    /// Rule 1: Remove loop nodes and add them to the solution. 
    /// Rule 2: Remove sources and sinks.
    /// Rule 3.1: Replace node v with one incoming neighbor u and one outgoing neighbor w by
    /// the edge (u,w).
    /// Rule 3.2: Merge node v with only one incoming neighbor u, back into u.
    /// Rule 3.3: Merge node v with only one outgoing neighbor w, into w.
    pub fn apply_simple_rules(&mut self) -> bool {
        let mut changed = true;
        let mut rounds = 0;
        while changed {
            rounds +=1;
            changed = false;
            let nodes = self.graph.nodes().collect::<Vec<usize>>();
            for node in nodes {
                // Rule 1: Add nodes with loops to the solution, remove it:
                if self.graph.in_neighbors(node).as_ref().expect("`node` is in graph.nodes()").contains(&node) {
                    self.add_to_solution(node);
                    changed = true;
                // Rule 2: Remove sources and sinks:
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 0 || self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 0 {
                    self.graph.remove_node(node);
                    changed = true;
                // Rule 3.1: Replace node v with one incoming neighbor u and one outgoing neighbor w by
                // the edge (u,w):
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 1 && self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let src = self.graph.in_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("Indegree of `node` > 0");
                    let trg = self.graph.out_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("Outdegree of `node` > 0");
                    self.graph.replace(node, (src, trg));
                    changed = true;
                // Rule 3.2: Merge node v with only one incoming neighbor u, back into u.
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let into = self.graph.in_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("`node` has indegree > 0");
                    self.graph.merge_into_back(node, into);
                    changed = true;
                // Rule 3.3: Merge node v with only one outgoing neighbor w, into w.
                } else if self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let into = self.graph.out_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("`node` has outdegree > 0");
                    self.graph.merge_into_front(node, into);
                    changed = true;
                }
            }
            // Rule 0: Remove duplicate edges is implicit. 
        }
        rounds>1
    }

    /// Applies the `LinkNode`-rule that contracts nodes of a strict strong degree of 2 (and no
    /// other incident edges) and merges its neighbors `neighs[0]` and `neighs[1]`
    /// (or adds both of them to the solution if they are strongly connected. 
    /// If there exists a strictly weak path between `neighs[0]` and `neighs[1]`, without there
    /// being an edge connecting those two, this rule is not applicable.
    /// "Small clique rule" only needs to be applied if clique rule wasn't applied before this
    /// rules.
    pub fn apply_link_node_rules(&mut self) -> bool {
        let mut changed = false;
        'outer: loop {
            let nodes = self.graph.nodes().collect::<Vec<usize>>();
            for node in nodes {
                if let Some(neighs) = self.graph.graph.is_link_node(node) {
                    if self.graph.graph.strongly_connected(neighs[0], neighs[1]) {
                        // Clique rule has not been applied yet. 
                        self.add_to_solution(neighs[0]);
                        self.add_to_solution(neighs[1]);
                        changed = true;
                        continue 'outer
                    } else {
                        // if there does not exist a edge from neighs[0] to neighs[1]:
                        if !self.graph.in_neighbors(neighs[1]).as_ref().expect("`neighs[1]` exists")
                            .contains(&neighs[0]) {
                            //  check if a weak path from neighs[1] to neighs[0] exist. 
                            let mut weak_in_0 = self.graph.weak_in_neighbors(neighs[0]).expect("`neighs[0]` exists");
                            weak_in_0.remove(&neighs[1]);
                            let mut weak_out_1 = self.graph.weak_out_neighbors(neighs[1]).expect("`neighs[1]` exists");
                            weak_out_1.remove(&neighs[0]);
                            if !(weak_in_0.is_empty() || weak_out_1.is_empty()) {
                                if self.graph.graph.weak_path_exists_between(&weak_out_1, &weak_in_0) {
                                    continue
                                }
                            }
                        }
                        // if there does not exist a edge from neighs[1] to neighs[0]:
                        if !self.graph.in_neighbors(neighs[0]).as_ref().expect("`neighs[0]` exists")
                            .contains(&neighs[1]) {
                            //  check if a weak path from neighs[1] to neighs[0] exist. 
                            let mut weak_in_1 = self.graph.weak_in_neighbors(neighs[1]).expect("`neighs[1]` exists");
                            weak_in_1.remove(&neighs[0]);
                            let mut weak_out_0 = self.graph.weak_out_neighbors(neighs[0]).expect("`neighs[0]` exists");
                            weak_out_0.remove(&neighs[1]);
                            if !(weak_in_1.is_empty() || weak_out_0.is_empty()) {
                                if self.graph.graph.weak_path_exists_between(&weak_out_0, &weak_in_1) {
                                    continue
                                }
                            }
                        }
                        self.contract_link_node(node, &neighs);
                        changed = true;
                        continue 'outer
                    }
                }
            }
            break 'outer
        }
        changed
    }

    /// Traverses over all nodes with a strict strong degree of 3 (and no weak degree) and stores
    /// the neighbors `neighs` of each such nodes either in `connects` if any pair in `neighs` is
    /// connected by a strong edge, or in `un_connects` otherwise. 
    /// If for any set of `neighs` there already exists an entry in `connects`, adds all nodes in
    /// `neighs` to the solution and returns. If `neighs` are already in `un_connects`, stores
    /// `neighs` also in `connects`.
    pub fn apply_twin_nodes_rule(&mut self) -> bool {
        let mut connects = HashSet::new();
        let mut un_connects = HashSet::new();
        let nodes = self.graph.nodes().collect::<Vec<usize>>();
        for node in nodes {
            if let Some(neighs) = self.graph.graph.is_trip_node(node) {
                if connects.contains(&neighs) {
                    self.add_to_solution(neighs[0]);
                    self.add_to_solution(neighs[1]);
                    self.add_to_solution(neighs[2]);
                    return true
                }
                else if un_connects.contains(&neighs) {
                    connects.insert(neighs);
                } else {
                    if self.graph.graph.has_strong_edge(&neighs.iter().copied().collect()) {
                        connects.insert(neighs);
                    } else {
                        un_connects.insert(neighs);
                    }
                }
            }
        }
        false
    }
                   
    /// Finds all strongly connected components and removes the edges between them and single node
    /// components.
    ///
    /// Returns true if at least one edge or node has been removed.
    pub fn apply_scc_rule(&mut self) -> bool {
        let mut change = false;
        let sccs = self.graph.find_strongly_connected_components_iter();
        let matter_sccs: Vec<FxHashSet<usize>> = sccs.into_iter()
            .filter(|scc| {
                if scc.len() == 1 {
                    let elem = scc.iter().next().expect("`scc` holds one element");
                    self.graph.remove_node(*elem);
                    change = true;
                    false
                } else {
                    true
                }
            }).collect();
        for i in 0..matter_sccs.len() {
            for j in 0..matter_sccs.len(){
                if i == j {
                    continue
                }
                let edges_from_to: Vec<_> = self.graph.edges_from_to(&matter_sccs[i], &matter_sccs[j]);
                if !edges_from_to.is_empty() {
                    self.graph.remove_edges(edges_from_to);
                    change = true;
                }
            }
        }
        change
    }

    /// Finds all strongly connected components ignoring strong edges and removes the edges
    /// between.
    ///
    /// Returns true if at least one edge has been removed.
    pub fn apply_advanced_scc_rule(&mut self) -> bool {
        let sccs = self.graph.find_weak_strongly_connected_components_iter();
        let mut change = false;
        for i in 0..sccs.len() {
            for j in 0..sccs.len(){
                if i == j {
                    continue
                }
                let edges_from_to: Vec<_> = self.graph.weak_edges_from_to(&sccs[i], &sccs[j]);
                if !edges_from_to.is_empty() {
                    self.graph.remove_edges(edges_from_to);
                    change = true;
                }
            }
        }
        change
    }

    /// An (hopefully) efficient implementation of the local k-flower rule. 
    /// 
    /// Considers every daisy `D` with at least effective upper bound - effective lower bound + 2 patels
    /// until some node was added to the solution.
    /// Checks if `D` has at least effective upper bound - effective lower bound of G/D + 1 patels,
    /// if so, adds the daisy core to the solution.
    ///
    /// Returns `true` if some node was added to the solution.
    /// For now also returns `false` if execution was interrupted by `self.interrupter`.
    pub fn apply_local_k_daisy(&mut self) -> bool {
        if self.upper_bound.is_some() && self.lower_bound.is_some() {
            for node in self.graph.nodes().collect::<Vec<_>>() {
                let daisy_patels = self.graph.strong_neighbors(node).expect("`node` exists");
                if daisy_patels.len() < self.effective_upper_bound().expect("`self.upper_bound` is some") - 
                    self.effective_lower_bound().expect("`self.lower_bound` is some") + 2 {
                        continue;
                }
                let mut left_over = self.graph.clone();
                left_over.remove_node(node);
                left_over.remove_nodes(daisy_patels.clone());
                let mut left_over_ins = DFVSInstance::new(left_over, None, None);
                if !left_over_ins.compute_and_set_lower(false){
                    return false
                }
                if daisy_patels.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") - 
                    left_over_ins.effective_lower_bound().expect("`self.lower_bound` is some") {
                        self.add_to_solution(node);
                        return true;
                }
            }
        }
        false 
    }

    /// Looks for a node with more than `self.upper_bound` node disjunct cycles and removes it. This algorithm
    /// stops as soon as a node was found, or all nodes have been tested. 
    ///
    /// Returns true if a node was found and false if no node was found.
    #[deprecated(since="1.1.1", note="Does not care about interrupts")]
    pub fn apply_k_flower(&mut self) -> bool {
        if self.upper_bound.is_some() {
            for node in self.graph.nodes().collect::<Vec<usize>>() {
                let (cycles, _) = self.graph.smallest_first_heuristic_node(node).expect("`node` exists");
                if cycles.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") {
                    self.add_to_solution(node);
                    return true
                }
                // If the heuristic was close to `k` we also try a different heuristic: 
                if cycles.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") - 2 {
                    let (cycles2, _) = self.graph.small_degree_heuristic_node(node).expect("`node` exists");
                    if cycles2.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") {
                        self.add_to_solution(node);
                        return true
                    }
                }
            }
        }
        false
    }

    /// Checks if `node` has more than `self.upper_bound` minus some lower bound node disjunct
    /// cycles adjacent to it. If so, `node` is added to the solution and removed from the graph. 
    ///
    /// TODO: does this even consider effective lower of left over?
    #[deprecated(since="1.1.1", note="please use `apply_local_k_daisy` instead")]
    pub fn apply_local_k_flower_node(&mut self, node: usize) -> bool {
        if self.upper_bound.is_some() {
            let (cycles, left_over) = self.graph.smallest_first_heuristic_node(node).expect("`node` exists");
            let mut left_over_ins = DFVSInstance::new(RebuildDigraph::new(left_over), None, None);
            if left_over_ins.compute_and_set_lower(false){
                if cycles.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") - left_over_ins.effective_lower_bound().expect("is some") {
                    self.add_to_solution(node);
                    return true
                }
            }
        }
        false
    }

    /// Checks if the node with the highest minimum directed degree has more than `self.upper_bound` minus some lower bound node disjunct
    /// cycles adjacent to it. If so, `node` is added to the solution and removed from the graph. 
    ///
    /// TODO: does this even consider effective lower of left over?
    #[deprecated(since="1.1.1", note="please use `apply_local_k_daisy` instead")]
    pub fn apply_local_k_flower(&mut self) -> bool {
        if self.upper_bound.is_some() {
            if let Some(node) = self.graph.get_max_min_direct_degree_node() {
                let (cycles, left_over) = self.graph.smallest_first_heuristic_node(node).expect("`node` exists");
                let mut left_over_ins = DFVSInstance::new(RebuildDigraph::new(left_over), None, None);
                if left_over_ins.compute_and_set_lower(false){
                    if cycles.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") - left_over_ins.effective_lower_bound().expect("is_some") {
                        self.add_to_solution(node);
                        return true
                    }
                }
            }
        }
        false
    }

    pub fn apply_crown_rule(&mut self) -> bool {
        // find max matching under nodes that are not connected to weak edges 
        let strong_ones: FxHashSet<usize> = self.graph.graph.only_strong_nodes().collect();
        // Get outsiders 
        let (_, outsiders) = self.graph.graph.strong_max_matching_between(&strong_ones, &None);
        // find max matching between outsiders and N(outsiders) M2
        // Get I0 of unmatched nodes in outsiders 
        let n_out = self.graph.graph.open_out_neighbors_of_set(&outsiders);
        let (m2,mut spikes) = self.graph.graph.strong_max_matching_between(&outsiders, &Some(n_out));
        let mut s = 0;
        let mut head = FxHashSet::default();
        while s != spikes.len() {
        // repeat until spikes does not change 
            s = spikes.len();
            //  let H = N(I) 
            head = self.graph.graph.open_out_neighbors_of_set(&spikes);
            // TODO: could we here just check one side of each edge?
            // TODO: can be made more efficient
            spikes.extend(m2.iter().filter_map(|(a, b)| {
                    if head.contains(a) {
                        Some(*b)
                    } else if head.contains(b) {
                        Some(*a)
                    } else {
                        None
                    }
                }));
        }
        if !spikes.is_empty() {
            // add head to the solution.
            self.add_all_to_solution(head);
            // remove spikes 
            self.graph.remove_nodes(spikes);
            return true
        }
        false 
    }

    /// Looks for an unconfined vertex and adds it to the solution if one was find.
    /// Returns `true` if a vertex was added to the solution and `false` otherwise.
    pub fn apply_dominion_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<usize>>();
        for node in nodes {
            let mut set: FxHashSet<usize> = vec![node].into_iter().collect();
            // Only the strong neighbors matter here.
            let mut set_closed_n = self.graph.strong_neighbors(node).expect("`node` exists");
            set_closed_n.insert(node);
            loop {
                let set_closed_n_clone = set_closed_n.clone();
                let opt = set_closed_n_clone.iter().filter_map(|neigh| {
                    if set.contains(neigh) {
                        return None
                    }
                    // TODO This can probably be relaxed:
                    if self.graph.min_weak_degree(*neigh).expect("`neigh` exists") != 0 {
                        return None
                    }
                    let nn = self.graph.strong_neighbors(*neigh).expect("`neigh` exists");
                    if nn.intersection(&set).count() != 1 {
                        return None
                    }
                    Some(nn.difference(&set_closed_n).copied().collect::<FxHashSet<usize>>())
                }).min_by_key(|diff| diff.len());
                if let Some(diff) = opt {
                    if diff.is_empty() {
                        self.add_to_solution(node);
                        return true
                    } else if diff.len() == 1 {
                        let s_prime = diff.into_iter().next().expect("`diff.len()` == 1");
                        set.insert(s_prime);
                        set_closed_n.extend(self.graph.strong_neighbors(s_prime).expect("`s_prime` exists"));
                        set_closed_n.insert(s_prime);
                        continue
                    }
                }
                // TODO: do diamond instead of break 
                break
            }
        }
        false
    }

    /// Applies the clique rule: 
    /// Let `C` be an induced strongly connected clique and `n` a node of `C` with n incoming-
    /// or no outgoing neighbors that are not in `C`. Then all nodes of `C`, that are not `n` can be added
    /// to the solution then all nodes in `C` can be removed.
    ///
    /// Better use `apply_exhaustive_clique_rule()`
    pub fn apply_clique_rule(&mut self) -> bool {
        let mut greedy_clique = self.graph.greedy_max_clique();
        let mut found_node = None;
        for node in &greedy_clique {
            // no in or out degree?
            if self.graph.min_direct_degree_outside(*node, &greedy_clique).expect("`node` is part of the cliqu") == 0 {
                found_node = Some(*node);
                break;
            }
        }
        if let Some(node) = found_node {
            greedy_clique.remove(&node);
            self.graph.remove_node(node);
            greedy_clique.iter().for_each(|other_node| {
                self.add_to_solution(*other_node);
            });
            return true
        }
        false
    }

    /// Applies the clique rule (see documentation) by first finding a greedy 
    /// clique with `self.graph.greedy_max_clique()`, and then going through all the nodes of the
    /// graph, checking whether a node is found which forms a clique with part if the greedy
    /// clique, while having only incoming or outgoing neighbors outisde that clique: 
    ///
    /// Better use `apply_exhaustive_clique_rule()`
    pub fn apply_advanced_clique_rule(&mut self) -> bool {
        let greedy_cluster = self.graph.greedy_max_clique();
        let mut addable_cluster = FxHashSet::default();
        let mut removable_node = None;
        for node in self.graph.nodes() {
            if self.graph.min_direct_degree_outside(node, &greedy_cluster).expect("`node` is in `.nodes()`") == 0 {
                let into_sol = self.graph.strong_neighbors(node).expect("`node` is in `.nodes()`").clone();
                if !into_sol.is_empty() && self.graph.min_direct_degree_outside(node, &into_sol).expect("`node` exists") == 0 {
                    addable_cluster = into_sol;
                    removable_node = Some(node);
                    break;
                }
            }
        }
        if let Some(node) = removable_node {
            self.graph.remove_node(node);
            addable_cluster.iter().for_each(|other_node| {
                self.add_to_solution(*other_node);
            });
            return true
        }
        false
    }

    /// Applies the clique rule (see documentation) exhaustively by checking for each node
    /// `v` in `self.graph` if `v` has a minimum weak degree of 0. If so, checks if the
    /// strong neighborhood of `v` is a clique. If this is given as well, adds all strong neighbors
    /// of `v` to the solution and removes the neighborhood as well as `v`.
    /// Returns `true` if anything happened.
    pub fn apply_exhaustive_clique_rule(&mut self) -> bool {
        // Go through all the nodes 
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for node in nodes {
            if self.graph.min_weak_degree(node) != Some(0) {
                continue
            }
            // Get strong neighborhood 
            let strong_neighborhood = self.graph.strong_neighbors(node).expect("`node` exists");
            // Check if strong neighborhood is a cluster 
            if self.graph.is_cluster(&strong_neighborhood) {
                for neigh in strong_neighborhood {
                    self.add_to_solution(neigh);
                }
                self.graph.remove_node(node);
                return true
            }
        }
        false 
    }

    /// Applies the core rule (see documentation) by greedily finding a near maximal clique
    /// `clique`, then looking for the node `v` with the least minimum amount of direct neighbors `
    /// neighs` outside of the clique and then checking, if the `neighs` have a common clique
    /// core in `clique`. This core is then added to the solution and removed from the graph.
    /// Returns `true` if a core was added to the solution.
    pub fn apply_core_clique_rule(&mut self) -> bool {
        // Use heuristic to find max cluster 
        let mut greedy_cluster = self.graph.greedy_max_clique();
        // Pick `node` with lowest min directed degree (outside)
        if !greedy_cluster.is_empty(){
            // TODO: why do we need the clone here?
            let greedy_clone = greedy_cluster.clone();
            let (min_node, its_min_neighbors) = greedy_clone.iter()
                .map(|node| (node, self.graph.min_direct_neighbors_outside(*node, &greedy_cluster).expect("`node` exists")))
                .min_by_key(|(_, neighs)| neighs.len()).expect("`greedy_cluster` is not empty");
            // Remove `node` from cluster
            greedy_cluster.remove(min_node);
            // Intesect the strong neighbors of all min directed neighbors of `node` with max cluster. 
            for node in its_min_neighbors {
                greedy_cluster = self.graph.strong_neighbors_in(node, &greedy_cluster).expect("`node` exists");
            }
        }
        // Add intersection to solution.
        for node in &greedy_cluster {
            self.add_to_solution(*node);
        }
        !greedy_cluster.is_empty() 
    }

    /// Applies the core rule (see documentation) on the maximal daisy.
    /// Returns `true` if daisy core was added to the solution.
    pub fn apply_daisy_core_rule(&mut self) -> bool {
        // Get max daisy
        if let Some((daisy_core, daisy_leaves)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            // Create union
            let mut daisy = daisy_leaves.clone();
            daisy.insert(daisy_core);
            // find node in daisy that has a min_direct_degree_outside union of 0. 
            for node in daisy_leaves {
                if self.graph.min_direct_degree_outside(node, &daisy) == Some(0) {
                    self.add_to_solution(daisy_core);
                    return true
                }
            }
        }
        false
    }

    /// Applies the core rule on all possible daisies.
    /// Returns `true` if a core was added to the solution.
    pub fn apply_exhaustive_daisy_core_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for daisy_core in nodes {
            if let Some(daisy_leaves) = self.graph.strong_neighbors(daisy_core) {
                // Create union
                let mut daisy = daisy_leaves.clone();
                daisy.insert(daisy_core);
                // find node in daisy that has a min_direct_degree_outside union of 0. 
                for node in daisy_leaves {
                    if self.graph.min_direct_degree_outside(node, &daisy) == Some(0) {
                        self.add_to_solution(daisy_core);
                        return true
                    }
                }
            }
        }
        false
    }

    /// Applies the core rule on the node with the minimal amount of minimal directed neighbors and
    /// looks for a core under those neighbors.
    /// Returns true if something happened.
    pub fn apply_min_direct_core_rule(&mut self) -> bool {
        // Get min direct degree node
        if let Some((node, min_neighbors)) = self.graph.get_min_min_direct_degree_node_and_neighbors() {
            // Get strong neighbors of `node`:
            if let Some(strong_neighbors) = self.graph.strong_neighbors(node) {
                // find core of `min_neighbors` in `strong_neighbors`: 
                let mut core = strong_neighbors;
                for nigh in min_neighbors {
                    core = self.graph.strong_closed_neighbors_in(nigh, &core).expect("`node` exists");
                    if core.is_empty() {
                        return false
                    }
                }
                for cs in core {
                    self.add_to_solution(cs);
                }
                return true
            }
        }
        false
    }

    /// Applies the core rule for all the nodes in the graph, where we look for a core under the
    /// minimal directed neighbors of the node.
    /// Returns true if something happened.
    pub fn apply_exhaustive_min_direct_core_rule(&mut self) -> bool {
        // Get min direct degree node
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        'outer: for node in nodes {
            if let Some(min_neighbors) = self.graph.min_direct_neighbors(node) {
                // Get strong neighbors of `node`:
                if let Some(strong_neighbors) = self.graph.strong_neighbors(node) {
                    // find core of `min_neighbors` in `strong_neighbors`: 
                    let mut core = strong_neighbors.clone();
                    for nigh in min_neighbors {
                        core = self.graph.strong_closed_neighbors_in(nigh, &core).expect("`node` exists");
                        if core.is_empty() {
                            continue 'outer
                        }
                    }
                    for cs in core {
                        self.add_to_solution(cs);
                    }
                    return true
                }
            }
        }
        false
    }

    /// Applies the core rule for all the nodes in the graph, where those nodes are the fix point
    /// of the core.
    ///
    /// Returns true if something happened.
    pub fn apply_exhaustive_core_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for node in nodes {
            let strong_neighbors = self.graph.strong_neighbors(node).expect("`node` exists");
            let in_neighbors = self.graph.in_neighbors(node).as_ref().expect("`node` exists");
            let out_neighbors = self.graph.out_neighbors(node).as_ref().expect("`node` exists");
            let mut core = strong_neighbors.clone();
            for nigh in in_neighbors {
                core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
                if core.is_empty() {
                    break;
                }
            }
            if !core.is_empty() {
                for cn in core {
                    self.add_to_solution(cn);
                }
                return true
            }
            let mut core = strong_neighbors.clone();
            for nigh in out_neighbors {
                core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
                if core.is_empty() {
                    break;
                }
            }
            if !core.is_empty() {
                for cn in core {
                    self.add_to_solution(cn);
                }
                return true
            }
        }
        false
    }

    /// Applies the dome rule that removes edges that are not part of a minimal cycle.
    /// Returns true if something changed.
    pub fn apply_dome_rule(&mut self) -> bool {
        let weak_edges = self.graph.weak_edges().collect::<Vec<_>>();
        let mut change = false;
        for (src, trg) in weak_edges {
            let weak_ins = self.graph.weak_in_neighbors(src).expect("`src` exists");
            let ins = self.graph.in_neighbors(trg).as_ref().expect("`trg` exists");
            if weak_ins.is_subset(ins) {
                self.graph.remove_edge(&(src, trg));
                change = true;
            }
            let weak_outs = self.graph.weak_out_neighbors(trg).expect("`src` exists");
            let outs = self.graph.out_neighbors(src).as_ref().expect("`trg` exists");
            if weak_outs.is_subset(outs) {
                self.graph.remove_edge(&(src, trg));
                change = true;
            }
        }
        change
    }


    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// The priority order should roughly be chosen by the time consumption of the respective rules. 
    ///
    /// Simple rules have to be the first rules applied.
    ///
    /// Returns `true` if reductions finished before an interrupt and `false` otherwise.
    ///
    /// TODO: short k flower (only on a few nodes with high strong degree).
    pub fn exhaustive_reductions(&mut self, priority_list: &Vec<Rule>) -> bool {
        assert_eq!(priority_list[0], Rule::SimpleRules);
        'outer: while !self.check_interrupt() {
            for rule in priority_list {
                match rule {
                    Rule::SimpleRules => {
                        self.apply_simple_rules();
                    },
                    Rule::SCC => {
                        if self.apply_advanced_scc_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Dome => {
                        if self.apply_dome_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Clique => {
                        if self.apply_exhaustive_clique_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Core => {
                        if self.apply_exhaustive_core_rule() {
                            continue 'outer
                        }
                    },
                    Rule::KFlower => {
                        self.compute_and_set_upper_lower(true);
                        if self.apply_local_k_daisy() {
                            continue 'outer
                        }
                    },
                    Rule::LinkNode => {
                        if self.apply_link_node_rules() {
                            continue 'outer
                        }
                    },
                    Rule::Crown => {
                        if self.apply_crown_rule() {
                            continue 'outer
                        }
                    },
                    Rule::TwinNodes => {
                        if self.apply_twin_nodes_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Dominion => {
                        if self.apply_dominion_rule() {
                            continue 'outer
                        }
                    },
                }
            }
            return true
        }
        return false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::digraph::{Digraph, RebuildDigraph};
    use crate::dfvs_instance::DFVSInstance;

    #[test]
    fn rebuild_graph_test() {
        let gr = Cursor::new("8 12 0\n2 3\n\n2 5\n1\n2 7\n1 4 7\n8 1\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let check_graph = graph.clone();
        let mut branch_instance = DFVSInstance::new(graph, None, None);
        assert!(branch_instance.apply_simple_rules());
        assert_eq!(branch_instance.graph.num_nodes(), 0);
        assert_eq!(branch_instance.solution, vec![0usize].into_iter().collect::<FxHashSet<usize>>());
        assert!(branch_instance.redo_changes().is_ok());
        assert_eq!(branch_instance.graph, check_graph);
        assert!(branch_instance.solution.is_empty());
    }

    #[test]
    fn merge_node_test() {
        let gr = Cursor::new("6 8 0\n3\n3\n4 5\n5 6\n1\n2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let check_graph = graph.clone();
        let mut branch_instance = DFVSInstance::new(graph, None, None);
        assert!(branch_instance.apply_simple_rules());
        assert_eq!(branch_instance.graph.num_nodes(), 0);
        assert_eq!(branch_instance.solution, vec![2usize].into_iter().collect::<FxHashSet<usize>>());
        assert!(branch_instance.redo_changes().is_ok());
        assert_eq!(branch_instance.graph, check_graph);
        assert!(branch_instance.solution.is_empty());
        let gr_rev = Cursor::new("6 8 0\n5\n6\n1 2\n3\n3 4\n4\n");
        let g = Digraph::read_graph(gr_rev);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let check_graph = graph.clone();
        let mut branch_instance = DFVSInstance::new(graph, None, None);
        assert!(branch_instance.apply_simple_rules());
        assert_eq!(branch_instance.graph.num_nodes(), 0);
        assert_eq!(branch_instance.solution, vec![2usize].into_iter().collect::<FxHashSet<usize>>());
        assert!(branch_instance.redo_changes().is_ok());
        assert_eq!(branch_instance.graph, check_graph);
        assert!(branch_instance.solution.is_empty());
    }
    #[test]
    fn scc_rule_and_rebuild_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7 8 9\n8 9\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let mut graph_check = graph.clone();
        let graph_full_check = graph.clone();
        let mut instance = DFVSInstance::new(graph, None, None);
        instance.apply_scc_rule();
        graph_check.remove_node(0);
        graph_check.remove_node(5);
        graph_check.remove_node(6);
        assert_eq!(instance.graph, graph_check);
        assert!(instance.redo_changes().is_ok());
        assert_eq!(instance.graph, graph_full_check);
        let gr = Cursor::new("9 14 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7\n8\n9\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let mut graph_check = graph.clone();
        let graph_full_check = graph.clone();
        let mut instance = DFVSInstance::new(graph, None, None);
        instance.apply_scc_rule();
        graph_check.remove_node(0);
        graph_check.remove_edge(&(3,5));
        graph_check.remove_edge(&(3,6));
        graph_check.remove_edge(&(4,6));
        graph_check.remove_edge(&(4,5));
        assert_eq!(instance.graph.graph, graph_check.graph);
        assert!(instance.redo_changes().is_ok());
        assert_eq!(instance.graph, graph_full_check);
    }
    #[test]
    fn flower_rule_and_rebuild_test() {
        let gr = Cursor::new("7 9 0\n2 4 6\n3\n1\n5\n1\n7\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let check_graph = graph.clone();
        let mut branch_instance = DFVSInstance::new(graph, Some(2), None);
        assert!(branch_instance.apply_k_flower());
        assert_eq!(branch_instance.graph.num_nodes(), 6);
        assert!(branch_instance.apply_simple_rules());
        assert_eq!(branch_instance.graph.num_nodes(), 0);
        assert_eq!(branch_instance.solution, vec![0usize].into_iter().collect::<FxHashSet<usize>>());
        assert!(branch_instance.redo_changes().is_ok());
        assert_eq!(branch_instance.graph, check_graph);
        assert!(branch_instance.solution.is_empty());
    }
    #[test]
    fn local_flower_rule_and_rebuild_test() {
        let gr = Cursor::new("11 16 0\n2 4 6\n3\n1 8\n\
                            5\n1\n7\n1 11\n2 9\n10\n8\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = RebuildDigraph::new(g);
        let check_graph = graph.clone();
        let mut branch_instance = DFVSInstance::new(graph, Some(3), None);
        assert!(branch_instance.apply_local_k_flower_node(0));
        assert_eq!(branch_instance.graph.num_nodes(), 10);
        branch_instance.apply_simple_rules();
        assert_eq!(branch_instance.graph.num_nodes(), 0);
        assert_eq!(branch_instance.solution, vec![0usize, 7, 10].into_iter().collect::<FxHashSet<usize>>());
        assert!(branch_instance.redo_changes().is_ok());
        assert_eq!(branch_instance.graph, check_graph);
        assert!(branch_instance.solution.is_empty());
    }

    #[test]
    fn reduction_test() {
        let gr = Cursor::new("8 12 0\n2 3\n\n2 5\n1\n2 7\n1 4 7\n8 1\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        instance.apply_simple_rules();
        assert_eq!(instance.graph.num_nodes(), 0);
        assert_eq!(instance.solution, vec![0usize].into_iter().collect::<FxHashSet<usize>>());
    }


    #[test]
    fn scc_reduction_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7 8 9\n8 9\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut g_check = RebuildDigraph::new(g.clone());
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        instance.apply_scc_rule();
        g_check.remove_node(0);
        g_check.remove_node(5);
        g_check.remove_node(6);
        let instance_check = DFVSInstance::new(g_check, None, None);
        assert_eq!(instance, instance_check);
        let gr = Cursor::new("9 14 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7\n8\n9\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut g_check = RebuildDigraph::new(g.clone());
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        instance.apply_scc_rule();
        g_check.remove_node(0);
        g_check.remove_edge(&(3,5));
        g_check.remove_edge(&(3,6));
        g_check.remove_edge(&(4,5));
        g_check.remove_edge(&(4,6));
        let instance_check = DFVSInstance::new(g_check, None, None);
        assert_eq!(instance.graph.graph, instance_check.graph.graph);
    }

    #[test]
    fn flower_rule_test() {
        let gr = Cursor::new("7 9 0\n2 4 6\n3\n1\n5\n1\n7\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), Some(2), None);
        assert!(instance.apply_k_flower());
        assert_eq!(instance.graph.num_nodes(), 6);
        instance.apply_simple_rules();
        assert_eq!(instance.graph.num_nodes(), 0);
        assert_eq!(instance.solution, vec![0usize].into_iter().collect::<FxHashSet<usize>>());
    }

    #[test]
    fn local_flower_rule_test() {
        let gr = Cursor::new("11 16 0\n2 4 6\n3\n1 8\n\
                            5\n1\n7\n1 11\n2 9\n10\n8\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = RebuildDigraph::new(g.unwrap());
        let mut instance = DFVSInstance::new(g.clone(), Some(4), None);
        assert!(!instance.apply_local_k_flower_node(0));
        let mut instance = DFVSInstance::new(g, Some(3), None);
        assert!(instance.apply_local_k_flower_node(0));
        assert_eq!(instance.graph.num_nodes(), 10);
        instance.apply_simple_rules();
        assert_eq!(instance.graph.num_nodes(), 0);
        assert_eq!(instance.solution, vec![0usize, 7, 10].into_iter().collect::<FxHashSet<usize>>());
    }

    #[test]
    fn clique_rule_test() {
        let gr = Cursor::new("7 22 0\n2\n1 3 4 5\n2 4 5\n2 3 5 6 7\n2 3 4 7\n4 7\n4 6 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_clique_rule());
        assert_eq!(instance.solution, vec![1,3,4].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 3);
    }

    #[test]
    fn adv_clique_rule_test() {
        let gr = Cursor::new("7 26 0\n2 3 4 6\n1 3 4 7\n1 2 4 5 7\n1 2 3 5 6\n3 4\n1 4 5 7\n2 3 5 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_advanced_clique_rule());
        assert_eq!(instance.solution, vec![2,3].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);
    }

    #[test]
    fn exh_clique_rule_test() {
        let gr = Cursor::new("7 27 0\n2 3 4 5 6\n1 3 4 7\n1 2 4 7\n1 2 3 6\n6 7\n1 4 5 7\n2 3 5 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_clique_rule());
        assert_eq!(instance.solution, vec![5,6].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);
    }

    #[test]
    fn adv_scc_rule_test() {
        let gr = Cursor::new("6 13 0\n2 3\n3 5\n1 2 6\n3 5\n4 6\n4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_advanced_scc_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 11);
    }

    #[test]
    fn adv_scc_plus_dome_rule_test() {
        let gr = Cursor::new("6 13 0\n3 5\n1 3\n1 4 6\n2 5\n2 6\n4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(!instance.apply_advanced_scc_rule());
        assert!(instance.apply_dome_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 10);
    }

    #[test]
    fn core_clique_rule_test() {
        let gr = Cursor::new("8 21 0\n2 3 8\n1 3 5 6\n1 2 7 8\n1\n1 2\n2 7\n3 6\n3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_core_clique_rule());
        assert_eq!(instance.solution, vec![2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 7);
    }

    #[test]
    fn daisy_core_rule_test() {
        let gr = Cursor::new("5 13 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_daisy_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);
    }

    #[test]
    fn exh_daisy_core_rule_test() {
        let gr = Cursor::new("11 29 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n8 9 10 11\n8 9 10 11\n6 7\n6 7\n6 7\n6 7\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(!instance.apply_daisy_core_rule());
        assert!(instance.apply_exhaustive_daisy_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 10);
    }

    #[test]
    fn min_direct_core_rule_test() {
        let gr = Cursor::new("5 13 0\n2 4 5\n1 3 5\n1 2\n1 3 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_min_direct_core_rule());
        assert_eq!(instance.solution, vec![1].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);
    }

    #[test]
    fn exh_min_direct_core_rule_test() {
        let gr = Cursor::new("8 25 0\n2 3 4 5\n1 3 4 5\n1 2 4 5\n2 3 6\n2 3 8\n1 7\n1 6 8\n1 7");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(!instance.apply_min_direct_core_rule());
        assert!(instance.apply_exhaustive_min_direct_core_rule());
        assert_eq!(instance.solution, vec![1, 2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 6);
    }

    #[test]
    fn exhaustive_core_rule_test() {
        let gr = Cursor::new("8 21 0\n2 3 8\n1 3 5 6\n1 2 7 8\n1\n1 2\n2 7\n3 6\n3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 7);

        let gr = Cursor::new("5 13 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);

        let gr = Cursor::new("11 29 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n8 9 10 11\n8 9 10 11\n6 7\n6 7\n6 7\n6 7\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 10);

        let gr = Cursor::new("5 13 0\n2 4 5\n1 3 5\n1 2\n1 3 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![1].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);

        let gr = Cursor::new("8 25 0\n2 3 4 5\n1 3 4 5\n1 2 4 5\n2 3 6\n2 3 8\n1 7\n1 6 8\n1 7");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![1, 2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 6);
    }

    #[test]
    fn only_adv_scc_rule_test() {
        let gr = Cursor::new("8 17 0\n2 6\n3 7\n1 4 8\n7 5\n8 6\n4 1\n2 4\n3 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(!instance.apply_dome_rule());
        assert!(instance.apply_advanced_scc_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 16);
    }

    #[test]
    fn local_k_daisy_rule_test() {
        let gr = Cursor::new("6 8 0\n2 3 4\n1\n1\n1\n6\n5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), Some(3), Some(2));
        assert!(instance.apply_local_k_daisy());
        assert_eq!(instance.solution, vec![0].into_iter().collect());
    }

    #[test]
    fn link_node_rule_test() {
        let gr = Cursor::new("7 15 0\n2 7\n1 3\n2 4\n1 3 5\n4 6\n\
                             5 7\n1 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_link_node_rules());
        assert!(instance.redo_changes().is_ok());
        assert!(instance.apply_link_node_rules());
        instance.finallize_solution(); 
        assert_eq!(instance.solution.len(), 4);
    }

    #[test]
    fn crown_rule_test() {
        //let gr = Cursor::new("12 30 0\n4 5\n5 6\n6 7\n1 5 8 12\n1 2 4 6\n\
        //                     2 3 5 7\n3 6 9 11 12\n4 7 9\n7 8\n8\n10\n4 7\n");
        let gr = Cursor::new("6 9 0\n4\n4\n4\n1 2 3 6\n4\n5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_crown_rule());
        assert_eq!(instance.graph.num_nodes(), 2);
        assert_eq!(instance.solution.len(), 1);
    }

    #[test]
    fn advances_link_node_test() {
        let gr = Cursor::new("13 20 0\n2 3\n1\n1 6 7 8\n2\n2\n7\n\
                              12\n10\n5 10\n9\n4 13\n11\n12\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        let src = vec![5, 6, 7].into_iter().collect();
        let trg = vec![3, 4].into_iter().collect();
        assert!(instance.graph.graph.weak_path_exists_between(&src, &trg));
        assert!(!instance.apply_link_node_rules());
        instance.graph.add_edge((1, 2));
        assert!(instance.apply_link_node_rules());
        instance.apply_simple_rules();
        instance.finallize_solution(); 
        assert_eq!(instance.solution.len(), 3);
        assert!(instance.solution.contains(&0));
    }

    #[test]
    fn twin_node_rule_test() {
        let gr = Cursor::new("5 15 0\n3 4 5\n3 4 5\n1 2 5\n1 2 3\n1 2 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(!instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 0);
        let gr = Cursor::new("6 18 0\n3 4 5\n3 4 5\n1 2 6\n1 2 6\n1 2 6\n\
                              3 4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 3);
        let gr = Cursor::new("5 14 0\n3 4 5\n3 4 5\n1 2\n1 2 5\n1 2 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 3);
    }

    #[test]
    fn dominion_rule_test() {
        let gr = Cursor::new("14 47 0\n5 6 7\n3 4 7\n2 5 8\n2 5 7\n1 3 4 6 14\n\
                             1 5 7 14\n1 2 4 6 11\n3 9 11 12 13\n10\n8\n7 8 12 13\n\
                             8 11 14\n8 11 14\n5 6 12 13\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(RebuildDigraph::new(g), None, None);
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
    }

}
