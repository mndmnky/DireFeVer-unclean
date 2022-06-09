# DireFeVerSS

## DIREcted FEedback VERtex Set Solver

A solver for the directed vertex feedback set problem. This library provides the following algorithms and data structures:

### Graph
* A directed, simple graph data structure that is implemented with two adjacency sets for each node.
* A extention of the graph datastructure that alows rebuilding after alteration.

### Kernelization
* Simple kernelization rules as partially described in [^fn2] that include:
	* Removing multi-edges.
	* Adding loop nodes to the solution and removing them from the graph.
	* Removing sink and source nodes.
	* Merging nodes with incomming, or outgoing degree of 1 into its sole in-(resp. out-)neighbor. 

* Strongly connected component rules, that include:
	* Removal of edges connecting strongly connected components.
	* Removal of edges that connect strongly connected components in a subgraph, where strong edges were removed. This rule is described in [^fn3].

* `k`-flower rules, that include:
	* Variations of the removal of nodes that are connected to `k`+1 node disjunct cycles. Where `k` denotes a upper bound of the DFVS instance. This rule is described in [^fn4].
	* Variations of the removal of nodes that are connected to `k`+1-`l` node disjunct cycles `C`. Where `k` denotes a upper bound of the DFVS instance, and `l` denotes a lower bound of the subgraph, after the deletion of all nodes in `C`.

* Core rules, that include:
	* Variations of the clique rule, as described in [^fn3], there called "core" rule.
	* Variations of our core rule, which is a generalized version of the clique rule.

* The dome rule, as described in [^fn3].

* Vertex Cover rules applied on subgraphs with reciprocal edges:
	* The link node rule, that contract nodes adjacent to and only to exactly two reciprocal edges, going to and coming from nodes `v` and `w`. As described for the vertex cover problem in [^fn5]. To adapt this rule for the directed feedback vertex set problem, one can not use cases where there exists a directed path from `v` to `w` since that might lead to a directed circle that was not present before.
	* The twin node rule, as described in [^fn6]. This rule looks for a pair of nodes with identical 3 reciprocal neighbors `v`, `w` and `u`, if there exists a reciprocal edge between any two nodes of `v`, `w` and `u`, those three nodes can be added to the cover.
	* The unconfined (here dubbed dominion) rule as described in [^fn7].
	* And the crown rule as described in [^fn8].


### Heuristics
* Lower bound heuristics:
	* A heuristic that counts and removes small cycles.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| - 1 to the lower bound. This heuristic is described in [^fn3].

* Upper bound heuristics:
	* Variations of the big degree heuristic.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| to the upper bound.

* Vertex-cover heuristic:
	* Transforms strongly connected components to a vertex cover instance by only regarding reciprocal edges as single undirected edges and discarding the rest. 
	* Then, we run a vertex cover solver over that instance. Here we use Duck and Cover [^fn9].
	* If the returned solution is also a solution for the original strongly connected component, then the solution is optimal.
	* Otherwise, we can use the solution as a lower bound, which can also be extended to an upper bound by solving the leftover graph.

### Exact Branching Algorithms
* A simple branching algorithm that should not currently be used.
* An advanced branching algorithm that first branches on cliques, then on daisies, and finally on small cycles.

### Parameterized Algorithm
* A solver that transforms the instance into a directed arc feedback set instance, then solves by iterative compression and transformation into skew edge multicut instances. As described by [^fn1]

## Changelog (at 1.0.0)

### 1.0.1 
* Speed up `apply_advanced_scc_rule()`.
* Added `apply_advanced_scc_rule()` to `BranchInstance`.
* Added `cai_weight()` and `lin_weight()` to `Digraph`.
* Added bottom-up and top-down weight heuristics.
* Added top-bottom-switch weight heuristics.
* Added local-search heuristics.
* Added clique upper bound branch heuristic.
* Simple statistics for the weight heuristics.
* Binary for statistics on the heuristics.

### 1.1.0
* Merge `DirectedFeedbackVertexSetInstance` and `BranchInstance`, only using `RebuildGraph`.
* Binary that records reduction statistics on multiple files.
* Iterative approach for scc- and advanced scc rule.
* Further speed up of the scc- and advanced scc rule.

### 1.2.0 
* Added `exhaustive_reductions()`.
* Added `apply_exhaustive_core_rule()`.
* Added `compute_and_set_lower()`.
* Added `apply_local_k_daisy()`.
* Added `compute_and_set_simple_upper_lower()`.
* Added `exhaustive_local_search()`.
* Fixed `mult_reduction_stats.rs`
* Struct interrupter:
	* Only handles time_outs
	* Implemented for exhaustive reductions
* Overhauled `advanced_branching()`.
* Put `exhaustive_reductions()` in each heuristic and branch and bound. 
* Added `current_best` to `dfvs_instance`.
* Proper bin for statistics on single instance heuristics.
* `advanced_clique_heuristic()` branches only to a maximal depth of 3, which could still be too much.
* Fixed a bug where `{RebuildGraph,DFVSInstance}::_reset_reductions()` did not reset `changes` (resp. `reductions`).
* Changes to `advanced_branching()`: 
	* now branches on highest cai weight (0.3) node if no more daisy or clique remains.
	* checks if the current lower bound (+ the size of the current solution) is greater or equal than the current best solution. If so, returns the best current solution. 
	* Handles sccs separately.

### 1.2.1 (before interrupt in recursive branching)
* Added branching on double paths in `advanced_branching()`.

### 1.2.2
* Complete interrupt for BST.
* binary for exact statistics

### 1.3.0
* Link node rule + restructuring of DFVSInstance.

### 1.3.1
* Use fxhash for all hashmaps and hashsets that use usize as a key.

### 1.3.2
* Crown rule

### 1.3.3 
* Twin rule
* Dominion rule

### 1.4.0
* Changes when computing and updating the upper and lower bound.
* Heuristic over vertex cover
* SCC split + vc heursitic in exact bin

### 1.4.1
* Via vertex_cover now also computes a upper bound

## TODO
* Final cleanup


[^fn1]: Chen, Jianer, et al. "A fixed-parameter algorithm for the directed feedback vertex set problem." Proceedings of the fortieth annual ACM symposium on Theory of computing. 2008.

[^fn2]: Levy, Hanoch, and David W. Low. "A contraction algorithm for finding small cycle cutsets." Journal of algorithms 9.4 (1988): 470-493.

[^fn3]: Lin, Hen-Ming, and Jing-Yang Jou. "On computing the minimum feedback vertex set of a directed graph by contraction operations." IEEE Transactions on computer-aided design of integrated circuits and systems 19.3 (2000): 295-307.

[^fn4]: Fleischer, Rudolf, Xi Wu, and Liwei Yuan. "Experimental study of FPT algorithms for the directed feedback vertex set problem." European Symposium on Algorithms. Springer, Berlin, Heidelberg, 2009.

[^fn5]: Chen, Jianer, Iyad A. Kanj, and Weijia Jia. "Vertex cover: further observations and further improvements." Journal of Algorithms 41.2 (2001): 280-301.

[^fn6]: Xiao, Mingyu, and Hiroshi Nagamochi. "Confining sets and avoiding bottleneck cases: A simple maximum independent set algorithm in degree-3 graphs." Theoretical Computer Science 469 (2013): 92-104.

[^fn7]: Akiba, Takuya, and Yoichi Iwata. "Branch-and-reduce exponential/FPT algorithms in practice: A case study of vertex cover." Theoretical Computer Science 609 (2016): 211-225.

[^fn8]: Abu-Khzam, Faisal N., et al. "Kernelization algorithms for the vertex cover problem." (2017).

[^fn9]: https://github.com/mndmnky/duck-and-cover
