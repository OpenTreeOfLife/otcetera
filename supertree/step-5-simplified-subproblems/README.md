# reduced problems
The move-subproblems-if-differing.py script should populate this directory with 
simplified version of the subproblems.

Some simplifications:
    1. if a nontaxonomic input is a polytomy, it can be removed.
    2. if there exists an ID, x, which is always on the same side of 
        every rooted split (phylogenetic statement) as another ID, y; then
        we can cull x and reattach it later.
