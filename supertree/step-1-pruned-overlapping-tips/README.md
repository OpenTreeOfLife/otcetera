# Pruned overlapping names.

This directory should contain:

  1. `tree-ranking.txt` - a file with each line being the name
        of a newick tree also found in this directory.
  2. each input tree in newick format (with filenames corresponding
        to the names listed in `tree-ranking.txt`). **TODO** internal
        node labels should probably be replaced with the nodeID in the 
        NexSON so that we can later tie support in the synthetic tree
        to particular edges in the input trees.
  3. some provenance in a JSON file for each tree describing what
        leaves were pruned and why.

It should be populated by some script that prunes 
tips such that to establish the guarantees listed below.

**TODO** currently MTH is just putting the newick's from treemachine 
    obtained from JWB in here to bootstrap the process of building the
    supertree.

## guarantees

  1. Each line of `tree-ranking.txt` corresponds to a tree file in this directory.

  2. For each tree:

    1. Unmapped tips have been pruned.

    2. No two tips are mapped to the same OTT id,

    3. No nesting - no tip is mapped to an OTT taxon that is a descendant of a taxon
        that is associated with another tip in this tree.