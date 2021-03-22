# otc-explain-phylo-diffs
This tool takes two filepaths to newick tree files as arguments, 
and an optional `-O` flag to specify a starting stem of an output
filepath name for the JSON
representation of the differences between the phylogenetic trees.

The motivating use case is that each treefile has a collection of
trees inferred using different methods.
The tool conducts a quick (greedy) analysis of the branches that
can be cut convert the two rooted trees to forests of rooted trees
that are isomorphic to each other.

## Algorithm
### Division of the trees into slices.
To improve the performance of the tree comparisons, the trees are 
divided into slices by identifying nodes that are parents to the
same sets of leaves.

If the children of a node `X` in tree 1 have the same topology
as the children of a node `Y` in tree 2, then the the children 
of node `X` will be merged into the same slice as node `X`.
Similarly node `Y` will be merged into a slice with its children.
In this way contiguous sections of the tree that are isomorphic
with each other will be placed in a single slice.

After each tree has been sliced, the slices are compared between
trees.
The leaf set of a slice can include leaves from the full trees, or
the labels of internal nodes in the original trees.
However the slicing guarantees that the leaves that represent
subtrees contain the same set of taxa.

The `compress_tree_slices.py` script currently at https://gist.github.com/mtholder/bb4dc9c2ed82c570993e26e4e0e40392
can be used to further compress the output `otc-explain-phylo-diffs`
into fewer slices that have isomorphic trees for internal portions
of the trees.

### Comparison of trees within a slice
Treating each slice of the trees as a tree, `otc-explain-phylo-diffs`
computes the rooted triple distance between the two trees.
If this is 0, then the slices are identical.

If the distance is non-zero, the algorithm determines the set of 
leaves that have the highest fraction of conflicting triplets vs 
total triplet comparisons.
In the current implementation (March 22, 2021), every triplet of
taxa is considered, so "the highest fraction of conflicting triplets"
is equivalent to find the highest number of conflicting triplets.
Compatible polytomies are counted as comparisons, but not conflicts.
Arguably one could count them as differences. Or as a lack of a comparison.

One of the taxa with the highest differing triplet fractions
(as of the March 22, 2021 we simply use the leaf with the lowest
index) is then pruned from the trees.

This calculation of the triplet distances and pruning constitutes a "round"
in the "comp_pruning_rounds" object that summarizes the differences.
The process is iterated until the trees have no conflicting triplets.

### Summarization the incompatibility scores by leaves.
The tool also reports a taxon-specific statistic of the extent to 
which a taxon may be a rogue taxon, or may be a part of "rogue" clades.
This statistic has not been studied much, and has some obvious issues 
(e.g the tie-breaking can result in different sets of taxa getting high scores).
So, treat it as an untested, rough statistic of instability.

Each pruning operation can be thought of as a bisection of the tree by 
cutting a branch.
While the operations are always cutting terminal branches in the slice-view
of the trees.
Some tips in a slice correspond to larger clades in the full trees.

For each pruning operation performed in the algorithm, the program
calculates a score for both "sides" of the branch that is 1 divided 
by the number of taxa on that side of the branch in the full tree
That score is added to each of taxon. 
Thus if a pruned branch is  a terminal branch in a 100 taxon tree, 
the leaf pruned will have 1 added to its score, and the other 99 leaves
will have 1/99 added to their scores.

The sum for all prunings is reported in "fraction_pruning_incompat_score".
A normalized form of the statistic is in "norm_pruning_incompat_score" where
1 indicates that the taxon had the highest value out of all leaves in the
tree, and 0 indicates the taxon had the lowest value. 
The normalization is just `(f - min)/(max - min)` for any unnormalized 
stat `f`.

## Handling of multiple trees for each file.
The tool compares the trees on line #1 of each file, then if there
are more lines to each file, compares the trees from line #2, etc.

Each pair of trees must have the same leaf label set.

Each tree comparison produces a object describing the differences
represented as JSON.
If there are multiple lines in the inputs and `-Otag.json` is used
as the output specifier, then the JSON for the first comparison 
will be written in `tag.json`, the JSON for the second pair of trees
will be in `tag.json-1`, the JSON for the third pair will be in
`tag.json-2`, *etc*.



## Output schema
The schema of the JSON object comparing trees is an object holding
slices of the tree:

  * `"incompat_scores_by_leaf"` holds an array of objects (one per leaf) summarizing
  the index for the leaf in the tree slices, the label, and the taxon-specific
  incompatibility scores mentioned above.
  * `"root_id"` -> a string key in tree_comp_slices_by_root valid label is a SimpleNewick tree.
  * `"tot_num_prunings_all_slices"` -> an integer. The number of node prunings to needed to make each retained backbone in a slice identical.
  * `"tree_comp_slices_by_root"` -> an object with keys that are valid label for a SimpleNewick tree. The label indicated by root_id must be present as  a key. The values are objects. Each object is one of 2 types: "same tree" or "differing trees" types.

 A same tree has a `"both_trees"` key that contains a subtree slice object.

 A differing trees object has the keys:
  * `"tree_1"` a subtree slice object for the tree in the first treefile
  * `"tree_2"` a subtree slice object for the tree in the second treefile
  * `"comparison"` -> an object with:
      * `"num_prunings"` -> an integer # of prunings needed to make this slice topology identical
      * `"comp_pruning_rounds"` -> an array of triplet comparison rounds object. Each object has:
       * `"num_triplets_comp"` -> # of triplets compared in this round
       * `"num_triplets_diff"` -> # of triplets differing in this round
       * `"prop_triplets_differing"` -> float num_triplets_diff/num_triplets_comp if num_triplets_comp > 0
       * `"pruned"` ->(if not first round) the label pruned to start this round

Subtree slice objects have the keys:
  * `"newick"` -> string SimpleNewick form of the tree for this slice 
  * `"node_data"` -> object with keys for nodel labels in newick. each key maps to an object with metadata for the node.  For example, the user-displayable label for the node will be in `"label"`
