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

If the children of a node `X` in tree 1 have the same descendant sets
as the children of a node `Y` in tree 2, then the the children 
of node `X` will be merged into the same slice as node `X`.
Similarly node `Y` will be merged into a slice with its children.
In this way contiguous sections of the tree that are isomorphic
with each other will be placed in a single slice.

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
