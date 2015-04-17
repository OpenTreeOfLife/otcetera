# Step 6 pruning of the input taxonomy to tips in phylogenetic statements

Currently produced by `otc-nonterminals-to-exemplars` used in step 5. See See [step_5/README.md](../step_5/README.md)

## Contents after `make` runs
This directory should hold `taxonomy.tre` which is the induced tree of `step_1/taxonomy.tre` when
restricted to the set of taxa that are found in the leaves of the trees in `step_5`. 
These are the "phylo-taxa".
All pruned taxa are only found in the taxonomy from `step_1`, so they can be added to the supertree
at the end of the pipeline.