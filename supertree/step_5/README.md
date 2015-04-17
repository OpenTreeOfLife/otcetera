# Expanding tips mapped to non-terminal taxa

Produced by running:

    otc-nonterminals-to-exemplars -estep_5 step_1/taxonomy.tre -fstep_4/step-5-phylo-args.txt

where `step_4/step-5-phylo-args.txt` is a file that make should create by prepending
`step_4/` to the filename of every tree file listed in `step_1/tree-ranking.txt`

## Contents after `make` runs
This directory should hold versions of each tree file. Each file should be:

  1. phylogenetically identical to the file of the same name in `step_4` if none of the 
tips were mapped to non-terminal, OR
  2. have each leaf mapped to a higher taxon replace by an expansion of that tip
to a polytomy of all of the terminal taxa that exemplify the higher taxon. These tips
will be attached to the parent node of the leaf that was replaced.