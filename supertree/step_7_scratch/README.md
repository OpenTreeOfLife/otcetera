# Step 7 scratch dir
Filled by make's use of the:

    otc-uncontested-decompose -estep_7_scratch/export-sub-temp step_6/taxonomy.tre -fstep_7/step-7-phylo-args.txt

command where `step_7/step-7-phylo-args.txt` is the filename of each tree with the `step_5/` prefix
added to it.
If that succeeds, make will also invoke:

    bash checksum-tree-files.sh step_7_scratch/export-sub-temp


## Expected outcome
Subproblems are written in the `export-sub-temp` directory here. 
`dumped-subproblem-ids.txt` will be written if the decompostition
succeeds. It contains the tree filenames of the subproblems trees (just their names, not their paths).

If that succeeds, then each subproblem is check summed with md5.
If that also succeeds, then `checksummed-subproblem-ids.txt` will be created.
It is just a copy of `dumped-subproblem-ids.txt` that serves as a sentinel for the successful 
execution of the checksum.

### `export-sub-temp` subdirectory contents
The `otc-uncontested-decompose` invocation produces a subproblem 
    file of newicks (called `ott123.tre` if the OTT Id was "123")
    for each uncontested taxon and a file of tree names.

The file of names will be called `ott123-tree-names.txt` if the 
    uncontested taxon was OTT 123.
This file will have a line for each newick in the `.tre` file.
Each line will hold either tree file name for each phylo tree that contributed
    to the subproblem or the word "TAXONOMY" if the tree is from the taxonomy.
The taxonomy is always last, and the order of the other trees is determined
    by their order of input.
That is determined ultimately, but the order of the tree identifiers in
    `step_1/tree-ranking.txt`, so the order in the subproblems reflects the tree rankings.

The `checksum-tree-files.sh` step create the `ott###.md5` files.
These will be used to quickly determine if a subproblem has changed 
from the previous invocation.

The idea is that the tree files will be moved to their their final location 
(the `step_7` subdirectory) only if they differ from content there
(based on diffing their checksums).
This should enable caching of solutions to subproblems.

