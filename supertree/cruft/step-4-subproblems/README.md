The contents are to be filled by running 
otcscaffoldedsupertree with the taxonomy, input trees and the -e 
option with the ARG for -e being the path to this directory.

Each contested taxon is collapsed, but rather than trying to 
find the solution for each subproblem, that invocation will
cause a ott###.tre file to be written with a tree per line 
for every tree that is relevant to this subproblem.

The trees have ott IDs that have been created by pruning each 
tree where it intersects with an uncontested taxon. Thus each tree
is pruned down to include only the set of child taxa that must be
present in the subproblem.

For each subproblem, an ott####-tree-names.txt file has also 
been written. For the tar archive, I move those files to a subdir
called tree-names.

After these subproblems are solved, we should be able to piece them
back together, collapse and unsupported edges, and then add 
all of the "taxonomy only" nodes to get a full tree.

