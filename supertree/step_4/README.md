# Step 4 cleaned input trees
Remove unmapped tips, and prune tips with repeted OTT IDs.

Currently obtained from Joseph Brown's repo (the tree comes from treemachine).
Assuming that you cloned this in [step_1](../step_1/README.md)

    $ cd synthesis_trees/Source_info
    $ tar xfvz Newicks_OTTID.tgz
    $ cd ../../step_4
    $ for f in $(cat ../step_1/tree-ranking.txt) ; do ln -s ../synthesis_trees/Source_info/Newicks_OTTID/$f . ; done


Quoting from the synthesis_tree/Source_info/Notes.txt:

"All 484 nexsons are in the archive 'Source_nexsons.tgz'. DO NOT download the studies anew; they may have been updated, and would not correspond exactly with other analyses.

Some source trees contain tips with duplicate taxa; these must be pruned to one exemplar each. Because alternative pruning approaches can lead to distinct trees (phylogenetic statements of relationship), we store here the versions used in all downstream analyses. The archive 'Newicks_OTTID.tgz' contain the 484 pruned trees in newick format. These trees have tips labelled with OTTIDs alone rather than taxon names. All OTTIDs correspond with the original curated tip labels; that is, none have been "mapped deeper" to more inclusive taxa in the context of the tree.

"
