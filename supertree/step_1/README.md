# Clean taxonomy and tree-ranking

The tree-ranking.txt file that MTH has been using can be obtained by

    $ wget http://phylo.bio.ku.edu/ot/tree-ranking.txt

**NOTE**: need to verify with the gcmdr folks that this the correct ranking!

Currently the taxonomy is obtained from Joseph Brown's repo (the tree comes from treemachine)

    $ git clone git@bitbucket.org:josephwb/synthesis_trees.git
    $ cd synthesis_trees/Source_info/
    $ tar xfvz Filtered_OTT_taxonomy.tre.tgz
    $ chmod -w Filtered_OTT_taxonomy.tre
    $ cd ../../step_1
    $ ln -s ../synthesis_trees/Source_info/Filtered_OTT_taxonomy.tre taxonomy.tre



Quoting from the synthesis_tree/Source_info/Notes.txt:

"The archive 'Filtered_OTT_taxonomy.tsv.tgz' contains the taxonomy that is used for all analyses; a newick version of the taxonomy is available in the archive 'Filtered_OTT_taxonomy.tre.tgz'. This corresponds to OTT version 2.8draft5. The 'Filtered' indicates that the taxonomy has been pruned of 'dubious' taxa, indicated by the following flags:

"major_rank_conflict","major_rank_conflict_direct","major_rank_conflict_inherited", "environmental","unclassified_inherited","unclassified_direct","viral","nootu","barren", "not_otu","incertae_sedis","incertae_sedis_direct","incertae_sedis_inherited", "extinct_inherited","extinct_direct","hidden","unclassified","tattered"
"
