# Pruned taxonomy

Just to make the size of the inputs more reasonable, this directory will contain
pruned-taxonomy.tre a version of the `taxonomy.tre` with pruned to remove any
taxon found only in the taxonomy.

Specifically, if *S* is the set of taxa that are mapped to the tips of the entire
set of inputs, then a taxon *t* will be pruned from the taxonomy if:

  1. *t* is not in *S*

  2. None of the ancestors of *t* are in *S*

  3. None of the descendants of *t* are in *S*

