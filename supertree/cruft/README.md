# Cruft docs from previous version of a pipeline

You need to prime the process by:

  1. getting the `taxonomy.tre` pruned by treemachine and placing it in `input`
  2. getting the trees with unmapped, dup, or nested tip mappings fixed by treemachine and placing them in `step-1-pruned-overlapping-tips`
  3. copying `step-1-pruned-overlapping-tips/tree-ranking.txt.example` to `step-1-pruned-overlapping-tips/tree-ranking.txt`

and then running steps 2 and 3 mentioned above in the "eventually..." section.