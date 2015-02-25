# make-based supertree pipeline

The initial thoughts for this are at https://github.com/OpenTreeOfLife/opentree/wiki/pseudo-code-supertree but that is likely to be a bit out-of-sync as this is a work in progress.


## usage

### eventually...

The tree will be buildable by:

  1. downloading the files described in [input/README.md](./input/README.md)

  2. and putting the otcetera tools on your path with something like:

    $ export PATH=$PATH:$PWD/../buildreleaseclang/installed/bin

  3. and then building with:

    $ make

### currently

You need to prime the process by:

  1. getting the `taxonomy.tre` pruned by treemachine and placing it in `input`
  2. getting the trees with unmapped, dup, or nested tip mappings fixed by treemachine and placing them in `step-1-pruned-overlapping-tips`
  3. copying `step-1-pruned-overlapping-tips/tree-ranking.txt.example` to `step-1-pruned-overlapping-tips/tree-ranking.txt`

and then running steps 2 and 3 mentioned above in the "eventually" section.