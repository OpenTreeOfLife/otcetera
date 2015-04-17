# make-based supertree pipeline

The initial thoughts for this are at https://github.com/OpenTreeOfLife/opentree/wiki/pseudo-code-supertree
but that is likely to be a bit out-of-sync as this is a work in progress.


# eventually...

The tree will be build-able by:

  1. downloading the files described in [input/README.md](./input/README.md)
  2. and putting the otcetera tools on your path with something like:
    `$ export PATH=$PATH:$PWD/../buildreleaseclang/installed/bin`
  3. and then building with:
    `$ make`

# Currently
This pipeline is not working yet. Running `make` will ultimately exit with an error, 
but it should run long enough to produce the decomposition into subproblems

There is some cruft from a previous stab at the pipeline (see below).
I'm in the process of refactoring to agree with the pipeline described in the doc subdir:
`../doc/summarizing-taxonomy-plus-trees.tex` (also occasionally posted to
    [this spot](phylo.bio.ku.edu/ot/summarizing-taxonomy-plus-trees.pdf) on the Holder lab site).
Each of the steps in section 2 of that doc are mapped to `step_#` directories in this new
system.  So `step_1` is described in section 2.1 of that doc; `step_2` is section 2.2; *etc*.

## Usage

  1. **Get the inputs** Currently you can get inputs from Joseph Brown's repos of treemachine output.
The instructions are in [./step_1/README.md](./step_1/README.md) and
[./step_4/README.md](./step_4/README.md)

  2. put the `otcetera` tools on your `PATH`. I like to:
    1. use the `--prefix="$PWD/installed"` as an argument to configure when I am building `otcetera`. 
    2. run `make && make check && make install` to put all of the tools in `./installed/bin`. 
    3. `cd ./installed/bin ; export PATH="${PATH}:${PWD}; cd -" to put the installed `bin` on my `PATH`.

  3. run `make` from the supertree directory.

## Expected outcome
It shold eventually crash with a message like:

    python move-subproblems-if-differing.py step_7_scratch/checksummed-subproblem-ids.txt step_7_scratch/export-sub-temp step_7 step_8 step_9
    Depends on some "bleeding edge" feature on the peyotl supertree branch
    make: *** [step_7/subproblem-ids.txt] Error 1

If you put that branch on peyotl on your PYTHONPATH, then you'll get a little further (but none of the difficult subproblems will be solved).
However, it should run long enough to produce:

  1. tree files in `step_5`. These should have the expanded tip version of each tree, so that there are no longer tips mapped to 
    non terminal taxa. See [step_5/README.md](./step_5/README.md)

  2. a pruned taxonomy in `step_6/taxonomy.tre`. See [step_6/README.md](./step_6/README.md)

  3. the subproblems in `step_7_scratch/export-sub-temp`. See [step_7_scratch/README.md](./step_7_scratch/README.md)

# Cruft docs from previous version of a pipeline

You need to prime the process by:

  1. getting the `taxonomy.tre` pruned by treemachine and placing it in `input`
  2. getting the trees with unmapped, dup, or nested tip mappings fixed by treemachine and placing them in `step-1-pruned-overlapping-tips`
  3. copying `step-1-pruned-overlapping-tips/tree-ranking.txt.example` to `step-1-pruned-overlapping-tips/tree-ranking.txt`

and then running steps 2 and 3 mentioned above in the "eventually..." section.