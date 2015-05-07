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
**Warning** the 2 Makfiles (`Makefile.exemplar` and `Makefile.synth-v3)
here clash because they reuse the same set of temporary
directories!.
You choose the one you want (via adding a `-f` flag to `make` calls
below to specify the Makefile or by renaming one to `Makefile`).
If you change to the other one, you need to clean the temporaries!

`Makefile.exemplar` exemplifies taxa mapped to non-terminal taxa.

`Makefile.synth-v3` is the version used to create the decompositions (on 29 Apr, 2014)
for the 3rd version of the synthetic tree.
That version of the decomposition does not replace higher taxa with exemplars.
The subproblems from that run 
are posted at http://phylo.bio.ku.edu/ot/subproblems-v3.tar.gz
for the time being.
That archive was created by renaming `step_7_scratch/export-sub-temp` and then tar+gzipping it.
The pruned taxonomy (copied from `step_6/taxonomy.tre`) is posted at http://phylo.bio.ku.edu/ot/synth-v3-pruned-taxonomy.tre

This pipeline is not working yet. Running `make` will ultimately exit with an error, 
but it should run long enough to produce the decomposition into subproblems

There is some cruft from a previous stab at the pipeline (see [here](./cruft/README.md)).
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
It should eventually crash with a message like:

    python move-subproblems-if-differing.py step_7_scratch/checksummed-subproblem-ids.txt step_7_scratch/export-sub-temp step_7 step_8 step_9
    Depends on some "bleeding edge" feature on the peyotl supertree branch
    make: *** [step_7/subproblem-ids.txt] Error 1

If you put that branch on peyotl on your PYTHONPATH, then you'll get a little further (but none of the difficult subproblems will be solved).
However, it should run long enough to produce:

  1. tree files in `step_5`. These should have the expanded tip version of each tree, so that there are no longer tips mapped to 
    non terminal taxa. See [step_5/README.md](./step_5/README.md)

  2. a pruned taxonomy in `step_6/taxonomy.tre`. See [step_6/README.md](./step_6/README.md)

  3. the subproblems in `step_7_scratch/export-sub-temp`. See [step_7_scratch/README.md](./step_7_scratch/README.md)

If you have the appropriate branch of `peyotl` installed, then you should get at least some of the 
subproblems moved from `step_7_scratch/export-sub-temp` to `step_7`, but that directory
is likely to be incomplete if the `move-subproblems-if-differing.py` step did not succeed.
Use `step_7_scratch/export-sub-temp` if you just want the subproblems!
