# otcetera - phylogenetic file format parser in C++
[![Build Status](https://secure.travis-ci.org/OpenTreeOfLife/otcetera.png)](http://travis-ci.org/OpenTreeOfLife/otcetera)

otcetera owes a lot of code and ideas to Paul Lewis' Nexus Class Library.
  See http://hydrodictyon.eeb.uconn.edu/ncl/ and
  https://github.com/mtholder/ncl

It also uses easyloggingpp which is distributed under an MIT License. See
  http://github.com/easylogging/ for info on that project. The file from
  that project is otc/easylogging++.h

Some set comparisons (in util.h) were based on
   http://stackoverflow.com/posts/1964252/revisions
by http://stackoverflow.com/users/127669/graphics-noob

The gitversion trick for the otc-version-reporter is from
  http://stackoverflow.com/questions/6526451/how-to-include-git-commit-number-into-a-c-executable

https://peerj.com/preprints/2538/ describes some of the tools that are a part of otcetera.

# Installation

The instructions below contain all of the gory detail. There are a few quirks with OS X installation. See [Short OSX instructions](#short-osx-instructions) for an overview of the process on OS X. 

## prerequisites

### autotools
You also need the a fairly recent version of the whole autotools stack
including libtool. MTH had problems with automake 1.10. If you can't install
these with something like apt, then you can grab the sources. The following
worked for MTH on Mac on 28-Feb-2015:

    wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
    tar xfvz autoconf-2.69.tar.gz
    cd autoconf-2.69
    ./configure
    make
    sudo make install
    cd ..
    wget http://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz
    tar xfvz automake-1.15.tar.gz
    cd automake-1.15
    ./configure
    make
    sudo make install


### BOOST C++ libraries
You also need the BOOST C++ source libraries.  You should install the BOOST libraries and
header files using your operating system's package manager. version 1.58 of boost works; earlier versions might.

Many BOOST modules are header-only.  However, some modules require linking to an installed
library archive.  You must install library archives for at least these BOOST libraries:

    cd "${BOOST_ROOT}"
    b2 program_options
    b2 system
    b2 filesystem

If you do not do a full installation of BOOST, then you will need to add the
libraries to your dynamic library loading path. On Mac:

    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$BOOST_ROOT/bin.v2/libs/program_options/build/darwin-4.2.1/release/threading-multi/:$BOOST_ROOT/bin.v2/libs/system/build/darwin-4.2.1/release/threading-multi:$BOOST_ROOT/bin.v2/libs/filesystem/build/darwin-4.2.1/release/threading-multi"

on most other unix variants the variable is called `LD_LIBRARY_PATH` (without the `DY`
prefix); also note that the path to the libraries in the build BOOST dir is platform-dependent.


# requests

The python requests package is need for running the `make check` target because it runs tests in the `ws` subdirectory.
# restbed

We are using the [Restbed framework](https://github.com/corvusoft/restbed) to implement web services for the tree of life. This is work in progress. By default, otcetera will NOT include restbed unless you run configure with the `--with-webservices=yes` option.

On debian or ubuntu:

    sudo apt-get install g++-5
    sudo apt-get install cmake
    sudo apt-get install libtool
    sudo apt-get install libboost-all-dev

On a mac:

    sudo brew install cmake
    sudo brew install openssl
    SSL=/usr/local/opt/openssl
    export CPPFLAGS="$CPPFLAGS -I${SSL}/include"
    export LDFLAGS="$LDFLAGS -L${SSL}/lib"

then:

    git clone --recursive https://github.com/corvusoft/restbed.git
    mkdir restbed/build
    cd restbed/build/
    cmake -DBUILD_SSL=NO -DCMAKE_INSTALL_PREFIX=$PWD/install ..
    make install
    export CPPFLAGS="$CPPFLAGS -I$PWD/install/include "
    export LDFLAGS="$LDFLAGS -I$PWD/install/library "

To build including restbed, you will also need to set:

    CPPFLAGS=-Ipath_to_restbed_include
    LDFLAGS=-Lpath_to_restbed_library

On a Mac, you will need to set:

    CPPFLAGS=-Ipath_to_homebrew_openssl_include -Ipath_to_restbed_include 
    LDFLAGS=-Lpath_to_homebrew_openssl_include -Lpath_to_restbed_library

These variables should be set by the lines above that begin with `export`.

## configuration + building

To run the whole autoreconf stuff in a manner that will add missing bits as needed,
run:

    $ sh bootstrap.sh

Then to configure and build with clang use:

    $ mkdir buildclang
    $ cd buildclang
    $ bash ../reconf-clang.sh
    $ make
    $ make check
    $ make install
    $ make installcheck

To use g++, substitute `reconf-gcc.sh` for `reconf-clang.sh` in that work flow.
g++ version 5.4 works earlier versions might work.

Python 2 (recent enough to have the subprocess module as part of the standard lib)
is required for the `make check` operation to succeed.

Those `reconf-...sh` scripts set the installation prefix to an `installed` sub-directory
of your build directory as the prefix for the installation. So, you will need to
add `$PWD/installed/bin` to your `PATH` environmental variable to use the version
of the otc tools that you just installed. Depending on your platform, you may have to
add the `$PWD/installed/lib` to your `LD_LIBRARY_PATH` variable.

# Short OSX instructions

Tested by kcranston on OS X Sierra 10.12.4. These instructions do not include building the web services. Assumes you have the [Homebrew package manager](https://brew.sh/) installed. 

    $ brew install autoconf
    $ brew install automake
    $ brew install boost
    $ brew install libtool
    $ brew install md5sha1sum
    
When you install libtool, you will get the following warning:

```
In order to prevent conflicts with Apple's own libtool we have prepended a "g"
so, you have instead: glibtool and glibtoolize.
```

so on OS X the bootstrap script calls the `g*` versions.

    $ sh bootstrap.sh
    $ mkdir build
    $ cd build
    $ bash ../reconf-clang.sh
    $ make
    $ make check
    $ make install
    $ make installcheck

Finally, add the `bin` directory to your $PATH:

    $ export PATH=$PATH:$PWD/installed/bin
    

# Documentation
A LaTeX documentation file is [./doc/summarizing-taxonomy-plus-trees.tex](./doc/summarizing-taxonomy-plus-trees.tex)
periodically, that is compiled and posted.
The currently URL for that compiled documentation is http://phylo.bio.ku.edu/ot/summarizing-taxonomy-plus-trees.pdf

# Usage
See the [supertree/README.md](./supertree/README.md) for instructions on using
`otcetera` to build a supertree (work in progress).

## Common command line flags
The tools use the same (OTCLI) class to process command line arguments.
This provides the following command line flags:
  * `-h` for help
  * `-fFILE` to treat every line of FILE as if it were a command line argument
      (useful for processing hundreds of filenames)
  * `-v` for verbose output
  * `-q` for quieter than normal output
  * `-t` for trace level (extremely verbose) output

Unless otherwise stated:

  1. the command line tools that need a tree take a filepath
to a newick tree file. The numeric suffix of each label in the tree is taken to
be the OTT id. This accommodates the name munging that some of the open tree of
life tools perform on taxonomic names with special characters (because only the
OTT id is used to associate labels in different trees)
  2. a full supertree tree and taxonomy tree have the same leaf set in terms of OTT ids

## Config file

### The `~/.opentree` file
You may optionally initialize the global config file `~/.opentree` to specify
the location of the OpenTree Taxonomy (OTT).  Currently the only use of this
file in otcetera is to avoid specifying the taxonomy argument on the command-line
to a few commands.

The config file should contain the [opentree] section with a definition for the variable
`ott`:

    [opentree]
    home = /home/USER/OpenTree
    ...
    ott = %(home)s/ott/ott2.9draft12/
    ...

You can optionally define a variable such as `home` to point to the parent directory.
Then you can reference that directory by writing `%(home)s` in other variables in the same section.

## Tools for checking a supertree against inputs

### Checking for nodes in a supertree that have no tree supporting them

    otc-check-supertree taxonomy.tre synth.tre inp1.tre inp2.tre ...

will report any nodes in `synth.tre` that are not named and which do not have
any [ITEB support](https://github.com/OpenTreeOfLife/treemachine/blob/nonsense-1/nonsense/iteb_support_theorem.md)

The taxonomy is just used for the ottID validation (on the assumption that
the nodes supported by the taxonomy and the the otcchecktaxonomicnodes tool
can help identify problems with those nodes).

Note that this check identifies a nodes that could be collapsed to produce
  a "minimal" tree (*sensu* [Semple, 2003](http://ir.canterbury.ac.nz/handle/10092/1652)).
As discussed in section 1.3 "Trees without unsupported groups" of
[the docs](http://phylo.bio.ku.edu/ot/summarizing-taxonomy-plus-trees.pdf), if you
want to get rid of such groups, then you should remove them one at a time and rerun
the check. Otherwise you may cause a supertree to display fewer of the input clusters.


### Checking for incorrect internal labels in a full tree

If you add a `-x` argument to the invocation, then the program will act like the taxonomy
is also a source of support for nodes.  Furthermore a report on problems with
taxonomic labels in synth.tre will be reported before the summary.

Using both `-x` and `-r` will create the taxonomic report, but then clear the taxonomic
support stats before analyzing the subsequent inputs. So the final summary should
be equivalent to what you get by dropping both the `-x` and `-r` args.

    otc-check-supertree -x -d taxonomy.tre synth.tre

will check every labelled internal node is correctly labelled. To do this, it
verifies that the set of OTT ids associated with tips that descend from the
node is identical to the set of OTT ids associated with terminal taxa below
the corresponding node in the taxonomic tree.

`otc-taxon-conflict-report` takes at least 2 newick file paths: a full tree, and some number of input trees.
It will write a summary of the difference in taxonomic inclusion for nodes that are in conflict:

    otc-taxon-conflict-report taxonomy.tre inp1.tre inp2.tre

#### Removal of `otc-find-unsupported-nodes` and `otc-check-taxonomic-nodes`
 The functionality that was previously in `otc-find-unsupported-nodes`
and `otc-check-taxonomic-nodes` is now implemented in `otc-check-supertree`.
This new tool has the same interface as `otc-find-unsupported-nodes` but the `-d`
option from `otc-check-taxonomic-nodes` was also added. Thus `otc-check-taxonomic-nodes`
is no longer necessary, and the name of the tool was changed to reflect its
broader set of checks.

### Counting input tree groups displayed by a full tree
The `otc-displayed-stats` analyzes the nodes of the inputs in the context of a summary tree.

    otc-displayed-stats -x taxonomy.tre synth.tre inp1.tre inp2.tre ...

writes tab-separated output. The first column is the number of non-redundant input
tree nodes that are displayed by the `synth.tre`. This vector of groupings displayed directly
corresponds to "Weighted Input Phylogenetic Statements Displayed" described in the documentation.
If you were to assign each tree a weight (based on its rank), you could then calculate
a score by multiplying the tree weight by the number displayed in the first column of the
output.

Or you could simply view the score of the `synth.tre` to be a vector of numbers that corresponds
to this first column of output. The goal of they Open Tree of Life supertree operation is
to maximize this score (in a lexicographic ordering with the first tree being the most significant)
while introducing no unsupported groups.

The `-x` flag above tells the tool to treat the taxonomy as the last input. (if this is lacking
the taxonomy is only used for validation of OTT identifiers in the other trees.

#### Explanation of the `otc-displayed-stats` output
The full output of the `otc-displayed-stats` is explained below.
Each row of the output reports the number of internal nodes of the input tree that
fall into each category. The two "axes" that the statistics explore are support and out-degree.

Columns starting with "F" are "forking" internal nodes with out-degree > 1.
Columns starting with "R" are "redundant" internal nodes with out-degree = 1.
A "D" suffix to a column header means that the node is displayed by the summary tree.
A "CR" suffix means that the node is could resolve a polytomy in the summary tree (so the
    summary tree is not unambiguously in conflict in the node).
An "I" suffix to a column header means that the node is incompatible with every resolution ofthe summary tree.

For the redundant nodes, the report indicates the conflict status of their closest non-redundant descendant.
A redundant node can also be marked "T" (for "trivial")if it is an ancestor of only 1 leaf or of the root.

The "F" and "R" column are just the sums for forking and redundant entries.

The "label" shows the tree name or "Total of # trees" for the global sum
The ordering of the rows is the input order. The final row shows the totals.
.
For columns, the order is: FD FCR FI F RD RCR RI RT R label.


### Checking for additional splits that could be added

    otc-find-resolution taxonomy.tre synth.tre tree1.tre tree2.tre ...

will look for groups in the input trees (`tree1.tre`, `tree2.tre`...) which could
resolve polytomies in `synth.tre`.  `taxonomy.tre` is used for label validation
and expanding any tips in input trees that are mapped to non-terminal taxa.

## Tools used in the supertree pipeline
### expanding tips mapped to higher taxa and pruning the taxonomy
`otc-nonterminals-to-exemplars` takes an -e flag specifying an export diretory and at least 2 newick file paths: a full taxonomy tree some number of input trees.
Any tip in non-taxonomic input that is mapped to non-terminal taoxn will be remapped such
that the parent of the non-terminal tip will hold all of the expanded exemplars.
The exemplars will be the union of tips that (a) occur below this non-terminal taxon in the taxonomy
and (b) occur, or are used as an exemplar, in another input tree.
The modified version of each input will be written in the export directory.
Trees with no non-terminal tips should be unaltered.
The taxonomy written out will be the taxonomy restricted to the set of leaves that are leaves of the exported trees:

    otc-nonterminals-to-exemplars -estep_5 taxonomy.tre inp1.tre inp2.tre ...

This is intended to perform steps 2.5 and 2.6 of the supertree pipeline mentioned in the `doc` subdirectory.

### pruning a taxonomy

    otc-prune-taxonomy taxonomy.tre inp1.tre inp2.tre ...

will write (to stdout) a newick version of the taxonomy that has been pruned to
not include subtrees that do not include any of the tips of the input trees.  See
[supertree/step-2-pruned-taxonomy/README.md](./supertree/step-2-pruned-taxonomy/README.md)
for a more precise description of the pruning rules. This is intended to be used in the
[ranked tree supertree pipeline](./supertree/README.md),

### decomposing set of  trees into subproblems based on uncontested taxa

    otc-uncontested-decompose -eEXPORT taxonomy.tre -ftree-list.txt

will create subproblems in the (existing) subdirectory EXPORT using the taxonomy.tre
as the taxonomy and every tree listed in tree-list.txt. (each line of that file)
should be an input tree filepath. Each output will have:
  * a name that corresponds to the OTT taxon,
  * the trees pruned down for each subproblem (in the) same order as the trees were
      provided in the invocation, and
  * a corresponding ott###-tree-names.txt file that list the input filenames for each
    tree (or "TAXONOMY" for taxonomy, which will always be the last tree).

*NOTE*: phylogenetic tips mapped to internal labels in the taxonomy will be pruned if
   the taxon is contested. This is probably not what one usually wants to do...

### supertree using the decomposition and a greedy subproblem solver
`otc-scaffolded-supertree` is incomplete. If completed it will produces a supertree
of the its inputs.

### Solving a subproblem

    otc-solve-subproblem subproblem.tre
    otc-solve-subproblem tree1.tre tree2.tre taxonomy.tre

This will construct a synthesis tree and write it out in newick format.  Here subproblem.tre
contains a list of newick trees ending in the taxonomy.  If more than one tree file is supplied,
the trees are concatenated to form a single subproblem.  Earlier trees are ranking higher.

The current solution algorithm attempts to add splits one-at-a-time, checking to see whether
the split set is consistent using the BUILD algorithm.

Non-terminal taxa in the input trees are allowed if they occur in the taxonomy.  Each terminal
taxon contained in the non-terminal taxon is attached to the parent of the non-terminal taxon.
The non-terminal taxon is them removed.  This behavior can be changed to reject non-terminal taxa
with
  * `-ifalse` for rejecting non-terminal taxa

Flags allow running the solver on non-standard input.
  * `-ofalse` for handling tree files without OTT ids
  * `-T` for handling subproblems without a taxonomy.
  * `-S` writes out a standardized subproblem instead of running a solver.

## Miscellaneous tree manipulations and tree statistics

### Calculating stats for the subproblems
This works on the outputs of `otc-uncontested-decompose`. Running:

    otc-subproblem-stats *.tre > stats.tsv

Will create a tab-separated file of stats for the subproblems.
As of 5, May 2015, the columns of the report are:
  *  Subproblem name
  * InSp = # of informative (nontrivial) splits
  * LSS = size of the leaf label set
  * ILSS = size of the set of labels included in at least one "ingroup"
  * NT = The number of trees.
  * TreeSummaryName = tree index or summary name where the summary name can be Phylo-only or Total.
  "Total" summarizes info all trees in the file (including the taxonomy).
  "Phylo-only" former summarizes all of the phylogenetic inputs.

Use the `-h` option to see an explanation of the columns if they differ from this list.

### Grafting subproblems back together
This works on the outputs of `otc-solve-subproblem`.  Running

```sh
	otc-graft-solutions ott*-solution.tre > grafted_solution.tre
```
or
```sh
	cat ott*-solution.tre > solutions.tre
	otc-graft-solutions solutions.tre > grafted_solution.tre
```
will produce a newick tree file containing the grafted solution.

If the sub-problems do not connect into a single component, the program will exit
with error code 1.  The program will write multiple trees, where each tree is a
connected component whose root is not found in the other trees.

The `-n` argument can be used to name the root if desired:

```sh
	otc-graft-solutions solutions.tre -nlife > grafted_solution.tre
```

### Unpruning the grafted solution

This tool takes the grafted solution and re-attaches leaves that were pruned

```sh
	otc-unprune-solution grafted_solution.tre cleaned_ott.tre > full_supertree.tre
```

The first argument is the grafted solution.  This is a solution on a reduced taxon set.

The second argument is a full (cleaned) taxonomy.  This contains leaves that have been pruned.

In order for this tool to work, the grafted solution must have internal nodes corresponding to the
taxonomy labelled with their OTT Ids.  Currently the generation of subproblems, solution of subproblems,
and grafting of solutions preserve these labels.

Typically, many leaves on the grafted solution are internal nodes in the full taxonomy.  In this
case, the leaves in the grafted solution are expanded to match the taxonomy.

Since many nodes in the taxonomy may have out-degree 1, unpruning involves re-inserting such nodes
into the grafted solution to form the full supertree.

In theory, one could use the sub-problem solver to unprune, if the sub-problem solver would handle
taxonomy nodes with out-degree 1.

### Naming unnamed nodes

This tool takes a series of trees, names the unnamed nodes, and writes out the resulting trees:

```sh
otc-name-unnamed-nodes tree1.tre > tree1-named.tre
```

It is assumed that monotypic nodes always have OTT Ids, and are therefore named.  Names for
unnamed nodes are of the form mrca-ottX-ottY.  To find X and Y in a unique, repeatable way,
each node in the tree is annotated with the OTT Id of the smallest leaf in the include group
for that node.  X and Y are then the annotations of the child nodes
with the smallest, and second-smallest annotations, respectively.

### Creating annotations for the synthesis tree

This tool takes a series of newick trees: a full supertree, and some number of input trees.

```sh
otc-annotate-synth super.tre inp1.tre inp2.tree ...
```

It outputs a JSON document with fields describing relationships between the input tree edges and the
supertree.  Relationships include conflict, support, etc and are described in the
[OpenTree v3 conflict API](https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3).

### Renaming tree nodes

This tool takes a Newick tree and writes out a relabelled tree.

```sh
otc-relabel-tree in.tre --format-tax="%N ott%I" --taxonomy=<ott-dir> --del-monotypic > out.tre
```

Format codes are given in `otc-relabel-tree -h`.  It is also possible to relabel
non-taxonomy nodes, but without refering to taxonomy fields.

It is possible to avoid specifying the taxonomy, if the the file `~/.opentree` contains
a config file specifying the location of OTT.

### getting the full distribution of out-degree counts for a tree

    otc-degree-distribution sometree.tre

will write out a tab-separated pair of columns of "out degree" and "count" that
shows how many nodes in the tree tree have each outdegree (0 are leaves. 1 are
redundant nodes. 2 are fully resolved internals...)

### counting the number of polytomies in a tree

    otc-polytomy-count sometree.tre

will write out the number of nodes with out degree greater than 2 to stdout. This
is just a summary of the info reported by `otcdegreedistribution`.

**Untested**

### counting the number of leaves in a tree
`otc-count-leaves` takes a filepath to a newick file and reports the number of leaves:

    otc-count-leaves sometree.tre

### Detecting contested taxa

`otc-detect-contested` takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees.
It will print out the OTT IDs of clades in the taxonomy whose monophyly is questioned by at least one input:

    otc-detect-contested taxonomy.tre inp1.tre inp2.tre

### Get an induced subtree

`otc-induced-subtree` takes at least 2 newick file paths: a full tree, and some number of input trees.
It will print a newick representation of the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes):

    otc-induced-subtree taxonomy.tre inp1.tre

**Untested**

### Extract a subtree from a larger tree
`otc-prune-to-subtree`: Reads a large tree and takes a set of OTT Ids.
It finds the MRCA of the OTT Ids, and writes the subtree for that MRCA as newick.
The flag preceding the comma-separated list of IDs indicates whether the user
want the subtree for the MRCA node (`-n` flag), its parent(`-p` flag), each
of its children (`-c` flag and writing one line per child), or each
of its siblings (`-s` flag and writing one line per sib):

    otc-prune-to-subtree -p5315,3512 some.tre
    otc-prune-to-subtree -n5315,3512 some.tre
    otc-prune-to-subtree -c5315,3512 some.tre
    otc-prune-to-subtree -s5315,3512 some.tre


**Untested**

### Find distance between a supertree and the input trees
`otc-disance` takes at least 2 newick file paths: a supertree, and some number of input trees.
It will print the Robinson-Foulds symmetric difference between the induced tree from the full tree to each
input tree (one RF distance per line), or the number of groupings in each input tree that are
either displayed or not displayed by the supertree

    otc-distance -r taxonomy.tre inp1.tre inp2.tre

Note the `otc-missing-splits` script reports just the splits in the induced tree that are
missing from the subsequent trees.
Comparing this number to the RF would reveal the number of groupings that are missing from the induced
tree but present in a subsequent tree.
Thus, one can calculate "missing" and "extra" grouping counts from the output of both tools.

**Untested**

### Suppress nodes of outdegree=1
`otc-suppress-monotypic` takes a filepath to a newick file and writes a newick
without any nodes that have just one child:

    otc-suppress-monotypic taxonomy.tre

### Suppress nodes of outdegree=1
`otc-suppress-monotypic` takes a filepath to a newick file and writes a newick
without any nodes that have just one child:

    otc-suppress-monotypic taxonomy.tre

### finding the union or intersection of the OTT Ids in a set of trees

    otc-set-of-ids tree1.tre tree2.tre

will print out the union of OTT Ids in tree1 and tree2.
The `-i` flag requests the intersection rather than the union.
The `-t` flag requests that only the tips be considered.
The `-n` flag requests that the output should be a newick tree (a polytomy) rather than a list.

# Testing
`otcetera` is still very much under development. You can trigger the running of the
tests by:

    $ make
    $ make check

(currently there are no tests in the `make installcheck` target).
The data for running these tests is in the `data` subdirectory (but the tests
are supposed to know how to find that data, so users do not need to know the location).

## unit tests

Some of the operations have unit-tests. These tests are found in the `test` subdir.
Successful execution of these tests results in a row of periods (one per test) appearing
when the make check enters the test directory.

## tools tests
Some of the executables in the `tools` subdirectory have tests. These
are executed as a part of the normal `make check` target.
The output of the tool can be check using text comparison or tree comparisons
  (to handle cases in which branch rotation might result in multiple valid outputs
  of the same operation).
Some of the tests just check the exit code.

The syntax used to describe a new test is described in [../expected/README.md](../expected/README.md)
and the directories that describe the expected behavior are in the `expected` subdirectory.

## ACKNOWLEDGEMENTS
See comments above about usage of [easyloggingpp](https://github.com/easylogging/)
and [rapidjson](https://github.com/miloyip/rapidjson)

To acknowledge the contributions of the NCL code and ideas, a snapshot of the
NCL credits taken from the version of NCL used to jump start otcetera is:

As of March 09, 2012, NCL is available under a Simplified BSD license (see
BSDLicense.txt) in addition to the GPL license.

NCL AUTHORS -- the author of the NEXUS Class Library (NCL) version 2.0 is

  Paul O. Lewis, Ph.D.
  Department of Ecology and Evolutionary Biology
  The University of Connecticut
  75 North Eagleville Road, Unit 3043
  Storrs, CT 06269-3043
  U.S.A.

  WWW: http://lewis.eeb.uconn.edu/lewishome
  Email: paul.lewis@uconn.edu


Versions after 2.0 contain changes primarily made by:
  Mark T. Holder  mholder@users.sourceforge.net

Other contributors to these versions include: Derrick Zwickl, Brian O'Meara, Brandon Chisham, François Michonneau, and Jeet Sukumaran

The code in examples/phylobase... was written by Brian O'Meara and Derrick Zwickl
for phylobase.

David Suárez Pascal contributed SWIG bindings which heavily influenced those
   found in branches/v2.2. Thanks to David for blazing the way on the swig binding,
    Google for funding, and NESCent (in particular Hilmar Lapp) for getting the
    NESCent GSoC program going.

The 2010 GSoC effort also led to enhancements in terms of annotation storage and
xml parsing which are currently on. Michael Elliot contributed some code to the branches/xml branch.
Thanks to NESCent and  Google for supporting that work.

Many of the files used for testing were provided by Arlin Stoltzfus (see
http://www.molevol.org/camel/projects/nexus/ for more information), the Mesquite
package, and from TreeBase (thanks, Bill Piel!).
