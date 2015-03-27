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

## Installation

### prerequisites

#### rapidjson
To facilitate parsing of NexSON, this version of requires rapidjson.
Download it from https://github.com/miloyip/rapidjson

**Note**: you'll need to put the full path to the `rapidjson/include` subdir 
after a "-I" flag in either the CPPFLAGS variable or the CXXFLAGS variable
when you run configure for otcetera



#### autotools
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


### configuration + building

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

Python 2 (recent enough to have the subprocess module as part of the standard lib)
is required for the `make check` operation to succeed.


## Usage
See the [supertree/README.md](./supertree/README.md) for instructions on using
`otcetera` to build a supertree (work in progress).

### Common command line flags
The tools use the same (OTCLI) class to process command line arguments. 
This provides the following command line flags:
  * `-h` for help
  * `-fFILE` to treat every line of FILE as if it were a command line argument
      (useful for processing hundreds of filenames)
  * `-v` for verbose output
  * `-q` for quieter than normal output
  * `-t` for trace level (extremely verbose) output

Unless otherwise stated, the command line tools that need a tree take a filepath 
to a newick tree file. The numeric suffix of each label in the tree is taken to
be the OTT id. This accommodates the name munging that some of the open tree of
life tools perform on taxonomic names with special characters (because only the
OTT id is used to associate labels in different trees)

### Checking for incorrect internal labels in a full tree

    otc-check-taxonomic-nodes synth.tre taxonomy.tre

will check every labelled internal node is correctly labelled. To do this, it 
verifies that the set of OTT ids associated with tips that descend from the 
node is identical to the set of OTT ids associated with terminal taxa below
the corresponding node in the taxonomic tree.

A report will be issued for every problematic labeling. 

Assumptions:
  1. synth tree and taxonomy tree have the same leaf set in terms of OTT ids
  2. each label has numeric suffix, which is treated as the OTT id.

### Checking for unnamed nodes in a full tree that have no tree supporting them

    otc-find-unsupported-nodes taxonomy.tre synth.tre inp1.tre inp2.tre ...

will report any nodes in `synth.tre` that are not named and which do not have
any [ITEB support](https://github.com/OpenTreeOfLife/treemachine/blob/nonsense-1/nonsense/iteb_support_theorem.md)

The taxonomy is just used for the ottID validation (on the assumption that 
the nodes supported by the taxonomy and the the otcchecktaxonomicnodes tool
can help identify problems with those nodes).

### pruning a taxonomy

    otc-prune-taxonomy taxonomy.tre inp1.tre inp2.tre ...

will write (to stdout) a newick version of the taxonomy that has been pruned to 
not include subtrees that do not include any of the tips of the input trees.  See
[supertree/step-2-pruned-taxonomy/README.md](./supertree/step-2-pruned-taxonomy/README.md)
for a more precise description of the pruning rules. This is intended to be used in the
[ranked tree supertree pipeline](./supertree/README.md),


### getting the full distribution of out degree counts for a tree

    otc-degree-distribution sometree.tre

will write out a tab-separated pair of columns of "out degree" and "count" that
shows how many nodes in the tree tree have each outdegree (0 are leaves. 1 are
redundant nodes. 2 are fully resolved internals...)

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

### decomposing a set of trees into subproblems by uncontested taxonomic group

    otc-polytomy-count sometree.tre

will write out the number of nodes with out degree greater than 2 to stdout. This
is just a summary of the info reported by `otcdegreedistribution`.

### debugging otcetera handling of trees

    otc-assert-invariants sometree.tre

will parse a tree and run through lots of traversals, asserting various
invariants. Should exit silently if there are no bugs.



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

Other contributors to these versions include:
  Derrick Zwickl
  Brian O'Meara
  Brandon Chisham
  François Michonneau
  Jeet Sukumaran

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

