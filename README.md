# otcetera - phylogenetic file format parser in C++ 
[![Build Status](https://secure.travis-ci.org/OpenTreeOfLife/otcetera.png)](http://travis-ci.org/OpenTreeOfLife/otcetera)

otcetera owes a lot of code and ideas to Paul Lewis' Nexus Class Library.
  See http://hydrodictyon.eeb.uconn.edu/ncl/ and 
  https://github.com/mtholder/ncl

It also uses easylogginpp which is distributed under an MIT License. See
  http://github.com/easylogging/ for info on that project. The file from
  that project is otc/easylogging++.h
  

## prerequisites
To facilitate parsing of NexSON, this version of requires rapidjson.
Download it from https://github.com/miloyip/rapidjson
and put the path to its include subdir in a "-I" CPPFLAGS or CXXFLAGS
when you run configure.

You also need the whole autotools stack including libtool.

# Installation

To run the whole autoreconf stuff in a manner that will add missing bits as needed,
run:

    $ sh bootstrap.sh

Then to configure and build with clang use:

    $ mkdir buildclang
    $ cd buildclang
    $ sh ../reconf-clang.sh
    $ make
    $ make check
    $ make install
    $ make installcheck


## NCL credits
As of March 09, 2012, NCL is available under a Simplified BSD license (see
BSDLicense.txt) in addition to the GPL license.

# ACKNOWLEDGEMENTS
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

