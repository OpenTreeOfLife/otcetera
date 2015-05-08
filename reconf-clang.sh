#!/bin/bash
if test -z $RAPID_JSON_INC
then
    echo RAPID_JSON_INC must be in your env
    exit 1
fi
set -x
CPPFLAGS="-I$RAPID_JSON_INC" CXX=$(which clang++) CC=$(which clang) \
    CXXFLAGS="-DDEBUGGING_PHYLO_STATEMENTS -Wno-c++98-compat -Weverything -Wpadded -pedantic -g -O0 -std=c++1y -cxx-isystem=/usr/lib/gcc/x86_64-linux-gnu/4.9" \
    ../configure --prefix=$PWD/installed
