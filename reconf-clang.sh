#!/bin/bash
set -x
CXX=$(which clang++) CC=$(which clang) \
    CXXFLAGS="-DDEBUGGING_PHYLO_STATEMENTS -Wno-c++98-compat -Wno-documentation -Wno-global-constructors -Wno-exit-time-destructors -Wno-c++98-compat-bind-to-temporary-copy -Wno-old-style-cast -Weverything -Wpadded -Wno-documentation-unknown-command -pedantic -pedantic -g -O0 -std=c++1y -cxx-isystem=/usr/lib/gcc/x86_64-linux-gnu/4.9 -D__extern_always_inline=inline" \
    ../configure --prefix=$PWD/installed
