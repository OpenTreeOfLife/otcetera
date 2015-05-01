#!/bin/bash
if test -z $RAPID_JSON_INC
then
    echo RAPID_JSON_INC must be in your env
    exit 1
fi
set -x
CPPFLAGS="-I$RAPID_JSON_INC" CXX=$(which clang++-3.5) CC=$(which clang-3.5) \
    CXXFLAGS="-Wno-c++98-compat -Wno-old-style-cast -Weverything -Wpadded -pedantic -O3 -std=c++1y -cxx-isystem=/usr/lib/gcc/x86_64-linux-gnu/4.9" \
    ../configure --prefix=$PWD/installed


