#!/bin/bash
CXX=$(which g++-5) CC=$(which gcc-5) \
    CXXFLAGS="-Waddress -Warray-bounds -Wc++11-compat -Wchar-subscripts -Wcomment -Wformat -Wmain -Wmaybe-uninitialized -Wmissing-braces -Wnonnull -Wopenmp-simd -Wparentheses -Wreorder -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow -Wswitch -Wtrigraphs -Wuninitialized -Wno-unknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -Wno-pragmas -pedantic -g -O0 -std=c++14" \
    ../configure --prefix=$PWD/installed


