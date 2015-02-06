if test -z $RAPID_JSON_INC
then
    echo RAPID_JSON_INC must be in your env
    exit 1
fi
CPPFLAGS="-I$RAPID_JSON_INC" CXX=/usr/bin/clang++ CC=/usr/bin/clang CXXFLAGS="-Weverything -pedantic -g -O0 -std=c++11" ../configure --prefix=$PWD/install


