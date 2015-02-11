if test -z $RAPID_JSON_INC
then
    echo RAPID_JSON_INC must be in your env
    exit 1
fi
CPPFLAGS="-I$RAPID_JSON_INC" CXX=/usr/bin/clang++ CC=/usr/bin/clang \
    CXXFLAGS="-Wno-c++98-compat -Weverything -Wpadded -pedantic -g -O0 -std=c++11" \
    ../configure --prefix=$PWD/installed


