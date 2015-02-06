#!/bin/sh
set -x
aclocal -I config || exit
libtoolize || exit
autoheader || exit
automake --add-missing || exit
autoconf
