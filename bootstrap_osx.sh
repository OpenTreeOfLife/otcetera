#!/bin/sh
set -x
# on mac I had to run libtoolize before aclocal.
#   if i recall correctly, the reverse was true on another
#   platform. So we'll try each and allow failure, and then repeat the
#   invocations and bailout if they still fail.
aclocal -I config
glibtoolize
aclocal -I config || exit
glibtoolize || exit
autoheader || exit
automake --add-missing || exit
autoconf
