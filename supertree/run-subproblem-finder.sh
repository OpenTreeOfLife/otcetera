#!/usr/bash
# Takes 2 script args  and then the args that will be passed to otc-uncontested-decompose
#   1. the export dir (any .tre files here will be removed, and new ones will be written here
#   2. a filepath to store the .tre filenames for all of the tree files created by the
#       the decomposition. This will only be created if the decomposition exits without error.
exportdir="$1"
shift
dumpedidfile="$1"
shift
set -x
if test -d "${exportdir}"
then
    rm -f "${exportdir}"/*.tre || exit
else
    mkdir "${exportdir}" || exit
fi
otc-uncontested-decompose -e"${exportdir}" $@ || exit
cd "${exportdir}"
ls *.tre | sort > "${dumpedidfile}"

