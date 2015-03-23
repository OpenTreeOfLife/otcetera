#!/usr/bash
set -x
scriptdir="$(dirname $0)"
exportdir="export-sub-temp"
destdir="$1"
shift
if ! test -d $destdir
then
    echo "expecting a destination directory as the first argument"
    exit 1
fi
if test -d "${exportdir}"
then
    if ! rmdir  "${exportdir}"
    then
        echo "${PWD}/${exportdir} is in the way (and not empty). empty it or remove it"
    fi
fi
mkdir "${exportdir}" || exit
otcuncontesteddecompose -e"${exportdir}" $@ || exit
python $"{scriptdir}"/move-subproblems-if-differing.py exportsubtemp "${destdir}"