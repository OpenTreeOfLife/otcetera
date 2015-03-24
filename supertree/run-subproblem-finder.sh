#!/usr/bash
set -x
scriptdir="$(dirname $0)"
exportdir="export-sub-temp"
rawsubproblems="$1"
shift
if ! test -d $rawsubproblems
then
    echo "expecting a raw subproblems output directory as the first argument"
    exit 1
fi

simplifiedsubproblems="$1"
shift
if ! test -d $simplifiedsubproblems
then
    echo "expecting a simplified subproblems output directory as the second argument"
    exit 1
fi

solutionsdir="$1"
shift
if ! test -d $solutionsdir
then
    echo "expecting a solutions output directory as the third argument"
    exit 1
fi

if false
then
if test -d "${exportdir}"
then
    if ! rmdir  "${exportdir}"
    then
        echo "${PWD}/${exportdir} is in the way (and not empty). empty it or remove it"
    fi
fi
mkdir "${exportdir}" || exit
otcuncontesteddecompose -e"${exportdir}" $@ || exit
fi
cd "${exportdir}"
ls *.tre | sort > ../dumped-subproblem-ids.txt

