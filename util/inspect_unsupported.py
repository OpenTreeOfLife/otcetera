#!/usr/bin/env python
import subprocess
import sys
import re
import os
PAT = re.compile('Unsupported.*ott([0-9]+)".*ott([0-9]+)"')
DIR = os.path.split(sys.argv[0])[0]
FIND_SUB = os.path.join(DIR, 'find_subproblem_by_mrca_designators.py')
assert os.path.exists(FIND_SUB)
try:
    SUBPROB_DIR = sys.argv[1]
    assert os.path.isdir(SUBPROB_DIR)
except:
    sys.exit('Expecting the first argument to be the directory holding the subproblems')
if len(sys.argv) > 2:
    inp = open(sys.argv[2], 'rU')
else:
    inp = sys.stdin
for line in inp:
    m = PAT.match(line)
    if m:
        f, s = m.group(1), m.group(2)
        print line.strip()
        subprocess.call([sys.executable, FIND_SUB, SUBPROB_DIR, f, s])