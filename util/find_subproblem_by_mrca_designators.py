#!/usr/bin/env python
import subprocess
import sys
import os
import re
try:
    PEYOTL_ROOT = os.environ.get('PEYOTL_ROOT')
    mrca_anc_lister = os.path.join(PEYOTL_ROOT, 'tutorials', 'ot-taxo-mrca-to-root.py')
    assert os.path.exists(mrca_anc_lister)
except:
    sys.exit('Expecting PEYOTL_ROOT to be in your env, and $PEYOTL_ROOT/tutorials/ot-taxo-mrca-to-root.py to exist.\n')
try:
    sub_prob_dir = sys.argv[1]
    assert os.path.isdir(sub_prob_dir)
except:
    sys.exit('Expecting the first argument to be a directory holding the subproblems.\n')
try:
    id_list = sys.argv[2:]
    assert len(id_list) > 1
    [int(i) for i in id_list]
except:
    sys.exit('Expecting the second, third... arguments to be OTT Ids (numeric only, no ott prefix).\n')
VERBOSE = False
INVOC = (sys.executable, mrca_anc_lister)
ID_PAT = re.compile('.*ott([0-9]+)\.tre')
def get_anc_ids(id_list):
    x = list(INVOC)
    x.extend(id_list)
    return[i.strip() for i in subprocess.check_output(x).split('\n')]

def find_mrca_subprob(anc_id_list):
    for anc_id in anc_id_list:
        fp = os.path.join(sub_prob_dir, 'ott' + anc_id + '.tre')
        if os.path.exists(fp):
            return fp
        if VERBOSE:
            sys.stderr.write('  {} not found.\n'.format(fp))

def get_ott_id_from_path(fp):
    m = ID_PAT.match(fp)
    assert m
    return m.group(1)

anc_id_list = get_anc_ids(id_list)
anc_id_set = set(anc_id_list)
rename_map = {}
for ott_id in id_list:
    a = get_anc_ids([ott_id])
    afp = find_mrca_subprob(a)
    afp_id = get_ott_id_from_path(afp)
    rename_map[ott_id] = ott_id if afp_id in anc_id_set else afp_id

subprob_fp = find_mrca_subprob(anc_id_list)
sys.stdout.write('{} contains the MRCA of the requested ids.\n'.format(subprob_fp))
for ott_id in id_list:
    sys.stdout.write(' {} appears as {} in this subproblem.\n'.format(ott_id, rename_map[ott_id]))
