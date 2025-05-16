#!/usr/bin/env python3
import sys
import re
import os
from peyutil import read_as_json, write_as_json

no_app_pat = re.compile(r"^amende?ment [#](\d+) not applied:")
homonym_pat = re.compile(r"^amende?ment [#](\d+) not applied: ([A-Za-z][-A-Za-z0-9 ]+[A-Za-z0-9]) is a homonym of (\d+)")

src_fn_pat = re.compile(r"^(additions\-\d+\-\d+):(\d+)$")
by_study_id = {}

def check_if_del_works(amend_num, name, ott_id, edott, amendments_repo):
    exp_amend_idx = amend_num - 1
    rel_amends = []
    for amend_idx, amend in enumerate(edott):
        if amend.get("action", "") != "add":
            continue
        taxon = amend.get("taxon", {})
        tn = taxon.get("name", "")
        if tn.strip() == name:
            rel_amends.append((amend_idx, amend))
    if len(rel_amends) < 2:
        return (False, "solo")
    if rel_amends[0][0] == exp_amend_idx:
        return False, "first"
    found = None
    for ra in rel_amends[1:]:
        if ra[0] == exp_amend_idx:
            found = ra[1]
            break
    if found is None:
        return False, "notfound"
    taxon = found["taxon"]
    src = taxon["sourceinfo"]
    m = src_fn_pat.match(src)
    if not m:
        raise ValueError(f"'sourceinfo' {src} does not fit pattern.")
    fn_frag = m.group(1)
    bogus_id = int(m.group(2))
    fn = f"{fn_frag}.json"
    fp = os.path.join(amendments_repo, "amendments", fn)
    if not os.path.isfile(fp):
        raise RuntimeError(f"amendments file {fp} does not exist")
    offending_amend = read_as_json(fp)
    study_id = offending_amend['study_id']
    name_set = by_study_id.setdefault(study_id, set())
    name_set.add(name)
    return True, ""

def main(edott_fp,
         logfile,
         amendments_repo):
    to_del = []
    edott = read_as_json(edott_fp)
    with open(logfile, "r") as logf:
        for line in logf:
            m = no_app_pat.match(line)
            if m:
                hm = homonym_pat.match(line)
                if hm:
                    amend_num = int(hm.group(1))
                    name = hm.group(2)
                    ott_id = int(hm.group(3))
                    rc = check_if_del_works(amend_num, name, ott_id, edott, amendments_repo)
                    if rc[0]:
                        to_del.append(amend_num)
                    else:
                        prob = rc[1]
                        if prob == "solo":
                            sys.stderr.write(f"Atypical homonym. Solo in amendments {amend_num}: {edott[amend_num -1]}\n")
                        elif prob == "notfound":
                            sys.stderr.write(f"PROBLEM {amend_num} does not match: {edott[amend_num -1]}\n")
                        else:
                            assert(prob == "first")
                            sys.stderr.write(f"Atypical homonym. First in amendments {amend_num} is bad: {edott[amend_num -1]}\n")
                else:
                    sys.stderr.write(m.group(1) + " not a homonym")
    sk = list(by_study_id.keys())
    sk.sort()
    for study_id in sk:
        name_set = by_study_id[study_id]
        if len(name_set) == 1:
            name = next(name_set)
            sys.stderr.write(f"In https://tree.opentreeoflife.org/curator/study/view/{study_id} need to remap 1 taxon: \"{name}\"\n") 
        else:
            sys.stderr.write(f"In https://tree.opentreeoflife.org/curator/study/view/{study_id} need to remap {len(name_set)} taxa:\n")
            nl = list(name_set)
            nl.sort()
            for name in nl:
                sys.stderr.write(f"    \"{name}\"\n")

    if to_del:
        tds = set([i-1 for i in to_del])
        new_edott = []
        for amend_idx, amend in enumerate(edott):
            if amend_idx not in tds:
                new_edott.append(amend)
        write_as_json(new_edott, sys.stdout, indent=2)
        
if __name__ == "__main__":
    try:
        _args = list(sys.argv[1:4])
        assert(len(_args) == 3)
    except:
        sys.exit("Expecting 3 arguments: the edott JSON, the logfile from the failed run, and the directory that holds the amendments-1 repo.")
    main(edott_fp=_args[0],
         logfile=_args[1],
         amendments_repo=_args[2])
