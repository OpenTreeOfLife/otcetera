#!/usr/bin/env python3
import sys
import re
from peyutil import read_as_json

no_app_pat = re.compile(r"^amende?ment [#](\d+) not applied:")
homonym_pat = re.compile(r"^amende?ment [#](\d+) not applied: ([A-Za-z][-A-Za-z0-9 ]+[A-Za-z0-9]) is a homonym of (\d+)")

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
    found = False
    for ra in rel_amends[1:]:
        if ra[0] == exp_amend_idx:
            found = True
            break
    if not found:
        return False, "notfound"
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
                    print(amend_num, name, ott_id)
                    rc = check_if_del_works(amend_num, name, ott_id, edott, amendments_repo)
                    if rc[0]:
                        to_del.append(amend_num)
                    else:
                        prob = rc[1]
                        if prob == "solo":
                            print("Atypical homonym. Solo in amendments {amend_num}:", edott[amend_num -1])
                        elif prob == "notfound":
                            print("PROBLEM {amend_num} does not match:", edott[amend_num -1])
                        else:
                            assert(prob == "first")
                            print("Atypical homonym. First in amendments {amend_num} is bad:", edott[amend_num -1])
                else:
                    print(m.group(1), "not a homonym")

if __name__ == "__main__":
    try:
        _args = list(sys.argv[1:4])
        assert(len(_args) == 3)
    except:
        sys.exit("Expecting 3 arguments: the edott JSON, the logfile from the failed run, and the directory that holds the amendments-1 repo.")
    main(edott_fp=_args[0],
         logfile=_args[1],
         amendments_repo=_args[2])
