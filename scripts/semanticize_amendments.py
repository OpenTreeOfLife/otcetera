#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ConstGitRepo is a base class containing methods that specify
repo-specific information, but no methods that alter the state
of the repository. This includes methods that just read from
the git db. They do not checkout, or switch branches...


"""
import json
import sys
try:
    from git_wrap_opentree import ConstGitRepo
except ImportError:
    sys.exit('Could not import git_wrap_opentree. See https://github.com/OpenTreeOfLife/git_wrap_opentree for installation instructions.')

_SCRIPT_NAME = 'semanticize_amendments'

def warn(msg):
    sys.stderr.write('{} WARNING: {}\n'.format(SCRIPT_NAME, msg))

def processs_addition(amend_id, taxon):
    ott_id = taxon["ott_id"]
    sourceinfo = '{}:{}'.format(amend_id, ott_id)
    name = taxon["name"]
    parent = taxon["parent"]
    rank = taxon.get("rank", '')
    t = {"ott_id": ott_id,
         "sourceinfo": sourceinfo,
         "name": name,
         "parent": parent,
         "rank": rank, 
        }
    return {"action": "add", "taxon": t}
 


def processs_amendment(fp, content):
    if not content:
        return
    try:
        blob = json.loads(content)
    except:
        warn('Could not parse "{}" as JSON'.format(fp))
        return
    amend_id = blob["id"]
    terse_blobs = []
    if amend_id.startswith("additions"):
        for taxon in blob["taxa"]:
            terse_blobs.append(processs_addition(amend_id, taxon))
    else:
        raise NotImplementedError("only addtions are implemented.")
    return terse_blobs

def main(repo_top, last_handled_sha):
    repo = ConstGitRepo(repo_top=repo_top)
    all_terse_blobs = []
    for c in repo.commits_after(after_sha=last_handled_sha):
        for fp, oid in repo.files_touched(c):
            if '/' in fp: # top level files are for notes, ID minting...
                sys.stderr.write(f"processing {fp} ...\n")
                content = repo.get_file_contents(oid)
                tb = processs_amendment(fp, content)
                all_terse_blobs.extend(tb)
            else:
                sys.stderr.write(f"skipping {fp} ...\n")
    print(json.dumps(all_terse_blobs, indent=2))


if __name__ == '__main__':
    try:
        amend_repo, last_handled_sha = sys.argv[1], sys.argv[2]
    except:
        sys.exit("""Expecting 2 arguments: 
  1. the path to the amendments repo and,
  2. the last SHA that has been dealt with in a previous edOTT.
""")
    main(amend_repo, last_handled_sha)
