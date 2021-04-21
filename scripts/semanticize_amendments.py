#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ConstGitRepo is a base class containing methods that specify
repo-specific information, but no methods that alter the state
of the repository. This includes methods that just read from
the git db. They do not checkout, or switch branches...


"""
import sys
try:
    from git_wrap_opentree import ConstGitRepo
except ImportError:
    sys.exit('Could not import git_wrap_opentree. See https://github.com/OpenTreeOfLife/git_wrap_opentree for installation instructions.')

def main(repo_top, last_handled_sha):
    repo = ConstGitRepo(repo_top=repo_top)
    for c in repo.commits_after(after_sha=last_handled_sha):
        for fp, oid in repo.files_touched(c):
            if '/' in fp:
                print(fp)
                print(repo.get_file_contents(oid))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])