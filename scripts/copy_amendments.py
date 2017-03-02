#!/usr/bin/env python
#
import subprocess
import sys
import os

def get_list_of_commits(amend_dir, last_ignored_sha=None):
    """Returns a list of SHAs (oldest to newest).

    If `last_ignored_sha` is provided and is a commit SHA, then
    the first sha returned witll be the first commit after that SHA.
    """
    log_out = subprocess.check_output(["git", "log", "--pretty=oneline"], cwd=amend_dir)
    shas = []
    for line in log_out.split("\n"):
        if line:
            sha = line.split()[0].strip()
            if last_ignored_sha and last_ignored_sha == sha:
                break
            shas.append(sha)
    shas.reverse()
    return shas

def get_amend_filename_touched(amend_dir, old_sha, new_sha):
    log_out = subprocess.check_output(["git", "diff", "--name-only", old_sha, new_sha], cwd=amend_dir)
    amend_file_subpaths = []
    for line in log_out.split("\n"):
        l = line.strip()
        if l.startswith("amendments/"):
            amend_file_subpaths.append(l)
    return amend_file_subpaths


def write_files_at_sha(amend_dir, sha, amend_file_paths, dest):
    assert(len(amend_file_paths) == len(dest))
    for ind, afp in enumerate(amend_file_paths):
        dfp = dest[ind]
        show_out = subprocess.check_output(["git", "show", "{}:{}".format(sha, afp)], cwd=amend_dir)
        with open(dfp, "w") as dfo:
            dfo.write(show_out)
    

def copy_amend_files_in_commits(amend_dir, last_ignored_sha, sha_list, out_dir):
    ind_fn = os.path.join(out_dir, "index.txt")
    with open(ind_fn, "w") as ind_fo:
        prev_sha = last_ignored_sha
        for sha in sha_list:
            amend_file_paths = get_amend_filename_touched(amend_dir, prev_sha, sha)
            if amend_file_paths:
                ind_fo.write("{}\n".format(sha))
                sd = os.path.join(out_dir, sha)
                if not os.path.isdir(sd):
                    os.mkdir(sd)
                dest = []
                for afp in amend_file_paths:
                    fn = os.path.split(afp)[-1]
                    dest.append(os.path.join(sd, fn))
                write_files_at_sha(amend_dir, sha, amend_file_paths, dest)
            prev_sha = sha


if __name__ == '__main__':
    try:
        amend_dir, sha, out_dir = sys.argv[1:]
    except:
        sys.exit("""Expecting 3 arguments:
    1. path of the amendments shard, 
    2. the last SHA that you have already processed, and
    3. the output directory

This script will copy all of the amendments that have been committed after
the last SHA into subdirectories of the output dir.

The subdirectory name will be the SHA.

An index.txt file will be written that lists the order of the SHAs (oldest
first). Processing the amendments in that order will recaptiulate the
series of commits made since last SHA.
""")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_dir = os.path.abspath(out_dir)
    sha_list = get_list_of_commits(amend_dir, sha)
    if not sha_list:
        sys.stderr.write("No commits since {}\n".format(sha))
        sys.exit(0)
    copy_amend_files_in_commits(amend_dir, sha, sha_list, out_dir)