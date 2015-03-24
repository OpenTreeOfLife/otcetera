#!/usr/bin/env python
import shutil

if __name__ == '__main__':
    import sys
    import os
    args = sys.argv[1:]
    for i in args[1:]:
        if not os.path.isdir(i):
            sys.exit('Expecting "{}" to be a directory'.format(i))
    subprob_list_fp, fresh_decomp, full_subprob, simple_subprob, solutions = args
    context = {'raw': fresh_decomp,
               'full': full_subprob,
               'simple': simple_subprob,
               'solution': solutions}
    subprob_list_fn = os.path.split(subprob_list_fp)[1]
    id_list = [i.strip() for i in open(subprob_list_fn, 'rU').readlines() if i.strip()]
    for fid in id_list:
        classify_raw_prob(context, fid)

