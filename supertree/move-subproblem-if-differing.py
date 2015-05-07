#!/usr/bin/env python
import sys
try:
    from peyotl.supertree import OtcPipelineContext, OtcArtifact
except ImportError:
    sys.exit('Depends on some "bleeding edge" feature on the peyotl supertree branch')
import os

SCRIPT_NAME = os.path.split(sys.argv[0])[1]
args = sys.argv[1:]
for i in args[2:]:
    if not os.path.isdir(i):
        sys.exit('{}: Expecting "{}" to be a directory'.format(SCRIPT_NAME, i))
subprob_tree_fn, fresh_decomp, full_subprob, simple_subprob, solutions = args
subprob_tree_fn = os.path.split(subprob_tree_fn)[1]
context = OtcPipelineContext(raw_output_dir=fresh_decomp,
                             stage_output_dir=full_subprob,
                             simplified_output_dir=simple_subprob,
                             solution_dir=solutions)
fid = context.subproblem_id_from_filename(subprob_tree_fn)
try:
    context.process_raw_subproblem_output(fid)
except Exception as x:
    if 'OTC_VERBOSE' in os.environ: # stacktrace for developers
        raise
    sys.exit('{}: Exiting due to an excetion:\n{}\n'.format(SCRIPT_NAME, x))

klist = context.summary_stats.keys()
klist.sort()
for k in klist:
    s = context.summary_stats[k]
    print 'summary_stats[', k, '] -> ', len(s)

