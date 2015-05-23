#!/usr/bin/env python
import subprocess

def execution_log_file(outdir):
    return os.path.join(outdir, 'execution_log.txt')
def get_cleaned_initial_tree(outdir):
    return os.path.join(outdir, 'cleaned-initial.tre')
def get_cleaned_nonredundant_tree(outdir):
    return os.path.join(outdir, 'cleaned-no-redundant-nodes.tre')
def get_resolved_treefile(outdir):
    return os.path.join(outdir, 'cleaned-resolved.tre')
def log_execution(invocation,
                  stdout_filename='/dev/null',
                  stderr_filename='/dev/null',
                  outdir='.'):
    ex_log = execution_log_file(outdir)
    with open(ex_log, 'a') as out_stream:
        out_stream.write('"{i}" >{o} 2>{e}\n'.format(i='" "'.join(invocation),
                                                     o=stdout_filename,
                                                     e=stderr_filename))
        out_stream.flush()

def execute(invoc, out_f='/dev/null', err_f='/dev/null', outdir='.'):
    log_execution(invoc, out_f, err_f, outdir=outdir)
    with open(out_f, 'w') as out_stream:
        with open(err_f, 'w') as err_stream:
            p = subprocess.Popen(invoc, stderr=err_stream, stdout=out_stream)
            if p.wait() != 0:
                log_execution(['Failure!'])
                sys.exit('Failure of {}\n'.format(invoc[0]))

def clean_tree(args, outdir):
    invoc = ['otc-find-unsupported-nodes', '-x', '-r', '-c', '-k']
    invoc.append(args.taxonomy)
    invoc.append(args.summary)
    invoc.append('-f' + args.input_files)
    out_tree = get_cleaned_initial_tree(outdir)
    err_f = os.path.join(outdir, 'cleaning-err.txt')
    execute(invoc, out_tree, err_f, outdir)

def suppress_outdegree_one(args, outdir):
    invoc = ['otc-suppress-monotypic', get_cleaned_initial_tree(outdir)]
    out_f = get_cleaned_nonredundant_tree(outdir)
    err_f = os.path.join(outdir, 'remove-redundant-err.txt')
    execute(invoc, out_f=out_f, err_f=err_f, outdir=outdir)

def check_cleaned_nonredundant(args, outdir):
    invoc = ['otc-find-unsupported-nodes']
    invoc.append(args.taxonomy)
    invoc.append(get_cleaned_nonredundant_tree(outdir))
    invoc.append('-f' + args.input_files)
    out_f = os.path.join(outdir, 'cleaned-no-redundant-check-out.txt')
    err_f = os.path.join(outdir, 'cleaned-no-redundant-check-err.txt')
    execute(invoc, out_f, err_f, outdir)

def resolve_cleaned(args, outdir):
    invoc = ['otc-find-resolution', '-r', '-u', '-v']
    invoc.append(args.taxonomy)
    invoc.append(get_cleaned_nonredundant_tree(outdir))
    invoc.append('-f' + args.input_files)
    out_tree = get_resolved_treefile(outdir)
    err_f = os.path.join(outdir, 'resolved-err.txt')
    execute(invoc, out_tree, err_f, outdir)

def check_resolved_cleaned(args, outdir):
    invoc = ['otc-find-unsupported-nodes', '-x', '-r']
    invoc.append(args.taxonomy)
    invoc.append(get_resolved_treefile(outdir))
    invoc.append('-f' + args.input_files)
    out_f = os.path.join(outdir, 'resolved-cleaned-check-out.txt')
    err_f = os.path.join(outdir, 'resolved-cleaned-check-err.txt')
    execute(invoc, out_f, err_f, outdir)

if __name__ == '__main__':
    import argparse
    import tempfile
    import sys
    import os
    description = 'Uses otc tools to clean unsupported nodes and misnamed taxa, ' \
                  'removes any nodes of outdegree-1, and then adds groupings from ' \
                  'input trees.'
    parser = argparse.ArgumentParser(prog='clean-and-resolve.py', description=description)
    parser.add_argument('--taxonomy', default=None, type=str, required=True)
    parser.add_argument('--summary', default=None, type=str, required=True)
    parser.add_argument('--input-files',
                        default=None,
                        type=str, 
                        help='file listing the input tree',
                        required=True)
    parser.add_argument('--output-dir',
                        default=None,
                        type=str, 
                        help='output dir',
                        required=True)
    args = parser.parse_args(sys.argv[1:])
    outdir = args.output_dir
    if outdir is None:
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        sys.stderr.write('work dir {} created as output dir.\n'.format(outdir))
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        sys.stderr.write('Using {} as output dir.\n'.format(outdir))
    from_top = False
    if from_top:
        if os.path.exists(execution_log_file(outdir)):
            os.unlink(execution_log_file(outdir))
        clean_tree(args, outdir=outdir)
        suppress_outdegree_one(args, outdir=outdir)
        check_cleaned_nonredundant(args, outdir=outdir)
    resolve_cleaned(args, outdir=outdir)
    check_resolved_cleaned(args, outdir=outdir)