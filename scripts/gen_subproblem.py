#!/usr/bin/env python3
from dendropy.simulate import treesim
from dendropy import TaxonNamespace, Node
from dendropy.utility.error import SeedNodeDeletionException
from random import Random
import subprocess
import tempfile
import sys
import os
import re
from collections import defaultdict

def tree_print(t):
    return t.as_string(schema='newick',
                       suppress_edge_lengths=True,
                       suppress_rooting=True)


def random_resolve(node, rng):
    if len(node._child_nodes) < 3:
        return
    # Adapted from DendroPy code
    to_attach = rng.sample(node._child_nodes, 1)
    for child in to_attach:
        node.remove_child(child)
    attachment_points = list(node._child_nodes) + [node]
    while len(to_attach) > 0:
        next_child = to_attach.pop()
        next_sib = rng.choice(attachment_points)
        next_attachment = Node()
        if next_sib is node:
            cc = list(node._child_nodes)
            node.add_child(next_attachment)
            for c in cc:
                node.remove_child(c)
                next_attachment.add_child(c)
            node.add_child(next_child)
        else:
            p = next_sib._parent_node
            p.add_child(next_attachment)
            p.remove_child(next_sib)
            next_attachment.add_child(next_sib)
            next_attachment.add_child(next_child)
        next_attachment.edge.length = 0.0
        attachment_points.append(next_attachment)
        attachment_points.append(next_child)

def get_internal_nodes(tree):
    is_internal_node = lambda x: bool((len(x.child_nodes()) > 0) and x is not tree.seed_node)
    return [i for i in tree.postorder_node_iter(is_internal_node)]
    
def do_ecr_move(tree, rng):
    """Modifies a tree in place."""
    internal_nodes = get_internal_nodes(tree)
    if not internal_nodes:
        return
    rand_internal = rng.choice(internal_nodes)
    int_edge = rand_internal.edge
    par = int_edge.tail_node
    assert par is not rand_internal
    assert int_edge.head_node is rand_internal
    int_edge.collapse()
    random_resolve(node=par, rng=rng)

def add_error(tree, num_ecr, rng):
    for j in range(num_ecr):
        do_ecr_move(tree=tree, rng=rng)

def collapse_some(tree, edge_collapse_prob, rng):
    internal_nodes = get_internal_nodes(tree)
    for nd in internal_nodes:
        if rng.random() < edge_collapse_prob:
            nd.edge.collapse()

def gen_subprob(taxon_namespace,
                num_phylos=1,
                num_ecr=1,
                otu_inclusion_prob=0.5,
                tax_edge_collapse_prob=1.0,
                out_stream=None,
                err_stream=None,
                rng=None):
    if out_stream is None:
        out_stream = sys.stdout
    if err_stream is None:
        err_stream = sys.stderr
    out_fn = out_stream.name    
    true_tree = treesim.birth_death_tree(birth_rate=1.0,
                                         death_rate=0,
                                         num_extant_tips=len(taxon_namespace),
                                         taxon_namespace=taxon_namespace,
                                         rng=rng)

    err_stream.write('{} true-tree: {}'.format(out_fn, tree_print(true_tree)))
    for i in range(num_phylos):
        t = None
        while t is None:
            try:
                t = true_tree.extract_tree(node_filter_fn=lambda n: rng.random() < otu_inclusion_prob,
                                           is_apply_filter_to_internal_nodes=False,
                                           suppress_unifurcations=True)
            except SeedNodeDeletionException:
                pass
        add_error(t, num_ecr=num_ecr, rng=rng)
        if t.seed_node.is_leaf():
            r = tree_print(t).strip()
            assert not r.startswith('(')
            assert r.endswith(';')
            out_stream.write('({});\n'.format(r[:-1]))
        else:
            out_stream.write(tree_print(t))
    tax = true_tree.extract_tree()
    add_error(tax, num_ecr=num_ecr, rng=rng)
    collapse_some(tax, edge_collapse_prob=tax_edge_collapse_prob, rng=rng)
    out_stream.write(tree_print(tax))
    

def gen_b_o_i():
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                yield (i, j, k)

def main():
    import argparse
    p = argparse.ArgumentParser(description="Simulator of OpenTree synth subproblems.")
    p.add_argument("--num-phylos",
                    type=int, 
                    default=0,
                    help="The number of phylogenies to simulate")
    p.add_argument("--num-otus",
                    type=int, 
                    default=0,
                    help="The number of OTUs in the complete solution")
    p.add_argument("--num-ecr",
                    type=int, 
                    default=1,
                    help="The number of edge contract refine moves to add noise to trees")
    p.add_argument("--num-sim-reps",
                    type=int, 
                    default=1,
                    help="The number of simulations to perform")
    p.add_argument("--otu-inclusion-prob",
                    type=float, 
                    default=0.5,
                    help="The probability of the inclusion of an OTU in each phylogeny")
    p.add_argument("--tax-edge-collapse-prob",
                    type=float, 
                    default=0.75,
                    help="The per-edge probability that a terminal edge in the taxonomy will be collapsed")
    p.add_argument("--output-prefix",
                    type=str, 
                    default=None,
                    help="If not specified, problems will be emitted to standard output. Otherwise #-sim.tre will be appended to each prefix to generate the output file.")
    p.add_argument("--output-suffix",
                    type=str, 
                    default=None,
                    help="If specified, it alters the behavior of --output-prefix, such that this suffix is added after the rep number. If the suffix ends in / then the output file will be sim.tre . Otherwise -sim.tre will be appended to the suffix")
    p.add_argument("--seed",
                    type=int, 
                    default=None,
                    help="RNG seed (clock used if omitted)")
    p.add_argument("--run-otc",
                    action='store_true', 
                    default=False,
                    help="Run otc-solve-subproblem with all of the relevant options on the simulated output")
    p.add_argument("--summary-out",
                    type=str, 
                    default=None,
                    help="Only used if --run-otc is specified. Appends runtime summaries to this file.")
    args = p.parse_args()
    if args.num_phylos < 1:
        sys.exit('--num-phylos must be greater than 0.\n')
    if args.num_otus < 3:
        sys.exit('--num-otus must be greater than 2.\n')
    if args.num_ecr < 0:
        sys.exit('--num-ecr must not be negative.\n')
    if args.num_sim_reps < 1:
        sys.exit('--num0-sim-reps must be positive.\n')
    if args.otu_inclusion_prob < 0.0 or args.otu_inclusion_prob > 1.0:
        sys.exit('--otu-inclusion-prob must be a probability.\n')
    if args.tax_edge_collapse_prob < 0.0 or args.tax_edge_collapse_prob > 1.0:
        sys.exit('--tax-edge-collapse-prob must be a probability.\n')
    
    if args.seed is not None:
        if args.seed < 1:
            sys.exit('--seed must be positive')
        seed = args.seed
    else:
        import time
        seed = time.monotonic_ns()
    
    out_pref = args.output_prefix
    run_otc = args.run_otc
    ope_sum = False
    sum_stream = sys.stderr
    if run_otc:    
        sum_existed = False
        if args.summary_out:
            sum_stream, sum_existed = open_fp(args.summary_out, 'a')
            ope_sum = True
        if not sum_existed:
            columns = ["num_otus", "num_phylos", "rep_num", "seed", "num_ecr", "inc_prob", "collapse_prob", "batch", "oracle", "incr", "time"]
            headers = '\t'.join(columns)
            sum_stream.write('{}\n'.format(headers))
    try:
        if run_otc and out_pref is None:
            with tempfile.TemporaryDirectory() as tmpdirname:
                do_sim(tmpdirname + '/', args, seed, sum_stream)
        else:
            do_sim(out_pref, args, seed, sum_stream)
    finally:
        if ope_sum:
            sum_stream.close()

def open_fp(out_fp, mode):
    existed = os.path.exists(out_fp)
    sys.stderr.write('Opening = {}\n'.format(out_fp))
    par_dir = os.path.split(out_fp)[0]
    if par_dir and not os.path.isdir(par_dir):
        os.makedirs(par_dir)
    return open(out_fp, mode=mode), existed

timing_pat = re.compile(r'^timing +(.+) +seconds.$')
def time_otc_runs(inp_fp):
    invoc_pref = ["otc-solve-subproblem", "-m", inp_fp]
    time_dict = defaultdict(lambda: 0.0)
    for b, o, i in gen_b_o_i():
        opts = ['--batching={}'.format(b),
                '--oracle={}'.format(o),
                '--incremental={}'.format(i)]
        invoc = invoc_pref + opts
        rp = subprocess.run(invoc, capture_output=True, check=True)
        timing_float = None
        for line in rp.stderr.decode('utf-8').split('\n'):
            m = timing_pat.match(line)
            if m:
                timing_float = m.group(1)
                break
            #else:
            #    print('line "{}" did not match'.format(line))
        if timing_float is None:
            raise RuntimeError('"timing ... seconds." not found in {} run.'.format(invoc))
        time_dict[(b,o,i)] = timing_float
    return time_dict

def do_sim(out_pref, args, seed, sum_stream):
    rng = Random()
    rng.seed(seed)
    sys.stderr.write("gen_subproblem.py: seed = {s}\n".format(s=seed))
    num_otus = args.num_otus
    num_ecr =args.num_ecr
    num_phylos =args.num_phylos
    otu_inclusion_prob = args.otu_inclusion_prob
    tax_edge_collapse_prob = args.tax_edge_collapse_prob
    tns = TaxonNamespace(["ott{}".format(i) for i in range(1, 1+num_otus)])
    run_otc = args.run_otc
    if run_otc:
        sum_row_temp = [str(num_otus),
                        str(num_phylos),
                        "",
                        str(seed),
                        str(num_ecr), 
                        str(otu_inclusion_prob),
                        str(tax_edge_collapse_prob)]
        rep_num_idx = 2
        seed_idx = 3
    out_suff = args.output_suffix
    if out_suff is not None:
        if out_pref is None:
            sys.exit('--output-suffix can only be used when --output-prefix is also used.\n')
        final_suff = 'sim.tre' if out_suff.endswith('/') else '-sim.tre'
    out_fp = None
    for rep_index in range(args.num_sim_reps):
        rep_num = 1 + rep_index
        ope_f = False
        if out_pref is not None:
            rs = str(rep_num)
            if out_suff:
                out_fp = out_pref + rs + out_suff + final_suff 
            else:
                out_fp = out_pref + rs + "-sim.tre"
            out_stream = open_fp(out_fp, 'w')[0]
            ope_f = True
        else:
            out_stream = sys.stdout
        try:
            gen_subprob(num_phylos=num_phylos,
                        num_ecr=num_ecr,
                        otu_inclusion_prob=otu_inclusion_prob,
                        tax_edge_collapse_prob=tax_edge_collapse_prob,
                        out_stream=out_stream, 
                        taxon_namespace=tns,
                        rng=rng)
        finally:
            if ope_f:
                out_stream.close()
        if run_otc:
            assert out_fp is not None
            times = time_otc_runs(out_fp)
            sum_row_temp[rep_num_idx] = str(rep_num)
            for (b,o,i) in times:
                tr = sum_row_temp + [str(b), str(o), str(i), str(times[(b,o,i)])]
                sum_stream.write('{}\n'.format('\t'.join(tr)))
#            sum_row_temp[seed_idx] = ''
        

if __name__ == '__main__':
    main()
