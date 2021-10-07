#!/usr/bin/env python3
from dendropy.simulate import treesim
from dendropy import TaxonNamespace
from dendropy.utility.error import SeedNodeDeletionException
import sys
from random import Random


def tree_print(t):
    return t.as_string(schema='newick',
                       suppress_edge_lengths=True,
                       suppress_rooting=True)

def do_ecr_move(tree, rng):
    is_internal_node = lambda x: bool((len(x.child_nodes()) > 0) and x is not tree.seed_node)
    internal_nodes = [i for i in tree.postorder_node_iter(is_internal_node)]
    print(internal_nodes)

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
                                         taxon_namespace=taxon_namespace)

    err_stream.write('{} true-tree: {}'.format(out_fn, tree_print(true_tree)))
    for i in range(num_phylos):
        t = None:
        while t is None:
            try:
                t = true_tree.extract_tree(node_filter_fn=lambda n: rng.random() < otu_inclusion_prob,
                                           is_apply_filter_to_internal_nodes=False,
                                           suppress_unifurcations=True)
            except SeedNodeDeletionException:
                pass
        for j in range(num_ecr):
            do_ecr_move(tree=t, rng=rng)
        if t.seed_node.is_leaf():
            r = tree_print(t).strip()
            assert not r.startswith('(')
            assert r.endswith(';')
            out_stream.write('({});\n'.format(r[:-1]))
        else:
            out_stream.write(tree_print(t))
    tax = true_tree.extract_tree()
    out_stream.write(tree_print(tax))
    

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
                    help="If not specified, problems will be emitted to standard output. Otherwise sim #.tre will be appended to each prefix to generate the output file.")
    p.add_argument("--seed",
                    type=int, 
                    default=None,
                    help="RNG seed (clock used if omitted)")
    args = p.parse_args()
    if args.num_phylos < 1:
        sys.exit('--num-phylos must be greater than 0.')
    if args.num_otus < 3:
        sys.exit('--num-otus must be greater than 2.')
    if args.num_ecr < 0:
        sys.exit('--num-ecr must not be negative.')
    if args.num_sim_reps < 1:
        sys.exit('--num0-sim-reps must be positive.')
    if args.otu_inclusion_prob < 0.0 or args.otu_inclusion_prob > 1.0:
        sys.exit('--otu-inclusion-prob must be a probability.')
    if args.tax_edge_collapse_prob < 0.0 or args.tax_edge_collapse_prob > 1.0:
        sys.exit('--tax-edge-collapse-prob must be a probability.')
    
    rng = Random()
    if args.seed is not None:
        if args.seed < 1:
            sys.exit('--seed must be positive')
        seed = args.seed
    else:
        import time
        seed = time.monotonic_ns()
    rng.seed(seed)
    sys.stderr.write("gen_subproblem.py: seed = {s}\n".format(s=seed))    

    num_otus = args.num_otus
    tns = TaxonNamespace(["ott{}".format(i) for i in range(1, 1+num_otus)])
    for rep_index in range(args.num_sim_reps):
        rep_num = 1 + rep_index
        ope_f = False
        if args.output_prefix is not None:
            out_fp = args.output_prefix + rep_num + ".tre"
            out_stream = open(out_fp, mode="w")
            ope_f = True
        else:
            out_stream = sys.stdout
        try:
            gen_subprob(num_phylos=args.num_phylos,
                        num_ecr=args.num_ecr,
                        otu_inclusion_prob=args.otu_inclusion_prob,
                        tax_edge_collapse_prob=args.tax_edge_collapse_prob,
                        out_stream=out_stream, 
                        taxon_namespace=tns,
                        rng=rng)
        finally:
            if ope_f:
                out_stream.close()



if __name__ == '__main__':
    main()
