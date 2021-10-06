#!/usr/bin/env python3
from dendropy.simulate import treesim
import sys
from random import Random

def gen_subprob(num_phylos=1,
                num_otus=4,
                num_ecr=1,
                otu_inclusion_prob=0.5,
                tax_edge_collapse_prob=1.0,
                out_stream=None,
                rng=None):
    if out_stream is None:
        out_stream = sys.stdout
    out_stream.write("hi\n")

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
                        num_otus=args.num_otus,
                        num_ecr=args.num_ecr,
                        otu_inclusion_prob=args.otu_inclusion_prob,
                        tax_edge_collapse_prob=args.tax_edge_collapse_prob,
                        out_stream=out_stream, 
                        rng=rng)
        finally:
            if ope_f:
                out_stream.close()



if __name__ == '__main__':
    main()
