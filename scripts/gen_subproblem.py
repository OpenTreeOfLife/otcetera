#!/usr/bin/env python3
from dendropy.simulate import treesim
import sys



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
    args = p.parse_args()
    if args.num_phylos < 1:
        sys.exit('--num-phylos must be greater than 0.')
    if args.num_otus < 3:
        sys.exit('--num-otus must be greater than 2.')
    if args.num_ecr < 0:
        sys.exit('--num-ecr must not be negative.')
    if args.otu_inclusion_prob < 0.0 or args.otu_inclusion_prob > 1.0:
        sys.exit('--args.otu-inclusion-prob must be a probability.')
    if args.tax_edge_collapse_prob < 0.0 or args.tax_edge_collapse_prob > 1.0:
        sys.exit('--args.tax-edge-collapse-prob must be a probability.')
    print(args)


if __name__ == '__main__':
    main()
