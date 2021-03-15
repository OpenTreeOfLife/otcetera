// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include <iostream>
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "json.hpp"
#include "otc/conflict.h"
#include "otc/otcli.h"
#include "otc/util.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/triplet_analysis.h"


using json=nlohmann::json;


using namespace otc;
using std::vector;
using std::string;

namespace po = boost::program_options;
using po::variables_map;

// TODO: Could we exemplify tips here, if we had access to the taxonomy?
variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("first", value<string>(),"First tree filepath.")
        ("second", value<string>(),"Second tree filepath.")
        ;

    options_description reporting("Reporting options");
    // reporting.add_options()
    //     ("each",value<bool>()->default_value(false),"Show separate results for each input tree")
    //     ("all",value<bool>()->default_value(true),"Show accumulated over all input trees")
    //     ("switch","Count synth nodes instead of input tree nodes")
    //     ("names,N","Write out node names instead of counts.")
    //     ;

    options_description visible;
    visible.add(reporting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("first", 1);
    p.add("second", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-triple-distances <first> <second> [OPTIONS]\n"
                                                    "Reports stats based on tree differences using rooted-triple distances.",
                                                    visible, invisible, p);

    return vm;
}

using Tree_t = ConflictTree;

void triplet_dist_analysis(const Tree_t & inp_tre1,
                           const Tree_t & inp_tre2) {
    TripletDistAnalysis<Tree_t> tda{inp_tre1, inp_tre2};
    std::ostream & out =std::cout;
    out << "ndiff\tncomp\tfdiff\tpruned\tnpruned" << std::endl;
    const auto n_rounds = tda.get_num_rounds();
    for (std::size_t round_i = 0; round_i < n_rounds; ++round_i) {
        const auto dc = tda.get_tot_diff_comp_for_round(round_i);
        out << dc.first << '\t' << dc.second << '\t' << frac_diff_from_pair(dc) << '\t';
        if (round_i > 0) {
            const auto pn = tda.get_nodes_paired_after_round(round_i - 1);
            const auto fnp = pn.first;
            out << fnp->get_name();
        }
        out << '\t' << round_i << std::endl; 
    }
}

int main(int argc, char *argv[]) {
    try {
        variables_map args = parse_cmd_line(argc, argv);
        string first = args["first"].as<string>();
        string second = args["second"].as<string>();
        ParsingRules rules;
        // 1. Load and process summary tree.
        auto fir_trees = get_trees<Tree_t>(first, rules);
        auto sec_trees = get_trees<Tree_t>(second, rules);
        std::cerr << fir_trees.size() << " trees in " << first << std::endl;
        std::cerr << sec_trees.size() << " trees in " << second << std::endl;
        if (fir_trees.size() != sec_trees.size()) {
            std::cerr << "otc-triplet-distances: Error tree files must have the same number of trees" << std::endl;
            exit(1);
        }
        for (std::size_t ind = 0 ; ind < fir_trees.size(); ++ind) {
            const Tree_t & fir_tree = *(fir_trees.at(ind));
            const Tree_t & sec_tree = *(sec_trees.at(ind));
            triplet_dist_analysis(fir_tree, sec_tree);
        }
    } catch (std::exception& e) {
        std::cerr << "otc-triplet-distances: Error! " << e.what() << std::endl;
        exit(1);
    }
}
