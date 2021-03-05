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
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/quartet_dist.h"


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
                                                    "Usage: otc-quartet-distances <first> <second> [OPTIONS]\n"
                                                    "Reports stats based on tree differences using quartet distances.",
                                                    visible, invisible, p);

    return vm;
}

using Tree_t = ConflictTree;
using node_t = Tree_t::node_type;
using str_set = std::set<std::string>;
using TreeAsUIntSplits = GenTreeAsUIntSplits<Tree_t>;
using AllConflictQuartets = AllQuartets<Tree_t>;

void quartet_dist_analysis(const Tree_t & inp_tre1,
                           const Tree_t & inp_tre2) {
    TreeAsUIntSplits tas_1{inp_tre1};
    TreeAsUIntSplits tas_2{inp_tre2};
    if (tas_1.leaf_label_to_ind != tas_2.leaf_label_to_ind) {
        throw OTCError() << "trees must have the same leaf label set.\n";
    }
    AllConflictQuartets t_1_q{tas_1};
    AllConflictQuartets t_2_q{tas_2};
    QuartDist qdist{t_1_q, t_2_q};
    const auto dc = qdist.get_diff_comp();
    std::cout << dc.first << "\t" << dc.second << "\t" << frac_diff_from_pair(dc) << std::endl;
    throw OTCError() << "Early exit.\n";    
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
            std::cerr << "otc-quartet-distances: Error tree files must have the same number of trees" << std::endl;
            exit(1);
        }
        for (std::size_t ind = 0 ; ind < fir_trees.size(); ++ind) {
            const Tree_t & fir_tree = *(fir_trees.at(ind));
            const Tree_t & sec_tree = *(sec_trees.at(ind));
            quartet_dist_analysis(fir_tree, sec_tree);
        }
    } catch (std::exception& e) {
        std::cerr << "otc-quartet-distances: Error! " << e.what() << std::endl;
        exit(1);
    }
}
