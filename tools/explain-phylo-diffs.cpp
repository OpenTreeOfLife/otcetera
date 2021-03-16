// See https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs-v3#synthetic-tree
#include <iostream>
#include <tuple>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include "json.hpp"
#include "otc/conflict.h"
#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/quartet_dist.h"
#include "otc/triplet_analysis.h"


using json = nlohmann::json;


using namespace otc;
using std::vector;
using std::to_string;
using std::string;
using std::stack;

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

    
    options_description output("Output options");
    output.add_options()
        ("out,O",value<string>(),"output JSON filepath");
    // reporting.add_options()
    //     ("each",value<bool>()->default_value(false),"Show separate results for each input tree")
    //     ("all",value<bool>()->default_value(true),"Show accumulated over all input trees")
    //     ("switch","Count synth nodes instead of input tree nodes")
    //     ("names,N","Write out node names instead of counts.")
    //     ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("first", 1);
    p.add("second", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-explain-phylo-diffs <first> <second> [OPTIONS]\n"
                                                    "Writes a JSON summary of the difference between 2 rooted trees ",
                                                    visible, invisible, p);

    return vm;
}


enum IniNodePhyloStatus {IDENTICAL_SUBTREE,    // this subtree rooted at this node is found in both trees
                         SAME_CLADE_DIFF_TOPO, // this cluster is found, but at least on descendant has different topology
                         NO_COMMON_CLADE,
                         UNKNOWN_STATUS
                        };

const std::size_t MAX_SIZE_T = std::numeric_limits<std::size_t>::max();

struct ExplTreeDiffNode
{
   /* int depth = 0; // depth = number of nodes to the root of the tree including the  endpoints (so depth of root = 1)
    int n_tips = 0;
    int n_include_tips = 0; */
    otc::RootedTreeNode<ExplTreeDiffNode> * partner;
    IniNodePhyloStatus ini_stat;
    std::size_t node_index;

    ExplTreeDiffNode()
        :partner(nullptr),
        ini_stat(IniNodePhyloStatus::UNKNOWN_STATUS),
        node_index(MAX_SIZE_T) {
        }

};


using Tree_t = otc::RootedTree<ExplTreeDiffNode, otc::RTreeNoData>;
using node_t = Tree_t::node_type;
using str_set = std::set<std::string>;
using TreeAsUIntSplits = GenTreeAsUIntSplits<Tree_t>;
template <typename T, typename U>
using map = std::unordered_map<T,U>;



using ini_stat_ptr = std::pair<IniNodePhyloStatus, const node_t *>;

inline std::ostream & write_ini_phy_stat_pair(std::ostream & out, const ini_stat_ptr & isp) {
    if (isp.first == IniNodePhyloStatus::IDENTICAL_SUBTREE) {
        out << "<IDENTICAL_SUBTREE," << isp.second << '>';
    } else if (isp.first == IniNodePhyloStatus::SAME_CLADE_DIFF_TOPO) {
        out << "<SAME_CLADE_DIFF_TOPO," << isp.second << '>';
    } else {
        out << "NO_COMMON_CLADE";
    }
    return out;
}

template<typename T>
inline void add_node_data(const T *nd,
                          TreeAsUIntSplits & ,
                          json & node_data) {
    json cd = json::object();
    const auto & ndo = nd->get_data();
    bool add_to_blob = false;
    if (nd->is_tip()) {
        add_to_blob = true;
        cd["label"] = nd->get_name();
    }

    if (add_to_blob) {
        node_data[to_string(ndo.node_index)] = cd;
    }
    
}

template<typename T>
inline void write_ind_labelled_closing_newick(std::ostream & out,
                                              const T *nd,
                                              const T * r,
                                              TreeAsUIntSplits & tas,
                                              json & node_data) {
    out << ')';
    auto n = nd->get_parent();
    out << n->get_data().node_index;
    add_node_data(n, tas, node_data);
    if (n == r) {
        return;
    }
    while (n->get_next_sib() == nullptr) {
        out << ')';
        add_node_data(n, tas, node_data);
        n = n->get_parent();
        assert(n != nullptr);
        out << n->get_data().node_index;
        if (n == r) {
            return;
        }
    }
    out << ',';
}

template<typename T>
inline void write_ind_labelled_newick_subtree(std::ostream & out,
                                                const T *nd,
                                                TreeAsUIntSplits & tas,
                                                std::set<std::size_t> & leaf_inds,
                                                const std::set<std::size_t> & boundaries,
                                                json & node_data) {
    assert(nd != nullptr);
    auto & ndata = nd->get_data();
    auto nind = ndata.node_index;
    if (nd->is_tip() || contains(boundaries, nind)) {
        leaf_inds.insert(ndata.node_index);
    } else {
        out << '(';
        int cnum = 0;
        for (auto c : iter_child_const(*nd)) {
            if (cnum++ > 0) {
                out << ',';
            }
            write_ind_labelled_newick_subtree(out, c, tas, leaf_inds, boundaries, node_data);
        }
        out << ')';
    }
    out << nind;
    add_node_data(nd, tas, node_data);
}

inline json get_tree_for_slice(node_t * nd,
                         TreeAsUIntSplits & tas,
                         std::set<std::size_t> & leaf_inds, 
                         const std::set<std::size_t> boundaries) {
    std::ostringstream s;
    json node_data = json::object();
    write_ind_labelled_newick_subtree(s, nd, tas, leaf_inds, boundaries, node_data);
    s << ';';
    json iltree = json::object();
    iltree["newick"] = s.str();
    if (!node_data.empty()) {
        iltree["node_data"] = node_data;
    }
    return iltree;   
}

inline json get_ident_tree_comp_slice(node_t * nd, TreeAsUIntSplits & tas) {
    std::set<std::size_t> boundaries;
    std::set<std::size_t> leaf_inds;
    return get_tree_for_slice(nd, tas, leaf_inds, boundaries);
}


inline json get_comparison(const std::string & t1newick,
                           node_t * ,
                           TreeAsUIntSplits & ,
                           const std::string & t2newick,
                           node_t * ,
                           TreeAsUIntSplits & ,
                           const std::set<std::size_t> & ) {
    json iltree = json::object();
    //iltree["t1newick"] = t1newick;
    //iltree["t2newick"] = t2newick;

    std::unique_ptr<ConflictTree> tree1 = tree_from_newick_string<ConflictTree>(t1newick);
    std::unique_ptr<ConflictTree> tree2 = tree_from_newick_string<ConflictTree>(t2newick);

    TripletDistAnalysis<ConflictTree> tda{*tree1, *tree2};
    int int_num_prunings = 0;
    json pruning_rounds = json::array();
    const auto n_rounds = tda.get_num_rounds();
    int_num_prunings = n_rounds - 1;
    for (std::size_t round_i = 0; round_i < n_rounds; ++round_i) {
        json curr = json::object();
        const auto dc = tda.get_tot_diff_comp_for_round(round_i);
        curr["num_triplets_diff"] = dc.first;
        curr["num_triplets_comp"] = dc.second;
        if (dc.second > 0) {
            curr["prop_triplets_differing"] = frac_diff_from_pair(dc);
        }
        if (round_i > 0) {
            const auto pn = tda.get_nodes_paired_after_round(round_i - 1);
            const auto fnp = pn.first;
            std::string nn = fnp->get_name();
            curr["pruned"] = nn;
        }
        pruning_rounds.push_back(curr);
    }
    iltree["comp_pruning_rounds"] = pruning_rounds;
    iltree["num_prunings"] = int_num_prunings;
    return iltree;
}

using json_num_edits_pair = std::pair<json, std::size_t>;
inline json_num_edits_pair get_tree_comp_slice(node_t * t1nd,
                                TreeAsUIntSplits & tas_1,
                                TreeAsUIntSplits & tas_2, 
                                stack<node_t *> & to_do) {
    assert(t1nd);
    const auto & t1nd_data = t1nd->get_data();
    node_t * t2nd = t1nd_data.partner;
    if (t2nd == nullptr) {
        throw OTCError() << "Node \"" << t1nd_data.node_index << "\" cannot start a slice because it has no paired node.\n";
    }
    json outer = json::object();
    if (t1nd_data.ini_stat == IniNodePhyloStatus::IDENTICAL_SUBTREE) {
        outer["both_trees"] = get_ident_tree_comp_slice(t1nd, tas_1);
        return json_num_edits_pair{outer, 0};
    }
    stack<node_t *> curr_slice_to_do;
    for (auto c : iter_child(*t1nd)) {
        curr_slice_to_do.push(c);
    }
    std::set<std::size_t> boundaries;
    while (!curr_slice_to_do.empty()) {
        node_t * curr = curr_slice_to_do.top();
        curr_slice_to_do.pop();
        const auto & currd = curr->get_data();
        if (currd.ini_stat == IniNodePhyloStatus::NO_COMMON_CLADE) {
            for (auto c : iter_child(*curr)) {
                curr_slice_to_do.push(c);
            }
        } else {
            if (currd.ini_stat == IniNodePhyloStatus::UNKNOWN_STATUS) {
                throw OTCError() << "Node \"" << currd.node_index << "\" in a slice but with uninitialized ini_stat field.\n";
            }
            boundaries.insert(currd.node_index);
            if (!curr->is_tip()) {
                to_do.push(curr);
            }
        }
    }
    std::set<std::size_t> leaf1_inds, leaf2_inds;
    outer["tree_1"] = get_tree_for_slice(t1nd, tas_1, leaf1_inds, boundaries);
    outer["tree_2"] = get_tree_for_slice(t2nd, tas_2, leaf2_inds, boundaries);
    assert(leaf1_inds == leaf2_inds);
    assert(leaf1_inds == boundaries);
    const std::string & t1n = outer["tree_1"]["newick"].get<std::string>();
    const std::string & t2n = outer["tree_2"]["newick"].get<std::string>();
    outer["comparison"] = get_comparison(t1n, t1nd, tas_1, t2n, t2nd, tas_2, boundaries);
    int npi = outer["comparison"]["num_prunings"].get<int>();
    assert(npi > -1);
    if (npi == 0) {
        json rewrite = json::object();
        rewrite["both_trees"] = outer["tree_1"];
        return json_num_edits_pair{rewrite, 0};
    }
    std::size_t npst = npi;
    return json_num_edits_pair{outer, npst};
}

void explain_phylo_diffs(std::ostream & out,
                         Tree_t & inp_tre1,
                         Tree_t & inp_tre2) {
    TreeAsUIntSplits tas_1{inp_tre1};
    TreeAsUIntSplits tas_2{inp_tre2};
    if (tas_1.leaf_label_to_ind != tas_2.leaf_label_to_ind) {
        throw OTCError() << "trees must have the same leaf label set.\n";
    }
    
    std::size_t num_tips = tas_1.ind_to_nd.size();
    
    std::size_t matched_node_index = tas_1.ind_to_nd.size();
    tas_1.ind_to_nd.reserve(3*matched_node_index);
    tas_2.ind_to_nd.reserve(3*matched_node_index);
    std::list<node_t *> unmatched;

    for (auto nd1 : all_nodes_postorder(inp_tre1)) {
        auto & data1 = nd1->get_data();
        data1.ini_stat = IniNodePhyloStatus::IDENTICAL_SUBTREE;
        const auto & t1_ind_set = tas_1.nd_to_taxset[nd1]; 
        if (nd1->is_tip()) {
            data1.node_index = *t1_ind_set.begin();
            std::cerr << data1.node_index << " " << nd1->get_name() << '\n';
            data1.partner = const_cast<node_t *>(tas_2.ind_to_nd[data1.node_index]);
        } else {
            if (nd1 == tas_1.root) {
                data1.partner = const_cast<node_t *>(tas_2.root);
            } else {
                const auto t2_ind_set_it = tas_2.inf_taxset_to_nd.find(t1_ind_set);
                if (t2_ind_set_it == tas_2.inf_taxset_to_nd.end()) {
                    data1.ini_stat = IniNodePhyloStatus::NO_COMMON_CLADE;
                    unmatched.push_back(nd1);
                } else {
                    data1.partner = const_cast<node_t *>(t2_ind_set_it->second);
                }
            }
            if (data1.partner != nullptr) {
                for (auto c : iter_child_const(*nd1)) {
                    const auto & cdata = c->get_data();
                    if (cdata.ini_stat != IniNodePhyloStatus::IDENTICAL_SUBTREE) {
                        data1.ini_stat = IniNodePhyloStatus::SAME_CLADE_DIFF_TOPO;
                        break;
                    }
                }
            }
            if (data1.ini_stat != IniNodePhyloStatus::NO_COMMON_CLADE) {
                data1.node_index = matched_node_index++;
            }
        }

        if (data1.ini_stat == IniNodePhyloStatus::IDENTICAL_SUBTREE
            || data1.ini_stat == IniNodePhyloStatus::SAME_CLADE_DIFF_TOPO) {
            assert(data1.partner);
            auto & data2 = data1.partner->get_data();
            data2.partner = nd1;
            data2.ini_stat = data1.ini_stat;
            data2.node_index = data1.node_index;
            if (data1.node_index >= num_tips) {
                assert(tas_1.ind_to_nd.size() == data1.node_index);
                assert(tas_2.ind_to_nd.size() == data1.node_index);
                tas_1.ind_to_nd.push_back(nd1);
                tas_2.ind_to_nd.push_back(data2.partner);
            }
        }
    }
    for (auto ndp : unmatched) {
        auto & data1 = ndp->get_data();
        data1.node_index = matched_node_index++;
        assert(tas_1.ind_to_nd.size() == data1.node_index);
        assert(tas_2.ind_to_nd.size() == data1.node_index);
        tas_1.ind_to_nd.push_back(ndp);
        tas_2.ind_to_nd.push_back(nullptr);
    }
    for (auto nd2 : all_nodes_postorder(inp_tre2)) {
        auto & data2 = nd2->get_data();
        if (data2.ini_stat == IniNodePhyloStatus::UNKNOWN_STATUS) {
            data2.ini_stat = IniNodePhyloStatus::NO_COMMON_CLADE;
            data2.node_index = matched_node_index++;
            assert(tas_1.ind_to_nd.size() == data2.node_index);
            assert(tas_2.ind_to_nd.size() == data2.node_index);
            tas_1.ind_to_nd.push_back(nullptr);
            tas_2.ind_to_nd.push_back(nd2);
        }
    }
    const auto & rd = tas_1.root->get_data();

    node_t * root1 = const_cast<node_t *>(tas_1.root);
    json document;
    document["root_id"] = to_string(root1->get_data().node_index);
    json tree_slices = json::object();
    std::stack<node_t *> to_do;
    to_do.push(root1);
    int tot_num_prunings = 0;
    while (!to_do.empty()) {
        auto next = to_do.top();
        to_do.pop();
        auto next_id_str = to_string(next->get_data().node_index);
        auto jnpp = get_tree_comp_slice(next, tas_1, tas_2, to_do);
        tree_slices[next_id_str] = jnpp.first;
        tot_num_prunings += jnpp.second;
    }
    document["tree_comp_slices_by_root"] = tree_slices;
    document["tot_num_prunings_all_slices"] = tot_num_prunings;
    out << document.dump() << std::endl;
}

int main(int argc, char *argv[]) {
    try {
        variables_map args = parse_cmd_line(argc, argv);
        string first = args["first"].as<string>();
        string second = args["second"].as<string>();
        string jsonfp;
        if (args.count("out")) {
            jsonfp = args["out"].as<string>();
        }
        std::ostream * outptr = &(std::cout);
        std::ofstream jsonoutstream;
        if (jsonfp.size() > 0) {
            jsonoutstream.open(jsonfp);
            if (!jsonoutstream.good()) {
                std::cerr << "Could not open \"" << jsonfp << "\"\n";
                return 1;
            }
            outptr = &jsonoutstream;
        }
        ParsingRules rules;
        // 1. Load and process summary tree.
        auto fir_trees = get_trees<Tree_t>(first, rules);
        auto sec_trees = get_trees<Tree_t>(second, rules);
        std::cerr << fir_trees.size() << " trees in " << first << std::endl;
        std::cerr << sec_trees.size() << " trees in " << second << std::endl;
        if (fir_trees.size() != sec_trees.size()) {
            std::cerr << "otc-explain-phylo-diffs: Error tree files must have the same number of trees" << std::endl;
            exit(1);
        }

        for (std::size_t ind = 0 ; ind < fir_trees.size(); ++ind) {
            if (ind > 0 && jsonfp.size() > 0) {
                string nfp = jsonfp;
                nfp += "-";
                nfp += std::to_string(ind);
                jsonoutstream.close();
                jsonoutstream.open(nfp);
                if (!jsonoutstream.good()) {
                    std::cerr << "Could not open \"" << nfp << "\"\n";
                    return 1;
                }
                outptr = &jsonoutstream;
            }
            Tree_t & fir_tree = *(fir_trees.at(ind));
            Tree_t & sec_tree = *(sec_trees.at(ind));
            explain_phylo_diffs(*outptr, fir_tree, sec_tree);
        }
    } catch (std::exception& e) {
        std::cerr << "otc-explain-phylo-diffs: Error! " << e.what() << std::endl;
        exit(1);
    }
}
