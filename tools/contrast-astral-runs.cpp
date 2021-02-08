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

using json=nlohmann::json;


using namespace otc;
using std::vector;
template <typename T, typename U>
using Map = std::unordered_multimap<T,U>;
template <typename T, typename U>
using map = std::unordered_map<T,U>;
template <typename T>
using set = std::unordered_set<T>;
using std::string;
//using std::map;
using std::pair;
using std::tuple;
using std::unique_ptr;
using std::string;

namespace po = boost::program_options;
using po::variables_map;

// TODO: Could we exemplify tips here, if we had access to the taxonomy?
variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("first", value<vector<string>>()->composing(),"First assembled astral tree file (output followed by inputs in newick).")
        ("second", value<vector<string>>()->composing(),"Second assembled astral tree file (output followed by inputs in newick).")
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
                                                    "Usage: otc-contrast-astral-runs <first> <second> [OPTIONS]\n"
                                                    "Calculates conflict stats for astral runs.\n\n"
                                                    "Input files should be newick. Astral output first, followed by the inputs.\n"
                                                    "all trees rooted, and the order of gene trees identical between inputs.",
                                                    visible, invisible, p);

    return vm;
}

namespace std {

template<typename T, typename U>
struct hash<pair<T,U>> {
    std::size_t operator()(const std::pair<T,U>& p) const noexcept {return std::hash<T>()(p.first) * std::hash<U>()(p.second);}
};

} // namespace std

using Tree_t = ConflictTree;
using node_t = Tree_t::node_type;

using str_set = std::set<std::string>;
    
class TreeAsSplits {
    public:
    std::map<string, const node_t *> leaf_label_to_nd;
    std::map<const node_t *, str_set> nd_to_taxset;
    std::map<str_set, const node_t *> inf_taxset_to_nd;
    str_set leaf_labels;

    TreeAsSplits(const Tree_t & tree) {
        auto nodes = all_nodes_postorder(tree);
        for (auto nd : nodes) {
            str_set ts;
            if (nd->is_tip()) {
                auto & nl = nd->get_name();
                //std::cerr << nl << '\n';
                ts.insert(nl);
                leaf_labels.insert(nl);
                leaf_label_to_nd[nl] = nd;
            } else {
                for (auto c : iter_child_const(*nd)) {
                    auto cts = nd_to_taxset[c];
                    ts.insert(cts.begin(), cts.end());
                }
                if (nd->get_parent() != nullptr) {
                    inf_taxset_to_nd[ts] = nd;
                }
            }
            nd_to_taxset[nd] = ts;
        }
    }
};


enum in_full_status {IN_NEITHER = 0,
                     IN_ONE = 1,
                     IN_TWO =2,
                     IN_BOTH = 3
                    };

enum CONFLICT_STATUS {IRRELEVANT = 0,
                      CONFLICTS = 1,
                      COMPATIBLE = 2,
                      SUPPORTS = 6    // COMPATIBLE + PARENT's unions bigger
                     };

class SummaryOfSplit {
    public:
    using conf_pair = std::pair<CONFLICT_STATUS, CONFLICT_STATUS>;
    using list_conf = std::list<conf_pair>;
    in_full_status in_full;
    list_conf in_inputs;
    
    SummaryOfSplit(const in_full_status & full_stat_arg) 
      :in_full(full_stat_arg) {
    }
};

using split_to_summary = std::map<str_set, SummaryOfSplit>;


CONFLICT_STATUS calc_conf_status(const TreeAsSplits &full_tas, 
                                 const TreeAsSplits &inp_tas,
                                 const str_set & full_tax_set,
                                 const str_set & inp_tax_set,
                                 const str_set & restricted,
                                 bool rooted) {
    if (contains(inp_tas.inf_taxset_to_nd, restricted)) {
        return CONFLICT_STATUS::SUPPORTS;
    }
    if (rooted) {
        for (const auto & ts2nd_el : inp_tas.inf_taxset_to_nd) {
            const auto & inp_nd_ts = ts2nd_el.first;
            if (!have_intersection(inp_nd_ts, restricted)) {
                continue;
            }
            if (inp_nd_ts == restricted) {
               return CONFLICT_STATUS::SUPPORTS; // should be caught above...
            }
            if (is_subset(inp_nd_ts, restricted)) {
                continue;
            }
            if (is_subset(restricted, inp_nd_ts)) {
                continue;
            }
            return CONFLICT_STATUS::CONFLICTS;
        }
    } else {
        auto other = set_difference_as_set(inp_tax_set, restricted);
        if (contains(inp_tas.inf_taxset_to_nd, other)) {
            return CONFLICT_STATUS::SUPPORTS;
        }
        for (const auto & ts2nd_el : inp_tas.inf_taxset_to_nd) {
            const auto & inp_nd_ts = ts2nd_el.first;
            if (!have_intersection(inp_nd_ts, restricted)) {
                continue;
            }
            if (!have_intersection(inp_nd_ts, other)) {
                continue;
            }
            if (is_subset(inp_nd_ts, restricted)) {
                continue;
            }
            if (is_subset(restricted, inp_nd_ts)) {
                continue;
            }
            if (is_subset(inp_nd_ts, other)) {
                continue;
            }
            if (is_subset(other, inp_nd_ts)) {
                continue;
            }
            return CONFLICT_STATUS::CONFLICTS;
        }
    }
    return CONFLICT_STATUS::IRRELEVANT;
}

void analyze_inp_tree_pair(const TreeAsSplits & tas_1,
                           const TreeAsSplits & tas_2,
                           const Tree_t & inp_tre1,
                           const Tree_t & inp_tre2, 
                           split_to_summary & s2ss,
                           unsigned int tree_ind,
                           bool rooted) {
    const TreeAsSplits tas_inp1(inp_tre1);
    const TreeAsSplits tas_inp2(inp_tre2);
    const auto & inp_labels = tas_inp1.leaf_labels;
    if (tas_inp2.leaf_labels != inp_labels) {
        auto d = set_sym_difference_as_set(tas_inp2.leaf_labels, inp_labels);
        std::cerr << "otc-contrast-astral-runs: Some labels uniq to one tree from index " << tree_ind << std::endl;
        std::cerr << " size of diff set = " << d.size() << "\n";
        for (auto item : d) {
            std::cerr << "  " << item << "\n";
        }
        exit(1);
    }
    unsigned int split_index = 0;
    for (auto & m_el: s2ss) {
        std::cerr << "tree " << tree_ind << " split " << split_index++ << std::endl;
        const auto & full_tax_set = m_el.first;
        auto & sos = m_el.second;
        auto restricted = set_intersection_as_set(full_tax_set, inp_labels);
        if (restricted.size() == inp_labels.size() || restricted.size() < 2) {
            sos.in_inputs.emplace_back(CONFLICT_STATUS::IRRELEVANT, CONFLICT_STATUS::IRRELEVANT);
        } else {
            auto cs1 = calc_conf_status(tas_1, tas_inp1, full_tax_set, inp_labels, restricted, rooted);
            auto cs2 = calc_conf_status(tas_2, tas_inp2, full_tax_set, inp_labels, restricted, rooted);
            sos.in_inputs.emplace_back(cs1, cs2);
        }
    }
}

inline char conf_stat_out(const CONFLICT_STATUS cs) {
    if (cs == CONFLICT_STATUS::IRRELEVANT) {
        return 'I';
    }
    if (cs == CONFLICT_STATUS::CONFLICTS) {
        return 'X';
    }
    if (cs == CONFLICT_STATUS::SUPPORTS) {
        return 'S';
    }
    return 'C';
}

void report_summaries(const split_to_summary &s2ss, 
                     const std::string & json_fp,
                     const std::string & summary_tsv) {
    auto split_key_content = json::array();
    unsigned int index = 0;
    std::ofstream summ_out;
    summ_out.open(summary_tsv);
    unsigned shared = 0;
    unsigned unshared = 0;
    for (const auto & m_el: s2ss) {
        const auto & tax_set = m_el.first;
        const auto & sos = m_el.second;
        auto curr_tax_set = json::array();
        for (auto & name : tax_set) {
            curr_tax_set.push_back(name);

        }
        summ_out << index ; 
        if (sos.in_full == in_full_status::IN_BOTH) {
            summ_out << "\t+\t+";
            shared++;
        } else if (sos.in_full == in_full_status::IN_ONE) {
            summ_out << "\t+\t-";
            unshared++;
        } else if (sos.in_full == in_full_status::IN_TWO) {
            summ_out << "\t-\t+";
            unshared++;
        } else {
            std::cerr << "split " << index << " in neither full tree!\n";
            exit(1);
        }
        for (const auto & lit : sos.in_inputs) {
            summ_out << '\t' << conf_stat_out(lit.first) << '/' << conf_stat_out(lit.second);
        }
        summ_out << '\n';
        index += 1;
        split_key_content.push_back(curr_tax_set);
    }
    std::cerr << "full tree symm diff = " << unshared << " / " << unshared + 2*shared << '\n';
    std::ofstream kout;
    kout.open(json_fp);
    kout << split_key_content.dump(1) << std::endl;
    kout.close();
}
int main(int argc, char *argv[]) {
    try {
        variables_map args = parse_cmd_line(argc, argv);
        vector<string> fir_astral = args["first"].as<vector<string> >();
        vector<string> sec_astral = args["second"].as<vector<string> >();
        ParsingRules rules;
        // 1. Load and process summary tree.
        auto fir_trees = get_trees<Tree_t>(fir_astral, rules);
        auto sec_trees = get_trees<Tree_t>(sec_astral, rules);
        std::cerr << fir_trees.size() << " trees in " << fir_astral[0] << std::endl;
        std::cerr << sec_trees.size() << " trees in " << sec_astral[0] << std::endl;
        if (fir_trees.size() != sec_trees.size()) {
            std::cerr << "otc-contrast-astral-runs: Error tree files must have the same number of trees" << std::endl;
            exit(1);
        }
        if (fir_trees.size() < 2) {
            std::cerr << "otc-contrast-astral-runs: Error tree files must have at least 2 trees" << std::endl;
            exit(1);
        }
        const Tree_t & fir_sp_tree = *(fir_trees[0]);
        const TreeAsSplits tas_1(fir_sp_tree);
        const Tree_t & sec_sp_tree = *(sec_trees[0]);
        const TreeAsSplits tas_2(sec_sp_tree);
        if (tas_2.leaf_labels != tas_1.leaf_labels) {
            auto d = set_sym_difference_as_set(tas_2.leaf_labels, tas_1.leaf_labels);
            std::cerr << "otc-contrast-astral-runs: Some labels uniq to one tree:" << std::endl;
            std::cerr << " size of diff set = " << d.size() << "\n";
            for (auto item : d) {
                std::cerr << "  " << item << "\n";
            }
            exit(1);
        }
        const auto & full_leaf_set = tas_1.leaf_labels;
        const auto & to_inf_1 = tas_1.inf_taxset_to_nd;
        const auto & to_inf_2 = tas_2.inf_taxset_to_nd;
        auto all_splits = set_union_as_set(keys(to_inf_1), keys(to_inf_2));
        split_to_summary s2ss;
        for (const auto & sp : all_splits) {
            if (contains(to_inf_1, sp)) {
                if (contains(to_inf_2, sp)) {
                    s2ss.emplace(sp, in_full_status::IN_BOTH);
                } else {
                    s2ss.emplace(sp, in_full_status::IN_ONE);
                }
            } else {
                s2ss.emplace(sp, in_full_status::IN_TWO);
            }
        }

        unsigned inp_tree_index = 1;
        while (inp_tree_index < fir_trees.size()) {
            const auto & inp_tre1 = *(fir_trees[inp_tree_index]);
            const auto & inp_tre2 = *(sec_trees[inp_tree_index]);
            analyze_inp_tree_pair(tas_1, tas_2, inp_tre1, inp_tre2, s2ss, inp_tree_index, false);
            inp_tree_index++;
        }
        report_summaries(s2ss, "splits-key.json", "split-summary.tsv");


    } catch (std::exception& e) {
        std::cerr << "otc-contrast-astral-runs: Error! " << e.what() << std::endl;
        exit(1);
    }
}
