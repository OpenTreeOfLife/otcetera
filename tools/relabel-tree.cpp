#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <regex>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/tokenizer.hpp>
#include <bitset>
#include <set>
#include <map>
#include "json.hpp"

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;
using std::map;
using std::ifstream;
using json = nlohmann::json;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Node_t = RootedTreeNode<RTNodeNoData>;
using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

using Cause2IDSetMap = map<string, OttIdSet>;
Cause2IDSetMap globalCauseToPrunedMap;
using Cause2StrSetMap = map<string, std::set<string> >;
Cause2StrSetMap globalCauseToPrunedStrMap;

void pruneHigherTaxonTips(Tree_t& tree, const Taxonomy& taxonomy);

variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("tree", value<string>(),"Filename for the tree")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<OttId>(), "OTT id of root node of subtree to keep")
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;

    options_description relabel("Relabelling options");
    relabel.add_options()
        ("format-tax",value<string>()->default_value("%L"),"Form of labels to write for taxonomy nodes.")
        ("format-unknown",value<string>()->default_value("%L"),"Form of labels to write for non-taxonomy nodes.")
        ("remap",value<string>(),"Filepath for a tab-separated mapping of label to ott id")
        ("replace,R",value<string>(),"Perform a regex replacement on all labels")
//        ("label-regex", value<string>(), "Return name of the given ID")
        ;

    options_description tree("Tree options");
    tree.add_options()
        ("del-monotypic","Remove monotypic nodes.")
        ("del-higher-taxon-tips","Prune the tree such that no tips mapped to taxa that have ranks above the species level are tips.")
        ("prune-flags",value<string>(),"Comma-separate list of flags to prune from tree")
        ("filter-flags",value<string>(),"Comma-separate list of flags to filter from tree")
        ;

    options_description output("Output options");
    output.add_options()
        ("json,j", value<string>(), "filepath to an output JSON log")
        ;

    options_description visible;
    visible.add(taxonomy).add(relabel).add(tree).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("tree", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-relabel-tree <newick-file> [OPTIONS]\n"
                                                    "Rewrite node labels for a tree.\n\n"
                                                    "Format strings have the following interpretation:\n"
                                                    "  %I=id %N=name %U=uniqname %R=rank %S=sourceinfo %%=% %L=original label",
                                                    visible, invisible, p);
    return vm;
}

std::size_t n_nodes(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    std::size_t count = 0;
    for(auto nd: iter_post_const(T)) {
        count++;
    }
    return count;
}

void pruneHigherTaxonTips(Tree_t& tree, const Taxonomy& taxonomy) {
    auto & prunedTipSet = globalCauseToPrunedMap[string("higher-taxon-tip")];
    // first we collect the tips to prune
    set<Tree_t::node_type*> tipsToPrune;
    for (auto nd: iter_leaf(tree)) {
        if (not nd->has_ott_id()) {
            throw OTCError() << "Tip \"" << nd->get_name() << "\" lacks an OTT ID.";
        }
        const auto ottId = nd->get_ott_id();
        const auto & taxonrecord = taxonomy.record_from_id(ottId);
        const auto & rank = taxonrecord.rank;
        if (rank.length() > 0 && (rank != "no rank"
                                  && rank != "no rank - terminal"
                                  && rank != "forma"
                                  && rank != "species"
                                  && rank != "subspecies"
                                  && rank != "subvariety"
                                  && rank != "varietas"
                                  && rank != "variety")) {
            tipsToPrune.insert(nd);
            prunedTipSet.insert(ottId);
        }
    }
    const auto root = tree.get_root();
    // now delete the tips and put their parents in a queue to check (to see if they've become tips)
    set<Tree_t::node_type*> toCheckNext;
    for (auto nd: tipsToPrune) {
        if (nd == root) {
            throw OTCError() << "Please do not call this program with a dot tree.";
        }
        toCheckNext.insert(nd->get_parent());
        nd->detach_this_node();
    }
    auto & prunedInternalSet = globalCauseToPrunedMap[string("empty-after-higher-taxon-tip-prune")];
    while (!toCheckNext.empty()) {
        set<Tree_t::node_type*> parSet;
        for (auto nd: toCheckNext) {
            if (!nd->has_children()) {
                if (nd == root) {
                    throw OTCError() << "The tree was pruned to non-existence!";
                }
                if (not nd->has_ott_id()) {
                    throw OTCError() << "Node \"" << nd->get_name() << "\" lacks an OTT ID.";
                }
                const auto ottId = nd->get_ott_id();
                prunedInternalSet.insert(ottId);
                parSet.insert(nd->get_parent());
                nd->detach_this_node();
            }
        }
        toCheckNext.swap(parSet);
    }
}

void filterTreeByFlags(Tree_t& tree, const Taxonomy& taxonomy, std::bitset<32> prune_flags) {
    vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    for(auto nd: nodes) {
        if (not nd->has_ott_id()) {
            continue;
        }
        auto id = nd->get_ott_id();
        if ((taxonomy.record_from_id(id).flags & prune_flags).any()) {
            if (nd == tree.get_root()) {
                if (nd->is_outdegree_one_node()) {
                    auto newroot = nd->get_first_child();
                    newroot->detach_this_node();
                    tree._set_root(newroot);
                } else {
                    throw OTCError() << "The root has flags set for pruning, but is not monotypic!";
                }
            }
            while (nd->has_children()) {
                auto c = nd->get_first_child();
                c->detach_this_node();
                nd->add_sib_on_left(c);
            }
            nd->detach_this_node();
        }
    }
}

void pruneTreeByFlags(Tree_t& tree, const Taxonomy& taxonomy, std::bitset<32> prune_flags) {
    vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    for(auto nd: nodes) {
        if (not nd->has_ott_id()) {
            continue;
        }
        auto id = nd->get_ott_id();
        if ((taxonomy.record_from_id(id).flags & prune_flags).any()) {
            if (nd == tree.get_root()) {
                throw OTCError() << "The root has flags set for pruning!";
            }
            nd->detach_this_node();
        }
    }
}

template<typename C>
void fillCause2SetMap(const C &csmap, json & document) {
    for (auto csp: csmap) {
        const auto & cause = csp.first;
        const auto & idSet = csp.second;
        if (idSet.empty()) {
            continue;
        }
        json id_array = json::array();
        for (auto oid: idSet) {
            id_array.push_back(oid);
        }
        document[cause] = id_array;
    }
}

void writeLog(std::ofstream & out, const Cause2IDSetMap &csmap, const Cause2StrSetMap &cs2strmap) {
    json document;
    fillCause2SetMap(csmap, document);
    fillCause2SetMap(cs2strmap, document);
    out << document.dump(1) << std::endl;
}

long long_ott_id_from_numeric_string(const string & s) {
    long conv = -2;
    auto r = char_ptr_to_long(s.c_str(), &conv);
    assert(r);
    return conv;
}

map<string, OttId> read_mapping_file(const string & mf) {
    ifstream inp(mf);
    if (!inp.good()) {
        throw OTCError() << "Could not open the mapping file \"" << mf << "\".";
    }
    map<string, OttId> ret;
    unsigned int line_num = 1;
    for (string next_line; getline(inp, next_line); ++line_num) {
        const auto words = split_string(next_line, '\t');
        if (words.size() == 0) {
            continue;
        }
        if (words.size() != 2) {
            throw OTCError() << "Expecting  one tab in each line. Problem with line " << line_num << ": \"" << next_line << "\".";
        }
        const auto & label = *words.begin();
        if (ret.count(label)) {
            throw OTCError() << "Repeated label in mapping file: \"" << label << "\" on line " << line_num;
        }
        auto raw_ott_id = long_ott_id_from_numeric_string(*words.rbegin());
        if (raw_ott_id < 0) {
            throw OTCError() << "Negative, missing, or out of range OTT ID in mapping file: \"" << *words.rbegin() << "\" on line " << line_num;
        }
        OttId ott_id = check_ott_id_size(raw_ott_id);
        ret[label] = ott_id;
    }
    return ret;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::ofstream jlogf;
    std::ofstream * json_log = nullptr;
    try {
        variables_map args = parse_cmd_line(argc,argv);
        if (not args.count("tree")) {
            throw OTCError() << "Please specify the newick tree to be relabelled!";
        }
        if (args.count("json")) {
            string jf = args["json"].as<string>();
            jlogf.open(jf);
            if (!jlogf.good()) {
                throw OTCError() << "Could not open JSON log file at \"" << jf << "\"";
            }
            json_log = &jlogf;
        }
        string format_tax = args["format-tax"].as<string>();
        string format_unknown = args["format-unknown"].as<string>();
        boost::optional<Taxonomy> taxonomy = boost::none;
        const bool do_prune_flags = args.count("prune-flags");
        const bool do_filter_flags = args.count("filter-flags");
        const bool do_prune_higher = args.count("del-higher-taxon-tips");
        const bool needs_taxonomy = (format_needs_taxonomy(format_tax)
                                     or do_prune_flags
                                     or do_filter_flags
                                     or do_prune_higher);
        if (needs_taxonomy) {
            taxonomy = load_taxonomy(args);
        }
        if (char c = format_needs_taxonomy(format_unknown)) {
            throw OTCError() << "Cannot use taxonomy-based specifier '%" << c << "' in non-taxonomy format string '" << format_unknown << "'";
        }
        auto tree = get_tree<Tree_t>(args["tree"].as<string>());
        if (args.count("remap")) {
            string mapping_file = args["remap"].as<string>();
            const auto label_to_ott_id = read_mapping_file(mapping_file);
            std::cerr << label_to_ott_id.size() << " records read from " << mapping_file << endl;
            std::set<Node_t *> to_detach;
            for (auto nd: iter_leaf(*tree)) {
                auto n = nd->get_name();
                auto rmit = label_to_ott_id.find(n);
                if (rmit == label_to_ott_id.end()) {
                    LOG(INFO) << "Tip \"" << n << "\" not found in mapping file, will be deleted.";
                    to_detach.insert(nd);
                } else {
                    string nn = "ott";
                    nn += std::to_string(rmit->second);
                    nd->set_name(nn);
                }
            }
            int round = 0;
            auto root = tree->get_root();
            auto & prunedTipSet = globalCauseToPrunedStrMap[string("unmapped-tips")];
            auto & barrenInternal = globalCauseToPrunedStrMap[string("barren-internal")];
            while (!to_detach.empty()) {
                std::set<Node_t *> next_to_detach;
                for (auto n : to_detach) {
                    auto p = n->get_parent();
                    if (round == 0) {
                        LOG(INFO) << "detaching leaf" << n->get_name();
                        prunedTipSet.insert(n->get_name());
                    } else {
                        LOG(INFO) << "detaching internal " << n->get_name();    
                        barrenInternal.insert(n->get_name());
                    }
                    assert(!n->has_children());
                    n->detach_this_node();
                    if (p && !p->has_children()) {
                        next_to_detach.insert(p);
                    }
                    tree->mark_as_detached(n);
                }
            round += 1;
            to_detach = next_to_detach;
            }
        } else {
            if (do_prune_flags) {
                auto flags = flags_from_string(args["prune-flags"].as<string>());
                pruneTreeByFlags(*tree, *taxonomy, flags);
            }
            if (do_filter_flags) {
                auto flags = flags_from_string(args["filter-flags"].as<string>());
                filterTreeByFlags(*tree, *taxonomy, flags);
            }
            if (args.count("del-monotypic")) {
                suppress_monotypic_fast(*tree);
            }
            if (do_prune_higher) {
                pruneHigherTaxonTips(*tree, *taxonomy);
            }
            for(auto nd: iter_pre(*tree)) {
                string name;
                if (nd->has_ott_id()) {
                    int id = nd->get_ott_id();
                    if (taxonomy) {
                        const auto& record = (*taxonomy).record_from_id(id);
                        name = format_with_taxonomy(nd->get_name(), format_tax, record, *taxonomy);
                    } else {
                        name = format_without_taxonomy(nd->get_name(), format_tax);
                    }
                } else {
                    name = format_without_taxonomy(nd->get_name(), format_unknown);
                }
                nd->set_name(std::move(name));
            }
            if (args.count("replace")) {
                string match_replace = args["replace"].as<string>();
                if (match_replace.empty()) {
                    throw OTCError() << "Empty pattern for argument 'replace'!";
                }
                char sep = match_replace[0];
                if (std::count(match_replace.begin(), match_replace.end(), sep) !=3) {
                    throw OTCError() << "Delimiter '" << sep << "' does not occur exactly three times in '" << match_replace << "'!";
                }
                if (match_replace.back() != sep) {
                    throw OTCError() << "Pattern '" << match_replace << "' does not end with delimiter '" << sep << "'";
                }
                int loc2 = match_replace.find(sep, 1);
                int loc3 = match_replace.find(sep, loc2 + 1);
                std::regex match (match_replace.substr(1, loc2 - 1));
                string replace = match_replace.substr(loc2 + 1, loc3 - loc2 - 1);
                for(auto nd: iter_pre(*tree)) {
                    string name = nd->get_name();
                    name = std::regex_replace(name,match,replace);
                    nd->set_name(name);
                }
            }
        }
        write_tree_as_newick(cout, *tree);
        std::cout << std::endl;
        if (json_log && (!globalCauseToPrunedMap.empty() || !globalCauseToPrunedStrMap.empty())) {
            writeLog(*json_log, globalCauseToPrunedMap, globalCauseToPrunedStrMap);
        }
    }
    catch (std::exception& e) {
        cerr << "otc-relabel-tree: Error! " << e.what() << std::endl;
        exit(1);
    }
}
