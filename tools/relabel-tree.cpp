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
#include <boost/utility/string_ref.hpp>
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
using json = nlohmann::json;

using boost::spirit::qi::symbols;
using boost::string_ref;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

typedef map<string, OttIdOSet> Cause2IDSetMap;
Cause2IDSetMap globalCauseToPrunedMap;

void pruneHigherTaxonTips(Tree_t& tree, const Taxonomy& taxonomy);

variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("tree", value<string>(),"Filename for the taxonomy")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<long>(), "OTT id of root node of subtree to keep")
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;

    options_description relabel("Relabelling options");
    relabel.add_options()
        ("format-tax",value<string>()->default_value("%L"),"Form of labels to write for taxonomy nodes.")
        ("format-unknown",value<string>()->default_value("%L"),"Form of labels to write for non-taxonomy nodes.")
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

long n_nodes(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    long count = 0;
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
        if (not nd->hasOttId()) {
            throw OTCError() << "Tip \"" << nd->getName() << "\" lacks an OTT ID.";
        }
        const auto ottId = nd->getOttId();
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
    const auto root = tree.getRoot();
    // now delete the tips and put their parents in a queue to check (to see if they've become tips)
    set<Tree_t::node_type*> toCheckNext;
    for (auto nd: tipsToPrune) {
        if (nd == root) {
            throw OTCError() << "Please do not call this program with a dot tree.";
        }
        toCheckNext.insert(nd->getParent());
        nd->detachThisNode();
    }
    auto & prunedInternalSet = globalCauseToPrunedMap[string("empty-after-higher-taxon-tip-prune")];
    while (!toCheckNext.empty()) {
        set<Tree_t::node_type*> parSet;
        for (auto nd: toCheckNext) {
            if (!nd->hasChildren()) {
                if (nd == root) {
                    throw OTCError() << "The tree was pruned to non-existence!";
                }
                if (not nd->hasOttId()) {
                    throw OTCError() << "Node \"" << nd->getName() << "\" lacks an OTT ID.";
                }
                const auto ottId = nd->getOttId();
                prunedInternalSet.insert(ottId);
                parSet.insert(nd->getParent());
                nd->detachThisNode();
            }
        }
        toCheckNext.swap(parSet);
    }
}

void filterTreeByFlags(Tree_t& tree, const Taxonomy& taxonomy, std::bitset<32> prune_flags)
{
    vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    for(auto nd: nodes) {
        if (not nd->hasOttId()) {
            continue;
        }
        auto id = nd->getOttId();
        if ((taxonomy.record_from_id(id).flags & prune_flags).any()) {
            if (nd == tree.getRoot()) {
                if (nd->isOutDegreeOneNode()) {
                    auto newroot = nd->getFirstChild();
                    newroot->detachThisNode();
                    tree._setRoot(newroot);
                } else {
                    throw OTCError()<<"The root has flags set for pruning, but is not monotypic!";
                }
            }
            while (nd->hasChildren()) {
                auto c = nd->getFirstChild();
                c->detachThisNode();
                nd->addSibOnLeft(c);
            }
            nd->detachThisNode();
        }
    }
}

void pruneTreeByFlags(Tree_t& tree, const Taxonomy& taxonomy, std::bitset<32> prune_flags)
{
    vector<Tree_t::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    for(auto nd: nodes) {
        if (not nd->hasOttId()) {
            continue;
        }
        auto id = nd->getOttId();
        if ((taxonomy.record_from_id(id).flags & prune_flags).any()) {
            if (nd == tree.getRoot()) {
                throw OTCError()<<"The root has flags set for pruning!";
            }
            nd->detachThisNode();
        }
    }
}

void writeLog(std::ofstream & out, const Cause2IDSetMap &csmap) {
    json document;
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
    out << document.dump(1) << std::endl;
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);
    std::ofstream jlogf;
    std::ofstream * json_log = nullptr;
    try {
        variables_map args = parse_cmd_line(argc,argv);
        if (not args.count("tree")) {
            throw OTCError()<<"Please specify the newick tree to be relabelled!";
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
            throw OTCError()<<"Cannot use taxonomy-based specifier '%"<<c<<"' in non-taxonomy format string '"<<format_unknown<<"'";
        }
        auto tree = get_tree<Tree_t>(args["tree"].as<string>());
        
        if (do_prune_flags) {
            auto flags = flags_from_string(args["prune-flags"].as<string>());
            pruneTreeByFlags(*tree, *taxonomy, flags);
        }
        if (do_filter_flags) {
            auto flags = flags_from_string(args["filter-flags"].as<string>());
            filterTreeByFlags(*tree, *taxonomy, flags);
        }
        if (args.count("del-monotypic")) {
            suppressMonotypicFast(*tree);
        }
        if (do_prune_higher) {
            pruneHigherTaxonTips(*tree, *taxonomy);
        }
        for(auto nd: iter_pre(*tree)) {
            string name;
            if (nd->hasOttId()) {
                int id = nd->getOttId();
                if (taxonomy) {
                    const auto& record = (*taxonomy).record_from_id(id);
                    name = format_with_taxonomy(nd->getName(), format_tax, record);
                } else {
                    name = format_without_taxonomy(nd->getName(), format_tax);
                }
            } else {
                name = format_without_taxonomy(nd->getName(), format_unknown);
            }
            nd->setName(std::move(name));
        }
        if (args.count("replace")) {
            string match_replace = args["replace"].as<string>();
            if (match_replace.empty()) {
                throw OTCError()<<"Empty pattern for argument 'replace'!";
            }
            char sep = match_replace[0];
            if (std::count(match_replace.begin(), match_replace.end(), sep) !=3) {
                throw OTCError()<<"Delimiter '"<<sep<<"' does not occur exactly three times in '"<<match_replace<<"'!";
            }
            if (match_replace.back() != sep) {
                throw OTCError()<<"Pattern '"<<match_replace<<"' does not end with delimiter '"<<sep<<"'";
            }
            int loc2 = match_replace.find(sep, 1);
            int loc3 = match_replace.find(sep, loc2 + 1);
            std::regex match (match_replace.substr(1, loc2 - 1));
            string replace = match_replace.substr(loc2 + 1, loc3 - loc2 - 1);
            for(auto nd: iter_pre(*tree)) {
                string name = nd->getName();
                name = std::regex_replace(name,match,replace);
                nd->setName(name);
            }
        }
        writeTreeAsNewick(cout, *tree);
        std::cout<<std::endl;
        if (json_log && !globalCauseToPrunedMap.empty()) {
            writeLog(*json_log, globalCauseToPrunedMap);
        }
    }
    catch (std::exception& e) {
        cerr<<"otc-relabel-tree: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
