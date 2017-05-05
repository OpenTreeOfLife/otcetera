#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <bitset>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::set;
using std::endl;
using std::bitset;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Node_t = RootedTreeNode<RTNodeIncludeBoolData>;
using Tree_t = RootedTree<RTNodeIncludeBoolData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("in-tree",value<string>(),"Tree of OTT ids in tree")
        ("json",value<string>(),"filepath for an output file of the pruning points")
        ;

    options_description visible;
    visible.add(taxonomy).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-regraft-taxonomy-generator <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy, filter all taxa according to the cleaning_flags or additional_regrafting_flags, but only if they are not ancestors of a tip in the input tree",
                                                    visible, invisible, p);

    return vm;
}

OttIdSet get_ids_from_tree(const string& filename) {
    OttIdSet ids;
    auto tree = get_tree<Tree_t>(filename);
    for(auto nd: iter_leaf(*tree)) {
        auto id = nd->get_ott_id();
        if (id != -1) {
            ids.insert(id);
        }
    }
    return ids;
}


OttIdSet prune_if_flagged_and_not_anc(const Taxonomy & taxonomy,
                                  Tree_t & taxo_tree,
                                  const OttIdSet & ids,
                                  const tax_flags regraft_filters) {

    for (auto nd : iter_post(taxo_tree)) {
        auto oid = nd->get_ott_id();
        auto & ndd = nd->get_data();
        if (ids.count(oid) > 0) {
            ndd.include = true;
            for (auto a : iter_anc(*nd)) {
                auto & ad = a->get_data();
                if (ad.include) {
                    break;
                }
                ad.include = true;
            }
        }
        if (!ndd.include) {
            // not an ancestor of any ID in ids, check flags
            const auto & tax_record = taxonomy.record_from_id(oid);
            const bool should_prune = (tax_record.flags & regraft_filters).any();
            if (!should_prune) {
                ndd.include = true;
            }
        }
    }
    set<Node_t *> prune_points;
    for (auto nd : iter_post(taxo_tree)) {
        auto & ndd = nd->get_data();
        if (!ndd.include) {
            auto pp = nd;
            for (auto a : iter_anc(*nd)) {
                auto & ad = a->get_data();
                if (ad.include) {
                    break;
                }
                pp = a;
            }
            prune_points.insert(pp);
        }
    }
    OttIdSet prune_ids;
    for (auto to_prune : prune_points) {
        prune_ids.insert(to_prune->get_ott_id());
        to_prune->detach_this_node();
        taxo_tree.mark_as_detached(to_prune);
    }
    return prune_ids;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc,argv);
        if (!args.count("in-tree")) {
            cerr << "otc-regraft-taxonomy-generator: Error! Expected an in-tree argument!" << std::endl;
            return 2;
        }
        tax_flags regraft_filters;
        if (args.count("config")) {
            regraft_filters = regrafting_flags_from_config_file(args["config"].as<string>());
        } else {
            cerr << "otc-regraft-taxonomy-generator: Error! Expected a confg argument!" << std::endl;
            return 2;
        }
        auto taxonomy = load_taxonomy(args);
        auto ids = get_ids_from_tree(args["in-tree"].as<string>());
        auto nodeNamer = [](const auto& record){return string(record.name)+"_ott"+std::to_string(record.id);};
        auto taxo_tree = taxonomy.get_tree<Tree_t>(nodeNamer);
        auto pruned_ids = prune_if_flagged_and_not_anc(taxonomy, *taxo_tree, ids, regraft_filters);
        write_tree_as_newick(cout, *taxo_tree);
        std::cout << std::endl;
        
        // hacky write of pruned_ids to std::cerr
        if (args.count("json")) {
            string jfp = args["json"].as<string>();
            std::ofstream jout(jfp);
            jout << "{\n  \"pruned\": [";
            bool first = true;
            for (auto i : pruned_ids) {
                if (first) {
                    first = false;
                } else {
                    jout << ", ";
                }
                jout << i;
            }
            jout << "]\n}" << std::endl;
        }
    } catch (std::exception& e) {
        cerr << "otc-taxonomy-parser: Error! " << e.what() << std::endl;
        return 1;
    }
}

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

