#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <bitset>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/diff_maker.h"
#include "otc/config_file.h"

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
using std::string_view;
using json = nlohmann::json;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("oldtaxonomy", value<string>(),"Filename for the old taxonomy")
        ("newtaxonomy", value<string>(),"Filename for the new taxonomy")
        ;

    options_description output("Output options");
    output.add_options()
        ("write-to-stdout","Primarily for debugging. Writes contents of taxonomy output to stdout. Only used if write-taxonomy is not used.")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("oldtaxonomy", 1);
    p.add("newtaxonomy", 1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-diff-maker <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy and edit JSON files",
                                                    visible, invisible, p);

    return vm;
}

using id2name_t = std::unordered_map<OttId, string_view>;
using name2id_t = std::unordered_map<string_view, OttId>;
using nd2idset_t = std::unordered_map<const RTRichTaxNode *, OttIdSet>;

template<typename T>
bool fill_name_id_maps(const T & tree_data, id2name_t& id2name) {
    for (auto id_nd_pair : tree_data.id_to_node) {
        const auto & nd_ptr = id_nd_pair.second;
        const auto tax_id = id_nd_pair.first;
        if (nd_ptr == nullptr) {
            throw OTCError() << "Unexpected nullptr in id_to_node";
        }
        id2name[tax_id] = nd_ptr->get_name();
    }
    return true;
}

OttIdSet find_ids_with_same_names(const id2name_t & old_id2name, const id2name_t & new_id2name) {
    OttIdSet same_id_name;
    for (auto id_name_pair : old_id2name) {
        OttId tax_id = id_name_pair.first;
        auto new_it = new_id2name.find(tax_id);
        if (new_it != new_id2name.end()) {
            if (new_it->second == id_name_pair.second) {
                same_id_name.insert(tax_id);
            }
        }
    }
    return same_id_name;
}

nd2idset_t fill_des_id_set(const RichTaxTree & tree, const OttIdSet & relevant_ids) {
    nd2idset_t nd2idset;
    for (auto nd : iter_post_const(tree)) {
        const OttId tax_id = nd->get_ott_id();
        const bool is_relevant = contains(relevant_ids, tax_id);
        OttIdSet & dest = nd2idset[nd];
        if (is_relevant) {
            dest.insert(tax_id);
        }
        if (!nd->is_tip()) {
            for (auto child : iter_child_const(*nd)) {
                const auto & child_set = nd2idset.at(child);
                dest.insert(child_set.begin(), child_set.end());
            }
        }
    }
    return nd2idset;
}

bool diff_from_taxonomies(std::ostream & out,
                          const TaxonomyDiffMaker & old_tax,
                          const TaxonomyDiffMaker & new_tax) {
    const auto & old_tree = old_tax.get_tax_tree();
    const auto & new_tree = new_tax.get_tax_tree();
    const auto & old_td = old_tree.get_data();
    const auto & new_td = new_tree.get_data();
    id2name_t old_id2name, new_id2name;
    LOG(DEBUG) << "old_td.name_to_node.size() = " << old_td.name_to_node.size() ;
    LOG(DEBUG) << "new_td.name_to_node.size() = " << new_td.name_to_node.size() ;
    fill_name_id_maps(old_td, old_id2name);
    fill_name_id_maps(new_td, new_id2name);
    LOG(DEBUG) << old_id2name.size() << " " << new_id2name.size() << std::endl;
    const OttIdSet same_id_name = find_ids_with_same_names(old_id2name, new_id2name);
    LOG(DEBUG) << same_id_name.size() << " IDs with the same name between versions";
    const nd2idset_t old_nd2ids = fill_des_id_set(old_tree, same_id_name);
    const nd2idset_t new_nd2ids = fill_des_id_set(new_tree, same_id_name);
    return true;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        std::ostream & out = std::cout;
        if (!args.count("oldtaxonomy")) {
            cerr << "oldtaxonomy expected as first unnamed argument\n";
            return 1;
        }
        if (!args.count("newtaxonomy")) {
            cerr << "newtaxonomy expected as second unnamed argument\n";
            return 1;
        }
        string otd = args["oldtaxonomy"].as<string>();
        string ntd = args["newtaxonomy"].as<string>();
        OttId keep_root = -1;
        bitset<32> cleaning_flags = 0;
        LOG(INFO) << "loading old taxonomy\n";
        Taxonomy::tolerate_synonyms_to_unknown_id = true;
        TaxonomyDiffMaker otaxonomy = {otd, cleaning_flags, keep_root};
        LOG(INFO) << "loading new taxonomy\n";
        TaxonomyDiffMaker ntaxonomy = {ntd, cleaning_flags, keep_root};
        diff_from_taxonomies(out, otaxonomy, ntaxonomy);
        
    } catch (std::exception& e) {
        cerr << "otc-taxonomy-diff-maker: Error! " << e.what() << std::endl;
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

