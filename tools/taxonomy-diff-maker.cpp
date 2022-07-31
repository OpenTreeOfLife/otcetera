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

template<typename T>
bool fill_name_id_maps(const T & tree_data, id2name_t& id2name, name2id_t & name2id) {
    LOG(INFO) << tree_data.id_to_record.size() << " tree_data.id_to_record";
    for (auto id_rec_pair : tree_data.id_to_record) {
        const auto & record = id_rec_pair.second;
        const auto & tax_id = id_rec_pair.first;
        const auto name = record->name;
        if (contains(name2id, name)) {
            LOG(ERROR) << "name \"" << name << "\" repeated";
        } else {
            if (contains(id2name, tax_id)) {
                LOG(ERROR) << "name \"" << tax_id << "\" repeated";
                return false;
            }
            id2name[tax_id] = name;
            name2id[name] = tax_id;
        }
    }
    return true;
}

bool diff_from_taxonomies(std::ostream & out,
                          const TaxonomyDiffMaker & old_tax,
                          const TaxonomyDiffMaker & new_tax) {
    const auto & old_tree = old_tax.get_tax_tree();
    const auto & new_tree = new_tax.get_tax_tree();
    const auto & old_td = old_tree.get_data();
    const auto & new_td = new_tree.get_data();
    id2name_t old_id2name, new_id2name;
    name2id_t old_name2id, new_name2id;
    LOG(DEBUG) << "old_td.name_to_node.size() = " << old_td.name_to_node.size() ;
    LOG(DEBUG) << "new_td.name_to_node.size() = " << new_td.name_to_node.size() ;
    
    if (!fill_name_id_maps(old_td, old_id2name, old_name2id)) {
        return false;
    }
    if (!fill_name_id_maps(new_td, new_id2name, new_name2id)) {
        return false;
    }
    out << old_id2name.size() << " " << new_id2name.size() << std::endl;
    LOG(ERROR) << "hi";
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

