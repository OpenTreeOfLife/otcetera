#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>

#include <regex>
#include <tuple>
#include <string>
#include <iomanip>
#include <algorithm>
#include <queue>
#include <map>
#include <stack>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"
#include "otc/ctrie/str_utils.h"
#include "otc/ctrie/str_utils.cpp"
#include "otc/ctrie/ctrie_db.h"


INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;


variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    options_description invisible("Invisible options");
    invisible.add_options()("taxonomy", value<string>(), "Filename for the taxonomy");
    options_description visible;
    visible.add(otc::standard_options());
    positional_options_description p;
    p.add("taxonomy", -1);
    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-tnrs-cli <taxonomy-dir> [OPTIONS]\n"
                                                    "Build data structures for name matching, and allow interactive testing.",
                                                    visible, invisible, p);

    return vm;
}

OttIdSet diagnose_name(const RTRichTaxTreeData & rt_data,
                   std::ostream * out,
                   const std::string &name) {
    OttIdSet ret;
    if (out != nullptr) {*out << "Diagnosing name \"" << name << "\":\n";}
    auto n2nit = rt_data.name_to_node.find(name);
    if (n2nit != rt_data.name_to_node.end()) {
        const RTRichTaxNode * nd = n2nit->second;
        if (nd != nullptr) {
            ret.insert(nd->get_ott_id());
            if (out != nullptr) {*out << "node in taxonomy: ott_id = " << nd->get_ott_id() << " name = \"" << nd->get_name() << "\"\n";}
        }
    }
    auto n2hit = rt_data.homonym_to_node.find(name);
    if (n2hit != rt_data.homonym_to_node.end()) {
        const auto & ndv = n2hit->second;
        for (auto nd : ndv) {
            if (nd != nullptr) {
                ret.insert(nd->get_ott_id());
                if (out != nullptr) {*out << "homnym to node in taxonomy: ott_id = " << nd->get_ott_id() << " name = \"" << nd->get_name() << "\"\n";}
            }
        }
    }
    auto n2snit = rt_data.name_to_record.find(name);
    if (n2snit != rt_data.name_to_record.end()) {
        const TaxonomyRecord * rec = n2snit->second;
        if (rec != nullptr) {
            ret.insert(rec->id);
            if (out != nullptr) {*out << "suppressed record: ott_id = " << rec->id << " name = \"" << rec->name << "\"\n";}
        }
    }
    auto n2shit = rt_data.homonym_to_record.find(name);
    if (n2shit != rt_data.homonym_to_record.end()) {
        const auto & ndv = n2shit->second;
        for (auto rec : ndv) {
            if (rec != nullptr) {
                ret.insert(rec->id);
                if (out != nullptr) {*out << "suppressed homnym to record: ott_id = " << rec->id << " name = \"" << rec->name << "\"\n";}
            }
        }
    }
    return ret;
}

void analyze_case_sensitivity(const RTRichTaxTreeData & rt_data,
                              const std::set<std::string_view> & all_names) {
    std::map<std::string, std::string> lc_names_set;
    int inc=1;
    for (auto n : all_names) {
        std::string ncn{std::begin(n), std::end(n)};
        std::string lcv = lower_case_version(ncn);
        if (contains(lc_names_set, lcv)) {
            const auto & pn = lc_names_set[lcv];
            std::cerr << "LC name clash \"" << ncn << "\" and \"" << pn << "\"\n";
            auto prev_id_set = diagnose_name(rt_data, &std::cerr, ncn);
            auto curr_id_set = diagnose_name(rt_data, &std::cerr, pn);
            if (prev_id_set != curr_id_set) {
                std::cerr << inc << "  \"" << ncn << "\" maps to ";
                write_ott_id_set(std::cerr, " ", prev_id_set, " ");
                std::cerr << '\n' <<  inc++ << "  \"" << pn << "\" maps to ";
                write_ott_id_set(std::cerr, " ", curr_id_set, " ");
                std::cerr << '\n';
            }
        } else {
            lc_names_set[lcv] = ncn;
        }
    }
    std::cout << all_names.size() << " names when honoring case.\n";
    std::cout << lc_names_set.size() << " names when not honoring case.\n";
    return;
}

void process_taxonomy(const RichTaxonomy & taxonomy) {
    const auto & rich_tax_tree = taxonomy.get_tax_tree();
    const auto & rt_data = rich_tax_tree.get_data();
    std::set<std::string_view> all_names;
    auto insert_hint = all_names.begin();
    for (auto const & name2nd : rt_data.name_to_node) {
        insert_hint = all_names.insert(insert_hint, name2nd.first);
    }
    insert_hint = all_names.begin();
    for (auto name2ndvec : rt_data.homonym_to_node) {
        insert_hint = all_names.insert(insert_hint, name2ndvec.first);
    }
    // filtered
    insert_hint = all_names.begin();
    for (auto name2rec : rt_data.name_to_record) {
        insert_hint = all_names.insert(insert_hint, name2rec.first);
    }
    insert_hint = all_names.begin();
    for (auto name2recvec : rt_data.homonym_to_record) {
        insert_hint = all_names.insert(insert_hint, name2recvec.first);
    }
    for (const auto & tjs : taxonomy.get_synonyms_list()) {
        all_names.insert(std::string_view{tjs.name});
    }
    


    CompressedTrieBasedDB ct{all_names};
    
    std::cout << "Enter a query and hit return:\n";
    std::string query;
    while (std::getline(std::cin, query)) {
        std::cout << "query =\"" << query << "\"\n";
        auto results = ct.fuzzy_query(query);
        std::cout << results.size() << " matches:\n";
        for (auto res : results) {
            std::cout << "res.match = \"" << res.match()
                      << "\", distance = " << res.distance 
                      << " score = " << res.score << "\n";
        }
        std::cout << "Enter a query and hit return:\n";
    }
    std::cerr << "EOF\n";
}


int main(int argc, char* argv[]) {
    if (set_global_conv_facet() != 0) {
        return 1;
    }
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        auto taxonomy = load_rich_taxonomy(args);
        process_taxonomy(taxonomy);
    } catch (std::exception& e) {
        cerr << "otc-tnrs-cli: Error! " << e.what() << std::endl;
        return 1;
    }
}
