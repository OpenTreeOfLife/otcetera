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
#include <chrono>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"
#include "otc/ctrie/str_utils.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/ctrie/context_ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/ws/tolws.h"

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
    auto n2hit = rt_data.homonym_to_nodes.find(name);
    if (n2hit != rt_data.homonym_to_nodes.end()) {
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

void interactive_tests() {
    const CTrie2_t testtrie{"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"};
    const std::string p1 = "Enter a query:\n";
    const std::string p2 = "Enter a trie:\n";
    const std::string p3 = "max distance:\n";
    auto & out = std::cerr;
    std::string trash;
    for (;;) {
        std::string query, trie;
        out << p1;
        if (!std::getline(std::cin, query)) {
            break;
        }
        out << p2;
        if (!std::getline(std::cin, trie)) {
            break;
        }
        out << p3;
        unsigned int dist_threshold;
        std::cin >> dist_threshold;
        std::getline(std::cin, trash);
        std::ofstream lastinteractive("lastinteractivetests.txt");
        lastinteractive << query << '\n';
        lastinteractive << trie << '\n';
        lastinteractive << dist_threshold << '\n';
        lastinteractive.close();
        auto wq = to_u32string(query);
        auto wqi = testtrie.using_letter_order_to_encode(wq);
        auto wt = to_u32string(trie);
        auto wti = testtrie.using_letter_order_to_encode(wt);
        out << "query \"" << query << "\"\n";
        out << "       "; for(auto q : wqi) {out << (unsigned int) q << ", ";}; out << "\n";
        out << "trie  \"" << trie << "\"\n";
        out << "       "; for (auto q : wti) {out << (unsigned int) q << ", ";}; out << "\n";
        out << "max_dist = " << dist_threshold << '\n';

        auto resdist = testtrie._calc_dist_prim_impl(NO_MATCHING_CHAR_CODE,
                                                     &(wqi[0]),
                                                     wqi.size(),
                                                     &(wti[0]),
                                                     wti.size(),
                                                     dist_threshold,
                                                     NO_MATCHING_CHAR_CODE);
        out << "result dist = " << resdist << "\n";
    }
    out << "EOF\n";
}

void process_taxonomy(const RichTaxonomy & taxonomy) {
    const Context * c = determine_context({});
    if (c == nullptr) {
        throw OTCError() << "no context found for entire taxonomy";
    }
    ContextAwareCTrieBasedDB ct{*c, taxonomy};

    using time_diff_t = std::chrono::duration<double, std::milli>;
    time_diff_t total_time;
    auto num_q = 0;
    
    std::cout << "Enter a query and hit return:\n";
    std::string query;
    const std::string context_name = "All life";
    const std::vector<std::string> empty_ids_vec;
    while (std::getline(std::cin, query)) {
        std::cout << "query =\"" << query << "\"\n";
        std::vector<std::string> qvec;
        qvec.push_back(query);
        auto start_time = std::chrono::steady_clock::now();
        auto tnrs_res = tnrs_match_names_ws_method(qvec,
                                                   context_name,
                                                   true,
                                                   empty_ids_vec,
                                                   false,
                                                   taxonomy);
        auto end_time = std::chrono::steady_clock::now();
        num_q++;
        time_diff_t time_diff = end_time - start_time;
        std::cerr << "took: " << time_diff.count() << " milliseconds.\n";
        total_time += time_diff;
        std::cout << tnrs_res;
        std::cout << "\nEnter a query and hit return:\n";
    }
    std::cerr << "EOF\n";
    std::cerr << num_q << " queries " << total_time.count() << " milliseconds\n";
}


int main(int argc, char* argv[]) {
    if (set_global_conv_facet() != 0) {
        return 1;
    }
    std::ios::sync_with_stdio(false);
    if (argc == 1) {
        interactive_tests();
        return 0;
    }
    try {
        auto args = parse_cmd_line(argc, argv);
        auto taxonomy = load_rich_taxonomy(args);
        process_taxonomy(taxonomy);
    } catch (std::exception& e) {
        cerr << "otc-tnrs-cli: Error! " << e.what() << std::endl;
        return 1;
    }
}
