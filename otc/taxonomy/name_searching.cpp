// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <bitset>
#include <fstream>
#include <regex>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
namespace fs = boost::filesystem;

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/ctrie/str_utils.h"
namespace otc
{
vec_tax_nodes_t RichTaxonomy::exact_name_search(const std::string & query_ref,
                                                const std::string * normalized_query,
                                                const RTRichTaxNode* context_root,
                                                std::function<bool(const RTRichTaxNode*)> ok) const {
    auto fuzzy = this->get_fuzzy_matcher();
    if (context_root == nullptr) {
        context_root = get_tax_tree().get_root();
    }
    vec_tax_nodes_t hits;
    if (fuzzy == nullptr) {
        std::string query{query_ref};
        for (auto& c: query) {
            c = std::tolower(c);
        }
        for(auto taxon: iter_post_n_const(*context_root)) {
            if (not ok(taxon)) {
                continue;
            }
            if (lcase_string_equals(query, taxon->get_data().get_nonuniqname())) {
                hits.push_back(taxon);
            }
        }
    } else {
        std::string scratch;
        if (normalized_query == nullptr) {
            scratch = normalize_query(query_ref);
            normalized_query = &scratch;
        }
        
    }
    return hits;
}



} //namespace otc
