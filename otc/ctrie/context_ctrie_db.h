#ifndef OTC_CONTEXT_CTRIE_DB_H
#define OTC_CONTEXT_CTRIE_DB_H

#include <vector>
#include <set>
#include "otc/ctrie/ctrie_db.h"
namespace otc {

class RichTaxonomy;
struct Context;

using const_rich_taxon_and_syn_ptr = std::pair<const RTRichTaxNode *, const void *>;
using vec_taxon_and_syn_ptrs = std::vector<const_rich_taxon_and_syn_ptr>;
class ContextAwareCTrieBasedDB {
    public:
    ContextAwareCTrieBasedDB(const Context &, const RichTaxonomy &);
    ContextAwareCTrieBasedDB(const Context &, const RichTaxonomy &, const std::set<std::string_view> & keys);
    sorted_q_res_set  fuzzy_query(const std::string & query_str) const;
    vec_q_res_w_taxon fuzzy_query_to_taxa(const std::string & query_str,
                                          const RTRichTaxNode * context_root,
                                          const RichTaxonomy & taxonomy,
                                          bool include_suppressed) const {
        const auto sorted = fuzzy_query(query_str);
        return tie_to_taxa(sorted, query_str, context_root, taxonomy, include_suppressed, nullptr);
    }
    
    std::optional<FuzzyQueryResult>  exact_query(const std::string & rqw_query, const std::string & norm_query) const;
    
    vec_q_res_w_taxon exact_query_to_taxa(const std::string & raw_query, 
                                          const std::string & norm_query,
                                          const RTRichTaxNode * context_root,
                                          const RichTaxonomy & taxonomy,
                                          bool include_suppressed) const {
        if (include_suppressed) {
            std::function<bool(const RTRichTaxNode*)> ok = [](const RTRichTaxNode*){return true;};
            return exact_query_to_taxa(raw_query, norm_query, context_root, taxonomy, ok);
        } else {
            std::function<bool(const RTRichTaxNode*)> ok = [&](const RTRichTaxNode* taxon) {
                return (taxon != nullptr) && (not taxonomy.node_is_suppressed_from_tnrs(taxon));
            };
            return exact_query_to_taxa(raw_query, norm_query, context_root, taxonomy, ok);
        }
    }

    vec_q_res_w_taxon exact_query_to_taxa(const std::string & raw_query, 
                                          const std::string & norm_query,
                                          const RTRichTaxNode * context_root,
                                          const RichTaxonomy & taxonomy,
                                          std::function<bool(const RTRichTaxNode*)> keep) const {
        const auto res = exact_query(raw_query, norm_query);
        sorted_q_res_set sorted;
        if (res) {
            sorted.insert(*res);
        }
        return tie_to_taxa(sorted, raw_query, context_root, taxonomy, keep, &raw_query);
    }
        
    private:

    vec_q_res_w_taxon tie_to_taxa(const sorted_q_res_set &sorted,
                                  const std::string & query_str,
                                  const RTRichTaxNode * context_root,
                                  const RichTaxonomy & taxonomy, 
                                  bool include_suppressed,
                                  const std::string * exact_string) const {
        if (include_suppressed) {
            std::function<bool(const RTRichTaxNode*)> ok =  [](const RTRichTaxNode*){return true;};
            return tie_to_taxa(sorted, query_str, context_root, taxonomy, ok, exact_string);
        } else {
            std::function<bool(const RTRichTaxNode*)> ok = [&](const RTRichTaxNode* taxon) {
                return (taxon != nullptr) && (not taxonomy.node_is_suppressed_from_tnrs(taxon));
            };
            return tie_to_taxa(sorted, query_str, context_root, taxonomy, ok, exact_string);
        }        
    }
    vec_q_res_w_taxon tie_to_taxa(const sorted_q_res_set &sorted,
                                  const std::string & query_str,
                                  const RTRichTaxNode * context_root,
                                  const RichTaxonomy & , 
                                  keep_taxon_pred_t keep,
                                  const std::string * exact_string) const;
    const Context & context;
    std::vector<const ContextAwareCTrieBasedDB *> children;
    CompressedTrieBasedDB trie;
    std::uint32_t filter_trav_enter = 0;
    std::uint32_t trav_exit = UINT32_MAX;
    std::map<std::string,  vec_taxon_and_syn_ptrs> match_name_to_taxon;

};



} // namespace otc
#endif
