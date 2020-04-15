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
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    std::set<FuzzyQueryResult, SortQueryResByNearness>  exact_query(const std::string & query_str) const;

    std::vector<FuzzyQueryResultWithTaxon> fuzzy_query_to_taxa(const std::string & query_str,
                                                               const RTRichTaxNode * context_root,
                                                               const RichTaxonomy & taxonomy,
                                                               bool include_suppressed) const;
    private:
    const Context & context;
    std::vector<const ContextAwareCTrieBasedDB *> children;
    CompressedTrieBasedDB trie;
    std::uint32_t filter_trav_enter = 0;
    std::uint32_t trav_exit = UINT32_MAX;
    std::map<std::string,  vec_taxon_and_syn_ptrs> match_name_to_taxon;

};



} // namespace otc
#endif
