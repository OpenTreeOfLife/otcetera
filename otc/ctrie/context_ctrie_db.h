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

    // What strings (for names or synonyms) match the normalized query string?
    std::set<FuzzyQueryResult, SortQueryResByNearness> fuzzy_query(const std::string & query_str) const;

    // Does anything match this normalized query string?
    std::optional<std::string> exact_query(const std::string & query_str) const;

    // Does anything match this normalized query string?
    std::vector<std::string> prefix_query(const std::string & query_str) const;

    std::vector<FuzzyQueryResultWithTaxon> fuzzy_query_to_taxa(const std::string & query_str,
                                                               const RTRichTaxNode * context_root,
                                                               const RichTaxonomy & taxonomy,
                                                               bool include_suppressed) const;

    std::vector<FuzzyQueryResultWithTaxon> to_taxa(const std::set<FuzzyQueryResult, SortQueryResByNearness>& sorted_results,
                                                   const RTRichTaxNode * context_root,
                                                   const RichTaxonomy & taxonomy,
                                                   bool include_suppressed) const;

    std::vector<TaxonResult> to_taxa(const std::optional<std::string>& result,
                                     const RTRichTaxNode * context_root,
                                     const RichTaxonomy & taxonomy,
                                     bool include_suppressed) const;

    std::vector<TaxonResult> to_taxa(const std::vector<std::string>& result,
                                     const RTRichTaxNode * context_root,
                                     const RichTaxonomy & taxonomy,
                                     bool include_suppressed) const;

    void add_key(const std::string& s, OttId id, const RichTaxonomy&);

private:
    const Context & context;
    CompressedTrieBasedDB trie;
    std::map<std::string,  vec_taxon_and_syn_ptrs> match_name_to_taxon;

};



} // namespace otc
#endif
