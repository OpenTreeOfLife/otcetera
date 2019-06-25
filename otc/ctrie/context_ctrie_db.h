#ifndef OTC_CONTEXT_CTRIE_DB_H
#define OTC_CONTEXT_CTRIE_DB_H

#include <vector>
#include <set>
#include "otc/ctrie/ctrie_db.h"

namespace otc {

class RichTaxonomy;
struct Context;

class ContextAwareCTrieBasedDB {
    public:
    ContextAwareCTrieBasedDB(const Context &, const RichTaxonomy &);
    ContextAwareCTrieBasedDB(const Context &, const RichTaxonomy &, const std::set<std::string_view> & keys);
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    private:
    const Context & context;
    std::vector<const ContextAwareCTrieBasedDB *> children;
};



} // namespace otc
#endif
