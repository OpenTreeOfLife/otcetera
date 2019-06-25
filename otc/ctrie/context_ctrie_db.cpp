#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

namespace otc {

ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                         const RichTaxonomy &)
    :context(context_arg) {
}


ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                         const RichTaxonomy &,
                         const std::set<std::string_view> &)
    :context(context_arg) {
}

std::set<FuzzyQueryResult, SortQueryResByNearness> ContextAwareCTrieBasedDB::fuzzy_query(const std::string & query_str) const {
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
    if (context.name_matcher != nullptr) {
        sorted = context.name_matcher->fuzzy_query(query_str);
    }
    for (auto c :children) {
        if (c->context.name_matcher) {
            auto csorted = c->context.name_matcher->fuzzy_query(query_str);
            sorted.insert(std::begin(csorted), std::end(csorted));
        }
    }
    return sorted;
}


} // namespace otc
