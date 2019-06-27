#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ctrie/ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

namespace otc {

ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                                                   const RichTaxonomy &taxonomy)
    :context(context_arg) {
    if (context_arg.name_matcher != nullptr) {
        return; // already initialized
    }
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
    
    trie.initialize(all_names);
    context_arg.name_matcher = &trie;
}


ContextAwareCTrieBasedDB::ContextAwareCTrieBasedDB(const Context &context_arg,
                         const RichTaxonomy &,
                         const std::set<std::string_view> &)
    :context(context_arg) {
    throw OTCError() << "partitioning by context is not implemented yet...";
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
