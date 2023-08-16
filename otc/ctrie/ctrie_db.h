#ifndef OTC_CTRIE_DB_H
#define OTC_CTRIE_DB_H


#include "otc/ctrie/ctrie.h"
#include <memory>

namespace otc {

class CompressedTrieBasedDB {
public:
    void initialize(const std::set<std::string> & keys);
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    std::set<FuzzyQueryResult, SortQueryResByNearness>  exact_query(const std::string & query_str) const;
    std::vector<std::string>                            prefix_query(const std::string & query_str) const;

    void add_key(const std::string& s);

    void rebuild_new_trie();

private:
    CompressedTrie wide_trie;
    CompressedTrie thin_trie;

    std::shared_ptr<CompressedTrie> new_trie;
    std::set<std::string> new_keys;
};


} // namespace otc
#endif
