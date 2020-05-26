#ifndef OTC_CTRIE_DB_H
#define OTC_CTRIE_DB_H


#include "otc/ctrie/ctrie.h"

namespace otc {

class CompressedTrieBasedDB {
    public:
    void initialize(const std::set<std::string> & keys);
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    std::set<FuzzyQueryResult, SortQueryResByNearness>  exact_query(const std::string & query_str) const;
    std::vector<std::string>                            prefix_query(const std::string & query_str) const;
    private:
    CompressedTrie wide_trie;
    CompressedTrie thin_trie;
};


} // namespace otc
#endif
