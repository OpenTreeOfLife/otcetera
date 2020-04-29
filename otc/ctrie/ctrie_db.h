#ifndef OTC_CTRIE_DB_H
#define OTC_CTRIE_DB_H


#include "otc/ctrie/ctrie.h"

namespace otc {
using CTrie3_t = CompressedTrie<CTrie3Node>;
using CTrie2_t = CompressedTrie<CTrie2Node>;

class CompressedTrieBasedDB {
    public:
    void initialize(const std::set<std::string> & keys);
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    std::set<FuzzyQueryResult, SortQueryResByNearness>  exact_query(const std::string & query_str) const;
    private:
    CTrie2_t wide_trie;
    CTrie2_t thin_trie;
};


} // namespace otc
#endif
