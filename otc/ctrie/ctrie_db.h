#ifndef OTC_CTRIE_DB_H
#define OTC_CTRIE_DB_H


#include "otc/ctrie/ctrie.h"

namespace otc {
using CTrie128_t = CompressedTrie<CTrie128Node>;
using CTrie80_t = CompressedTrie<CTrie80Node>;

class CompressedTrieBasedDB {
    public:
    void initialize(const std::set<std::string> & keys);
    std::set<FuzzyQueryResult, SortQueryResByNearness>  fuzzy_query(const std::string & query_str) const;
    private:
    // CTrie3_t wide_trie;
    CTrie128_t thin_trie;
};


inline std::set<FuzzyQueryResult, SortQueryResByNearness> CompressedTrieBasedDB::fuzzy_query(const std::string & query_str) const {
#   if defined(U32_TRIE_QUERIES)
        auto conv_query = to_query_str(query_str);
#   else
        const std::string & conv_query = lower_case_version(query_str);
#   endif
    unsigned int max_dist;
    // defaults taken from taxomachine...
    const unsigned int SHORT_NAME_LENGTH = 9;
    const unsigned int MEDIUM_NAME_LENGTH = 14;
    const unsigned int LONG_NAME_LENGTH = 19;
    std::size_t iql = conv_query.length();
    if (iql < SHORT_NAME_LENGTH) {
        max_dist = 1;
    } else if (iql < MEDIUM_NAME_LENGTH) {
        max_dist = 2;
    } else {
        max_dist = (iql < LONG_NAME_LENGTH ? 3 : 4);
    }
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
    
    auto from_thin = thin_trie.fuzzy_matches(conv_query, max_dist);
    sorted.insert(std::begin(from_thin), std::end(from_thin));

    //auto from_full = wide_trie.fuzzy_matches(conv_query, max_dist);
    //sorted.insert(std::begin(from_full), std::end(from_full));
    return sorted;
}


inline void CompressedTrieBasedDB::initialize(const std::set<std::string> & keys) {
    ctrie_init_set_t for_wide;
    ctrie_init_set_t for_thin;
    // could fit a couple more non-funky, if we want <- 76, I think...
    bool trimming_by_funky = false;
    auto nonfunky = " \'()-.0123456789:,_aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ/?";
    std::ostream & out = std::cout;
    std::map<stored_char_t, unsigned int> letter_counts;
    std::set<stored_char_t> thin_letter_set;
    unsigned mem_str = 0;
    out << keys.size() << " keys\n";
    for (auto i : keys) {
        mem_str += i.length();
        auto widestr = to_stored_str_type(i);
        bool has_funky = false;
        if (trimming_by_funky) {
            for (auto c : i) {
                if (std::strchr(nonfunky, c) == nullptr) {
                    has_funky = true;
                    break;
                }
            }
        }
        if (has_funky) {
            for_wide.insert(widestr);
            for (auto letter : widestr) {
                letter_counts[letter] += 1;
            }
        } else {
            for_thin.insert(widestr);
            for (auto letter : widestr) {
               thin_letter_set.insert(letter);
            }
        }
        if (contains(thin_letter_set, '\0')) {
            LOG(DEBUG) << "\\0 in \"" << i << "\"";
        }
        //std::cerr << glob_conv8.to_bytes(widestr) << '\n';
    }
    stored_str_t wide_letters;
    stored_str_t thin_letters;
    thin_letters.insert(std::begin(thin_letters), std::begin(thin_letter_set), std::end(thin_letter_set));
    std::map<unsigned int, stored_str_t> by_count;
    for (auto lcp : letter_counts) {
        wide_letters.push_back(lcp.first);
        by_count[lcp.second].push_back(lcp.first);
    }
    //std::cerr << "set size = " << (sizeof(std::string *) + sizeof(char *) + 8)*keys.size() + mem_str << "bytes\n";
    //wide_trie.init(for_wide, wide_letters);
    thin_trie.init(for_thin, thin_letters);

    std::cerr << "writing input words:\n";
    thin_trie.db_write_words(std::cerr);

}


} // namespace otc
#endif
