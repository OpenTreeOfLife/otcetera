#ifndef OTC_CTRIE_DB_H
#define OTC_CTRIE_DB_H


#include "otc/ctrie/ctrie.h"

namespace otc {
using CTrie3_t = CompressedTrie<CTrie3Node>;
using CTrie2_t = CompressedTrie<CTrie2Node>;


class CompressedTrieBasedDB {
    public:
    CompressedTrieBasedDB(const std::set<std::string_view> & keys);
    std::list<FuzzyQueryResult> fuzzy_query(const std::string & query_str);
    private:
    CTrie3_t fat_trie;
    CTrie2_t thin_trie;
};



std::list<FuzzyQueryResult> CompressedTrieBasedDB::fuzzy_query(const std::string & query_str) {
    auto conv_query = to_u32string(query_str);
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
    auto from_thin = thin_trie.fuzzy_matches(conv_query, max_dist);
    auto from_full = fat_trie.fuzzy_matches(conv_query, max_dist);
    from_thin.insert(std::end(from_thin), std::begin(from_full), std::end(from_full));
    return from_thin;
}




CompressedTrieBasedDB::CompressedTrieBasedDB(const std::set<std::string_view> & keys) {
    ctrie_init_set_t for_fat;
    ctrie_init_set_t for_thin;
    // could fit a couple more non-funky, if we want <- 76, I think...
    auto nonfunky = " \'()-.0123456789:,_aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ/?";
    std::ostream & out = std::cout;
    std::map<stored_char_t, unsigned int> letter_counts;
    std::set<stored_char_t> thin_letter_set;
    unsigned mem_str = 0;
    out << keys.size() << " keys\n";
    for (auto i : keys) {
        mem_str += i.length();
        auto widestr = to_u32string(i);
        bool has_funky = false;
        for (auto c : i) {
            if (std::strchr(nonfunky, c) == nullptr) {
                has_funky = true;
                break;
            }
        }
        if (has_funky) {
            for_fat.insert(widestr);
            for (auto letter : widestr) {
                letter_counts[letter] += 1;
            }
        } else {
            for_thin.insert(widestr);
            for (auto letter : widestr) {
               thin_letter_set.insert(letter);
            }
        }
        //std::cerr << glob_conv8.to_bytes(widestr) << '\n';
    }
    stored_str_t fat_letters;
    stored_str_t thin_letters;
    thin_letters.insert(std::begin(thin_letters), std::begin(thin_letter_set), std::end(thin_letter_set));
    std::map<unsigned int, stored_str_t> by_count;
    for (auto lcp : letter_counts) {
        fat_letters.push_back(lcp.first);
        by_count[lcp.second].push_back(lcp.first);
    }
    /* 
    int i = 0;
    for (auto bcit = by_count.rbegin(); bcit != by_count.rend(); ++bcit) {
        for (auto curr_let : bcit->second) {
            out << i++ << " \"" << to_char_str(curr_let) <<  "\" " << bcit->first << '\n';
        }
    }
    */
    std::cerr << "set size = " << (sizeof(std::string *) + sizeof(char *) + 8)*keys.size() + mem_str << "bytes\n";
    fat_trie.init(for_fat, fat_letters);
    thin_trie.init(for_thin, thin_letters);
}



} // namespace otc
#endif
