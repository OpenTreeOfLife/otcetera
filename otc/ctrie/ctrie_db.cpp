#include "otc/ctrie/ctrie_db.h"

namespace otc {

std::set<FuzzyQueryResult, SortQueryResByNearness> CompressedTrieBasedDB::fuzzy_query(const std::string & query_str) const {
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
    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;
    
    auto from_thin = thin_trie.fuzzy_matches(conv_query, max_dist);
    sorted.insert(std::begin(from_thin), std::end(from_thin));

    auto from_full = wide_trie.fuzzy_matches(conv_query, max_dist);
    sorted.insert(std::begin(from_full), std::end(from_full));
    return sorted;
}

std::set<FuzzyQueryResult, SortQueryResByNearness> CompressedTrieBasedDB::exact_query(const std::string & query_str) const
{
    auto conv_query = to_u32string(query_str);

    std::set<FuzzyQueryResult, SortQueryResByNearness> sorted;

    auto from_thin = thin_trie.fuzzy_matches(conv_query, 0);
    sorted.insert(std::begin(from_thin), std::end(from_thin));

    auto from_full = wide_trie.fuzzy_matches(conv_query, 0);
    sorted.insert(std::begin(from_full), std::end(from_full));

    return sorted;
}

void CompressedTrieBasedDB::initialize(const std::set<std::string> & keys) {
    ctrie_init_set_t for_wide;
    ctrie_init_set_t for_thin;

    // We don't need capital letters here, since we lcase queries when we normalize them.
    auto nonfunky = " \"\'()[]+-%.&0123456789:<=>,^_abcdefghijklmnopqrstuvwxyz/?#*!";

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
        //std::cerr << glob_conv8.to_bytes(widestr) << '\n';
    }
    std::cerr<<for_thin.size()<<" keys for thin ctrie\n";
    std::cerr<<for_wide.size()<<" keys for wide ctrie\n";
    stored_str_t wide_letters;
    stored_str_t thin_letters;
    thin_letters.insert(std::begin(thin_letters), std::begin(thin_letter_set), std::end(thin_letter_set));
    std::map<unsigned int, stored_str_t,std::greater<int>> by_count;


    std::cerr<<"thin letters: "<<thin_letters.size()<<"\n";
    for (auto letter : thin_letters)
        std::cerr<<"  "<<to_char_str(letter)<<"\n";
    std::cerr<<"done\n";

    for (auto [letter,count] : letter_counts)
        by_count[count].push_back(letter);

    std::cerr<<"wide letters: "<<letter_counts.size()<<" letters\n";
    int n_dropped_letters = 0;
    for(auto [count,letters]: by_count)
    {
        std::cerr<<"  "<<count<<" : \n";
        for(auto letter : letters)
        {
            std::cerr<<"    "<<to_char_str(letter)<<" ("<<letter<<")";
            if (wide_letters.size() < 64)
                wide_letters.push_back(letter);
            else
            {
                n_dropped_letters++;
                std::cerr<<"  DROPPED!";
            }
            std::cerr<<"\n";
        }
    }
    std::cerr<<"done\n";

    if (n_dropped_letters)
        std::cerr<<"dropped "<<n_dropped_letters<<" letters.\n";

    //std::cerr << "set size = " << (sizeof(std::string *) + sizeof(char *) + 8)*keys.size() + mem_str << "bytes\n";
    wide_trie.init(for_wide, wide_letters);
    thin_trie.init(for_thin, thin_letters);

    /*
    wide_trie.db_write_words(std::cerr);
    thin_trie.db_write_words(std::cerr);
    */
}

}
