#include "otc/ctrie/ctrie_db.h"

using std::vector;
using std::string;

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

    auto from_new = new_trie->fuzzy_matches(conv_query, max_dist);
    sorted.insert(std::begin(from_new), std::end(from_new));

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

    auto from_new = new_trie->fuzzy_matches(conv_query, 0);
    sorted.insert(std::begin(from_new), std::end(from_new));

    return sorted;
}

vector<string> CompressedTrieBasedDB::prefix_query(const std::string & query_str) const
{
    auto conv_query = to_u32string(query_str);

    auto sorted = thin_trie.prefix_query(conv_query);

    auto from_full = wide_trie.prefix_query(conv_query);
    sorted.insert(sorted.end(), std::begin(from_full), std::end(from_full));

    auto from_new = new_trie->prefix_query(conv_query);
    sorted.insert(sorted.end(), std::begin(from_new), std::end(from_new));

    // I'm not sure this is a good idea...
    std::sort(sorted.begin(), sorted.end());

    return sorted;
}

void CompressedTrieBasedDB::add_key(const std::string& s)
{
    new_keys.insert(s);
    rebuild_new_trie();
}

void CompressedTrieBasedDB::rebuild_new_trie()
{
    LOG(INFO)<<"Rebuilding ctree for additions";
    LOG(INFO)<< new_keys.size() << " keys\n";

    ctrie_init_set_t for_new;
    std::map<stored_char_t, unsigned int> letter_counts;
    std::set<stored_char_t> new_letter_set;
    // unsigned mem_str = 0;
    for (auto i : new_keys)
    {
        // mem_str += i.length();
        auto widestr = to_u32string(i);

	for_new.insert(widestr);

	for (auto letter : widestr) {
	    new_letter_set.insert(letter);
	}
        //LOG(INFO) << glob_conv8.to_bytes(widestr) << '\n';
    }
    LOG(INFO)<<for_new.size()<<" keys for new ctrie\n";

    stored_str_t new_letters;
    new_letters.insert(std::begin(new_letters), std::begin(new_letter_set), std::end(new_letter_set));
    std::map<unsigned int, stored_str_t,std::greater<int>> by_count;

    LOG(INFO)<<"new letters: "<<new_letters.size()<<"\n";
    for (auto letter : new_letters)
        LOG(INFO)<<"  "<<to_char_str(letter)<<"\n";
    LOG(INFO)<<"done\n";

    for (auto [letter,count] : letter_counts)
        by_count[count].push_back(letter);

    int n_dropped_letters = 0;
    for(auto [count,letters]: by_count)
    {
        LOG(INFO)<<"  "<<count<<" : \n";
        for(auto letter : letters)
        {
            LOG(INFO)<<"    "<<to_char_str(letter)<<" ("<<(uint32_t)(letter)<<")";
            if (new_letters.size() < 64)
                new_letters.push_back(letter);
            else
            {
                n_dropped_letters++;
                LOG(INFO)<<"  DROPPED!";
            }
            LOG(INFO)<<"\n";
        }
    }
    LOG(INFO)<<"done\n";

    if (n_dropped_letters)
    {
        LOG(INFO)<<"dropped "<<n_dropped_letters<<" letters.\n";
    }

    //LOG(INFO) << "set size = " << (sizeof(std::string *) + sizeof(char *) + 8)*keys.size() + mem_str << "bytes\n";

    new_trie = std::make_shared<CompressedTrie>();
    new_trie->init(for_new, new_letters);
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
            std::cerr<<"    "<<to_char_str(letter)<<" ("<<(uint32_t)(letter)<<")";
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

    new_trie = std::make_shared<CompressedTrie>();
    new_trie->init({},thin_letters);

    /*
    wide_trie.db_write_words(std::cerr);
    thin_trie.db_write_words(std::cerr);
    */
}

}
