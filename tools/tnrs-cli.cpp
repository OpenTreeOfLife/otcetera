#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <bitset>
#include <regex>
#include <tuple>
#include <string>
#include <locale>
#include <iomanip>
#include <codecvt>
#include <algorithm>
#include <queue>
#include <map>
#include <stack>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;


/* Compressed Trie
  based on, but not identical to structure by Maly 1976
*/


const std::ctype<char> * glob_facet;
using stored_index_t = unsigned char;

#define ENCODE_AS_CHAR 0
#if ENCODE_AS_CHAR
    using stored_char_t = char;
    using stored_str_t = std::string;
    inline std::string_view to_u32string(const std::string_view & undecoded) {
        return undecoded;
    }

    std::string to_char_str(const stored_str_t & undecoded);
    std::string to_char_str(const stored_char_t & undecoded);

    inline std::string to_char_str(const stored_str_t & undecoded) {
        return undecoded;
    }

    inline std::string to_char_str(const stored_char_t & undecoded) {
        return std::string{undecoded};
    }
#else
    //conversion from https://en.cppreference.com/w/cpp/locale/codecvt
    // utility wrapper to adapt locale-bound facets for wstring/wbuffer convert
    template<class Facet>
    struct deletable_facet : Facet {
        template<class ...Args>
        deletable_facet(Args&& ...args) : Facet(std::forward<Args>(args)...) {}
        ~deletable_facet() {}
    };
    std::wstring_convert<deletable_facet<std::codecvt<char32_t, char, std::mbstate_t>>, char32_t> glob_conv32;
    std::wstring_convert<std::codecvt_utf8_utf16<char32_t>, char32_t> glob_conv8;
    using stored_str_t = std::u32string;
    using stored_char_t = char32_t;
    inline std::u32string to_u32string(const std::string_view & undecoded) {
        return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    }

    inline std::u32string to_u32string(const std::string & undecoded) {
        return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    }

    std::string to_char_str(const stored_str_t & undecoded);
    std::string to_char_str(const stored_char_t & undecoded);

    inline std::string to_char_str(const stored_str_t & undecoded) {
        return glob_conv8.to_bytes(undecoded);
    }

    inline std::string to_char_str(const stored_char_t & undecoded) {
        return glob_conv8.to_bytes(undecoded);
    }

#endif


inline std::u32string to_u32string_ci(const std::string_view & uncap_mod) {
    std::string undecoded{uncap_mod};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    std::u32string ret = glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    return ret;
}

inline std::string lower_case_version(const std::string & arg) {
    std::string undecoded{arg};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}

inline std::string upper_case_version(const std::string & arg) {
    std::string undecoded{arg};
    glob_facet->toupper(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}


constexpr std::size_t num_index_bits = 50;

// top:
//      highest bit: is terminal node
//      If that is 0:
//          second highest bit: has key that terminates with this node
//          bits 0 - 61 letter codes

constexpr uint64_t ZERO_64 = 0;
constexpr uint64_t ONE_64 = 1;
constexpr uint64_t HIGHEST_BIT = ONE_64 << 63;
constexpr uint64_t SECOND_HIGHEST_BIT = ONE_64 << 62;
constexpr uint64_t INDEX_MASK = (ONE_64 << num_index_bits) - 1;
constexpr stored_index_t NO_MATCHING_CHAR_CODE = 255;

//std::size_t max_node_index = 0;

template<typename T> void set_first_child_index(T& node, std::size_t index);
template<typename T> void flag_as_key_terminating(T & node);
template<typename T> void flag_as_terminal(T & node);
template<typename T> void flag_as_suffix(T & node, std::size_t pos);
template<typename T> void set_index(T& node, std::size_t index);

template<typename T>
inline void flag_as_key_terminating(T & node) {
    node.top |= SECOND_HIGHEST_BIT;
}

template<typename T>
inline void flag_as_terminal(T & node) {
    node.top |= HIGHEST_BIT;
}

template<typename T>
inline bool is_terminal(const T & node) {
    return node.top & HIGHEST_BIT;
}


template<typename T>
uint64_t get_index(const T & node) {
    uint64_t masked = node.bot & INDEX_MASK;
    return masked;
}


template<typename T>
inline void set_index(T& node, std::size_t index) {
    uint64_t ind = index;
    if ((ind & INDEX_MASK) != ind) {
        throw OTCError() << "not enough index field to hold pos = " << index;
    }
    // if (index > max_node_index) {
    //     max_node_index = index;
    // }
    node.bot = ind;
}

template<typename T>
inline void flag_as_suffix(T & node, std::size_t pos) {
    flag_as_terminal(node);
    set_index(node, pos);
}

template<typename T>
inline void set_first_child_index(T& node, std::size_t index) {
    assert((node.top & HIGHEST_BIT) == 0);
    set_index(node, index);
}

using ind_pair_t = std::pair<stored_index_t, uint64_t>;
using vec_ind_pair_t = std::vector<ind_pair_t>;

inline void fill_letter_and_node_indices(unsigned char curr_byte, stored_index_t offset, vec_ind_pair_t & ret, uint64_t node_index) {
    unsigned char curr_bit = 1 << 7;
    for (unsigned char i = 0; i < 8; ++i) {
        if (curr_byte & curr_bit) {
            ret.push_back(ind_pair_t{i + offset, node_index + i + offset});
        }
        curr_bit >>= 1;
    }
}


constexpr unsigned char top_first_byte = 0x03F; // 64 and 128 bit used for top 2 flags
constexpr unsigned char full_byte = 0x0FF;
constexpr unsigned char bot_last_byte = 0x0FC; // bits for 1 and 2 used for upper bits of index;

class CTrie3Node {

    
    public:
    uint64_t top, mid, bot;
    CTrie3Node() :top{ZERO_64}, mid{ZERO_64}, bot{ZERO_64} {
        //log_state();
    }
    void log_state() {
       LOG(TRACE) << " CTrie3Node( top = " << bitset<64>{top} << " mid = " << bitset<64>{mid} << " bot = " << bitset<64>{bot} << ")\n";
    }
    void flag_letter(unsigned int i) {
        uint64_t bit = 1;
        //log_state();
        if (i < 62) {
            const uint64_t shifted = (bit << (61 - i));
            top |= shifted;
            //std::cerr << " flag_letter( " << i << ") top shifted = " << bitset<64>{shifted} << " top = " << bitset<64>{top} << '\n';
            //log_state();
        } else if (i < 126) {
            bit <<= (125 - i);
            mid |= bit;
        } else {
            assert(i < 140);
            bit <<= (139 - i);
            bot |= bit;
        }
    }

    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const {
        assert(!is_terminal(*this));
        vec_ind_pair_t ret;
        ret.reserve(62 + 64 + 14);
        auto node_index = get_index(*this);
        unsigned char curr_byte = (top >> 24) & top_first_byte;
        fill_letter_and_node_indices(curr_byte, 0, ret, node_index);
        curr_byte = (top >> 16) & full_byte;
        fill_letter_and_node_indices(curr_byte, 6, ret, node_index);
        curr_byte = (top >> 8) & full_byte;
        fill_letter_and_node_indices(curr_byte, 14, ret, node_index);
        curr_byte = top & full_byte;
        fill_letter_and_node_indices(curr_byte, 22, ret, node_index);

        curr_byte = (mid >> 24) & full_byte;
        fill_letter_and_node_indices(curr_byte, 30, ret, node_index);
        curr_byte = (mid >> 16) & full_byte;
        fill_letter_and_node_indices(curr_byte, 38, ret, node_index);
        curr_byte = (mid >> 8) & full_byte;
        fill_letter_and_node_indices(curr_byte, 46, ret, node_index);
        curr_byte = mid & full_byte;
        fill_letter_and_node_indices(curr_byte, 54, ret, node_index);

        curr_byte = (bot >> 24) & full_byte;
        fill_letter_and_node_indices(curr_byte, 62, ret, node_index);
        curr_byte = (bot >> 16) & bot_last_byte;
        fill_letter_and_node_indices(curr_byte, 70, ret, node_index);
        return ret;
    }

    const static std::size_t max_num_letters = 3*64 - num_index_bits;
};


class CTrie2Node {
    public:
    uint64_t top, bot;
    CTrie2Node() :top(0),  bot(0) {
    }
    void log_state() {
       std::cerr << " CTrie2Node( top = " << bitset<64>{top} << " bot = " << bitset<64>{bot} << ")\n";
    }
    
    void flag_letter(unsigned int i) {
        uint64_t bit = 1;
        if (i < 62) {
            bit <<= (61 - i);
            top |= bit;
        }  else {
            assert(i < 76);
            bit <<= (75 - i);
            bot |= bit;
        }
    } 

    vec_ind_pair_t get_letter_and_node_indices_for_on_bits() const {
        assert(!is_terminal(*this));
        vec_ind_pair_t ret;
        ret.reserve(62 + 64 + 14);
        auto node_index = get_index(*this);
        unsigned char curr_byte = (top >> 24) & top_first_byte;
        fill_letter_and_node_indices(curr_byte, 0, ret, node_index);
        curr_byte = (top >> 16) & full_byte;
        fill_letter_and_node_indices(curr_byte, 6, ret, node_index);
        curr_byte = (top >> 8) & full_byte;
        fill_letter_and_node_indices(curr_byte, 14, ret, node_index);
        curr_byte = top & full_byte;
        fill_letter_and_node_indices(curr_byte, 22, ret, node_index);

        curr_byte = (bot >> 24) & full_byte;
        fill_letter_and_node_indices(curr_byte, 30, ret, node_index);
        curr_byte = (bot >> 16) & bot_last_byte;
        fill_letter_and_node_indices(curr_byte, 38, ret, node_index);
        return ret;
    }

    const static std::size_t max_num_letters = 2*64 - num_index_bits; 
};

using ctrie_init_set_t = std::set<stored_str_t>;
template <typename T>
class CTrieCtorHelperTemp {
    public:
    stored_str_t prefix;
    T * node_ptr;
    ctrie_init_set_t::const_iterator lower;
};

class FuzzyQueryResult {
    public:
    std::string match;
    unsigned int distance;
    float score;
    FuzzyQueryResult(): distance(0), score(0.0) {
    }
};

template <typename T>
class PartialMatch {
    public:
    PartialMatch(const stored_str_t &q, std::size_t pos, const T *nextn)
        :query(q),
        qpos(pos),
        distance(0),
        next_node(nextn),
        prev_mismatch(NO_MATCHING_CHAR_CODE){
    }
    const stored_str_t & query;
    
    std::size_t qpos;
    const stored_str_t growing_match;
    unsigned int distance;
    const T * next_node;
    stored_index_t prev_mismatch{'\0'};
};

template<typename T>
std::list<FuzzyQueryResult> fuzzy_matches(const T & trie, const stored_str_t & query_str, unsigned int max_dist);

template <typename T>
class CompressedTrie {
    public:
    CompressedTrie() {}
    private:
    using CTrieCtorHelper = CTrieCtorHelperTemp<T>;
    void init(const ctrie_init_set_t & keys, const stored_str_t & letter_var);
    void _process_prefix(const stored_str_t & curr_pref,
                         std::stack<CTrieCtorHelper> & todo_q,
                         const stored_str_t & rev_letters,
                         const ctrie_init_set_t & keys,
                         T & par_node,
                         std::map<std::string, std::size_t> & suffix2index);
    
    void _store_suffix_node(T & curr_node,
                            const stored_str_t & curr_str,
                            const stored_str_t & handled,
                            std::map<std::string, std::size_t> & suffix2index);

    T & append_node() {
        T empty;
        node_list.push_back(empty);
        return *(node_list.rbegin());
        

    }

    
    void clear() {
        letters.clear();
        node_list.clear();
        concat_suff.clear();
        node_vec.clear();
        letter_to_ind.clear();
        equivalent_letter.clear();
    }
    std::list<FuzzyQueryResult> fuzzy_matches(const stored_str_t & query_str, unsigned int max_dist) const;

    stored_index_t get_index_for_letter(const stored_char_t & c) const {
        auto ltiit = letter_to_ind.find(c);
        if (ltiit == letter_to_ind.end()) {
            return NO_MATCHING_CHAR_CODE;
        }
        return ltiit->second;
    }

    std::vector<stored_index_t> encode_as_indices(const stored_str_t & query_str) const {
        std::vector<stored_index_t> ret;
        ret.reserve(query_str.length());
        for (auto c: query_str) {
            ret.push_back(get_index_for_letter(c));
        }
        return ret;
    }
    std::unordered_map<stored_char_t, stored_index_t> letter_to_ind;
    std::vector<stored_index_t> equivalent_letter;
    stored_str_t letters;
    std::list<T> node_list;
    std::vector<char> concat_suff;
    std::vector<T> node_vec;

    friend class CompressedTrieBasedDB;
};

template<typename T>
std::list<FuzzyQueryResult> CompressedTrie<T>::fuzzy_matches(const stored_str_t & query_str,
                                                             unsigned int max_dist) const {
    std::list<FuzzyQueryResult> results;
    if (query_str.length() == 0) {
        return results;
    }
    std::list<PartialMatch<T> > alive;
    const std::vector<stored_index_t> query_as_indices = encode_as_indices(query_str);
    const T * nd = &(node_vec.at(0));
    std::size_t max_skip = std::min(std::size_t {max_dist}, query_as_indices.size());
    for (std::size_t i = 0; i < max_skip; ++i) {
        auto inds_on = nd->get_letter_and_node_indices_for_on_bits();
        auto qi = query_as_indices[i];
        auto altqi = equivalent_letter[qi];
        for (auto x : inds_on) {
            auto prefchar = x.first;
            auto next_ind = x.second;
            const T * next_nd = &(node_vec[next_ind]); 
            if (prefchar == qi || prefchar == altqi) {
                alive.push_back(PartialMatch<T>{query_str, i, next_nd});
            } else if (i + 1 < max_dist) {
                alive.push_back(PartialMatch<T>{query_str, i + 1, next_nd});
            }
        }
    }
    return results;
}


inline bool starts_with(const stored_str_t & full, const stored_str_t & pref) {
    if (full.length() < pref.length()) {
        return false;
    }
    //std::cerr << "starts_with(" << to_char_str(full) << ", " << to_char_str(pref) << ")\n";
    return 0 == full.compare(0, pref.length(), pref);
}

template <typename T>
void CompressedTrie<T>::_process_prefix(const stored_str_t & curr_pref,
                                        std::stack<CTrieCtorHelper> & todo_q,
                                        const stored_str_t & rev_letters,
                                        const ctrie_init_set_t & keys,
                                        T & par_node,
                                        std::map<std::string, std::size_t> & suffix2index) {
    stored_str_t next_pref;
    CTrieCtorHelper ctch;
    unsigned int curr_letter_index = 0;
    std::list<CTrieCtorHelper> to_queue;
    ctrie_init_set_t::const_iterator lb;
    bool has_indexed_par = false;
    for (auto letter : rev_letters) {
        next_pref = curr_pref;
        next_pref.push_back(letter);
        lb = keys.lower_bound(next_pref);
        if (lb == keys.end()) {
            break;
        }
        if (starts_with(*lb, next_pref)) {
            auto advit = lb;
            T & next_node = append_node();
            ctch.node_ptr = &next_node;
            //std::cerr << "next_node: "; next_node.log_state();
            advit++;
            if (advit != keys.end() && starts_with(*advit, next_pref)) {
                // std::cerr << " pref \"" << to_char_str(next_pref) 
                //          << "\" found in key \"" << to_char_str(*lb) << "\"\n";
                ctch.prefix = next_pref;
                ctch.lower = lb;
                todo_q.push(ctch);
            } else {
                _store_suffix_node(next_node, *lb, curr_pref, suffix2index);
            }
            par_node.flag_letter(curr_letter_index);
            if (!has_indexed_par) {
                set_first_child_index(par_node, node_list.size() - 1);
                has_indexed_par = true;
            }
        }
        curr_letter_index++;
    }
    assert(has_indexed_par);
}


template <typename T>
void CompressedTrie<T>::_store_suffix_node(T & curr_node,
                        const stored_str_t & curr_str,
                        const stored_str_t & handled,
                        std::map<std::string, std::size_t> & suffix2index) {
    const stored_str_t suffix = curr_str.substr(handled.length());
    const std::string suff_as_char = to_char_str(suffix);
    // std::cerr << " handled \"" << to_char_str(handled) << "\" suffix = \"" << suff_as_char << "\"\n";
    auto mit = suffix2index.find(suff_as_char);
    if (mit != suffix2index.end()) {
        flag_as_suffix(curr_node, mit->second);
    } else {
        std::size_t pos = concat_suff.size();
        concat_suff.insert(std::end(concat_suff), std::begin(suff_as_char), std::end(suff_as_char));
        concat_suff.push_back('\0');
        flag_as_suffix(curr_node, pos);
        suffix2index[suff_as_char] = pos;
        std::size_t suff_pref = 1;
        while (suff_pref < suff_as_char.length()) {
            std::string tmp = suff_as_char.substr(suff_pref);
            if (suffix2index.find(tmp) != suffix2index.end()) {
                // std::cerr << " found tmp \"" << tmp << "\"\n";
                break;
            }
            // std::cerr << " adding tmp \"" << tmp << "\" pos + suff_pref = " << pos + suff_pref << "\n";
            suffix2index[tmp] = pos + suff_pref;
            suff_pref++;
        }
    }
}

template <typename T>
void CompressedTrie<T>::init(const ctrie_init_set_t & keys, const stored_str_t & letter_var) {
    clear();
    // max_node_index = 0;
    if (keys.empty()) {
        return;
    }
    // sort the letters in strings, and make sure they are uniq
    std::set<stored_char_t> let_set{letter_var.begin(), letter_var.end()};
    letters = stored_str_t{let_set.begin(), let_set.end()};
    stored_str_t rev_letters = stored_str_t{letters.rbegin(), letters.rend()};
    if (letters.length() > T::max_num_letters) {
        throw OTCError() << "# of letters (" << letters.length() << ") exceeds size of CompressedTrie node type";
    }
    if (letters.length() > 253) {
        throw OTCError() << "# of letters (" << letters.length() << ") exceeds 253, so letter_to_ind value type needs to be changed.";
    }
    stored_index_t curr_ind = 0;
    for (auto nl : letters) {
        letter_to_ind[nl] = curr_ind++;
    }
    equivalent_letter.reserve(letters.length());
    equivalent_letter.clear();
    for (auto nl : letters) {
        std::string uncov = to_char_str(nl);
        std::string lccov = lower_case_version(uncov);
        stored_index_t char_ind = NO_MATCHING_CHAR_CODE;
        if (lccov != uncov) {
            auto alt = to_u32string(lccov);
            if (alt.length() != 1) {
                throw OTCError() << "lower case version of \"" << uncov << "\" was not one character: \"" << lccov << "\"\n";
            }
            char_ind = get_index_for_letter(alt[0]);
        } else {
            std::string uccov = upper_case_version(uncov);
            if (uccov != uncov) {
                auto alt = to_u32string(uccov);
                if (alt.length() != 1) {
                    throw OTCError() << "lower case version of \"" << uncov << "\" was not one character: \"" << uccov << "\"\n";
                }
                char_ind = get_index_for_letter(alt[0]);
            }
        }
        equivalent_letter.push_back(char_ind);
    }
    
    for (unsigned int eli = 0; eli < equivalent_letter.size(); ++eli) {
        if (equivalent_letter[eli] == NO_MATCHING_CHAR_CODE) {
            std::cerr << to_char_str(letters[eli]) << " = <nothing>\n";
        } else {
            std::cerr << to_char_str(letters[eli]) << " = " << to_char_str(letters[equivalent_letter[eli]]) << "\n";
        }
    }


    std::stack<CTrieCtorHelper> todo_q;
    stored_str_t curr_pref;
    std::map<std::string, std::size_t> suffix2index;
    std::string mt; 
    concat_suff.push_back('\0');
    suffix2index[mt] = 0;
    T & root_node = append_node();
    _process_prefix(curr_pref, todo_q, letters, keys, root_node, suffix2index);
    CTrieCtorHelper curr_ctch;
    while (!todo_q.empty()) {
        curr_ctch = todo_q.top();
        todo_q.pop();
        curr_pref = curr_ctch.prefix;
        T & curr_node = *(curr_ctch.node_ptr);
        bool done_with_curr = false;
        if (*curr_ctch.lower == curr_pref) {
            curr_ctch.lower++;
            if (curr_ctch.lower != keys.end() && starts_with(*curr_ctch.lower, curr_pref)) {
                flag_as_key_terminating(curr_node);
            } else {
                done_with_curr = true;
                flag_as_terminal(curr_node);
            }
        }
        if (!done_with_curr) {
            _process_prefix(curr_pref, todo_q, letters, keys, curr_node, suffix2index);
        }
    }
    // move to vector...
    node_vec.clear();
    node_vec.insert(node_vec.begin(), node_list.begin(), node_list.end());
    node_list.clear();
    auto nvs = sizeof(T)*node_vec.size();
    auto suffs = concat_suff.size();
    std::cerr << "vecsize = " << nvs << " bytes\n";
    std::cerr << "concat_suff length = " << suffs << " bytes\n";
    std::cerr << "compressed tree size = " << 4*letters.size() + nvs + suffs << " bytes\n";
    // std::cerr << "max_node_index = " << max_node_index << "\n";
    /* 
    std::cerr << "concat_suff = \"";
    for (auto c : concat_suff) {
        if (c == '\0') {
            std::cerr << "\\0";
        } else {
            std::cerr << c;
        }
    }
    std::cerr <<  "\"\n";
    */
}

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
    cout << keys.size() << " keys\n";
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


/* End compressed tree
*/
variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    options_description invisible("Invisible options");
    invisible.add_options()("taxonomy", value<string>(), "Filename for the taxonomy");
    options_description visible;
    visible.add(otc::standard_options());
    positional_options_description p;
    p.add("taxonomy", -1);
    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-tnrs-cli <taxonomy-dir> [OPTIONS]\n"
                                                    "Build data structures for name matching, and allow interactive testing.",
                                                    visible, invisible, p);

    return vm;
}


void process_taxonomy(const RichTaxonomy & taxonomy) {
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
    
    CompressedTrieBasedDB ct{all_names};


    std::cout << "Enter a query and hit return:\n";
    std::string query;
    while (std::getline(std::cin, query)) {
        std::cout << "query =\"" << query << "\"\n";
        auto results = ct.fuzzy_query(query);
        std::cout << results.size() << " matches:\n:";
        for (auto res : results) {
            std::cout << "res.match = \"" << res.match << "\", distance = " << res.distance << " score = " << res.score << "\n";
        }
        std::cout << "Enter a query and hit return:\n";
    }
    std::cerr << "EOF\n";
}


int main(int argc, char* argv[]) {
    std::locale global_locale;
    try {
        global_locale = std::locale("en_US.utf8");
    } catch (const std::exception &) {
        try {
            global_locale = std::locale("en");
        } catch (const std::exception & x) {
             std::cerr << "locale \"en_US.utf8\" or \"en\" must be supported on the system to run this program.\n";
             std::cerr << x.what() << "\n";
             return 1;
        }
    }
    auto & f = std::use_facet<std::ctype<char> >(global_locale);
    glob_facet = &f;

    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        auto taxonomy = load_rich_taxonomy(args);
        process_taxonomy(taxonomy);
    } catch (std::exception& e) {
        cerr << "otc-tnrs-cli: Error! " << e.what() << std::endl;
        return 1;
    }
}
