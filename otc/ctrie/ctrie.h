#ifndef OTC_CTRIE_H
#define OTC_CTRIE_H

#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <stack>
#include <deque>
#include <climits>
#include "otc/otc_base_includes.h"
#include "otc/ctrie/search_data_models.h"
#include "otc/ctrie/ctrie_node.h"

namespace otc {
constexpr bool DB_FUZZY_MATCH = false;
/* Compressed Trie
  based on, but not identical to structure by Maly 1976
*/

using ctrie_init_set_t = std::set<stored_str_t>;

class CTrieCtorHelper
{
public:
    stored_str_t prefix;
    CTrieNode* node_ptr;
    ctrie_init_set_t::const_iterator lower;
};


using suff_map_t = std::map<std::vector<stored_index_t> , std::size_t>;

class CompressedTrie {
    public:
     
    CompressedTrie() {
    }

    // test
    CompressedTrie(const std::string &inp_letters) {
        std::set<char> ils{std::begin(inp_letters), std::end(inp_letters)};
        std::string uniq{std::begin(ils), std::end(ils)};
        letters = to_u32string(uniq);
        null_char_index = letters.size();
    }

    std::vector<stored_index_t> using_letter_order_to_encode(const stored_str_t & query_str,
                                                  bool null_terminate=false) const {
        std::vector<stored_index_t> ret;
        ret.reserve(query_str.length() + (null_terminate ? 1 : 0));
        for (auto c: query_str) {
            auto li = letters.find(c);
            if (li == std::string::npos) {
                ret.push_back(NO_MATCHING_CHAR_CODE);
            } else {
                ret.push_back((stored_char_t) li);
            }
        }
        if (null_terminate) {
            ret.push_back(null_char_index);
        }
        return ret;
    }

    std::vector<stored_index_t> encode_as_indices(const stored_str_t & query_str,
                                                  bool null_terminate=false) const {
        std::vector<stored_index_t> ret;
        ret.reserve(query_str.length() + (null_terminate ? 1 : 0));
        for (auto c: query_str) {
            ret.push_back(ctrie_get_index_for_letter(c));
        }
        if (null_terminate) {
            ret.push_back(null_char_index);
        }
        return ret;
    }
    
    std::vector<FuzzyQueryResult> fuzzy_matches(const stored_str_t & query_str,
                                                unsigned int max_dist) const;
    
    void db_write(std::ostream & out) const;
    
    void db_write_words(std::ostream & out) const;
    
    void db_write_node(std::ostream & out, const CTrieNode & nd) const;

    std::string to_char_from_inds(const stored_index_t * p, std::size_t len) const {
        std::string ret;
        for (std::size_t i = 0; i < len; ++i) {
            auto let_ind = p[i];
            if (let_ind == NO_MATCHING_CHAR_CODE) {
                ret += "?";
            } else if (let_ind == null_char_index) {
            } else {
                auto nl = letters[let_ind];
                try {
                    ret += to_char_str(nl);
                } catch (...) {
                    LOG(ERROR) << "error translating p[" << i << "] = "<< int(p[i]) << " where letters = \"" << to_char_str(letters) << "\"\n";
                    throw;
                }
            }
        }
        return ret;
    }
    
    stored_str_t get_suffix(std::size_t suff_ind) const {
        auto sip = get_suffix_as_indices(suff_ind);
        stored_str_t ret;
        while (*sip != null_char_index) {
            ret.push_back(letters[*sip++]);
        }
        return ret;
    }

    std::size_t count_suff_len(const stored_index_t *c) const {
        std::size_t i = 0;
        while (*c != null_char_index) {
            assert(*c != NO_MATCHING_CHAR_CODE);
            c++;
            i++;
        }
        return i;
    }
    
    const stored_index_t * get_suffix_as_indices(std::size_t suff_ind) const {
        return &(concat_suff.at(suff_ind));
    }
    
    const stored_index_t * get_suffix_ptr(const CTrieNode& node) const
    {
        return get_suffix_as_indices( node.get_index() );
    }

    std::size_t get_suffix_length(const CTrieNode& node) const
    {
        assert(node.is_terminal());
        auto suffix_index = node.get_index();
        auto ptr = get_suffix_as_indices(suffix_index);
        return count_suff_len(ptr);
    }

    private:
    
    void init(const ctrie_init_set_t & keys, const stored_str_t & letter_var);
    void _process_prefix(const stored_str_t & curr_pref,
                         std::stack<CTrieCtorHelper> & todo_q,
                         const stored_str_t & rev_letters,
                         const ctrie_init_set_t & keys,
                         CTrieNode & par_node,
                         suff_map_t& suffix2index);
    
    void _store_suffix_node(CTrieNode & curr_node,
                            const stored_str_t & curr_str,
                            const stored_str_t & handled,
                            suff_map_t & suffix2index);
   
    void _finish_query_result(FuzzyQueryResult & res) const;
   
    CTrieNode & append_node() {
        CTrieNode empty;
        node_list.push_back(empty);
        return *(node_list.rbegin());
    }

    void clear() {
        letters.clear();
        node_list.clear();
        concat_suff.clear();
        node_vec.clear();
        letter_to_ind.clear();
    }

    void extend_partial_match(const std::vector<stored_index_t>& query,
                              const unsigned int max_dist,
                              const CTrieNode* curr_node,
                              class dp_matrix& score,
                              std::vector<stored_index_t>& match_coded,
                              std::vector<FuzzyQueryResult> & results) const;

    stored_index_t ctrie_get_index_for_letter(const stored_char_t & c) const {
        auto ltiit = letter_to_ind.find(c);
        if (ltiit == letter_to_ind.end()) {
            return NO_MATCHING_CHAR_CODE;
        }
        return ltiit->second;
    }

    void all_descendants(stored_str_t& prefix, uint64_t index, std::vector<std::string>& results) const;
    std::vector<std::string> prefix_query(const stored_str_t& uquery) const;

    std::unordered_map<stored_char_t, stored_index_t> letter_to_ind;
    stored_str_t letters;
    std::list<CTrieNode> node_list;
    std::vector<stored_index_t> concat_suff;
    std::vector<CTrieNode> node_vec;
    stored_index_t null_char_index;

    friend class CompressedTrieBasedDB;
};

} // namespace otc

#endif
