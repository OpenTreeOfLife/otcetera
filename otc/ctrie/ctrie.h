#ifndef OTC_CTRIE_H
#define OTC_CTRIE_H

#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <stack>
#include <deque>
#include "otc/otc_base_includes.h"
#include "otc/ctrie/ctrie_node.h"
#include "otc/ctrie/str_utils.h"

namespace otc {

/* Compressed Trie
  based on, but not identical to structure by Maly 1976
*/

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
    const std::vector<stored_index_t> match_coded;
    const stored_str_t match_suffix;
    FuzzyQueryResult(const std::vector<stored_index_t> & mc, unsigned int d):
        distance(d),
        score(0.0), 
        match_coded(mc),
        match_suffix() {
    }
    FuzzyQueryResult(const std::vector<stored_index_t> & mc,
                     const stored_str_t & ms,
                     unsigned int d):
        distance(d),
        score(0.0), 
        match_coded(mc),
        match_suffix(ms) {
    }
};

template <typename T>
class PartialMatch {
    public:
    PartialMatch(std::size_t max_match_len, std::size_t pos, unsigned int dist, stored_index_t match_char_code, const T *nextn)
        :qpos(pos),
         distance(dist),
         next_node(nextn),
         prev_mismatched_q(NO_MATCHING_CHAR_CODE) {
        match_coded.reserve(max_match_len);
        if (match_char_code != NO_MATCHING_CHAR_CODE) {
            match_coded.push_back(match_char_code);
        }
    }
    
    std::size_t qpos;
    const stored_str_t growing_match;
    unsigned int distance;
    const T * next_node;
    stored_index_t prev_mismatched_q;
    std::vector<stored_index_t> match_coded;
};

class FQuery {
    public: 
    FQuery(const stored_str_t & q, const std::vector<stored_index_t> as_inds)
        :fat_str(q),
        as_indices(as_inds) {
        suffixes.resize(q.length() + 1);
    }

    const stored_str_t & fat_str;
    const std::vector<stored_index_t> as_indices;
    std::vector<stored_str_t> suffixes;
};

using suff_map_t = std::map<std::vector<stored_index_t> , std::size_t>;
template <typename T>
class CompressedTrie {
    public:
     
    CompressedTrie() {}
    std::list<FuzzyQueryResult> fuzzy_matches(const stored_str_t & query_str,
                                              unsigned int max_dist) const;
    void db_write(std::ostream & out) const;
    void db_write_words(std::ostream & out) const;
    void db_write_node(std::ostream & out, const T & nd) const;
    stored_str_t get_suffix(std::size_t suff_ind) const {
        auto sip = get_suffix_as_indices(suff_ind);
        stored_str_t ret;
        while (*sip != null_char_index) {
            ret.push_back(letters[*sip++]);
        }
        return ret;
    } 
    const stored_index_t * get_suffix_as_indices(std::size_t suff_ind) const {
        return &(concat_suff.at(suff_ind));
    }
    private:
    using CTrieCtorHelper = CTrieCtorHelperTemp<T>;
    void init(const ctrie_init_set_t & keys, const stored_str_t & letter_var);
    void _process_prefix(const stored_str_t & curr_pref,
                         std::stack<CTrieCtorHelper> & todo_q,
                         const stored_str_t & rev_letters,
                         const ctrie_init_set_t & keys,
                         T & par_node,
                         suff_map_t& suffix2index);
    
    void _store_suffix_node(T & curr_node,
                            const stored_str_t & curr_str,
                            const stored_str_t & handled,
                            suff_map_t & suffix2index);

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
    void extend_partial_match(const PartialMatch<T> &pm,
                              const FQuery & query,
                              unsigned int max_dist,
                              std::list<FuzzyQueryResult> & results,
                              std::list<PartialMatch<T> > & next_alive) const;

    stored_index_t ctrien_get_index_for_letter(const stored_char_t & c) const {
        auto ltiit = letter_to_ind.find(c);
        if (ltiit == letter_to_ind.end()) {
            return NO_MATCHING_CHAR_CODE;
        }
        return ltiit->second;
    }

    std::vector<stored_index_t> encode_as_indices(const stored_str_t & query_str,
                                                  bool null_terminate=false) const {
        std::vector<stored_index_t> ret;
        ret.reserve(query_str.length() + (null_terminate ? 1 : 0));
        for (auto c: query_str) {
            ret.push_back(ctrien_get_index_for_letter(c));
        }
        if (null_terminate) {
            ret.push_back(null_char_index);
        }
        return ret;
    }
    std::unordered_map<stored_char_t, stored_index_t> letter_to_ind;
    std::vector<stored_index_t> equivalent_letter;
    stored_str_t letters;
    std::list<T> node_list;
    std::vector<stored_index_t> concat_suff;
    std::vector<T> node_vec;
    stored_index_t null_char_index;

    friend class CompressedTrieBasedDB;
};



unsigned int calc_damerau_levenshtein(const stored_char_t * s1,
                                      std::size_t s1len,
                                      const stored_char_t * s2,
                                      std::size_t s2len,
                                      unsigned int max_dist) {
    // return distance if <= max_dist or max_dist + 1;
    throw OTCError() << "calc_dist not implemented";
    /* 
    code here, this is a mess...
    assert(s1len > 0);
    assert(s2len > 0);
    if (s1[0] == s2[0]) {
        std::size_t ns1len = s1len - 1;
        std::size_t ns2len = s2len - 1;
        if (ns1len == 0) {
            return ns2len;
        }
        if (ns2len == 0) {
            return ns1len;
        }
        return calc_damerau_levenshtein(&s1[1], ns1len, &s2[1], ns2len, max_dist)
    }
    // first pos mismatch
    if (s1len == 1 || s2len == 1) {
        if (s1len != s2len) {
            return 2; // mismatch and length diff
        }
        if (s1[1] == s2[1] || s1[1] == s2[0] || s1[0] == s2[1]) {
            return 1; // total of 1 subst or transposition
        }
        return 2;
    }
    */

}

template<typename T>
void CompressedTrie<T>::extend_partial_match(const PartialMatch<T> & pm,
                                             const FQuery & query,
                                             unsigned int max_dist,
                                             std::list<FuzzyQueryResult> & results,
                                             std::list<PartialMatch<T> > & next_alive) const {
    const T * trienode = pm.next_node;
    std::size_t num_unchecked = query.as_indices.size() - pm.qpos;
    if (ctrien_is_terminal(*trienode)) {
        auto suffix_index = ctrien_get_index(*trienode);
        /* 
        const stored_char_t * trie_suff = this->concat_suff[suffix_index];
        std::size_t triesuff_len = strlen(trie_suff);
        // check if the first char of the suffix is part of a 2-letter inversion
        bool is_inv = false;
        if (triesuff_len > 0 && pm.prev_mismatched_q != NO_MATCHING_CHAR_CODE) {
            if (trie_suff[0] == letters[pm.prev_mismatched_q]) {
                is_inv = true;
            } else {
                auto eli = equivalent_letter[pm.prev_mismatched_q];
                if (eli != NO_MATCHING_CHAR_CODE && trie_suff[0] == letters[eli]) {
                    is_inv = true;
                }
            }
        }
        if (is_inv) {
            trie_suff++;
            triesuff_len -= 1;
        }

        const stored_char_t * qsuff = &(query.fat_str[pm.qpos]);
        std::size_t abs_len_diff = (num_unchecked > triesuff_len ? num_unchecked - triesuff_len : triesuff_len - num_unchecked);
        if (pm.distance + abs_len_diff <= max_dist) {
            if (triefuff_len > 0 && num_unchecked > 0) {
                auto suff_dist = calc_damerau_levenshtein(qsuff,
                                                        num_unchecked,
                                                        trie_suff,
                                                        triesuff_len,
                                                        max_dist - pm.distance);
                if (suff_dist + pm.distance <= max_dist) {
                    results.push_back(FuzzyQueryResult{pm.match_coded,
                                                    stored_str_t{trie_suff, triesuff_len},
                                                    pm.distance + num_unchecked});
                }
            } 
        }
        */
        throw OTCError() << "terminal ctrie node node imple";
        return;
    }
    if (pm.distance + num_unchecked <= max_dist && ctrien_is_key_terminating(*pm.next_node)) {
        results.push_back(FuzzyQueryResult{pm.match_coded, pm.distance + num_unchecked});
    }

}

template<typename T>
std::list<FuzzyQueryResult> CompressedTrie<T>::fuzzy_matches(const stored_str_t & query_str,
                                                             unsigned int max_dist) const {
    std::list<FuzzyQueryResult> results;
    if (query_str.length() == 0) {
        return results;
    }
    FQuery query{query_str, encode_as_indices(query_str)};
    std::list<PartialMatch<T> > alive;
    const std::size_t max_match_len = query.as_indices.size() + 5;
    const T * root_nd = &(node_vec.at(0));
    std::size_t max_skip = std::min(std::size_t {max_dist}, query.as_indices.size());
    for (std::size_t i = 0; i < max_skip; ++i) {
        auto inds_on = root_nd->get_letter_and_node_indices_for_on_bits();
        auto qi = query.as_indices[i];
        auto altqi = equivalent_letter[qi];
        for (auto x : inds_on) {
            auto prefchar = x.first;
            auto next_ind = x.second;
            const T * next_nd = &(node_vec[next_ind]); 
            if (prefchar == qi || prefchar == altqi) {
                alive.push_back(PartialMatch<T>{max_match_len, 1, i, x.first, next_nd});
            } else if (i + 1 < max_dist) {
                PartialMatch<T> pm{max_match_len, 1, i + 1, x.first, next_nd};
                pm.prev_mismatched_q = qi;
                alive.push_back(pm);
            }
        }
    }
    // might be sufficient for full alive list?
    alive.push_back(PartialMatch<T>{max_match_len, 0, 0, NO_MATCHING_CHAR_CODE, root_nd});
    while (!alive.empty()) {
        std::list<PartialMatch<T> > next_alive;
        for (const auto & pm : alive) {
            extend_partial_match(pm, query, max_dist, results, next_alive);
        }
        std::swap(alive, next_alive);
    }    
    return results;
}




template <typename T>
void CompressedTrie<T>::_process_prefix(const stored_str_t & curr_pref,
                                        std::stack<CTrieCtorHelper> & todo_q,
                                        const stored_str_t & rev_letters,
                                        const ctrie_init_set_t & keys,
                                        T & par_node,
                                        suff_map_t & suffix2index) {
    stored_str_t next_pref;
    CTrieCtorHelper ctch;
    unsigned int curr_letter_index = 0;
    std::list<CTrieCtorHelper> to_queue;
    ctrie_init_set_t::const_iterator lb;
    bool has_indexed_par = false;
    for (auto letter : rev_letters) {
        if (letter == '\0') {
            assert(curr_letter_index == rev_letters.length() - 1);
            break;
        }
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
                ctrien_set_first_child_index(par_node, node_list.size() - 1);
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
                        suff_map_t & suffix2index) {
    const stored_str_t suffix = curr_str.substr(handled.length() + 1);
    auto suff_as_inds = encode_as_indices(suffix, true);
    //const std::string suff_as_char = to_char_str(suffix);
    // std::cerr << " handled \"" << to_char_str(handled) << "\" suffix = \"" << suff_as_char << "\"\n";
    auto mit = suffix2index.find(suff_as_inds);
    if (mit != suffix2index.end()) {
        ctrien_flag_as_suffix(curr_node, mit->second);
    } else {
        std::size_t pos = concat_suff.size();
        concat_suff.insert(std::end(concat_suff), std::begin(suff_as_inds), std::end(suff_as_inds));
        concat_suff.push_back(null_char_index);
        ctrien_flag_as_suffix(curr_node, pos);
        suffix2index[suff_as_inds] = pos;
        std::size_t suff_pref = 1;
        auto sai_it = suff_as_inds.begin();
        sai_it++;
        while (suff_pref < suff_as_inds.size() - 1) {
            std::vector<stored_index_t> tmp{sai_it, suff_as_inds.end()};
            if (suffix2index.find(tmp) != suffix2index.end()) {
                // std::cerr << " found tmp \"" << tmp << "\"\n";
                break;
            }
            // std::cerr << " adding tmp \"" << tmp << "\" pos + suff_pref = " << pos + suff_pref << "\n";
            suffix2index[tmp] = pos + suff_pref;
            suff_pref++;
            sai_it++;
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
            char_ind = ctrien_get_index_for_letter(alt[0]);
        } else {
            std::string uccov = upper_case_version(uncov);
            if (uccov != uncov) {
                auto alt = to_u32string(uccov);
                if (alt.length() != 1) {
                    throw OTCError() << "lower case version of \"" << uncov << "\" was not one character: \"" << uccov << "\"\n";
                }
                char_ind = ctrien_get_index_for_letter(alt[0]);
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
    null_char_index = letters.length();
    letters.append(1, '\0');

    std::stack<CTrieCtorHelper> todo_q;
    stored_str_t curr_pref;
    suff_map_t suffix2index;
    concat_suff.push_back(null_char_index);
    std::vector<stored_index_t> mt{1, null_char_index};
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
                ctrien_flag_as_key_terminating(curr_node);
            } else {
                done_with_curr = true;
                ctrien_flag_as_terminal(curr_node);
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



template <typename T>
void CompressedTrie<T>::db_write_node(std::ostream & out, const T & nd) const {
    if (ctrien_is_terminal(nd)) {
        auto suff_index = ctrien_get_index(nd);
        auto suff = get_suffix(suff_index);
        auto suff_str = to_char_str(suff);
        out << "TerminalNode suffix_ind=" << suff_index 
            << " suffix=" << suff 
            << " char_str=\"" << suff_str << "\"\n";
    } else {
        out << "InternalNode" << (ctrien_is_key_terminating(nd) ? "* " : " ");
        out << "  offset = " << ctrien_get_index(nd) << "\n";
        //out << "  letterbits = ";
        //nd.db_write_state(out);
        //out << "\n";
        auto vipt = nd.get_letter_and_node_indices_for_on_bits();
        for (auto & ind_pair : vipt) {
            out << "  " << to_char_str(letters[ind_pair.first]);
            out << " => node[" << std::dec << ind_pair.second << "]\n";
        }
    }
}

template <typename T>
void CompressedTrie<T>::db_write(std::ostream & out) const {
    out << "CompressedTrie<with " << sizeof(T) << " byte> nodes. Letters = \"" << to_char_str(letters) << "\"\n";
    out << "  " << node_vec.size() << " nodes:\n";
    std::size_t i = 0;
    for (auto nd : node_vec) {
        out << "node_vec[" << i++ << "] = ";
        db_write_node(out, nd); 

    }

}
template <typename T>
void CompressedTrie<T>::db_write_words(std::ostream & out) const {
    using nd_pref_pair = std::pair<const T *, stored_str_t>;
    std::deque<nd_pref_pair> todo;
    stored_str_t mt;
    todo.push_back(nd_pref_pair{&(node_vec[0]), mt});
    std::size_t i = 0;
    while (!todo.empty()) {
        auto curr_nd_pref = todo.front();
        todo.pop_front();
        auto nd_ptr = curr_nd_pref.first;
        if (ctrien_is_terminal(*nd_ptr)) {
            auto suff_index = ctrien_get_index(*nd_ptr);
            auto suff = get_suffix(suff_index);
            auto full = curr_nd_pref.second + suff;
            out << i++ << " = " << to_char_str(full) << '\n';
        } else {
            auto vipt = nd_ptr->get_letter_and_node_indices_for_on_bits();
            for (auto vipirit = vipt.rbegin(); vipirit != vipt.rend(); vipirit++) {
                const T * nn = &(node_vec[vipirit->second]);
                stored_str_t np = curr_nd_pref.second + letters[vipirit->first];
                todo.push_front(nd_pref_pair{nn, np});
            }
        }
    }
}


} // namespace otc
#endif
