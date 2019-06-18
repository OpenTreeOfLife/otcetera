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
    std::vector<stored_index_t> match_coded;
    const stored_str_t match_suffix;
    FuzzyQueryResult(const std::vector<stored_index_t> & mc,
                     const stored_index_t * trie_suff,
                     std::size_t suff_len,
                     unsigned int d)
      :distance(d),
      score(0.0), 
      match_coded(mc),
      match_suffix() {
        match_coded.insert(std::begin(match_coded), trie_suff, trie_suff + suff_len);
    }
    /* 
    FuzzyQueryResult(const std::vector<stored_index_t> & mc,
                     const stored_str_t & ms,
                     unsigned int d):
        distance(d),
        score(0.0), 
        match_coded(mc),
        match_suffix(ms) {
        }*/
};

class FQuery {
    public: 
    FQuery(const stored_str_t & q,
           const std::vector<stored_index_t> as_inds,
           unsigned int max_d)
        :fat_str(q),
        as_indices(as_inds),
        max_dist(max_d) {
        suffixes.resize(q.length() + 1);
    }

    const stored_str_t & fat_str;
    const std::vector<stored_index_t> as_indices;
    unsigned int max_dist;
    std::vector<stored_str_t> suffixes;
};

template <typename T>
class PartialMatch {
    public:
    PartialMatch(const FQuery & q,
                 const T *nextn)
        :query(q),
         qpos(0),
         distance(0),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE) {
        auto max_match_len = q.max_dist + q.as_indices.size();
        match_coded.reserve(max_match_len);
    }
    // create a partial match previous match and a char match
    PartialMatch(const PartialMatch & prevpm,
                 stored_index_t match_char,
                 unsigned int start_dist,
                 const T *nextn,
                 bool was_match)
        :query(prevpm.query),
         qpos(prevpm.qpos + 1),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE) {
        match_coded.reserve(prevpm.match_coded.capacity());
        match_coded = prevpm.match_coded;
        match_coded.push_back(match_char);
        if (!was_match) {
            prev_mismatched_trie = match_char;
        }
        assert(nextn != prevpm.next_node);
    }
    // create a partial match from a gap, moving through query but not trie
    PartialMatch(const PartialMatch & prevpm,
                 unsigned int start_dist,
                 const T *nextn)
        :query(prevpm.query),
         qpos(prevpm.qpos + 1),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE) {
        match_coded.reserve(prevpm.match_coded.capacity());
        match_coded = prevpm.match_coded;
        assert(nextn == prevpm.next_node);

    }
    // create a partial match from a gap, moving through trie but not query
    PartialMatch(const PartialMatch & prevpm,
                 unsigned int start_dist,
                 const T *nextn, 
                 stored_index_t match_char)
        :query(prevpm.query),
         qpos(prevpm.qpos),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE) {
        match_coded.reserve(prevpm.match_coded.capacity());
        match_coded = prevpm.match_coded;
        match_coded.push_back(match_char);
        assert(nextn != prevpm.next_node);
    }

    
    const T * get_next_node() const {
        return next_node;
    }
    
    std::size_t num_q_char_left() const {
        return query.as_indices.size() - qpos;
    }

    unsigned int curr_distance() const {
        return distance;
    }

    unsigned int max_distance() const {
        return query.max_dist;
    }

    const std::vector<stored_index_t> & get_prev_match_coded() const {
        return match_coded;
    }

    const stored_index_t * query_data() const {
        return &(query.as_indices[0]);
    }

    stored_index_t query_char() const {
        return query.as_indices[qpos];
    }

    std::size_t query_len() const {
        return query.as_indices.size();
    }

    std::size_t query_pos() const {
        return qpos;
    }

    stored_index_t get_prev_mismatched_trie() const {
        return prev_mismatched_trie;
    }

    private:
    const FQuery & query;
    std::size_t qpos;
    stored_str_t growing_match;
    unsigned int distance;
    const T * next_node;
    stored_index_t prev_mismatched_trie;
    std::vector<stored_index_t> match_coded;
};


using suff_map_t = std::map<std::vector<stored_index_t> , std::size_t>;
template <typename T>
class CompressedTrie {
    public:
     
    CompressedTrie() {
    }

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
                              std::list<FuzzyQueryResult> & results,
                              std::list<PartialMatch<T> > & next_alive) const;

    void _check_suffix_for_match(const PartialMatch<T> &pm,
                                 const stored_index_t * suffix,
                                 std::list<FuzzyQueryResult> & results) const;

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

template<typename T>
void CompressedTrie<T>::_check_suffix_for_match(const PartialMatch<T> &pm,
                                 const stored_index_t * trie_suff,
                                 std::list<FuzzyQueryResult> & results) const {
    auto num_tr_left = count_suff_len(trie_suff);
    auto num_q_left = pm.num_q_char_left();
    std::size_t abs_len_diff = (num_tr_left > num_q_left ? num_tr_left - num_q_left : num_q_left - num_tr_left);
    if (abs_len_diff + pm.curr_distance() > pm.max_distance()) {
        return;
    }
    if (num_q_left == 0 || num_tr_left == 0) {
        results.push_back(FuzzyQueryResult{pm.get_prev_match_coded(), trie_suff, num_tr_left, abs_len_diff + pm.curr_distance()});
        return;
    }
    const unsigned int md = pm.max_distance();
    std::vector<unsigned int> prev_row;
    std::pair<std::size_t, std::size_t> prev_lt_coord{pm.query_pos(), 0};
    prev_row.reserve(1 + 2*md);
    unsigned int cd = pm.curr_distance();
    while (cd < md) {
        prev_row.push_back(cd++);
    }
    const stored_index_t * q_suff = pm.query_data();
    
    std::vector<unsigned int> curr_row;
    std::pair<std::size_t, std::size_t> curr_lt_coord{pm.query_pos(), 1};
    curr_row.reserve(1 + 2*md);
    auto q_size = pm.query_len();
    unsigned int leftside_cost = cd + 1;
    stored_index_t prev_trie_match_char = pm.get_prev_mismatched_trie();
    for (;;) {
        curr_row.clear();
        curr_lt_coord = prev_lt_coord;
        curr_lt_coord.second += 1; // moving 1 down in the trie_suff
        if (curr_lt_coord.second > num_tr_left) {
            // ran out of trie characters. Add as many gap costs as needed to each el in prev_row
            unsigned int gd = num_q_left - curr_lt_coord.first;
            unsigned int d = md;
            for (auto psc : prev_row) {
                unsigned int elc = psc + gd;
                assert(gd > 0);
                gd--;
                if (elc < d) {
                    d = elc;
                }
            }
            if (d < md) {
                results.push_back(FuzzyQueryResult{pm.get_prev_match_coded(), trie_suff, num_tr_left, d});
            }
            return;
        }
        for (std::size_t qp_ind = 0;; ++qp_ind) {
            if (leftside_cost <= md) {
                curr_row.push_back(leftside_cost++);
            } else {
                std::size_t x_shift = 0;
                for (;;) {
                    if (x_shift >= prev_row.size()) {
                        return ; // ran out of possible matches
                    }
                    if (prev_row[x_shift] < md) {
                        curr_row.push_back(prev_row[0] + 1);
                        break;
                    } else {
                        // bump over the coordinate of the leftmost-topmost point we'll consider
                        curr_lt_coord.first += 1;
                    }
                    ++x_shift;
                }
            }
            assert(curr_row.size() == 1);
            // we have now filled in the first cell of costs in this row (or exited w/o a match)
            // grow curr_row while costs are still under the min and we have characters left in the query...
            std::size_t match_prev_index = curr_lt_coord.first - prev_lt_coord.first;
            std::size_t match_q_suff_pos = curr_lt_coord.first;
            stored_index_t trie_match_char = trie_suff[prev_lt_coord.second];
            if (match_q_suff_pos >= q_size) {
                // ran out of query characters. Add as many gap costs as needed to each el in prev_row
                assert(prev_row.size() == 1);
                unsigned int gd = num_tr_left - prev_lt_coord.second;
                unsigned int d = gd + prev_row[0];
                if (d < md) {
                    results.push_back(FuzzyQueryResult{pm.get_prev_match_coded(), trie_suff, num_tr_left, d});
                }
                return;
            }
            for (;;) {
                unsigned int cell_left_cost = *curr_row.rbegin();
                unsigned int cell_match_cost, cell_top_cost;
                stored_index_t q_match_char = q_suff[match_q_suff_pos];
                cell_match_cost = (match_prev_index >= prev_row.size() ? md : prev_row[match_prev_index]);
                if (q_match_char != trie_match_char) {
                    if (prev_trie_match_char == NO_MATCHING_CHAR_CODE) {
                        cell_match_cost += 1;
                    } else if (prev_trie_match_char == q_match_char && trie_match_char == q_suff[match_q_suff_pos - 1]) {
                        // transposition of 2 characters. Already penalized, don't add another ...
                    } else {
                        cell_match_cost += 1;
                    }
                }
                cell_top_cost = 1 + (match_prev_index + 1 >= prev_row.size() ? md : prev_row[match_prev_index + 1]);
                auto min_cost = std::min(cell_left_cost, std::min(cell_match_cost, cell_top_cost));
                if (min_cost > md) {
                    break;
                }
                match_q_suff_pos++;
                match_prev_index++;
            }
            prev_trie_match_char = trie_match_char;
        }
        std::swap(prev_row, curr_row);
        prev_lt_coord = curr_lt_coord;
    }
}

template<typename T>
void CompressedTrie<T>::extend_partial_match(const PartialMatch<T> & pm,
                                             std::list<FuzzyQueryResult> & results,
                                             std::list<PartialMatch<T> > & next_alive) const {
    const T * trienode = pm.get_next_node();
    if (ctrien_is_terminal(*trienode)) {
        auto suffix_index = ctrien_get_index(*trienode);
        _check_suffix_for_match(pm, get_suffix_as_indices(suffix_index), results);
        return;
    }
    const unsigned int max_dist = pm.max_distance();
    auto cd = pm.curr_distance();
    auto qc = pm.query_char();
    auto altqc = equivalent_letter[qc];
    auto inds_on = trienode->get_letter_and_node_indices_for_on_bits();
    for (auto x : inds_on) {
        auto trie_char = x.first;
        auto next_ind = x.second;
        const T * next_nd = &(node_vec[next_ind]);
        if (trie_char == qc || trie_char == altqc) {
            next_alive.push_back(PartialMatch<T>{pm, trie_char, cd, next_nd, false});
        } else if (cd + 1 < max_dist) {
            next_alive.push_back(PartialMatch<T>{pm, trie_char, cd + 1, next_nd, true});
            next_alive.push_back(PartialMatch<T>{pm, cd + 1, next_nd, trie_char}); // rightshift
        }
    }
    // frameshift
    if (cd + 1 < max_dist) {
        next_alive.push_back(PartialMatch<T>{pm, cd + 1, trienode}); //downshift
    }
 
}

template<typename T>
std::list<FuzzyQueryResult> CompressedTrie<T>::fuzzy_matches(const stored_str_t & query_str,
                                                             unsigned int max_dist) const {
    std::list<FuzzyQueryResult> results;
    if (query_str.length() == 0) {
        return results;
    }
    const FQuery query{query_str, encode_as_indices(query_str), max_dist};
    const T * root_nd = &(node_vec.at(0));
    std::list<PartialMatch<T> > alive;
    alive.push_back(PartialMatch<T>{query, root_nd});
    while (!alive.empty()) {
        std::list<PartialMatch<T> > next_alive;
        for (const auto & pm : alive) {
            extend_partial_match(pm, results, next_alive);
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
