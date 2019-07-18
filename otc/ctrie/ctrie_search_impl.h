#ifndef OTC_CTRIE_SEARCH_IMPL_H
#define OTC_CTRIE_SEARCH_IMPL_H

#include "otc/ctrie/ctrie.h"

namespace otc {

inline std::vector<unsigned int> _init_prev_row(unsigned int dist_threshold) {
    std::vector<unsigned int> prev_row;
    prev_row.reserve(2 + 2*dist_threshold);
    unsigned int cd = 0;
    while (cd <= dist_threshold) {
        prev_row.push_back(cd++);
    }
    return prev_row;
}

template<typename T>
inline unsigned CompressedTrie<T>::_calc_dist_impl(const PartialMatch<T> &pm,
                                                   const stored_index_t * trie_suff,
                                                   const std::size_t trie_len) const {
    const stored_index_t * quer_suff = pm.query_data();
    const unsigned int dist_threshold = pm.max_distance() - pm.curr_distance();
    stored_index_t prev_trie_match_char = pm.get_prev_mismatched_trie();
    stored_index_t prev_query_char = NO_MATCHING_CHAR_CODE;
    if (prev_trie_match_char != NO_MATCHING_CHAR_CODE) {
        assert(pm.query_pos() > 0);
        prev_query_char = *(pm.get_prev_match_coded().rbegin());
    }
    auto d = _calc_dist_prim_impl(prev_query_char,
                                  quer_suff + pm.query_pos(),
                                  pm.query_len() - pm.query_pos(),
                                  trie_suff,
                                  trie_len,
                                  dist_threshold,
                                  prev_trie_match_char);
    return pm.curr_distance() + d;
}

inline unsigned int _ran_out_of_trie_score(const std::vector<unsigned int> & prev_row,
                                           std::size_t first_quer_ind,
                                           std::size_t quer_len) {
    int num_q_left = quer_len - first_quer_ind;
    if (num_q_left < 0) {
        assert(prev_row.size() == 1);
        return prev_row[0];
    }
    // add 1 because we'll decrement in the loop below, before we use it.
    unsigned int gd = 1 + num_q_left;
    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "no more trie at gd=" << gd << '\n';}
    unsigned int d = UINT_MAX;
    for (auto psc : prev_row) {
        assert(gd > 0);
        gd--;
        unsigned int elc = psc + gd;
        if (elc < d) {
            d = elc;
        }
    }
    return d;
}

template<typename T>
inline unsigned int CompressedTrie<T>::_match_cost(stored_char_t prev_q_match_char,
                                                   stored_char_t q_match_char,
                                                   stored_char_t prev_trie_match_char,
                                                   stored_char_t trie_match_char) const {
    if (q_match_char == NO_MATCHING_CHAR_CODE || trie_match_char == NO_MATCHING_CHAR_CODE) {
        return 1;
    }
    if (q_match_char == trie_match_char) {
        return 0;
    }
#   if defined(USING_EQUIL_LETTER_ARRAY)
        if (q_match_char == equivalent_letter[trie_match_char]) {
            return 0;
        }
#   endif
    if (prev_trie_match_char == NO_MATCHING_CHAR_CODE) {
        // transposition is not possible
        return 1;
    }

#   if defined(USING_EQUIL_LETTER_ARRAY)
        if ((prev_q_match_char == trie_match_char || prev_q_match_char == equivalent_letter[trie_match_char])
            && (q_match_char == prev_trie_match_char || q_match_char == equivalent_letter[prev_trie_match_char])) {
            // transpostion, don't double penalize
            return 0;
        }
#   else
        if (prev_q_match_char == trie_match_char && (q_match_char == prev_trie_match_char)) {
            // transpostion, don't double penalize
            return 0;
        }
#   endif
    return 1;
}

template<typename T>
inline unsigned int CompressedTrie<T>::_match_cost_no_transp(stored_char_t q_match_char,
                                                             stored_char_t trie_match_char) const {
    if (q_match_char == NO_MATCHING_CHAR_CODE || trie_match_char == NO_MATCHING_CHAR_CODE) {
        return 1;
    }
    if (q_match_char == trie_match_char) {
        return 0;
    }
#   if defined(USING_EQUIL_LETTER_ARRAY)
        if (q_match_char == equivalent_letter[trie_match_char]) {
            return 0;
        }
#   endif
    return 1;
}



template<typename T>
inline bool CompressedTrie<T>::_are_equivalent(stored_char_t prev_q,
                                               const stored_index_t * quer_suff,
                                               const std::size_t quer_len,
                                               const stored_index_t * trie_suff,
                                               const std::size_t trie_len,
                                               stored_index_t prev_t) const {
    if (quer_len != trie_len) {
        return false;
    }
    if (trie_len == 0) {
        return true;
    }
    if (_match_cost(prev_q, quer_suff[0], prev_t, trie_suff[0]) > 0) {
        return false;
    }
    for (std::size_t i = 1; i < trie_len; ++i) {
        if (_match_cost_no_transp(quer_suff[i], trie_suff[i]) > 0) {
            return false;
        }
    }
    return true;
}

// checks for some easy optimizations, and calls dynamic programming version if needed.
template<typename T>
inline unsigned int CompressedTrie<T>::_calc_dist_prim_impl(stored_char_t prev_quer_char,
                                                            const stored_index_t * quer_suff,
                                                            const std::size_t quer_len,
                                                            const stored_index_t * trie_suff,
                                                            const std::size_t trie_len,
                                                            const unsigned int dist_threshold,
                                                            stored_index_t prev_trie_match_char) const {
    if (DB_FUZZY_MATCH) {
        LOG(DEBUG) << "_calc_dist_prim_impl(pqc=\"" << (prev_quer_char == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[prev_quer_char])) << ", \"";
        for (auto i=0U; i < quer_len; ++i) {LOG(DEBUG) << (quer_suff[i] == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[quer_suff[i]]));}
        LOG(DEBUG) << "\", " << quer_len << ", \"";
        for (auto i=0U; i < trie_len; ++i) {LOG(DEBUG) << (trie_suff[i] == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[trie_suff[i]]));}
        LOG(DEBUG) << "\", " << trie_len << ", " << dist_threshold << ", ptc=\"" << (prev_trie_match_char == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[prev_trie_match_char])) << "\")";
    }
    if (dist_threshold == 0) {
        if (_are_equivalent(prev_quer_char, quer_suff, quer_len, trie_suff, trie_len, prev_trie_match_char)) {
            return 0;
        }
        return 1;
    }
    if (trie_len == 0) {
        return quer_len;
    }
    if (quer_len == 0) {
        return trie_len;
    }
    std::size_t eff_quer_len = quer_len;
    std::size_t eff_trie_len = trie_len;
    // strip identical pref
    if (0 == _match_cost(prev_quer_char, quer_suff[0], prev_trie_match_char, trie_suff[0])) {
        --eff_quer_len;
        --eff_trie_len;
        quer_suff++;
        trie_suff++;
        prev_trie_match_char = NO_MATCHING_CHAR_CODE;
        prev_quer_char = NO_MATCHING_CHAR_CODE;
        for (;;) {
            if (eff_quer_len == 0) {
                return eff_trie_len;
            }
            if (eff_trie_len == 0) {
                return eff_quer_len;
            }
            if (0 < _match_cost_no_transp(*quer_suff, *trie_suff)) {
                break;
            }
            --eff_quer_len;
            --eff_trie_len;
            quer_suff++;
            trie_suff++;
        }
    }
    // trim off matches at the end...
    while (0 == _match_cost_no_transp(quer_suff[eff_quer_len - 1], trie_suff[eff_trie_len - 1])) {
        --eff_quer_len;
        --eff_trie_len;
        if (eff_quer_len == 0) {
            return eff_trie_len;
        }
        if (eff_trie_len == 0) {
            return eff_quer_len;
        }
    }
    if (eff_quer_len == 1 || eff_trie_len == 1) {
        const unsigned int ldc = (eff_quer_len > eff_trie_len ? eff_quer_len - eff_trie_len : eff_quer_len - eff_trie_len);
        if (0 == _match_cost(prev_quer_char, quer_suff[0], prev_trie_match_char, trie_suff[0])) {
            return ldc;
        }
        if (eff_quer_len == 1) {
            if (eff_trie_len == 1) {
                return 1; // mismatch in only place to check...
            }
            for (auto tp = 1U; tp < eff_trie_len; ++tp) {
                if (0 == _match_cost_no_transp(quer_suff[0], trie_suff[tp])) {
                    return ldc;
                }
            }
        } else {
            for (auto tp = 1U; tp < eff_quer_len; ++tp) {
                if (0 == _match_cost_no_transp(quer_suff[tp], trie_suff[0])) {
                    return ldc;
                }
            }
        }
        return 1 + ldc;
    }
    return _dp_calc_dist_prim_impl(prev_quer_char, quer_suff, eff_quer_len, 
                                                   trie_suff, eff_trie_len, dist_threshold, prev_trie_match_char);
}

// called after preprocessing by _calc_dist_prim_impl
template<typename T>
inline unsigned int CompressedTrie<T>::_dp_calc_dist_prim_impl(stored_char_t prev_quer_char,
                                                               const stored_index_t * quer_suff,
                                                               const std::size_t quer_len,
                                                               const stored_index_t * trie_suff,
                                                               const std::size_t trie_len,
                                                               const unsigned int dist_threshold,
                                                               stored_index_t prev_trie_match_char) const {
    std::size_t prev_quer_ind = 0;
    std::size_t trie_ind = 0;
    std::vector<unsigned int> prev_row = _init_prev_row(dist_threshold);
    std::vector<unsigned int> next_row;
    next_row.reserve(prev_row.capacity());

    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "    quer_len = " << quer_len << "\"\n";}
    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "    query suffix =\"" << to_char_from_inds(quer_suff, quer_len) << "\"\n";}
    
    unsigned int leftside_cost = 1;
    
    for (;;) {
        if (trie_ind >= trie_len) {
            return _ran_out_of_trie_score(prev_row, prev_quer_ind, quer_len);
        }
        if (DB_FUZZY_MATCH) { 
            LOG(DEBUG) << "rowchar = " << to_char_str(letters[trie_suff[trie_ind]]) ;
            LOG(DEBUG) << " prev_row  (" << prev_quer_ind << ", " << trie_ind << ") " ; 
            for (auto pr : prev_row) {LOG(DEBUG) << pr << ' ';}   LOG(DEBUG) << "\"\n";
        }
        
        next_row.clear();
        std::size_t next_quer_ind = prev_quer_ind;
        // initialize first el in curr_row and update next_quer_ind if needed
        if (leftside_cost <= dist_threshold) {
            if (DB_FUZZY_MATCH) {LOG(DEBUG) << "leftside_cost = " << leftside_cost << '\n';}
            next_row.push_back(leftside_cost++);
        } else {
            next_row.push_back(prev_row[0] + 1);
        }
        assert(next_row.size() == 1);
        assert(next_quer_ind < quer_len);
        const stored_index_t trie_match_char = trie_suff[trie_ind];
        auto min_in_next_row = next_row[0];
        // we have now filled in the first cell of costs in this row.
        // Now, we grow next_row
        std::size_t match_quer_pos = next_quer_ind;
        std::size_t match_prev_index = 0;
        stored_char_t prev_q_match_char = (match_quer_pos == 0 ? prev_quer_char : quer_suff[match_quer_pos - 1]);
        for (;;) {
            unsigned int cell_left_cost = 1 + *next_row.rbegin();
            
            const stored_index_t q_match_char = quer_suff[match_quer_pos];
            
            if (DB_FUZZY_MATCH) {LOG(DEBUG) << " q_match_char = " << to_char_str(letters[q_match_char]) << ' ';}
        
            assert(match_prev_index < prev_row.size());
            unsigned int cell_match_cost = prev_row[match_prev_index];
            cell_match_cost += _match_cost(prev_q_match_char, q_match_char, prev_trie_match_char, trie_match_char);

            unsigned int cell_top_cost = 1 + (match_prev_index + 1 >= prev_row.size() ? dist_threshold : prev_row[match_prev_index + 1]);
            
            if (DB_FUZZY_MATCH) {LOG(DEBUG) << "cell costs: (l = " << cell_left_cost << ", m = " << cell_match_cost << ", t = "  << cell_top_cost << ")\n";}
            
            auto min_cost = std::min(cell_left_cost, std::min(cell_match_cost, cell_top_cost));
            next_row.push_back(min_cost);
            if (min_cost < min_in_next_row) {
                min_in_next_row = min_cost;
            }
            match_quer_pos++;
            if (match_quer_pos >= quer_len) {
                break;
            }
            match_prev_index++;
            if (match_prev_index >= prev_row.size()) {
                break;
            }
            prev_q_match_char = q_match_char;
        }
        if (min_in_next_row > dist_threshold) {
            if (DB_FUZZY_MATCH) {LOG(DEBUG) << "exceeded match threshold dist\n";}
            return dist_threshold + 1;
        }
        prev_trie_match_char = trie_match_char;

        if (DB_FUZZY_MATCH) {LOG(DEBUG) << "next_row = "; for (auto pr : next_row) {LOG(DEBUG) << pr << ' ';}   LOG(DEBUG) << "\n";}
        
        prev_quer_ind = next_quer_ind;
        trie_ind++;
        while (*next_row.rbegin() > dist_threshold ) {
            next_row.pop_back();
        }
        if (next_row[0] > dist_threshold) {
            prev_row.clear();
            auto s = next_row.begin();
            s++;
            prev_row.assign(s, next_row.end());
            prev_quer_ind += 1;
        } else {
            std::swap(prev_row, next_row);
        }
    }
}


template<typename T>
bool CompressedTrie<T>::_check_suffix_for_match(const PartialMatch<T> &pm,
                                 const stored_index_t * trie_suff,
                                 std::list<FuzzyQueryResult> & results) const {
    if (pm.has_matched_suffix(trie_suff)) {
        return false;
    }
    if (DB_FUZZY_MATCH) {db_write_pm("_check_suffix", pm);}

    const auto trie_len = count_suff_len(trie_suff);
    const auto num_q_left_ini = pm.num_q_char_left();
    
    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "    trie suffix =\"" << to_char_from_inds(trie_suff, trie_len) << "\"\n";}
    
    std::size_t abs_len_diff = (trie_len > num_q_left_ini ? trie_len - num_q_left_ini : num_q_left_ini - trie_len);
    if (abs_len_diff + pm.curr_distance() > pm.max_distance()) {
        if (DB_FUZZY_MATCH) {LOG(DEBUG) << "    bailing out because abs_len_diff = " << abs_len_diff << '\n';}
        return false;
    }
    unsigned int d;
    if (num_q_left_ini == 0 || trie_len == 0) {
        if (DB_FUZZY_MATCH) {LOG(DEBUG) << "    Match via running out of trie \n";}
        d = abs_len_diff + pm.curr_distance();
    } else {
        d = _calc_dist_impl(pm, trie_suff, trie_len);
    }
    if (d <= pm.max_distance()) {
        return pm.store_result(results, trie_suff, trie_len, d);
    }
    return false;
}


template<typename T>
void CompressedTrie<T>::db_write_pm(const char * context, const PartialMatch<T> &pm) const {
    auto & out = LOG(DEBUG);
    if (context != nullptr) {
        out << context << " ";
    }
    out << "PM(query=\"" << to_char_from_inds(pm.query_data(), pm.query_len()) << "\"";
    out << ", qpos=" << pm.query_pos();
    const auto & mc = pm.get_prev_match_coded();
    if (mc.empty()) {
        out << ", no matches to trie";
    } else {
        out << ", matched to tree=\"" << to_char_from_inds(&(mc[0]), mc.size()) << "\"";
    }
    out << ", dist=" << pm.curr_distance();
    out << ")\n";
}

template<typename T>
void CompressedTrie<T>::extend_partial_match(const PartialMatch<T> & pm,
                                             const db_qu_let_tuple pm_coords,
                                             std::list<FuzzyQueryResult> & results,
                                             CompressedTrie<T>::partial_match_queue_t & next_alive) const {
    if (DB_FUZZY_MATCH) {db_write_pm("extend", pm);}
    const T * trienode = pm.get_next_node();
    if (trienode->is_terminal()) {
        auto suffix_index = trienode->get_index();
        _check_suffix_for_match(pm, get_suffix_as_indices(suffix_index), results);
        return;
    }
    unsigned int max_dist = pm.max_distance();
    const unsigned int max_offset = std::max(std::get<0>(pm_coords), std::get<1>(pm_coords));
    if (max_offset > 2) {
        max_dist = std::min(max_dist, max_offset - 1);
    }
    if (max_offset > 4) {
        max_dist = std::min(max_dist, max_offset - 2);
    }
    auto cd = pm.curr_distance();
    auto qc = pm.query_char();
#   if defined(USING_EQUIL_LETTER_ARRAY)
        auto altqc = equivalent_letter[qc];
#   endif
    if (DB_FUZZY_MATCH) {trienode->log_state(std::cerr);}
    auto inds_on = trienode->get_letter_and_node_indices_for_on_bits();
    for (auto & x : inds_on) {
        auto trie_char = x.first;
        auto next_ind = x.second;
        const T * next_nd = &(node_vec[next_ind]);
#   if defined(USING_EQUIL_LETTER_ARRAY)
            const bool matched = (trie_char == qc || trie_char == altqc);
#       else
            const bool matched = trie_char == qc;
#       endif
        const db_qu_let_tuple nmc{std::get<0>(pm_coords) + 1, std::get<1>(pm_coords) + 1, next_nd};
        auto nait = next_alive.find(nmc);
        stored_char_t trie_letter = letters[trie_char];
        if (matched) {
            if (NEW_DB_FUZZY_MATCH) {std::cerr << "matched  " << to_char_str(trie_letter) << " in pre adding extended pm.\n";}
            bool add_el = false;
            if (nait == next_alive.end()) {
                add_el = true;
            } else if (nait->second.curr_distance() > cd) {
                add_el = true;
                next_alive.erase(nait);
            } else if (nait->second.curr_distance() == cd) {
                nait->second.add_creation_mode(PartialMatch<T>::creation_modes::MATCH);
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "updated pm creation mode: "; emit_partial(std::cerr, nmc, nait->second);}
            } else {
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "no mod because cd= " << cd << " and pm->curr_distance() = " << nait->second.curr_distance() << '\n';}
            }
            if (add_el) {
                next_alive.insert(std::make_pair(nmc, PartialMatch<T>{pm, trie_char, cd, next_nd, true}));
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "added new pm m.";  emit_partial(std::cerr, nmc, next_alive.at(nmc));
                    std::cerr << "from "; emit_partial(std::cerr, pm_coords, pm);
                }
            }
        } else if (cd + 1 <= max_dist) {
            if (NEW_DB_FUZZY_MATCH) {std::cerr << "mismatched " << to_char_str(trie_letter) << " in pre adding extended pm.\n";}
            // to do need some || logic here for transposition
            bool add_el = false;
            if (nait == next_alive.end()) {
                add_el = true;
            } else if (nait->second.curr_distance() > cd + 1) {
                add_el = true;
                next_alive.erase(nait);
            }  else if (nait->second.curr_distance() == cd + 1) {
                nait->second.add_creation_mode(PartialMatch<T>::creation_modes::MATCH);
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "updated mm pm creation mode: "; emit_partial(std::cerr, nmc, nait->second);}
            } else {
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "no mod because cd= " << cd << " and pm->curr_distance() = " << nait->second.curr_distance() << '\n';}
            }
            if (add_el) {
                next_alive.insert(std::make_pair(nmc, PartialMatch<T>{pm, trie_char, cd + 1, next_nd, false}));
                if (NEW_DB_FUZZY_MATCH) {std::cerr << "added new pm mm.";  emit_partial(std::cerr, nmc, next_alive.at(nmc));}
            }
            if (pm.can_rightshift()) {
                const db_qu_let_tuple rsc{std::get<0>(pm_coords) + 1, std::get<1>(pm_coords), next_nd};
                nait = next_alive.find(rsc);
                bool add_el = false;
                if (nait == next_alive.end()) {
                    add_el = true;
                } else if (nait->second.curr_distance() > cd + 1) {
                    add_el = true;
                    next_alive.erase(nait);
                }   else if (nait->second.curr_distance() == cd + 1) {
                    if (NEW_DB_FUZZY_MATCH) {std::cerr << "updated rsc pm creation mode: "; emit_partial(std::cerr, rsc, nait->second);}
                    nait->second.add_creation_mode(PartialMatch<T>::creation_modes::RIGHT);
                } else {
                    if (NEW_DB_FUZZY_MATCH) {std::cerr << "no mod because cd= " << cd << " and pm->curr_distance() = " << nait->second.curr_distance()<< '\n';}
                }
                if (add_el) {
                    next_alive.insert(std::make_pair(rsc,  PartialMatch<T>{pm, cd + 1, next_nd, trie_char})); // rightshift
                    if (NEW_DB_FUZZY_MATCH) {std::cerr << "added new pm r.";  emit_partial(std::cerr, rsc, next_alive.at(rsc));}
                }
            }
        }
    }
    // frameshift
    if (cd + 1 <= max_dist && pm.can_downshift()) {
        const db_qu_let_tuple dsc{std::get<0>(pm_coords), std::get<1>(pm_coords) + 1, std::get<2>(pm_coords)};
        auto nait = next_alive.find(dsc);
        bool add_el = false;
        if (nait == next_alive.end() ) {
            add_el = true;
        } else if (nait->second.curr_distance() > cd + 1) {
            add_el = true;
            next_alive.erase(nait);
        }  else if (nait->second.curr_distance() == cd + 1) {
            if (NEW_DB_FUZZY_MATCH) {std::cerr << "updated dsc pm creation mode: "; emit_partial(std::cerr, dsc, nait->second);}
            nait->second.add_creation_mode(PartialMatch<T>::creation_modes::DOWN);
        } else {
            if (NEW_DB_FUZZY_MATCH) {std::cerr << "no mod because cd= " << cd << " and pm->curr_distance() = " << nait->second.curr_distance()<< '\n';}
        }
        if (add_el) {
            next_alive.insert(std::make_pair(dsc, PartialMatch<T>{pm, cd + 1, trienode})); //downshift
            if (NEW_DB_FUZZY_MATCH) {std::cerr << "added new pm d.";  emit_partial(std::cerr, dsc, next_alive.at(dsc));}
        }
    }
    if (trienode->is_key_terminating()) {
        auto d = pm.num_q_char_left() + pm.curr_distance();
        if (d <= max_dist) {
            pm.store_result(results, nullptr, 0, d);
        }
    }
}


template<typename T>
inline void CompressedTrie<T>::_finish_query_result(FuzzyQueryResult & res) const {
    res.match_wide_char.clear();
    for (auto ind : res.match_coded) {
        assert(ind != NO_MATCHING_CHAR_CODE);
        res.match_wide_char.push_back(letters[ind]);
    }
    float sl = res.match_wide_char.length();
    res.score = ((sl - (float)res.distance)/sl);
    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "res.score = " << res.score << " res.match: " << res.match() << '\n';}
}


template<typename T>
std::list<FuzzyQueryResult> CompressedTrie<T>::fuzzy_matches(const stored_str_t & query_str,
                                                             unsigned int max_dist) const {
    if (DB_FUZZY_MATCH) {LOG(DEBUG) << "fuzzy_matches (within " << max_dist << " edits) of \"" << to_char_str(query_str) << "\"\n";}
    if (query_str.length() == 0) {
        return std::list<FuzzyQueryResult>{};
    }
    const FQuery query{query_str, encode_as_indices(query_str), max_dist};
    unsigned int num_missing_in_letters = 0;
    for (auto qai : query.as_indices) {
        if (qai == NO_MATCHING_CHAR_CODE) {
            num_missing_in_letters++;
            if (num_missing_in_letters > max_dist) {
                if (DB_FUZZY_MATCH) {LOG(DEBUG) << "match infeasible because >= " << num_missing_in_letters << " positions in the query were not in the trie.\n";}
                return std::list<FuzzyQueryResult>{};
            }
        }
    }
    // non-trivial case
    std::list<FuzzyQueryResult> results;
    const T * root_nd = &(node_vec.at(0));
    partial_match_queue_t alive;
    alive.emplace(db_qu_let_tuple{0, 0, root_nd}, PartialMatch<T>{query, root_nd});
    while (!alive.empty()) {
        if (NEW_DB_FUZZY_MATCH) {
            std::cerr << "  " << alive.size() << " alive partial matches and " << results.size() << " hits.\nkeys:";
            //for (auto it : alive) {
            //    LOG(DEBUG) << (int) std::get<0>(it.first) << ' ' << (int) std::get<1>(it.first) << ' ' << (int) std::get<2>(it.first); 
            //}
        }
        const auto alive_it = alive.begin();
        const db_qu_let_tuple curr_coord = alive_it->first;
        const auto curr_pm = alive_it->second;
        alive.erase(alive_it);
        const auto prevnalen = alive.size();
        if (NEW_DB_FUZZY_MATCH) {std::cerr << "  " << alive.size() << " alive partial matches. "; emit_partial(std::cerr, curr_coord, curr_pm); }
        extend_partial_match(curr_pm, curr_coord, results, alive);
        if (NEW_DB_FUZZY_MATCH) {if (alive.size() != prevnalen) {std::cerr << "added " << alive.size() - prevnalen << " PMs.\n"; }}
    }
    for (auto & r : results) {
        _finish_query_result(r);
    }
    return results;
}

template<typename T>
void CompressedTrie<T>::emit_partial(std::ostream & out,
                                     const db_qu_let_tuple & curr_coord,
                                     const PartialMatch<T> & curr_pm,
                                     bool newline) const {
    out << "@(" << std::get<0>(curr_coord) << ", " << std::get<1>(curr_coord) << ", " << (long) std::get<2>(curr_coord) <<  ") ";
    out << "for PM(match_coded=[";
    for (auto i : curr_pm.match_coded) {
        out << (int)i  << ", ";
    }
    out << "], match=\"" << curr_pm.get_growing_match(letters) << "\"";
    out << ", q=\"" << curr_pm.get_growing_query() << "\"";
    out << ", dist=" << curr_pm.curr_distance();
    out << ")";
    if (newline) {
        out << '\n';
    }
}



} // namespace otc
#endif
