#include "otc/ctrie/ctrie.h"

namespace otc {

std::vector<unsigned int> _init_prev_row(unsigned int dist_threshold) {
    std::vector<unsigned int> prev_row;
    prev_row.reserve(2 + 2*dist_threshold);
    unsigned int cd = 0;
    while (cd <= dist_threshold) {
        prev_row.push_back(cd++);
    }
    return prev_row;
}

unsigned CompressedTrie::_calc_dist_impl(const PartialMatch &pm,
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

unsigned int _ran_out_of_trie_score(const std::vector<unsigned int> & prev_row,
                                    std::size_t first_quer_ind,
                                    std::size_t quer_len) {
    int num_q_left = quer_len - first_quer_ind;
    if (num_q_left < 0) {
        assert(prev_row.size() == 1);
        return prev_row[0];
    }
    // add 1 because we'll decrement in the loop below, before we use it.
    unsigned int gd = 1 + num_q_left;
    if (DB_FUZZY_MATCH) {std::cerr << "no more trie at gd=" << gd << '\n';}
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

unsigned int CompressedTrie::_match_cost(stored_char_t prev_q_match_char,
                                         stored_char_t q_match_char,
                                         stored_char_t prev_trie_match_char,
                                         stored_char_t trie_match_char) const {
    if (q_match_char == NO_MATCHING_CHAR_CODE || trie_match_char == NO_MATCHING_CHAR_CODE) {
        return 1;
    }
    if (q_match_char == trie_match_char || q_match_char == equivalent_letter[trie_match_char]) {
        return 0;
    }
    if (prev_trie_match_char == NO_MATCHING_CHAR_CODE) {
        // transposition is not possible
        return 1;
    }
    if ((prev_q_match_char == trie_match_char || prev_q_match_char == equivalent_letter[trie_match_char])
        && (q_match_char == prev_trie_match_char || q_match_char == equivalent_letter[prev_trie_match_char])) {
        // transpostion, don't double penalize
        return 0;
    }
    return 1;
}

unsigned int CompressedTrie::_match_cost_no_transp(stored_char_t q_match_char,
                                                             stored_char_t trie_match_char) const {
    if (q_match_char == NO_MATCHING_CHAR_CODE || trie_match_char == NO_MATCHING_CHAR_CODE) {
        return 1;
    }
    if (q_match_char == trie_match_char || q_match_char == equivalent_letter[trie_match_char]) {
        return 0;
    }
    return 1;
}



bool CompressedTrie::_are_equivalent(stored_char_t prev_q,
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
unsigned int CompressedTrie::_calc_dist_prim_impl(stored_char_t prev_quer_char,
                                                  const stored_index_t * quer_suff,
                                                  const std::size_t quer_len,
                                                  const stored_index_t * trie_suff,
                                                  const std::size_t trie_len,
                                                  const unsigned int dist_threshold,
                                                  stored_index_t prev_trie_match_char) const {
    if (DB_FUZZY_MATCH) {
        std::cerr << "_calc_dist_prim_impl(pqc=\"" << (prev_quer_char == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[prev_quer_char])) << ", \"";
        for (auto i=0U; i < quer_len; ++i) {std::cerr << (quer_suff[i] == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[quer_suff[i]]));}
        std::cerr << "\", " << quer_len << ", \"";
        for (auto i=0U; i < trie_len; ++i) {std::cerr << (trie_suff[i] == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[trie_suff[i]]));}
        std::cerr << "\", " << trie_len << ", " << dist_threshold << ", ptc=\"" << (prev_trie_match_char == NO_MATCHING_CHAR_CODE ? "NO_MATCHING_CHAR_CODE" : to_char_str(letters[prev_trie_match_char])) << "\")\n";
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
unsigned int CompressedTrie::_dp_calc_dist_prim_impl(stored_char_t prev_quer_char,
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

    if (DB_FUZZY_MATCH) {std::cerr << "    quer_len = " << quer_len << "\"\n";}
    if (DB_FUZZY_MATCH) {std::cerr << "    query suffix =\"" << to_char_from_inds(quer_suff, quer_len) << "\"\n";}
    
    unsigned int leftside_cost = 1;
    
    for (;;) {
        if (trie_ind >= trie_len) {
            return _ran_out_of_trie_score(prev_row, prev_quer_ind, quer_len);
        }
        if (DB_FUZZY_MATCH) { 
            std::cerr << "rowchar = " << to_char_str(letters[trie_suff[trie_ind]]) ;
            std::cerr << " prev_row  (" << prev_quer_ind << ", " << trie_ind << ") " ; 
            for (auto pr : prev_row) {std::cerr << pr << ' ';}   std::cerr << "\"\n";
        }
        
        next_row.clear();
        std::size_t next_quer_ind = prev_quer_ind;
        // initialize first el in curr_row and update next_quer_ind if needed
        if (leftside_cost <= dist_threshold) {
            if (DB_FUZZY_MATCH) {std::cerr << "leftside_cost = " << leftside_cost << '\n';}
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
            
            if (DB_FUZZY_MATCH) {std::cerr << " q_match_char = " << to_char_str(letters[q_match_char]) << ' ';}
        
            assert(match_prev_index < prev_row.size());
            unsigned int cell_match_cost = prev_row[match_prev_index];
            cell_match_cost += _match_cost(prev_q_match_char, q_match_char, prev_trie_match_char, trie_match_char);

            unsigned int cell_top_cost = 1 + (match_prev_index + 1 >= prev_row.size() ? dist_threshold : prev_row[match_prev_index + 1]);
            
            if (DB_FUZZY_MATCH) {std::cerr << "cell costs: (l = " << cell_left_cost << ", m = " << cell_match_cost << ", t = "  << cell_top_cost << ")\n";}
            
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
            if (DB_FUZZY_MATCH) {std::cerr << "exceeded match threshold dist\n";}
            return dist_threshold + 1;
        }
        prev_trie_match_char = trie_match_char;

        if (DB_FUZZY_MATCH) {std::cerr << "next_row = "; for (auto pr : next_row) {std::cerr << pr << ' ';}   std::cerr << "\n";}
        
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


bool CompressedTrie::_check_suffix_for_match(const PartialMatch &pm,
                                             const stored_index_t * trie_suff,
                                             std::vector<FuzzyQueryResult> & results) const {
    if (pm.has_matched_suffix(trie_suff)) {
        return false;
    }
    if (DB_FUZZY_MATCH) {db_write_pm("_check_suffix", pm);}

    const auto trie_len = count_suff_len(trie_suff);
    const auto num_q_left_ini = pm.num_q_char_left();
    
    if (DB_FUZZY_MATCH) {std::cerr << "    trie suffix =\"" << to_char_from_inds(trie_suff, trie_len) << "\"\n";}
    
    std::size_t abs_len_diff = (trie_len > num_q_left_ini ? trie_len - num_q_left_ini : num_q_left_ini - trie_len);
    if (abs_len_diff + pm.curr_distance() > pm.max_distance()) {
        if (DB_FUZZY_MATCH) {std::cerr << "    bailing out because abs_len_diff = " << abs_len_diff << '\n';}
        return false;
    }
    unsigned int d;
    if (num_q_left_ini == 0 || trie_len == 0) {
        if (DB_FUZZY_MATCH) {std::cerr << "    Match via running out of trie \n";}
        d = abs_len_diff + pm.curr_distance();
    } else {
        d = _calc_dist_impl(pm, trie_suff, trie_len);
    }
    if (d <= pm.max_distance()) {
        return pm.store_result(results, trie_suff, trie_len, d);
    }
    return false;
}


void CompressedTrie::db_write_pm(const char * context, const PartialMatch &pm) const {
    auto & out = std::cerr;
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

void CompressedTrie::extend_partial_match(const PartialMatch & pm,
                                          std::vector<FuzzyQueryResult> & results,
                                          std::list<PartialMatch> & next_alive) const {
    if (DB_FUZZY_MATCH) {db_write_pm("extend", pm);}
    const CTrieNode * trienode = pm.get_next_node();
    if (trienode->is_terminal()) {
        auto suffix_index = trienode->get_index();
        _check_suffix_for_match(pm, get_suffix_as_indices(suffix_index), results);
        return;
    }
    const unsigned int max_dist = pm.max_distance();
    auto cd = pm.curr_distance();
    auto qc = pm.query_char();
    auto altqc = equivalent_letter[qc];
    if (DB_FUZZY_MATCH) {trienode->log_state();}

    for (auto [trie_char, next_ind] : trienode->children())
    {
        const CTrieNode * next_nd = &(node_vec[next_ind]);
        if (trie_char == qc || trie_char == altqc) {
            if (DB_FUZZY_MATCH) {std::cerr << "matched " << to_char_str(letters[trie_char]) << " in pre adding extended pm.\n";}
            next_alive.push_back(PartialMatch{pm, trie_char, cd, next_nd, false});
        } else if (cd + 1 <= max_dist) {
            if (DB_FUZZY_MATCH) {std::cerr << "mismatched " << to_char_str(letters[trie_char]) << " in pre adding extended pm.\n";}
            next_alive.push_back(PartialMatch{pm, trie_char, cd + 1, next_nd, true});
            if (pm.can_rightshift()) {
                next_alive.push_back(PartialMatch{pm, cd + 1, next_nd, trie_char}); // rightshift
            }
        }
    }
    // frameshift
    if (cd + 1 <= max_dist && pm.can_downshift()) {
        next_alive.push_back(PartialMatch{pm, cd + 1, trienode}); //downshift
    }
    if (trienode->is_key_terminating()) {
        auto d = pm.num_q_char_left() + pm.curr_distance();
        if (d <= max_dist) {
            pm.store_result(results, nullptr, 0, d);
        }
    }
}


void CompressedTrie::_finish_query_result(FuzzyQueryResult & res) const {
    res.match_wide_char.clear();
    for (auto ind : res.match_coded) {
        assert(ind != NO_MATCHING_CHAR_CODE);
        res.match_wide_char.push_back(letters[ind]);
    }
    float sl = res.match_wide_char.length();
    res.score = ((sl - (float)res.distance)/sl);
    if (DB_FUZZY_MATCH) {std::cerr << "res.score = " << res.score << " res.match: " << res.match() << '\n';}
}


std::vector<FuzzyQueryResult> CompressedTrie::fuzzy_matches(const stored_str_t & query_str,
                                                            unsigned int max_dist) const {
    if (DB_FUZZY_MATCH) {std::cerr << "fuzzy_matches (within " << max_dist << " edits) of \"" << to_char_str(query_str) << "\"\n";}
    if (query_str.length() == 0) {
        return std::vector<FuzzyQueryResult>{};
    }
    const FQuery query{query_str, encode_as_indices(query_str), max_dist};
    unsigned int num_missing_in_letters = 0;
    for (auto qai : query.as_indices) {
        if (qai == NO_MATCHING_CHAR_CODE) {
            num_missing_in_letters++;
            if (num_missing_in_letters > max_dist) {
                if (DB_FUZZY_MATCH) {std::cerr << "match infeasible because >= " << num_missing_in_letters << " positions in the query were not in the trie.\n";}
                return std::vector<FuzzyQueryResult>{};
            }
        }
    }
    // non-trivial case
    std::vector<FuzzyQueryResult> results;
    results.reserve(20);
    const CTrieNode * root_nd = &(node_vec.at(0));
    std::list<PartialMatch> alive;
    alive.push_back(PartialMatch{query, root_nd});
    while (!alive.empty()) {
        if (DB_FUZZY_MATCH) {std::cerr << "  " << alive.size() << " alive partial matches and " << results.size() << " hits.\n";}
        std::list<PartialMatch> next_alive;
        for (const auto & pm : alive) {
            auto prevnalen = next_alive.size();
            extend_partial_match(pm, results, next_alive);
            if (DB_FUZZY_MATCH) {if (next_alive.size() != prevnalen) {std::cerr << "added " << next_alive.size() - prevnalen << " PMs.\n"; }}
        }
        std::swap(alive, next_alive);
    }
    for (auto & r : results) {
        _finish_query_result(r);
    }
    return results;
}


} // namespace otc

