#include "otc/ctrie/ctrie.h"

namespace otc {

void CompressedTrie::_process_prefix(const stored_str_t & curr_pref,
                                     std::stack<CTrieCtorHelper> & todo_q,
                                     const stored_str_t & rev_letters,
                                     const ctrie_init_set_t & keys,
                                     CTrieNode & par_node,
                                     suff_map_t & suffix2index) {
    stored_str_t next_pref;
    CTrieCtorHelper ctch;
    unsigned int curr_letter_index = 0;
    std::list<CTrieCtorHelper> to_queue;
    ctrie_init_set_t::const_iterator lb;
    bool has_indexed_par = false;
    static const std::string TARGET_THIN_STR{"A"};
    static const stored_str_t TARGET_STR = to_u32string(TARGET_THIN_STR);
    bool had_target_pref = false;
    for (auto letter : rev_letters) {
        if (letter == '\0') {
            assert(curr_letter_index == rev_letters.length() - 1);
            break;
        }
        next_pref = curr_pref;
        next_pref.push_back(letter);
        lb = keys.lower_bound(next_pref);
        if (lb != keys.end()) {
            if (starts_with(*lb, next_pref)) {
                auto advit = lb;
                CTrieNode & next_node = append_node();
                ctch.node_ptr = &next_node;
                advit++;
                if (advit != keys.end() && starts_with(*advit, next_pref)) {
                    ctch.prefix = next_pref;
                    ctch.lower = lb;
                    todo_q.push(ctch);
                } else {
                    _store_suffix_node(next_node, *lb, curr_pref, suffix2index);
                }
                par_node.flag_letter(curr_letter_index);
                if (!has_indexed_par) {
                    par_node.set_first_child_index(node_list.size() - 1);
                    has_indexed_par = true;
                }
            }
        }
        curr_letter_index++;
    }
    assert(has_indexed_par);

}



void CompressedTrie::_store_suffix_node(CTrieNode & curr_node,
                                        const stored_str_t & curr_str,
                                        const stored_str_t & handled,
                                        suff_map_t & suffix2index) {
    const stored_str_t suffix = curr_str.substr(handled.length() + 1);
    auto suff_as_inds = encode_as_indices(suffix, true);
    //const std::string suff_as_char = to_char_str(suffix);
    // std::cerr << " handled \"" << to_char_str(handled) << "\" suffix = \"" << suff_as_char << "\"\n";
    auto mit = suffix2index.find(suff_as_inds);
    if (mit != suffix2index.end()) {
        curr_node.flag_as_suffix(mit->second);
    } else {
        std::size_t pos = concat_suff.size();
        concat_suff.insert(std::end(concat_suff), std::begin(suff_as_inds), std::end(suff_as_inds));
        concat_suff.push_back(null_char_index);
        curr_node.flag_as_suffix(pos);
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


void CompressedTrie::fill_equivalent_letter_array() {
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
            char_ind = ctrie_get_index_for_letter(alt[0]);
        } else {
            std::string uccov = upper_case_version(uncov);
            if (uccov != uncov) {
                auto alt = to_u32string(uccov);
                if (alt.length() != 1) {
                    throw OTCError() << "lower case version of \"" << uncov << "\" was not one character: \"" << uccov << "\"\n";
                }
                char_ind = ctrie_get_index_for_letter(alt[0]);
            }
        }
        equivalent_letter.push_back(char_ind);
    }
}

void CompressedTrie::init(const ctrie_init_set_t & keys, const stored_str_t & letter_var) {
    clear();
    // max_node_index = 0;
    if (keys.empty()) {
        return;
    }
    // sort the letters in strings, and make sure they are uniq
    std::set<stored_char_t> let_set{letter_var.begin(), letter_var.end()};
    letters = stored_str_t{let_set.begin(), let_set.end()};
    stored_str_t rev_letters = stored_str_t{letters.rbegin(), letters.rend()};
    if (letters.length() >= CTrieNode::END_LETTER_INDEX) {
        throw OTCError() << "# of letters (" << letters.length() << ") exceeds size of CompressedTrie node type";
    }
    if (letters.length() > 253) {
        throw OTCError() << "# of letters (" << letters.length() << ") exceeds 253, so letter_to_ind value type needs to be changed.";
    }
    stored_index_t curr_ind = 0;
    for (auto nl : letters) {
        letter_to_ind[nl] = curr_ind++;
    }
    fill_equivalent_letter_array();
    null_char_index = letters.length();
    letters.append(1, '\0');

    std::stack<CTrieCtorHelper> todo_q;
    stored_str_t curr_pref;
    suff_map_t suffix2index;
    concat_suff.push_back(null_char_index);
    std::vector<stored_index_t> mt{1, null_char_index};
    suffix2index[mt] = 0;
    CTrieNode & root_node = append_node();
    static const std::string TARGET_THIN_STR{"A"};
    static const stored_str_t TARGET_STR = to_u32string(TARGET_THIN_STR);
    unsigned int target_ind = UINT_MAX;
    assert(node_list.size() == 1);
    // std::cerr << "ROOT before any children:"; root_node.log_state();
    _process_prefix(curr_pref, todo_q, letters, keys, root_node, suffix2index);
    // std::cerr << "ROOT after first _process_prefix:"; root_node.log_state();
    CTrieCtorHelper curr_ctch;
    while (!todo_q.empty()) {
        curr_ctch = todo_q.top();
        todo_q.pop();
        curr_pref = curr_ctch.prefix;
        CTrieNode & curr_node = *(curr_ctch.node_ptr);
        bool done_with_curr = false;
        if (*curr_ctch.lower == curr_pref) {
            curr_ctch.lower++;
            if (curr_ctch.lower != keys.end() && starts_with(*curr_ctch.lower, curr_pref)) {
                curr_node.flag_as_key_terminating();
            } else {
                done_with_curr = true;
                curr_node.flag_as_terminal();
            }
        }
        if (!done_with_curr) {
            _process_prefix(curr_pref, todo_q, letters, keys, curr_node, suffix2index);
        }
        if (curr_pref == TARGET_STR) {
            // std::cerr << "MATCH TARGET: "; curr_node.log_state();
            std::size_t i = 0;
            for (const auto & nd : node_list) {
                if (&(nd) == &curr_node) {
                    target_ind = i;
                    break;
                }
                i++;
            }
            // std::cerr << "MATCH TARGET at node " << target_ind << "\n";
        }
    }

    if (target_ind != UINT_MAX) {
        std::size_t i = 0;
        for (const auto & nd : node_list) {
            if (i++ == target_ind) {
                // std::cerr << "MATCH TARGET from node list spot " << target_ind << " = ";
                nd.log_state();
            }
        }
    }
            
    // move to vector...
    node_vec.clear();
    node_vec.insert(node_vec.begin(), node_list.begin(), node_list.end());
    node_list.clear();

    if (target_ind != UINT_MAX) {
        // std::cerr << "MATCH TARGET from node vector spot " << target_ind << " = ";
        node_vec[target_ind].log_state();
    }
    
    if (DB_FUZZY_MATCH) {node_vec[0].log_state();}
    // std::cerr << "ROOT:"; node_vec[0].log_state();
    
    for (auto [trie_char, next_ind] : node_vec[0].children()) {
        const CTrieNode * next_nd = &(node_vec[next_ind]);
        // std::cerr << "ROOT child for \"" << to_char_str(letters[trie_char]) <<  "\" "; next_nd->log_state();
    }
    
    for (unsigned int eli = 0; eli < equivalent_letter.size(); ++eli) {
        if (equivalent_letter[eli] == NO_MATCHING_CHAR_CODE) {
            std::cerr << to_char_str(letters[eli]) << " = <nothing>\n";
        } else {
            std::cerr << to_char_str(letters[eli]) << " = " << to_char_str(letters[equivalent_letter[eli]]) << "\n";
        }
    }
    auto nvs = sizeof(CTrieNode)*node_vec.size();
    auto suffs = concat_suff.size();
    std::cerr << "vecsize = " << nvs << " bytes\n";
    std::cerr << "concat_suff length = " << suffs << " bytes\n";
    std::cerr << "compressed tree size = " << 4*letters.size() + nvs + suffs << " bytes\n";
    std::cerr << "max_node_index = " << node_vec.size() << "\n";
    
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



void CompressedTrie::db_write_node(std::ostream & out, const CTrieNode & nd) const {
    if (nd.is_terminal()) {
        auto suff_index = nd.get_index();
        auto suff = get_suffix(suff_index);
        auto suff_str = to_char_str(suff);
        out << "TerminalNode suffix_ind=" << suff_index 
            << " suffix=" << to_char_str(suff)
            << " char_str=\"" << suff_str << "\"\n";
    } else {
        out << "InternalNode" << (nd.is_key_terminating() ? "* " : " ");
        out << "  offset = " << nd.get_index() << "\n";
        //out << "  letterbits = ";
        //nd.db_write_state(out);
        //out << "\n";
        for (auto [trie_char,index] : nd.children()) {
            out << "  " << to_char_str(letters[trie_char]);
            out << " => node[" << std::dec << index << "]\n";
        }
    }
}

void CompressedTrie::db_write(std::ostream & out) const {
    out << "CompressedTrie<with " << sizeof(CTrieNode) << " byte> nodes. Letters = \"" << to_char_str(letters) << "\"\n";
    out << "  " << node_vec.size() << " nodes:\n";
    std::size_t i = 0;
    for (auto nd : node_vec) {
        out << "node_vec[" << i++ << "] = ";
        db_write_node(out, nd); 

    }

}

void CompressedTrie::db_write_words(std::ostream & out) const {
    using nd_pref_pair = std::pair<const CTrieNode *, stored_str_t>;
    std::deque<nd_pref_pair> todo;
    stored_str_t mt;
    todo.push_back(nd_pref_pair{&(node_vec[0]), mt});
    std::size_t i = 0;
    while (!todo.empty()) {
        auto curr_nd_pref = todo.front();
        todo.pop_front();
        auto nd_ptr = curr_nd_pref.first;
        if (nd_ptr->is_terminal()) {
            auto suff_index = nd_ptr->get_index();
            auto suff = get_suffix(suff_index);
            auto full = curr_nd_pref.second + suff;
            out << i++ << " = " << to_char_str(full) << '\n';
        } else {
            vec_ind_pair_t vipt;
            for(auto x : nd_ptr->children())
                vipt.push_back(x);
            for (auto vipirit = vipt.rbegin(); vipirit != vipt.rend(); vipirit++) {
                const CTrieNode * nn = &(node_vec[vipirit->second]);
                stored_str_t np = curr_nd_pref.second + letters[vipirit->first];
                todo.push_front(nd_pref_pair{nn, np});
            }
        }
    }
}

} // namespace otc

// search impl in different file just to separate init from search.
