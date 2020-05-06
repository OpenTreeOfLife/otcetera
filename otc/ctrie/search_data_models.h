#ifndef OTC_CTRIE_SEARCH_DATA_MODELS_H
#define OTC_CTRIE_SEARCH_DATA_MODELS_H

#include <set>
#include <vector>
#include <algorithm>
#include <optional>

#include "otc/otc_base_includes.h"
#include "otc/ctrie/str_utils.h"
#include "otc/ctrie/ctrie_node.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/util.h"

namespace otc {

class FuzzyQueryResult {
    public:
    std::string match() const {
        return to_char_str(match_wide_char);
    }

    stored_str_t match_wide_char;
    unsigned int distance;
    float score;
    std::vector<stored_index_t> match_coded;
    FuzzyQueryResult(const std::vector<stored_index_t> & mc,
                     const stored_index_t * trie_suff,
                     std::size_t suff_len,
                     unsigned int d)
      :distance(d),
      score(0.0), 
      match_coded() {
        match_coded.reserve(mc.size() + suff_len);
        match_coded = mc;
        auto ip = std::end(match_coded);
        if (suff_len > 0) {
            match_coded.insert(ip, trie_suff, trie_suff + suff_len);
        }
    }

};

class TaxonResult
{
    const RTRichTaxNode * taxon = nullptr;
    const TaxonomyRecord * record = nullptr;
    bool matched_to_synonym;
    const std::string matched_name;

public:

    bool is_synonym() const {
        return matched_to_synonym;
    }

    std::string get_matched_name() const {
        return matched_name;
    }

    const RTRichTaxNode * get_taxon() const {
        return taxon;
    }

    const TaxonomyRecord * get_record() const {
        return record;
    }

    TaxonResult(const RTRichTaxNode * tax_arg)
        :taxon(tax_arg),
         matched_to_synonym(false),
         matched_name(tax_arg->get_data().get_nonuniqname())
    { }

    TaxonResult(const TaxonomyRecord * tax_rec)
        :record(tax_rec),
         matched_to_synonym(false),
         matched_name(tax_rec->name)
    { }

    TaxonResult(const RTRichTaxNode * tax_arg,
                const TaxonomicJuniorSynonym *syn)
        :taxon(tax_arg),
         matched_to_synonym(true),
         matched_name(syn->get_name())
    { }
};
        

class FuzzyQueryResultWithTaxon: public TaxonResult
{
    const FuzzyQueryResult query_result;
public:
    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const RTRichTaxNode * tax_arg)
        :TaxonResult(tax_arg),
         query_result(fqr)
        { }

    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const TaxonomyRecord * tax_rec)
        :TaxonResult(tax_rec),
         query_result(fqr)
        { }

    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const RTRichTaxNode * tax_arg,
                              const TaxonomicJuniorSynonym *syn)
        :TaxonResult(tax_arg,syn),
         query_result(fqr)
        { }

    float get_score() const {
        return query_result.score;
    }
};

struct SortQueryResByNearness {
    bool operator() (const FuzzyQueryResult & lhs,
                     const FuzzyQueryResult & rhs) const {
        if (lhs.score < rhs.score) {
            return false;
        } else if (rhs.score < lhs.score) {
            return true;
        }
        return rhs.match_wide_char < rhs.match_wide_char;
    }
};

class FQuery {
    using ptr_pair = std::pair<const void *, const stored_index_t *>;
    public: 
    FQuery(const stored_str_t & q,
           const std::vector<stored_index_t> as_inds,
           unsigned int max_d)
        :wide_str(q),
        as_indices(as_inds),
        max_dist(max_d) {
        suffixes.resize(q.length() + 1);
    }
    const stored_str_t & wide_str;
    const std::vector<stored_index_t> as_indices;
    unsigned int max_dist;
    std::vector<stored_str_t> suffixes;

    mutable std::set<ptr_pair> already_matched; // threading issue if we parallelize traversal of trie on same query.
    
    void store_result_ptrs(const void * node_ptr, const stored_index_t * trie_suff) const {
        already_matched.insert(ptr_pair{node_ptr, trie_suff});
    }

    bool has_matched_suffix(const void * node_ptr, const stored_index_t * trie_suff) const {
        ptr_pair check{node_ptr, trie_suff};
        return contains(already_matched, check);
    }
};

class PartialMatch {
    public:
    enum creation_modes {MATCH, DOWN, RIGHT};
    
    PartialMatch(const FQuery & q,
                 const CTrieNode *nextn)
        :query(q),
         qpos(0),
         distance(0),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::MATCH)
        {
        }

    // create a partial match previous match and a char match
    PartialMatch(const PartialMatch & prevpm,
                 stored_index_t match_char,
                 unsigned int start_dist,
                 const CTrieNode *nextn,
                 bool was_match)
        :prev_match(&prevpm),
         match_letter(match_char),
         query(prevpm.query),
         qpos(prevpm.qpos + 1),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::MATCH)
        {
            if (not was_match)
                prev_mismatched_trie = match_char;
            assert(nextn != prevpm.next_node);
        }

    // create a partial match from a gap, moving through query but not trie
    PartialMatch(const PartialMatch & prevpm,
                 unsigned int start_dist,
                 const CTrieNode *nextn)
        :prev_match(&prevpm),
         query(prevpm.query),
         qpos(prevpm.qpos + 1),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::DOWN)
        {
            assert(nextn == prevpm.next_node);
        }

    // create a partial match from a gap, moving through trie but not query
    PartialMatch(const PartialMatch & prevpm,
                 unsigned int start_dist,
                 const CTrieNode *nextn, 
                 stored_index_t match_char)
        :prev_match(&prevpm),
         match_letter(match_char),
         query(prevpm.query),
         qpos(prevpm.qpos),
         distance(start_dist),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::RIGHT)
        {
            assert(nextn != prevpm.next_node);
        }

    bool can_downshift() const {
        return create_mode != creation_modes::RIGHT;
    }

    bool can_rightshift() const {
        return create_mode != creation_modes::DOWN;
    }

    bool has_matched_suffix(const stored_index_t * trie_suff) const {
        return query.has_matched_suffix(next_node, trie_suff);
    }
    bool store_result(std::vector<FuzzyQueryResult> & results,
                      const stored_index_t * trie_suff,
                      std::size_t suff_len,
                      unsigned int distance) const {
        query.store_result_ptrs(next_node, trie_suff);
        results.push_back(FuzzyQueryResult{get_prev_match_coded(), trie_suff, suff_len, distance});
        return true;
    }
    
    const CTrieNode * get_next_node() const {
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

    void make_match_coded(std::vector<stored_index_t>& s) const
    {
        if (match_letter)
            s.push_back(*match_letter);
        if (prev_match)
            prev_match->make_match_coded(s);
    }

    std::vector<stored_index_t> get_prev_match_coded() const
    {
        std::vector<stored_index_t> s;
        s.reserve(query.max_dist + query.as_indices.size());
        make_match_coded(s);
        std::reverse(s.begin(), s.end());
        return s;
    }

    std::optional<stored_index_t> prev_match_letter() const
    {
        return match_letter;
    }
private:
    const PartialMatch* prev_match = nullptr;
    std::optional<stored_index_t> match_letter;

    const FQuery & query;
    std::size_t qpos;
    unsigned int distance;
    const CTrieNode * next_node;
    stored_index_t prev_mismatched_trie;
    const creation_modes create_mode;
};


} // namespace otc
#endif
