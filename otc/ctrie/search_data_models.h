#ifndef OTC_CTRIE_SEARCH_DATA_MODELS_H
#define OTC_CTRIE_SEARCH_DATA_MODELS_H

#include <set>
#include <vector>
#include <algorithm>

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
    
    FuzzyQueryResult(const stored_str_t & exact_match)
        :match_wide_char(exact_match),
        distance(0),
        score(1.0) {
    }
    

};

class FuzzyQueryResultWithTaxon {
    const FuzzyQueryResult query_result;
    const RTRichTaxNode * taxon = nullptr;
    const TaxonomyRecord * record = nullptr;
    bool matched_to_synonym;
    const std::string matched_name;
    float rescored  = -1;
    public:
    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const RTRichTaxNode * tax_arg,
                              const std::u32string * wide_query=nullptr
                              )
        :query_result(fqr),
        taxon(tax_arg),
        record(nullptr),
        matched_to_synonym(false),
        matched_name(tax_arg->get_data().get_nonuniqname()) {
        if (wide_query != nullptr) {
            do_rescore(*wide_query);
        }
    }

    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const TaxonomyRecord * tax_rec,
                              const std::u32string * wide_query=nullptr)
        :query_result(fqr),
        taxon(nullptr),
        record(tax_rec),
        matched_to_synonym(false),
        matched_name(tax_rec->name) {
        if (wide_query != nullptr) {
            do_rescore(*wide_query);
        }
    }

    FuzzyQueryResultWithTaxon(const FuzzyQueryResult & fqr,
                              const RTRichTaxNode * tax_arg,
                              const TaxonomicJuniorSynonym *syn,
                              const std::u32string * wide_query=nullptr)
        :query_result(fqr),
        taxon(tax_arg),
        record(nullptr),
        matched_to_synonym(true),
        matched_name(syn->get_name()) {
        if (wide_query != nullptr) {
            do_rescore(*wide_query);
        }
    }

    void do_rescore(const std::u32string &wide_query) {
        const auto lmn = lower_case_version(matched_name);
        auto wide_matched = to_u32string(lmn);
        auto dist = calc_damerau_levenshtein_dist(wide_query, wide_matched);
        float x;
        if (dist == 0) {
            x = 1.0;
        } else {
            const float denom = (float)wide_matched.length();
            x = (denom - (float)(dist))/denom;
        }
        set_rescore(x);
    }

    void set_rescore(float rescore_arg) {
        rescored = rescore_arg;
    }

    float get_score() const {
        return (rescored < -0.5 ? query_result.score : rescored);
    }

    bool is_synonym() const {
        return matched_to_synonym;
    }
    const std::string & get_matched_name() const {
        return matched_name;
    }

    const RTRichTaxNode * get_taxon() const {
        return taxon;
    }

    const TaxonomyRecord * get_record() const {
        return record;
    }

};

using vec_q_res_w_taxon = std::vector<FuzzyQueryResultWithTaxon>;

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

using sorted_q_res_set = std::set<FuzzyQueryResult, SortQueryResByNearness>;

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

template<typename T> class CompressedTrie;
template <typename T>
class PartialMatch {
    public:
    enum creation_modes {MATCH, DOWN, RIGHT, MULTIPLE};
    
    PartialMatch(const FQuery & q,
                 const T *nextn)
        :query(q),
         qpos(0),
         distance(0),
         next_node(nextn),
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::MATCH) {
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
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::MATCH) {
        match_coded.reserve(1 + prevpm.match_coded.capacity());
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
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::DOWN) {
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
         prev_mismatched_trie(NO_MATCHING_CHAR_CODE),
         create_mode(creation_modes::RIGHT) {
        match_coded.reserve(prevpm.match_coded.capacity());
        match_coded = prevpm.match_coded;
        match_coded.push_back(match_char);
        assert(nextn != prevpm.next_node);
    }
    void add_creation_mode(creation_modes ncm) {
        if (ncm != create_mode) {
            create_mode = creation_modes::MULTIPLE;
        }
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
    bool store_result(std::list<FuzzyQueryResult> & results,
                      const stored_index_t * trie_suff,
                      std::size_t suff_len,
                      unsigned int distance) const {
        query.store_result_ptrs(next_node, trie_suff);
        results.push_back(FuzzyQueryResult{get_prev_match_coded(), trie_suff, suff_len, distance});
        return true;
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
    
    const stored_str_t get_growing_match(const stored_str_t & letters) const {
        stored_str_t x;
        x.reserve(match_coded.size());
        for (auto i : match_coded) {
            if (i < letters.size()) {
                x.push_back(letters[i]);
            } else {
                x.push_back('#');
            }
        }
        return x;
    }

    stored_str_t get_growing_query() const {
        return query.wide_str.substr(0, qpos);
    }
    private:
    const FQuery & query;
    std::size_t qpos;
    unsigned int distance;
    const T * next_node;
    stored_index_t prev_mismatched_trie;
    std::vector<stored_index_t> match_coded;
    creation_modes create_mode;
    friend class CompressedTrie<T>;
};


} // namespace otc
#endif
