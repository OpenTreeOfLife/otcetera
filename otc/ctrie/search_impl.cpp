#include "otc/ctrie/ctrie.h"

using std::vector;
using std::string;

namespace otc {

// Dynamic programming matrix.
//   This implementation doesn't separate gap-opening from gap-extension.
class dp_matrix
{
    int* data;
    int x_width;
    int y_width;
    int* xmin;
    int* xmax;
public:
    int& operator()(int x,int y)       {assert(0 <= x and x < x_width); assert(0 <= y and y < y_width) ; return data[y*x_width + x];}
    int  operator()(int x,int y) const {assert(0 <= x and x < x_width); assert(0 <= y and y < y_width) ; return data[y*x_width + x];}

    int size1() const {return x_width;}
    int size2() const {return y_width;}

    int calc_row(int y, stored_index_t target_char, const vector<stored_index_t>& query)
    {
        int best = INT_MAX/2;
        for(int x=xmin[y]; x<=xmax[y]; x++)
        {
            // We have just matched the x-th char, which is query[x-1]
            int match_cost       = (target_char == query[x-1]) ? 0 : 1;
            // FIXME: if we have the PREVIOUS target_char, we could look at letter transpositions here.

            int match_score      = (*this)(x-1,y-1) + match_cost;
            int del_target_score = (*this)(x-1,y  ) + 1;
            int del_query_score  = (*this)(x  ,y-1) + 1;
            (*this)(x,y) = std::min(match_score,std::min(del_query_score, del_target_score));
            best = std::min((*this)(x,y), best);
        }
        return best;
    }

    int score_for_row(int y) const
    {
        return (*this)(x_width-1, y);
    }

    dp_matrix(int xw, int yw, int max_dist)
        :data(new int[xw*yw]),
         x_width(xw),
         y_width(yw),
         xmin(new int[yw]),
         xmax(new int[yw])
    {
        // Initialize cells to INT_MAX/2 so that they don't affect the calculation.
        for(int x=0;x<x_width;x++)
            for(int y=0;y<y_width;y++)
                (*this)(x,y) = INT_MAX/2;

        // Initialize the first row and first column.
        // This avoids handling special cases later.
        for(int x=0;x<x_width;x++)
            (*this)(x,0) = x;
        for(int y=0;y<y_width;y++)
            (*this)(0,y) = y;

        // The score cannot be higher than the score you would get if all letters matched.
        // And that score is |x-y| for cell (x,y).
        // So we can impose the constraint that |x-y| <= max_dist.
        // This means that x \in (y-max_dist,y+max_dist) \cap (0,x_width-1)
        // Since we are computed (x=0,y) previously, we set xmin to at least 1
        for(int y=0;y<y_width;y++)
        {
            xmin[y] = std::max(1,y-max_dist);
            xmax[y] = std::min(x_width-1,y+max_dist);
        }
    };

    ~dp_matrix()
     {
         delete[] data;
         delete[] xmin;
         delete[] xmax;
     }
};

void CompressedTrie::extend_partial_match(const vector<stored_index_t>& query,
                                          const unsigned int max_dist,
                                          const CTrieNode* curr_node,
                                          dp_matrix& score,
                                          vector<stored_index_t>& match_coded,
                                          std::vector<FuzzyQueryResult> & results) const
{
    // 2. Handle case where there only one target string with this prefix
    if (curr_node->is_terminal())
    {
        auto prev_length = match_coded.size();
        auto suffix_length = get_suffix_length(*curr_node);
        auto total_length = prev_length + suffix_length;

        // If the terminal suffix is too long, we can't match.
        if (total_length >= score.size2()) return;

        auto suffix_char_ptr = get_suffix_ptr(*curr_node);
        int next_best = 0;
        for(int y = prev_length+1; y <= total_length; y++, suffix_char_ptr++)
        {
            if (next_best > max_dist) return;
            next_best = score.calc_row(y, *suffix_char_ptr, query);
        }

        unsigned int dist = score.score_for_row(total_length);

        // FIXME: maybe stop passing in L just to do a reserve?
        if (dist <= max_dist)
            results.push_back( {match_coded, get_suffix_ptr(*curr_node), suffix_length, dist} );

        return;
    }

    // 3. Handle case where there are multiple target strings with this prefix.
    for (auto [letter, index] : curr_node->children())
    {
        match_coded.push_back(letter);
        {

            int best = score.calc_row(match_coded.size(), letter, query);

            const CTrieNode * next_nd = &(node_vec[index]);

            // Don't search branches of the tree that can never achieve dist <= max_dist
            if (best <= max_dist)
                extend_partial_match(query, max_dist, next_nd, score, match_coded, results);

        }
        match_coded.pop_back();
    }

    // 1. Handle case where one target string terminates here, but this is a prefix of other target strings.
    if (curr_node->is_key_terminating())
    {
        unsigned int dist = score.score_for_row( match_coded.size() );

        if (dist <= max_dist)
            results.push_back( {match_coded, nullptr, 0, dist} );
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


std::vector<FuzzyQueryResult> CompressedTrie::fuzzy_matches(const stored_str_t & query_str, unsigned int max_dist) const
{
    if (DB_FUZZY_MATCH) {std::cerr << "fuzzy_matches (within " << max_dist << " edits) of \"" << to_char_str(query_str) << "\"\n";}
    if (query_str.length() == 0) {
        return std::vector<FuzzyQueryResult>{};
    }

    // 1. Trivial case: unmatchable query.
    auto query = encode_as_indices(query_str);
    unsigned int num_missing_in_letters = 0;
    for (auto qai : query) {
        if (qai == NO_MATCHING_CHAR_CODE) {
            num_missing_in_letters++;
            if (num_missing_in_letters > max_dist) {
                if (DB_FUZZY_MATCH) {std::cerr << "match infeasible because >= " << num_missing_in_letters << " positions in the query were not in the trie.\n";}
                return {};
            }
        }
    }

    // non-trivial case
    std::vector<FuzzyQueryResult> results;
    results.reserve(20);

    // 2. Allocate a dynamic programming matrix with size (|query|+1, |target|+1)
    //    We can ignore target strings longer than query_str().size+max_dist.
    const int XW = query.size()+1;
    const int YW = query.size()+1+max_dist;

    //    The cell (x,y) indicates the score after having seen x letters of the query and y letters of the target.
    dp_matrix score(XW, YW, max_dist);

    // 4. Keep track of the path through the prefix ctrie as we walk it.
    vector<stored_index_t> match_coded;
    match_coded.reserve(YW-1);

    auto root_node = &(node_vec.at(0));

    // Do a depth-first search using the stack.
    extend_partial_match(query, max_dist, root_node, score, match_coded, results);

    for (auto & r : results)
        _finish_query_result(r);

    return results;
}

void CompressedTrie::all_descendants(stored_str_t& prefix, uint64_t index, vector<string>& results) const
{
    auto& node = node_vec[index];

    if (node.is_key_terminating())
        results.push_back(to_char_str(prefix));

    if (node.is_terminal()) {
        auto suffix_index = node.get_index();
        auto suffix = get_suffix(suffix_index);
        results.push_back(to_char_str(prefix + suffix));
    }
    else
    {
        for (auto [letter, next_index] : node.children())
        {
            prefix.push_back(letters[letter]);
            all_descendants(prefix, next_index, results);
            prefix.pop_back();
        }
    }
}

// can we have both is_terminal() and key_terminating() on the same node?
// if we could have them both set and have an empty suffix, then we would match the same string twice.
// therefore perhaps, we can have them both set, but we only have is_terminal() set if there is a non-empty suffix.
vector<string> CompressedTrie::prefix_query(const stored_str_t& uquery) const
{
    if (node_vec.empty()) return {};

    auto query_letters = encode_as_indices(uquery);

    std::size_t index = 0;
    int letters_matched = 0;
    for(int i=0;i<query_letters.size() and not node_vec[index].is_terminal();i++)
    {
        auto letter = query_letters[i];
        auto next_index = node_vec[index].child_index_for_letter(letter);
        if (next_index)
            index = *next_index;
        else
            return {};

        letters_matched = i+1;
    }

    // If we have not matched all the prefix letters, check the suffix
    if (letters_matched<uquery.size())
    {
        assert(node_vec[index].is_terminal());
        auto suffix_index = node_vec[index].get_index();
        auto suffix = get_suffix(suffix_index);

        // We can't match if there aren't enough letters in the suffix
        if (letters_matched + suffix.size() < uquery.size())
            return {};

        // Check that the remaining query letters match the suffix
        for(int i=0;i<uquery.size() - letters_matched;i++)
            if (suffix[i] != uquery[letters_matched+i])
                return {};

        // Compute the matched string.
        auto matched_string = uquery.substr(0,letters_matched)+suffix;

        // Return a vector with 1 entry.
        return {to_char_str(matched_string)};
    }
    // Otherwise find all descendants of the node we have found.
    else
    {
        vector<string> results;
        auto prefix = uquery;
        all_descendants(prefix, index, results);

        return results;
    }
}


} // namespace otc

