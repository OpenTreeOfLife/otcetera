#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
#include <memory> // for shared_ptr and unique_ptr

#include "rsplit.h"
#include "component.h"
#include "tree.h"

struct Solution
{
    std::vector<int> taxa;

    std::vector<ConstRSplit> implied_splits;

    std::vector< component_ref > component_for_index;
    std::vector< std::shared_ptr<component_t> > components;

    // Counter for determining if this is a new Solution or a pre-existing one.
    int visited = 0;

    bool all_taxa_in_one_component() const;

    void initialize_taxon_index_map() const;
    void clear_taxon_index_map() const;

    std::vector<ConstRSplit> non_implied_splits_from_components() const;
    std::vector<ConstRSplit> splits_from_components() const;

    int n_splits_from_components() const;

    std::unique_ptr<Tree_t> get_tree() const;

    bool valid() const;

    Solution& operator=(const Solution&) = default;
    Solution& operator=(Solution&&) = default;

    Solution(const Solution&) = default;
    Solution(Solution&&) = default;

    Solution(const std::vector<int>& t);
    Solution(const component_t& c, const std::vector<int> other_taxa);
};

#endif
