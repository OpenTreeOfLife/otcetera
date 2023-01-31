#ifndef ROLLBACK_H
#define ROLLBACK_H

#include <vector>
#include <list>
#include <memory> // for shared_ptr
#include <optional>

#include "solution.h"
#include "component.h"

struct MergeRollbackInfo
{
    component_ref c1;
    component_ref c2;
    std::list<int>::const_iterator it;

    std::shared_ptr<Solution> old_solution;

    void unmerge(Solution& S);
};

struct SolutionRollbackInfo
{
    std::shared_ptr<Solution> S;
    int n_old_implied_splits;
    std::vector<MergeRollbackInfo> merge_rollback_info;
    std::optional<int> n_orig_components;
    std::optional<std::vector< std::shared_ptr<component_t> >> old_components;

    void rollback();

    explicit SolutionRollbackInfo(const std::shared_ptr<Solution>& s)
        : S(s), n_old_implied_splits( s->implied_splits.size() )
    {
    }
};


#endif
