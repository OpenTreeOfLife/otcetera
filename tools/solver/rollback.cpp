#include "rollback.h"

using std::vector;
using std::list;
using std::shared_ptr;
using std::optional;

void MergeRollbackInfo::unmerge(Solution& S)
{
    // Undo merging two components.
    if (c2)
    {
        c2->elements.splice(c2->elements.end(), c1->elements, it, c1->elements.end());
        for(auto& x: c2->elements)
            S.component_for_index[x] = c2;
    }
    // Undo merge with trivial
    else
    {
        S.component_for_index[c1->elements.back()] = nullptr;
        c1->elements.pop_back();
    }

    if (old_solution)
        c1->solution = old_solution;

    c1->old_solutions.clear();
    assert(c1->new_splits.empty());
}

void SolutionRollbackInfo::rollback()
{
    assert(n_old_implied_splits <= S->implied_splits.size());
    S->implied_splits.resize(n_old_implied_splits);

    // If n_orig_component is set, then we didn't quit early before the Merge( ) step.
    if (n_orig_components and *n_orig_components == 0)
    {
        S->components.clear();

        // If we are just going to _delete_ S later on, then this is a waste of time.
        for(auto& c: S->component_for_index)
            c = nullptr;

        return;
    }
    else
        for(int i=(int)merge_rollback_info.size()-1; i >= 0; i--)
            merge_rollback_info[i].unmerge(*S);

    // NOTE: Some components are created during merging that are
    // (i) are not original components, and also
    // (ii) end up being empty.  So they are not final components.

    // * We need these components to survive (i.e not be destructed)
    //   so that we can temporarily add elements to them during rollback.

    // * We will then move these elements OUT of them into original
    //   components.

    // * These temporary components should all end up at the end of the
    //   original components vector<> (the non-packed vector<>) because they
    //   were added during merging by push_back( ).

    // * We record the original size of the components vector BEFORE merging
    //   and truncate to that length.  But before we truncate, we check that
    //   rolling back merges leaves all the dropped components with no elements.
    if (old_components)
    {
        assert(n_orig_components);
        std::swap(S->components, *old_components);

        for(int i=0;i<S->components.size();i++)
            assert(S->components[i]);

        for(int i=*n_orig_components;i<S->components.size();i++)
            assert(S->components[i]->elements.empty());

        assert(*n_orig_components <= S->components.size());
        S->components.resize(*n_orig_components);

        for(int i=0;i<S->components.size();i++)
            assert(not S->components[i]->elements.empty());
    }
}


