#include "build.h"
#include "rollback.h"

using std::vector;
using std::list;
using std::optional;
using std::shared_ptr;

vector<int> indices;

template <typename T>
bool sort_cmp(const vector<T>& v1,  const vector<T>& v2)
{
    if (v1.size() != v2.size()) return false;

    auto w1 = v1; std::sort(w1.begin(), w1.end());
    auto w2 = v2; std::sort(w2.begin(), w2.end());

    return (w1 == w2);
}

void Solution::initialize_taxon_index_map() const
{
    for(int k=0;k<indices.size();k++)
        assert(indices[k] == -1);
    for (int i=0;i<taxa.size();i++)
        indices[taxa[i]] = i;
}

void Solution::clear_taxon_index_map() const
{
    for(int id: taxa)
        indices[id] = -1;
}


bool exclude_group_intersects_component(const ConstRSplit& split, const component_t* component, const vector<component_ref>& component_for_index)
{
    for(int taxon: split->out)
    {
        int index = indices[taxon];

        if (index == -1) continue;

        if (component_for_index[index] == component) return true;
    }
    return false;
}

bool exclude_group_intersects_taxon_set(const ConstRSplit& split)
{
    for(int taxon: split->out)
    {
        int index = indices[taxon];

        if (index != -1) return true;
    }
    return false;
}


template <typename T>
void append(vector<T>& v1, const vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
}

/// Merge components c1 and c2 and return the component name that survived
component_ref merge_components(component_ref c1, component_ref c2, vector<component_ref>& component_for_index, vector<MergeRollbackInfo>& merge_rollback_info, bool record_component_mergers)
{
    if (c2->elements.size() > c1->elements.size())
        std::swap(c1, c2);

    for(int i: c2->elements)
        component_for_index[i] = c1;

    if (c1->solution)
    {
        // We should be able to revert the merge by doing {c1->solution = c1->old_solutions[0]; c1->old_solutions.clear();}
        assert(c1->old_solutions.empty());
        c1->old_solutions.push_back(c1->solution);
    }

    if (c2->solution)
        c1->old_solutions.push_back(c2->solution);

    if (record_component_mergers)
        merge_rollback_info.push_back({c1, c2, c2->elements.begin(), c1->solution});

    // One of these components could be new -- that is, composed only of previously-trivial components.
    append(c1->old_solutions, c2->old_solutions);

    c1->solution = {};

    c1->elements.splice(c1->elements.end(), c2->elements);

    return c1;
}

/// Merge components c1 and c2 and return the component name that survived
void merge_component_with_trivial(component_ref c1, int index2, vector<component_ref>& component_for_index, vector<MergeRollbackInfo>& merge_rollback_info, bool record_component_mergers)
{
    component_for_index[index2] = c1;

    if (c1->solution)
    {
        // We should be able to revert the merge by doing {c1->solution = c1->old_solutions[0]; c1->old_solutions.clear();}
        assert(c1->old_solutions.empty());
        c1->old_solutions.push_back(c1->solution);
    }

    if (record_component_mergers)
        merge_rollback_info.push_back({c1, nullptr, {}, c1->solution});

    c1->solution = {};

    c1->elements.push_back(index2);
}

void MaybeReuseSolution(shared_ptr<Solution>& solution, vector<shared_ptr<Solution>>& sub_solutions)
{
    // 0. If we found a solution to THIS exact problem then we can just re-use it.
    if (sub_solutions.size() == 1 and sub_solutions[0]->taxa.size() == solution->taxa.size())
    {
        auto prev_solution = sub_solutions[0];
        assert(sort_cmp(prev_solution->taxa, solution->taxa));

        // It only makes sense to switch to an old solution if we don't already have an old solution,
        // so check that that is the case.  This problem should have a new/empty solution object.
        assert(solution->non_implied_splits_from_components().empty());
        assert(solution->implied_splits.empty());

        // Move the components from the previous solution to this one!
        solution = prev_solution;
        sub_solutions.clear();

        // We are not done yet: we may need to add the `new_splits` to the (partial) solution that we just found.
        // So do NOT return yet.
    }
}



// TODO: If BUILD fails, can we rebuild the solution that we have modified in-place?
//       * we need to avoid modifying the old solutions (for a merged component) in-place.
//       * we need to restore the component->elements list.

/// Construct a tree with all the splits mentioned, and return false if this is not possible
///   You can get the resulting tree from it with solution.get_tree().
///   New splits are in both `new_splits` and `sub_solution`.
void RemoveImpliedSplits(const shared_ptr<Solution>& solution, vector<ConstRSplit>& new_splits, vector<shared_ptr<Solution>>& sub_solutions)
{
#pragma clang diagnostic ignored  "-Wsign-conversion"
#pragma clang diagnostic ignored  "-Wsign-compare"
#pragma clang diagnostic ignored  "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored  "-Wsign-compare"

    // 1. If there are no splits to add, then we are consistent.
    if (new_splits.empty() and sub_solutions.empty()) return;

    // 2. Initialize the mapping from taxa to indices.
    solution->initialize_taxon_index_map();

    // 3. Determine the new splits that go into each component (both satisfied AND unsatisfied)
    for(int k = new_splits.size()-1; k >= 0; k--)
    {
        auto& split = new_splits[k];
        bool implied = not exclude_group_intersects_taxon_set(split);
        if (implied)
        {
            // Copy the split to the implied_splits set.
            solution->implied_splits.push_back(split);

            // Remove it from the new_splits set.
            if (k < new_splits.size()-1)
                std::swap(split, new_splits.back());
            new_splits.pop_back();
        }
    }

    // 4. Check sub_solutions to see if they are punctured.
    for(int k = sub_solutions.size()-1; k >= 0; k--)
    {
        auto& sub_solution = sub_solutions[k];

        assert(solution != sub_solution);

        // I. Check if sub_solution is punctured.
        //    If so, then copy splits to solution->{implied,non_implied}_splits.
        bool punctured = false;
        for(int i=0; i < sub_solution->implied_splits.size(); i++)
        {
            auto& split = sub_solution->implied_splits[i];

            bool implied = not exclude_group_intersects_taxon_set(split);

            // If we just realized that this sub_solution is punctured, then copy the previously seen splits to new_splits
            if (implied and not punctured)
            {
                punctured = true;
                for(int j=0;j<i;j++)
                {
                    auto& split_prev = sub_solution->implied_splits[j];
                    new_splits.push_back(split_prev);
                }
            }

            // Copy the split to {implied,non_implied}_splits if the sub_solution is punctured.
            if (punctured)
            {
                if (implied)
                    solution->implied_splits.push_back(split);
                else
                    new_splits.push_back(split);
            }
        }

        // II. Replace punctured sub-solutions with their sub-component colutions.
        if (punctured)
        {
            // IIa. Add the sub-component solutions of the top-level sub_solution.
            for(auto& fragment: sub_solution->components)
                sub_solutions.push_back(fragment->solution);

            // IIb. Remove the punctured sub-solution.
            if (k != sub_solutions.size()-1)
                std::swap(sub_solutions[k], sub_solutions.back());
            sub_solutions.pop_back();
        }
    }

    // 5. Determine the new splits that go into each component (both satisfied AND unsatisfied)
    solution->clear_taxon_index_map();
}


void Merge(shared_ptr<Solution>& solution, const vector<ConstRSplit>& new_splits, vector<shared_ptr<Solution>>& sub_solutions, SolutionRollbackInfo& rollback_info)
{
    auto& component_for_index = solution->component_for_index;
    auto& components = solution->components;

    rollback_info.n_orig_components = solution->components.size();

    bool has_initial_components = not solution->components.empty();

    auto merge = [&](auto& group)
        {
            component_ref split_comp = nullptr;
            for(int taxon: group)
            {
                int index = indices[taxon];
                assert(index != -1);
                auto taxon_comp = component_for_index[index];
                if (not split_comp)
                {
                    if (not taxon_comp)
                    {
                        components.push_back(std::make_unique<component_t>());
                        taxon_comp = components.back().get();
                        merge_component_with_trivial(taxon_comp, index, component_for_index, rollback_info.merge_rollback_info, has_initial_components);
                    }
                    split_comp = taxon_comp;
                }
                else if (not taxon_comp)
                    merge_component_with_trivial(split_comp, index, component_for_index, rollback_info.merge_rollback_info, has_initial_components);
                else if (split_comp != taxon_comp)
                    split_comp = merge_components(split_comp,taxon_comp,component_for_index, rollback_info.merge_rollback_info, has_initial_components);
            }
        };

    // 3a. For each new split, all the leaves in the include group must be in the same component
    for(const auto& split: new_splits)
        merge(split->in);
    // 3b. For each sub_solution, all the leaves in the taxon set must be in the same component
    for(const auto& sub_solution: sub_solutions)
    {
        assert(sub_solution->taxa.size() < solution->taxa.size());
        merge(sub_solution->taxa);
    }

    // 4. Pack the components
    rollback_info.old_components = vector<shared_ptr<component_t>>();
    auto& packed_components = *rollback_info.old_components;
    assert(packed_components.empty());
    for(auto& component: components)
        if (not component->elements.empty())
            packed_components.push_back( component );
    std::swap(components, packed_components);

}

bool MaybeFail(shared_ptr<Solution>& solution)
{
    auto& components = solution->components;

    // SUCCESS!
    if (not solution->all_taxa_in_one_component())
        return false;

    // FAILURE!
    else
    {
        assert(components.size() == 1);
        solution->clear_taxon_index_map();

        // we failed!
        return true;
    }
}

void Assign(shared_ptr<Solution>& solution, vector<ConstRSplit>& new_splits, vector<shared_ptr<Solution>>& sub_solutions)
{
    auto& component_for_index = solution->component_for_index;

    // 1. Determine the new splits that go into each component.
    //    We will check if they are implied or unimplied when we call RemoveImpliedSplits( ) on the component.
    for(auto& split: new_splits)
    {
        int first = indices[*split->in.begin()];
        assert(first >= 0);
        auto component = component_for_index[first];

        component->new_splits.push_back(split);
    }

    // 2. Pass down sub_solutions into the correct component.
    //    They basically are bundles of splits to work on.
    //    All splits in the same bundle always go into the same component because we merged
    //       any intersecting components in 5b.
    //    We will check if they are punctured when we call RemoveImpliedSplits( ) on the component.
    for(auto& sub_solution: sub_solutions)
    {
        int first_taxon = sub_solution->taxa[0];
        int first_index = indices[first_taxon];
        auto component = component_for_index[first_index];
        component->old_solutions.push_back(sub_solution);
    }
}



bool BuildIncA(shared_ptr<Solution>& solution, vector<ConstRSplit>& new_splits, vector<shared_ptr<Solution>>& sub_solutions,
               vector<SolutionRollbackInfo>& all_rollback_info, bool top = false);

bool SolveSubproblems(shared_ptr<Solution>& solution, vector<SolutionRollbackInfo>& all_rollback_info)
{
    auto& taxa = solution->taxa;
    auto& components = solution->components;

    optional<int> failing_component;
    for(int i=0;i<components.size();i++)
    {
        auto& component = components[i];
        assert(component->elements.size() >= 2);

        vector<ConstRSplit> comp_new_splits;
        std::swap(component->new_splits, comp_new_splits);

        vector<shared_ptr<Solution>> comp_sub_solutions;
        std::swap(component->old_solutions, comp_sub_solutions);

        // If a previous component failed, we just want to clean up component->new_splits and component->old_solutions.
        if (failing_component) continue;

        // If we've invalidated the solution for this component because the component's taxon set increased,
        // then create an empty solution to use here.
        if (not component->solution)
            component->solution = std::make_shared<Solution>(*component, taxa);

        if (not BuildIncA(component->solution, comp_new_splits, comp_sub_solutions, all_rollback_info))
            failing_component = i;

        assert(component->old_solutions.empty());
        assert(component->solution);
    }

    return (not failing_component);
}

bool BuildIncA(shared_ptr<Solution>& solution, vector<ConstRSplit>& new_splits, vector<shared_ptr<Solution>>& sub_solutions,
               vector<SolutionRollbackInfo>& all_rollback_info, bool top)
{
    // 1. MaybeReuseSolution
    MaybeReuseSolution(solution, sub_solutions);

    // 2. Check if the solution is new.
    bool solution_is_new = (solution->visited == 0);
    solution->visited++;

    SolutionRollbackInfo sol_rollback_info(solution);

    // 3. Remove implied splits
    if (not top)
        RemoveImpliedSplits(solution, new_splits, sub_solutions);

    // A. If there are no splits to add, then we are consistent.
    if (new_splits.empty() and sub_solutions.empty())
    {
        if (not solution_is_new)
            all_rollback_info.push_back(sol_rollback_info);

        return true;
    }

    // B. Initialize the mapping from taxa to indices.
    solution->initialize_taxon_index_map();

    // 4. Merge components
    Merge(solution, new_splits, sub_solutions, sol_rollback_info);

    if (not solution_is_new)
        all_rollback_info.push_back(sol_rollback_info);

    // 5. Fail if there is only one component
    bool fail = MaybeFail(solution);
    if (fail)
        return false;

    // 6. Assign splits and sub_solutions to components
    Assign(solution, new_splits, sub_solutions);

    // C. Clear the taxon index map
    solution->clear_taxon_index_map();

    // 7. Recurse into sub-problems
    bool success = SolveSubproblems(solution, all_rollback_info);

    return success;
}


bool BUILDINC(shared_ptr<Solution>& solution, const vector<ConstRSplit>& new_splits)
{
    auto new_splits2 = new_splits;
    vector<shared_ptr<Solution>> sub_solutions;

    vector<SolutionRollbackInfo> all_rollback_info;
    bool ok =  BuildIncA(solution, new_splits2, sub_solutions, all_rollback_info, true);

    if (not ok)
    {
        for(int i = (int)all_rollback_info.size()-1; i>=0; i--)
            all_rollback_info[i].rollback();
    }

    return ok;
}

bool BUILD_check(const std::vector<int> all_leaves_indices, const std::vector<ConstRSplit>& splits)
{
    auto solution = std::make_shared<Solution>(all_leaves_indices);
    return BUILDINC(solution, splits);
}

std::unique_ptr<Tree_t> BUILD(const std::vector<int> all_leaves_indices, const std::vector<ConstRSplit>& splits)
{
    auto solution = std::make_shared<Solution>(all_leaves_indices);
    bool compatible = BUILDINC(solution, splits);
    if (compatible)
        return solution->get_tree();
    else
        return {};
}
