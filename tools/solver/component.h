#ifndef COMPONENT_H
#define COMPONENT_H

#include <list>
#include <vector>
#include <memory>  // for shared_ptr

struct Solution;


// A non-trivial component
struct component_t
{
    std::list<int> elements;

    std::shared_ptr<Solution> solution;

    std::vector<ConstRSplit> new_splits;
    std::vector<std::shared_ptr<Solution>> old_solutions;

    std::vector<int> get_taxa(const std::vector<int>& other_taxa) const
    {
        std::vector<int> taxa;
        taxa.reserve(elements.size());
        for(auto index: elements)
            taxa.push_back(other_taxa[index]);
        return taxa;
    }
};

typedef component_t* component_ref;

#endif
