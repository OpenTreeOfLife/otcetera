#include "solution.h"

using std::vector;
using std::unique_ptr;

bool Solution::all_taxa_in_one_component() const
{
    return component_for_index[0] and component_for_index[0]->elements.size() == taxa.size();
}

bool Solution::valid() const
{
    for(auto& component: components)
    {
        if (not component->solution) return false;

        if (not component->solution->valid()) return false;
    }

    return true;
}

unique_ptr<Tree_t> Solution::get_tree() const
{
    assert(taxa.size() > 1);

    // 1. Make a tree with just a root node
    std::unique_ptr<Tree_t> tree(new Tree_t());
    tree->create_root();

    // 2. Add children for non-trivial components
    for(auto& component: components)
        add_subtree(tree->get_root(), *component->solution->get_tree());

    // 3. Add children for trivial components
    for(int index=0;index<taxa.size();index++)
    {
        if (not component_for_index[index])
        {
            auto taxon = taxa[index];
            auto node = tree->create_child(tree->get_root());
            node->set_ott_id(taxon);
        }
    }

    return tree;
}

Solution::Solution(const vector<int>& t)
    :taxa(t), component_for_index(taxa.size())
{}

Solution::Solution(const component_t& c, const vector<int> other_taxa)
    :Solution(c.get_taxa(other_taxa))
{}
