#ifndef ORACLE_H
#define ORACLE_H

#include <vector>
#include <memory>
#include "tree.h"

void remove_conflicting_splits_from_tree(std::vector<std::unique_ptr<Tree_t>>& trees, int k);

template<typename N>
inline void collapse_node_(N* nd)
{
    assert(nd);
    assert(nd->has_children());

    while(nd->has_children())
    {
        auto child = nd->get_first_child();
        child->detach_this_node();
        nd->add_sib_on_left(child);
    }
    nd->detach_this_node();
}


#endif
