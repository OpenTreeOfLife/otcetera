#ifndef PRUNE_H
#define PRUNE_H

#include "otc/tree.h"
#include "otc/taxonomy/taxonomy.h"

namespace otc {

template <typename Tree>
std::pair<int,int> prune_unmapped_leaves(Tree& tree, const BaseTaxonomy& tax)
{
    int mapped_leaves = 0;
    int unmapped_leaves = 0;

    std::vector<typename Tree::node_type*> leaves;
    for(auto leaf: iter_leaf(tree))
    {
        if (leaf->has_ott_id())
        {
            auto id1 = leaf->get_ott_id();
            auto id2 = tax.get_unforwarded_id(id1);
            if (id2)
            {
                // Handle forwards
                if (*id2 != id1)
                    leaf->set_ott_id(*id2);

                // Count as mapped
                mapped_leaves++;
                continue;
            }
        }
        // Mark leaf for deletion
        leaves.push_back(leaf);
        unmapped_leaves++;
    }
    for(auto leaf: leaves) {
        while (leaf and leaf->is_tip()) {
            auto parent = leaf->get_parent();
            if (parent) {
                leaf->detach_this_node();
                delete leaf;
                leaf = parent;
            } else {
                delete leaf;
                tree._set_root(nullptr);
            }
        }
        assert(tree.get_root());
    }
    return {mapped_leaves, unmapped_leaves};
}

template <typename Tree>
void prune_duplicate_ottids(Tree& tree) {
    std::vector<typename Tree::node_type*> leaves;
    for(auto leaf: iter_leaf(tree)) {
        leaves.push_back(leaf);
    }
    std::map<OttId, typename Tree::node_type*> node_ptrs;
    for(auto leaf: leaves) {
        if (not leaf->has_ott_id()) {
            continue;
        }
        auto id = leaf->get_ott_id();
        // If the OTT id is new, then add the node as canonical representative of the OTT id
        if (not node_ptrs.count(id)) {
            node_ptrs.insert({id, leaf});
        } else {
            // Otherwise delete the non-canonical OTT id and its ancestors
            delete_tip_and_monotypic_ancestors(tree, leaf);
        }
    }
}

}

#endif
