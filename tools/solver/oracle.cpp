#include "oracle.h"
#include "otc/conflict.h"

using namespace otc;

using std::set;
using std::unique_ptr;
using std::vector;

set<Tree_t::node_type*> find_conflicting_nodes(unique_ptr<Tree_t>& ok_tree, unique_ptr<Tree_t>& tree_to_clean)
{
    set<node_t*> conflicting_nodes;

    typedef Tree_t::node_type node_t;
    typedef ConflictTree::node_type cnode_t;
    std::function<node_t*(node_t*,node_t*)> mrca_of_pair = [](node_t* n1, node_t* n2) {return mrca_from_depth(n1,n2);};


    auto tree_to_clean_ottid_to_node = get_ottid_to_node_map(*tree_to_clean);

    auto ok_tree_ottid_to_node = get_ottid_to_node_map(*ok_tree);

    compute_depth(*ok_tree);

    compute_depth(*tree_to_clean);

    auto [induced_tree_to_clean,to_induced] = get_induced_tree_and_node_map<ConflictTree>(*tree_to_clean,
                                                                                          tree_to_clean_ottid_to_node,
                                                                                          mrca_of_pair,
                                                                                          *ok_tree,
                                                                                          ok_tree_ottid_to_node);

    auto induced_ok_tree = get_induced_tree<ConflictTree>(*ok_tree,
                                                          ok_tree_ottid_to_node,
                                                          mrca_of_pair,
                                                          *tree_to_clean,
                                                          tree_to_clean_ottid_to_node);

    if (induced_tree_to_clean->get_root() and induced_ok_tree->get_root())
    {
        std::unordered_map<const cnode_t*, node_t*> from_induced;
        for(auto& [non_induced,induced]: to_induced)
        {
            from_induced[induced] = non_induced;
        }

        auto log_supported_by    = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_partial_path_of = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_resolved_by     = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_terminal        = [&](const cnode_t* /* node2 */, const cnode_t* /* node1 */) {};
        auto log_conflicts_with  = [&](const cnode_t* /* node2 */, const cnode_t* node1 )
        {
            assert(node1->has_children());
            // for each of node1 and all its monotypic ancestors.
            do
            {
                // mark the corresponding node in the original (not induced) tree as conflicting.
                auto node2 = from_induced.at(node1);
                assert(node2->has_children());
                conflicting_nodes.insert(node2);

                node1 = node1->get_parent();
            }
            while(node1->is_outdegree_one_node());
        };

        perform_conflict_analysis(*induced_tree_to_clean, *induced_ok_tree, log_supported_by, log_partial_path_of, log_conflicts_with, log_resolved_by, log_terminal);
    }

    return conflicting_nodes;
}

void remove_conflicting_splits_from_tree(unique_ptr<Tree_t>& ok_tree, unique_ptr<Tree_t>& tree_to_clean)
{
    for(auto& node: find_conflicting_nodes(ok_tree, tree_to_clean))
        collapse_node_(node);
}

void remove_conflicting_splits_from_tree(vector<unique_ptr<Tree_t>>& trees, int k)
{
    for(int i=0;i<k;i++)
        remove_conflicting_splits_from_tree(trees[i],trees[k]);
}

