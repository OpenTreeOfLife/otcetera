#ifndef OTCETERA_CONFLICT_H
#define OTCETERA_CONFLICT_H
#include "otc/induced_tree.h"
#include "otc/tree_operations.h"
#include <unordered_map>
#include <vector>

template <typename node_t>
int n_include_tips(const node_t* node) {
    return node->get_data().n_include_tips;
}

template <typename node_t>
int n_tips(const node_t* node) {
    return node->get_data().n_tips;
}

template <typename node_t>
int& n_include_tips(node_t* node) {
    return node->get_data().n_include_tips;
}

template <typename node_t>
int& n_tips(node_t* node) {
    return node->get_data().n_tips;
}

// Compute number of tips at this node or below it
template <typename Tree_t>
void compute_tips(Tree_t& tree) {
    // Iterate over nodes, leaves first
    for(auto nd: iter_post(tree)) {
        if (nd->isTip()) {
            nd->get_data().n_tips = 1;
        }
        auto p = nd->getParent();
        if (p) {
            p->get_data().n_tips += nd->get_data().n_tips;
        }
    }
}

template <typename N>
std::vector<N*> leaf_nodes_below(N* node) {
    std::vector<N*> nodes;
    for(auto nd: iter_leaf_n_const(*node)) {
        nodes.push_back(nd);
    }
    return nodes;
}

template <typename N>
std::vector<N*> map_to_summary(const std::vector<const N*>& nodes) {
    std::vector<N*> nodes2(nodes.size(),nullptr);
    for(int i=0;i<(int)nodes2.size();i++) {
        nodes2[i] = summary_node(nodes[i]);
    }
    return nodes2;
}

template <typename N>
auto find_induced_nodes(const std::vector<N*>& leaves) {
    std::unordered_set<N*> nodes;
    N* MRCA = nullptr;
    for(auto leaf: leaves) {
        nodes.insert(leaf);
        MRCA = trace_find_mrca(MRCA, leaf, nodes);
    }
    std::vector<N*> vnodes;
    for(auto nd: nodes) {
        vnodes.push_back(nd);
    }
    return vnodes;
}

template <typename Tree_t>
using node_logger_t = std::function<void(const typename Tree_t::node_type* node2, const typename Tree_t::node_type* node1)>;

// Map nodes of T1 onto T2.
// * Monotypic nodes of T1 are ignored.
// * Leaves of T1 should correspond to leaves of T2, or nothing in T2. (?)
// * Each non-monotypic node of T1 will equal (supported_by/partial_path_of), conflict with, or be a terminal of nodes in T2.
// * A node of T1 can equal one node of T2 (supported_by) or several nodes of T2 (partial_path_of).
// * A node of T1 can conflict with several nodes of T2.
template <typename Tree_t>
void perform_conflict_analysis(const Tree_t& tree1,
                               const std::unordered_map<long, const typename Tree_t::node_type*>& ottid_to_node1,
                               const Tree_t& tree2,
                               const std::unordered_map<long, const typename Tree_t::node_type*>& ottid_to_node2,
                               node_logger_t<Tree_t> log_supported_by,
                               node_logger_t<Tree_t> log_partial_path_of,
                               node_logger_t<Tree_t> log_conflicts_with,
                               node_logger_t<Tree_t> log_resolved_by,
                               node_logger_t<Tree_t> log_terminal) {
    // Handle non-leaf correspondence?
    auto induced_tree1 = get_induced_tree<Tree_t>(tree1, ottid_to_node1, tree2, ottid_to_node2);
    auto induced_tree2 = get_induced_tree<Tree_t>(tree2, ottid_to_node2, tree1, ottid_to_node1);

    computeDepth(*induced_tree1);
    compute_tips(*induced_tree1);

    computeDepth(*induced_tree2);
    compute_tips(*induced_tree2);

    std::vector<typename Tree_t::node_type*> conflicts;

    auto map1 = get_ottid_to_node_map(*induced_tree1);
    auto map2 = get_ottid_to_node_map(*induced_tree2);
        
    // make summary_node field of induced_tree1 leaves point to leaves of induced_tree2.
    for(auto leaf: iter_leaf(*induced_tree1)) {
        auto leaf2 = map2.at(leaf->get_ott_id());
        summary_node(leaf) = leaf2;
    }
        
    auto L = countLeaves(*induced_tree1);
    assert(L == countLeaves(*induced_tree2));
        
    std::vector<typename Tree_t::node_type*> tree_nodes;
    for(auto nd: iter_post(*induced_tree1)) {
        tree_nodes.push_back(nd);
    }

    for(auto nd: tree_nodes) {
        if (not nd->getParent()) {
            continue;
        }
        // Ignore knuckles in input trees.
        // (Note that in general, if we've pruned this tree down to match the shared taxon set
        //  then this could produce knuckles that were not originally there.)
        if (nd->isOutDegreeOneNode()) {
            continue;
        }
        // If this node contains all tips under it, then it doesn't correspond to a split.
        if (nd->get_data().n_tips == (int)L) {
            continue;
        }
        // If this node is a tip, the mark the corresponding nodes
        if (nd->isTip()) {
            auto nd2 = summary_node(nd);
            log_terminal(nd2, nd);
            nd2 = nd2->getParent();
            for(; nd2 and nd2->isOutDegreeOneNode(); nd2 = nd2->getParent()) {
                log_terminal(nd2, nd);
            }
            continue;
        }

        // Find the list of nodes in the input tree that are below nd.
        auto leaves1 = leaf_nodes_below(const_cast<const typename Tree_t::node_type*>(nd));

        // Since nd is not a tip, and not monotypic, it should have at least 2 leaves below it.
        assert(leaves1.size() >= 2);

        int L2 = 0;
        for(auto nd: leaves1) {
            L2 += n_tips(nd);
        }

        // Find the corresponding list of nodes in the summary tree
        auto leaves2 = map_to_summary(leaves1);

        // Find the nodes in the induced tree of those nodes
        auto nodes = find_induced_nodes(leaves2);

        // Sort the nodes by depth to ensure all children occur before their parents -- unfortunately n*log(n) !
        std::sort(nodes.begin(), nodes.end(), [](auto x, auto y){return depth(x) > depth(y);});

        // The MRCA should be the last node in the vector.
        auto MRCA = nodes.back();

        // The n_include_tips for a parent node should count the n_include_tips for this node
        if (nodes.size() > 1) {
            for(std::size_t i=0;i<nodes.size()-1;i++) {
                auto nd = nodes[i];
                if (nd->isTip()) {
                    n_include_tips(nd) = n_tips(nd);
                }
                auto p = nd->getParent();
                assert(p);
                assert(nd != MRCA);
                n_include_tips(p) += n_include_tips(nd);
                assert(n_include_tips(nd) <= n_tips(nd));
           }
        }
            
        // If MRCA includes all and only the tips under nd, then MRCA is supporting or partial_path_of
        bool conflicts_or_resolved_by = n_include_tips(MRCA) < n_tips(MRCA);

        // Supported_by or partial_path_of
        if (not conflicts_or_resolved_by) {
            assert(MRCA->getParent());
            if (MRCA->getParent()->get_data().n_tips > MRCA->get_data().n_tips) {
                log_supported_by(MRCA, nd);
            } else {
                for(auto nd2 = MRCA; nd2 and nd2->get_data().n_tips == MRCA->get_data().n_tips; nd2 = nd2->getParent()) {
                    log_partial_path_of(nd2, nd);
                }
            }
        }

        conflicts.clear();
        for(auto nd: nodes) {
            // If we have (a) some, but not all of the include group
            //            (b) any of the exclude group
            if (n_include_tips(nd) < n_tips(nd) and n_include_tips(nd) < L2) {
                conflicts.push_back(nd);
            }
        }
        for(auto nd: nodes) {
            n_include_tips(nd) = 0;
        }
            
        for(auto conflicting_node: conflicts) {
            log_conflicts_with(conflicting_node, nd);
        }

        if (conflicts.empty() and conflicts_or_resolved_by) {
            log_resolved_by(MRCA, nd);
        }

#ifdef CHECK_MARKS
        for(const auto nd2: iter_post_const(*induced_tree2)){
            assert(n_include_tips(nd2) == 0);
        }
#endif

        // nd -> MRCA
        if (not conflicts_or_resolved_by) {
            summary_node(nd) = MRCA;
            destroy_children(nd);
            destroy_children(MRCA);
        }
    }
}

#endif
