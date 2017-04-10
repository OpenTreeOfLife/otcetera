#ifndef OTCETERA_INDUCED_TREE_H
#define OTCETERA_INDUCED_TREE_H
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/otc_base_includes.h"

namespace otc {
template <typename N>
inline int depth(const N* node) {
    assert(node->get_data().depth > 0);
    return node->get_data().depth;
}


template <typename N>
inline int& depth(N* node) {
    assert(node->get_data().depth > 0);
    return node->get_data().depth;
}

template <typename node_t>
node_t* trace_to_parent(node_t* node, std::unordered_set<node_t*>& nodes) {
    assert(nodes.count(node));
    node = node->get_parent(); // move to parent and insert it.
    nodes.insert(node);
    return node;
}

template <typename N>
N* mrca_from_depth(N* node1, N* node2) {
    assert(node1 or node2);
    if (not node1) {
        return node2;
    }
    if (not node2) {
        return node1;
    }
    assert(node1 and node2);
    assert(get_root(node1) == get_root(node2));
    while (depth(node1) > depth(node2)) {
        node1 = node1->get_parent();
    }
    while (depth(node1) < depth(node2)) {
        node2 = node2->get_parent();
    }
    assert(depth(node1) == depth(node2));
    while (node1 != node2) {
        assert(node1->get_parent());
        assert(node2->get_parent());
        node1 = node1->get_parent();
        node2 = node2->get_parent();
    }
    assert(node1 == node2);
    return node1;
}


template <typename N>
N* MRCA_of_group(const std::vector<N*>& leaves,
                 std::function<N*(N*,N*)> MRCA_of_pair) {
    N* MRCA = nullptr;
    bool first = true;
    for(auto leaf: leaves) {
        if (first) {
            first = false;
            MRCA = leaf;
        }
        else {
            MRCA = MRCA_of_pair(MRCA, leaf);
        }
        if (not MRCA) {
            return nullptr;
        }
    }
    return MRCA;
}


template <typename N>
std::set<N*> find_induced_nodes(const std::vector<N*>& leaves, N* MRCA) {
    assert(MRCA);
    std::set<N*> nodes;
    nodes.insert(MRCA);
    for (auto leaf : leaves) {
        auto node = leaf;
        while (nodes.count(node) == 0) {
            nodes.insert(node);
            node = node->get_parent(); 
        } 
    }
    return nodes;
}

template <typename Tree_In_t, typename Tree_Out_t=Tree_In_t>
std::unique_ptr<Tree_Out_t>
get_induced_tree_from_leaves_and_MRCA(const std::vector<const typename Tree_In_t::node_type*>& leaves,
                                      const typename Tree_In_t::node_type* MRCA) {
    std::unique_ptr<Tree_Out_t> induced_tree(new Tree_Out_t());
    // 1. Find all nodes in the tree
    auto nodes = find_induced_nodes(leaves, MRCA);;
    // 2. Construct duplicate nodes for the induced tree, recording correspondence
    std::unordered_map<const typename Tree_In_t::node_type*, typename Tree_Out_t::node_type*> to_induced_tree;
    for(auto nd: nodes) {
        auto nd2 = induced_tree->create_node(nullptr);
        if (nd->has_ott_id()) {
            nd2->set_ott_id(nd->get_ott_id());
        }
        if (nd->get_name().size()) {
            nd2->set_name(nd->get_name());
        }
        to_induced_tree[nd] = nd2;
    }
    // 3. Link corresponding nodes to their corresponding parents
    for(auto nd: nodes) {
        auto p = nd->get_parent();
        auto nd2 = to_induced_tree.find(nd)->second;
        auto p2_it = to_induced_tree.find(p);
        if (p2_it == to_induced_tree.end()) {
            assert(nd == MRCA);
        } else {
            auto p2 = p2_it->second;
            p2->add_child(nd2);
        }
    }
    // 4. Set the root of the induced tree to node corresponding to the MRCA
    induced_tree->_set_root( to_induced_tree.at(MRCA) );
    return induced_tree;
}

template <typename Tree_In_t, typename Tree_Out_t=Tree_In_t>
std::unique_ptr<Tree_Out_t>
get_induced_tree(const std::vector<const typename Tree_In_t::node_type*>& leaves,
                 std::function<const typename Tree_In_t::node_type*(const typename Tree_In_t::node_type*,const typename Tree_In_t::node_type*)> MRCA_of_pair)
{
    return get_induced_tree_from_leaves_and_MRCA<Tree_In_t,Tree_Out_t>(leaves, MRCA_of_group(leaves, MRCA_of_pair));
}

// Get a list of leaves of tree 1 that are also in tree 2.
template <typename Tree1_t, typename Tree2_t>
std::vector<const typename Tree1_t::node_type*> get_induced_leaves(
                    const Tree1_t& T1,
                    const std::unordered_map<OttId, const typename Tree1_t::node_type*>& nodes1,
                    const Tree2_t& T2,
                    const std::unordered_map<OttId, const typename Tree2_t::node_type*>& nodes2) {
  std::vector<const typename Tree1_t::node_type*> leaves;
    if (nodes2.size() < nodes1.size()) {
        for(auto leaf: iter_leaf_const(T2)) {
            auto id = leaf->get_ott_id();
            auto it = nodes1.find(id);
            if (it != nodes1.end()) {
                leaves.push_back(it->second);
            }
        }
    } else {
        for(auto leaf: iter_leaf_const(T1)) {
            auto id = leaf->get_ott_id();
            if (nodes2.find(id) != nodes2.end()) {
                leaves.push_back(leaf);
            }
        }
    }
    return leaves;
}

// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree_In1_t, typename Tree_In2_t, typename Tree_Out_t>
std::unique_ptr<Tree_Out_t> get_induced_tree(const Tree_In1_t& T1,
                                             const std::unordered_map<OttId, const typename Tree_In1_t::node_type*>& nodes1,
                                             std::function<const typename Tree_In1_t::node_type*(const typename Tree_In1_t::node_type*,const typename Tree_In1_t::node_type*)> MRCA_of_pair,
                                             const Tree_In2_t& T2,
                                             const std::unordered_map<OttId, const typename Tree_In2_t::node_type*>& nodes2)
{
    auto induced_leaves = get_induced_leaves(T1, nodes1, T2, nodes2);
    auto induced_tree = get_induced_tree<Tree_In1_t, Tree_Out_t>(induced_leaves, MRCA_of_pair);
    induced_tree->set_name(T1.get_name());
    return induced_tree;
}


}// namespace otc
#endif
