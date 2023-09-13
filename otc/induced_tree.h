#ifndef OTCETERA_INDUCED_TREE_H
#define OTCETERA_INDUCED_TREE_H
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/otc_base_includes.h"

namespace otc {

template <typename node_t>
node_t* trace_to_parent(node_t* node, std::unordered_set<node_t*>& nodes) {
    assert(nodes.count(node));
    node = node->get_parent(); // move to parent and insert it.
    nodes.insert(node);
    return node;
}

template <typename N, typename F>
N* MRCA_of_group(const std::vector<N*>& leaves,
                 const F& MRCA_of_pair) {
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

template <typename Tree_Out_t, typename Node_In_t>
std::pair<std::unique_ptr<Tree_Out_t>,std::unordered_map<Node_In_t*, non_const_node_type<Tree_Out_t>*>>
get_induced_tree_from_leaves_and_MRCA(const std::vector<Node_In_t*>& leaves,
                                      Node_In_t* MRCA)
{

    std::unique_ptr<Tree_Out_t> induced_tree(new Tree_Out_t());
    std::unordered_map<Node_In_t*, non_const_node_type<Tree_Out_t>*> to_induced_tree;

    // 0. If there are no leaves, return an empty tree.
    if (leaves.empty()) return {std::move(induced_tree), to_induced_tree};

    // 1. Find all nodes in the tree
    auto nodes = find_induced_nodes(leaves, MRCA);;
    // 2. Construct duplicate nodes for the induced tree, recording correspondence
    for(auto nd: nodes)
    {
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
    return {std::move(induced_tree), to_induced_tree};
}

template <typename Tree_Out_t, typename Node_In_t, typename F>
auto get_induced_tree_and_node_map(const std::vector<Node_In_t*>& leaves,
                                   const F& MRCA_of_pair)
{
    return get_induced_tree_from_leaves_and_MRCA<Tree_Out_t>(leaves, MRCA_of_group(leaves, MRCA_of_pair));
}

template <typename Tree_Out_t, typename Node_In_t, typename F>
std::unique_ptr<Tree_Out_t>
get_induced_tree(const std::vector<Node_In_t*>& leaves,
                 const F& MRCA_of_pair)
{
    return get_induced_tree_and_node_map<Tree_Out_t>(leaves, MRCA_of_pair).first;
}

// Get a list of leaves of tree 1 that are also in tree 2.
template <typename Tree1_t, typename Tree2_t=Tree1_t>
std::vector<node_type<Tree1_t>*> get_induced_leaves(
                    Tree1_t& T1,
                    const std::unordered_map<OttId, node_type<Tree1_t>*>& nodes1,
                    Tree2_t& T2,
                    const std::unordered_map<OttId, node_type<Tree2_t>*>& nodes2)
{
  std::vector<node_type<Tree1_t>*> leaves;
    if (nodes2.size() < nodes1.size()) {
        for(auto leaf: iter_leaf(T2)) {
            auto id = leaf->get_ott_id();
            auto it = nodes1.find(id);
            if (it != nodes1.end()) {
                leaves.push_back(it->second);
            }
        }
    } else {
        for(auto leaf: iter_leaf(T1)) {
            auto id = leaf->get_ott_id();
            if (nodes2.find(id) != nodes2.end()) {
                leaves.push_back(leaf);
            }
        }
    }
    return leaves;
}

// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree_Out_t, typename Tree_In1_t, typename Tree_In2_t, typename F>
auto
get_induced_tree_and_node_map(Tree_In1_t& T1,
                                   const std::unordered_map<OttId, node_type<Tree_In1_t>*>& nodes1,
                                   const F& MRCA_of_pair,
                                   Tree_In2_t& T2,
                                   const std::unordered_map<OttId, node_type<Tree_In2_t>*>& nodes2)
{
    auto induced_leaves = get_induced_leaves(T1, nodes1, T2, nodes2);

    auto tree_and_node_map = get_induced_tree_and_node_map<Tree_Out_t>(induced_leaves, MRCA_of_pair);
    auto& induced_tree = tree_and_node_map.first;
    induced_tree->set_name(T1.get_name());
    return tree_and_node_map;
}

// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree_Out_t, typename Tree_In1_t, typename Tree_In2_t, typename F>
std::unique_ptr<Tree_Out_t> get_induced_tree(Tree_In1_t& T1,
                                             const std::unordered_map<OttId, node_type<Tree_In1_t>*>& nodes1,
                                             const F& MRCA_of_pair,
                                             Tree_In2_t& T2,
                                             const std::unordered_map<OttId, node_type<Tree_In2_t>*>& nodes2)
{
    return get_induced_tree_and_node_map<Tree_Out_t,Tree_In1_t,Tree_In2_t>(T1, nodes1, MRCA_of_pair, T2, nodes2).first;
}

template <typename T>
class has_id_to_node {
private:
    typedef char YesType[1];
    typedef char NoType[2];

    template <typename C> static YesType& test(decltype(&C::id_to_node));
    template <typename C> static NoType& test(...);

public:
    enum { value = sizeof(test<T>(0)) == sizeof(YesType) };
};

// Get a list of nodes in T2 that are leaves in T1.
// The nodes in T2 do NOT need to be leaves of T2.

template <typename T>
using data_type = typename T::data_type;

template <typename Tree1_t, typename Tree2_t>
typename std::enable_if<has_id_to_node<data_type<Tree2_t>>::value, std::vector<node_type<Tree2_t>*>>::type
get_induced_nodes(const Tree1_t& T1, Tree2_t& T2)
{
    auto& ott_to_nodes2 = T2.get_data().id_to_node;
    std::vector<node_type<Tree2_t>*> nodes;
    for(auto leaf: iter_leaf(T1))
    {
        auto id = leaf->get_ott_id();
        auto it = ott_to_nodes2.find(id);
        if (it != ott_to_nodes2.end()) {
            nodes.push_back(const_cast<node_type<Tree2_t>*>(it->second));
        }
    }
    return nodes;
}


// Get a list of nodes in T2 that are leaves in T1.
// The nodes in T2 do NOT need to be leaves of T2.

template <typename Tree1_t, typename Tree2_t>
typename std::enable_if<not has_id_to_node<data_type<Tree2_t>>::value, std::vector<node_type<Tree2_t>*>>::type
get_induced_nodes(const Tree1_t& T1, Tree2_t& T2)
{
    auto ott_to_nodes2 = otc::get_ottid_to_const_node_map(T2);

    std::vector<node_type<Tree2_t>*> nodes;
    for(auto leaf: iter_leaf(T1))
    {
        auto id = leaf->get_ott_id();
        auto it = ott_to_nodes2.find(id);
        if (it != ott_to_nodes2.end()) {
            nodes.push_back(const_cast<node_type<Tree2_t>*>(it->second));
        }
    }
    return nodes;
}


}// namespace otc
#endif
