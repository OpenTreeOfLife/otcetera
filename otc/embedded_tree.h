#ifndef OTCETERA_EMBEDDED_TREE_H
#define OTCETERA_EMBEDDED_TREE_H

#include <map>
#include <list>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_data.h"
#include "otc/node_embedding.h"

namespace otc {
class EmbeddedTree {
    protected:
    std::list<NodePairingWithSplits> nodePairings;
    std::list<PathPairingWithSplits> pathPairings;
    std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> scaffoldNdToNodeEmbedding;
    public:
    EmbeddedTree() {
    }
    void embed_new_tree(TreeMappedWithSplits & scaffold_tree,
                      TreeMappedWithSplits & tree,
                      std::size_t treeIndex);
    void embed_scaffold_clone(TreeMappedWithSplits & scaffold_tree,
                            TreeMappedWithSplits & tree,
                            std::size_t treeIndex);
    void write_dot_export(std::ostream & out,
                        const NodeEmbedding<NodeWithSplits, NodeWithSplits> & thr,
                        const NodeWithSplits * nd,
                        const std::vector<TreeMappedWithSplits *> &,
                        bool entireSubtree,
                        bool includeLastTree) const;
    // for testing...
    std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> &_get_scaffold_nd_to_node_embedding() {
        return scaffoldNdToNodeEmbedding;
    }
    protected:
    NodePairingWithSplits * _add_node_mapping(NodeWithSplits *taxo,
                                            NodeWithSplits *nd,
                                            std::size_t treeIndex);
    PathPairingWithSplits * _add_path_mapping(NodePairingWithSplits * parentPairing,
                                            NodePairingWithSplits * childPairing,
                                            std::size_t treeIndex);
    NodeEmbeddingWithSplits & _get_embedding_for_node(NodeWithSplits * nd);
    const NodeEmbeddingWithSplits & _get_embedding_for_node(const NodeWithSplits * nd) const {
        return scaffoldNdToNodeEmbedding.at(nd);
    }
    void embedTree(TreeMappedWithSplits & scaffold_tree,
                   TreeMappedWithSplits & tree,
                   std::size_t treeIndex,
                   bool isScaffoldClone);
};

inline NodeEmbeddingWithSplits & EmbeddedTree::_get_embedding_for_node(NodeWithSplits * nd) {
    const auto sIt = scaffoldNdToNodeEmbedding.find(nd);
    if (sIt == scaffoldNdToNodeEmbedding.end()) {
        auto eres = scaffoldNdToNodeEmbedding.emplace(std::piecewise_construct,
                                                      std::forward_as_tuple(nd),
                                                      std::forward_as_tuple(nd));
        assert(eres.second);
        return eres.first->second;
    }
    return sIt->second;
}

inline void EmbeddedTree::embed_new_tree(TreeMappedWithSplits & scaffold_tree,
                                TreeMappedWithSplits & tree,
                                std::size_t treeIndex) {
    embedTree(scaffold_tree, tree, treeIndex, false);
}

inline void EmbeddedTree::embed_scaffold_clone(TreeMappedWithSplits & scaffold_tree,
                                      TreeMappedWithSplits & tree,
                                      std::size_t treeIndex) {
    embedTree(scaffold_tree, tree, treeIndex, true);
}

} // namespace
#endif
