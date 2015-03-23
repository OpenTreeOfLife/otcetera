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
    void embedNewTree(TreeMappedWithSplits & scaffoldTree,
                      TreeMappedWithSplits & tree,
                      std::size_t treeIndex);
    void embedScaffoldClone(TreeMappedWithSplits & scaffoldTree,
                            TreeMappedWithSplits & tree,
                            std::size_t treeIndex);
    void writeDOTExport(std::ostream & out,
                        const NodeEmbedding<NodeWithSplits, NodeWithSplits> & thr,
                        const NodeWithSplits * nd,
                        const std::vector<TreeMappedWithSplits *> &,
                        bool entireSubtree,
                        bool includeLastTree) const;
    // for testing...
    std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> &_getScaffoldNdToNodeEmbedding() {
        return scaffoldNdToNodeEmbedding;
    }
    protected:
    NodePairingWithSplits * _addNodeMapping(NodeWithSplits *taxo,
                                            NodeWithSplits *nd,
                                            std::size_t treeIndex);
    PathPairingWithSplits * _addPathMapping(NodePairingWithSplits * parentPairing,
                                            NodePairingWithSplits * childPairing,
                                            std::size_t treeIndex);
    NodeEmbeddingWithSplits & _getEmbeddingForNode(NodeWithSplits * nd);
    const NodeEmbeddingWithSplits & _getEmbeddingForNode(const NodeWithSplits * nd) const {
        return scaffoldNdToNodeEmbedding.at(nd);
    }
    void embedTree(TreeMappedWithSplits & scaffoldTree,
                   TreeMappedWithSplits & tree,
                   std::size_t treeIndex,
                   bool isScaffoldClone);
};

inline NodeEmbeddingWithSplits & EmbeddedTree::_getEmbeddingForNode(NodeWithSplits * nd) {
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

inline void EmbeddedTree::embedNewTree(TreeMappedWithSplits & scaffoldTree,
                                TreeMappedWithSplits & tree,
                                std::size_t treeIndex) {
    embedTree(scaffoldTree, tree, treeIndex, false);
}

inline void EmbeddedTree::embedScaffoldClone(TreeMappedWithSplits & scaffoldTree,
                                      TreeMappedWithSplits & tree,
                                      std::size_t treeIndex) {
    embedTree(scaffoldTree, tree, treeIndex, true);
}

} // namespace
#endif
