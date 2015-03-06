#ifndef OTCETERA_EMBEDDED_TREE_H
#define OTCETERA_EMBEDDED_TREE_H

#include <map>
#include <list>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_data.h"

namespace otc {
class EmbeddedTree {
    protected:
    std::list<NodePairingWithSplits> nodePairings;
    std::list<PathPairingWithSplits> pathPairings;
    std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> taxoToEmbedding;
    
    public:
    EmbeddedTree() {
    }

    NodePairingWithSplits * _addNodeMapping(NodeWithSplits *taxo,
                                            NodeWithSplits *nd,
                                            std::size_t treeIndex);
    PathPairingWithSplits * _addPathMapping(NodePairingWithSplits * parentPairing,
                                            NodePairingWithSplits * childPairing,
                                            std::size_t treeIndex);
    void embedNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex);
    void threadTaxonomyClone(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex);
    void writeDOTExport(std::ostream & out,
                           const NodeEmbedding<NodeWithSplits, NodeWithSplits> & thr,
                           const NodeWithSplits * nd,
                           const std::vector<TreeMappedWithSplits *> &,
                           bool entireSubtree) const;
};



} // namespace
#endif
