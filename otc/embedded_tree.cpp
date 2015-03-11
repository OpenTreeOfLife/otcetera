#include <algorithm>
#include "otc/embedded_tree.h"
#include "otc/embedding.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/write_dot.h"
namespace otc {

NodeEmbeddingWithSplits & EmbeddedTree::_getEmdeddingForNode(NodeWithSplits * nd) {
    std::map<const NodeWithSplits *, NodeEmbeddingWithSplits>::iterator sIt = taxoToEmbedding.find(nd);
    if (sIt == taxoToEmbedding.end()) {
        auto eres = taxoToEmbedding.emplace(std::piecewise_construct,
                                            std::forward_as_tuple(nd),
                                            std::forward_as_tuple(nd));
        assert(eres.second);
        return eres.first->second;
    }
    return sIt->second;
}

NodePairingWithSplits * EmbeddedTree::_addNodeMapping(NodeWithSplits *taxo, NodeWithSplits *nd, std::size_t treeIndex) {
    assert(taxo != nullptr);
    assert(nd != nullptr);
    nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
    auto ndPairPtr = &(*nodePairings.rbegin());
    _getEmdeddingForNode(taxo).nodeEmbeddings[treeIndex].insert(ndPairPtr);
    return ndPairPtr;
}
PathPairingWithSplits * EmbeddedTree::_addPathMapping(NodePairingWithSplits * parentPairing,
                                        NodePairingWithSplits * childPairing,
                                        std::size_t treeIndex) {
    pathPairings.emplace_back(*parentPairing, *childPairing);
    auto pathPairPtr = &(*pathPairings.rbegin());
    // register a pointer to the path at each traversed...
    auto currTaxo = pathPairPtr->scaffoldDes;
    auto ancTaxo = pathPairPtr->scaffoldAnc;
    if (currTaxo != ancTaxo) {
        while (currTaxo != ancTaxo) {
            _getEmdeddingForNode(currTaxo).edgeBelowEmbeddings[treeIndex].insert(pathPairPtr);
            currTaxo = currTaxo->getParent();
            if (currTaxo == nullptr) {
                break;
            }
        }
    } else {
        _getEmdeddingForNode(currTaxo).loopEmbeddings[treeIndex].insert(pathPairPtr);
    }
    return pathPairPtr;
}
void EmbeddedTree::embedNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
    // do embedding
    std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
    std::set<NodePairingWithSplits *> tipPairings;
    for (auto nd : iter_post(tree)) {
        auto par = nd->getParent();
        if (par == nullptr) {
            continue;
        }
        NodePairingWithSplits * ndPairPtr = nullptr;
        NodeWithSplits * taxoDes = nullptr;
        if (nd->isTip()) {
            assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
            assert(nd->hasOttId());
            auto ottId = nd->getOttId();
            taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
            assert(taxoDes != nullptr);
            ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
            for (auto former : tipPairings) {
                if (areLinearlyRelated(taxoDes, former->scaffoldNode)) {
                    std::string m = "Repeated or nested OTT ID in tip mapping of an input tree: \"";
                    m += nd->getName();
                    m += "\" and \"";
                    m += former->phyloNode->getName();
                    m += "\" found.";
                    throw OTCError(m);
                }
            }
            tipPairings.insert(ndPairPtr);
            currTreeNodePairings[nd] = ndPairPtr;
        } else {
            auto reuseNodePairingIt = currTreeNodePairings.find(nd);
            assert(reuseNodePairingIt != currTreeNodePairings.end());
            ndPairPtr = reuseNodePairingIt->second;
            taxoDes = ndPairPtr->scaffoldNode;
            assert(taxoDes != nullptr);
        }
        NodePairingWithSplits * parPairPtr = nullptr;
        auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
        if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
            const auto & parDesIds = par->getData().desIds;
            auto taxoAnc = searchAncForMRCAOfDesIds(taxoDes, parDesIds);
            assert(taxoAnc != nullptr);
            parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
            currTreeNodePairings[par] = parPairPtr;
        } else {
            parPairPtr = prevAddedNodePairingIt->second;
        }
        _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
    }
}

void EmbeddedTree::embedTaxonomyClone(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
    // do embedding
    std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
    for (auto nd : iter_post(tree)) {
        auto par = nd->getParent();
        if (par == nullptr) {
            continue;
        }
        NodePairingWithSplits * ndPairPtr = nullptr;
        NodeWithSplits * taxoDes = nullptr;
        if (nd->isTip()) {
            assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
            assert(nd->hasOttId());
            auto ottId = nd->getOttId();
            taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
            assert(taxoDes != nullptr);
            ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
            currTreeNodePairings[nd] = ndPairPtr;
        } else {
            auto reuseNodePairingIt = currTreeNodePairings.find(nd);
            assert(reuseNodePairingIt != currTreeNodePairings.end());
            ndPairPtr = reuseNodePairingIt->second;
            taxoDes = ndPairPtr->scaffoldNode;
            assert(taxoDes != nullptr);
        }
        NodePairingWithSplits * parPairPtr = nullptr;
        auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
        if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
            auto pottId = par->getOttId(); // since it is a taxonomy, it will have internal node labels
            auto taxoAnc = scaffoldTree.getData().getNodeForOttId(pottId);
            assert(taxoAnc != nullptr);
            parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
            currTreeNodePairings[par] = parPairPtr;
        } else {
            parPairPtr = prevAddedNodePairingIt->second;
        }
        _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
    }
}

void EmbeddedTree::writeDOTExport(std::ostream & out,
                       const NodeEmbedding<NodeWithSplits, NodeWithSplits> & ,
                       const NodeWithSplits * nd,
                       const std::vector<TreeMappedWithSplits *> &t,
                       bool entireSubtree,
                       bool includeLastTree) const {
    writeDOTForEmbedding(out, nd, t, taxoToEmbedding, entireSubtree, includeLastTree);
}

} // namespace otc

