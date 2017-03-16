#include <algorithm>
#include "otc/embedded_tree.h"
#include "otc/node_embedding.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/write_dot.h"
namespace otc {

NodePairingWithSplits * EmbeddedTree::_add_node_mapping(
                NodeWithSplits *taxo,
                NodeWithSplits *nd, std::size_t treeIndex) {
    assert(taxo != nullptr);
    assert(nd != nullptr);
    nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
    auto ndPairPtr = &(*nodePairings.rbegin());
    _get_embedding_for_node(taxo).add_node_embeddings(treeIndex, ndPairPtr);
    return ndPairPtr;
}

PathPairingWithSplits * EmbeddedTree::_add_path_mapping(
                NodePairingWithSplits * parentPairing,
                NodePairingWithSplits * childPairing,
                std::size_t treeIndex) {
    pathPairings.emplace_back(*parentPairing, *childPairing);
    auto pathPairPtr = &(*pathPairings.rbegin());
    // register a pointer to the path at each traversed...
    auto currTaxo = pathPairPtr->scaffold_des;
    auto ancTaxo = pathPairPtr->scaffold_anc;
    auto & ne = _get_embedding_for_node(currTaxo);
    if (currTaxo != ancTaxo) {
        while (currTaxo != ancTaxo) {
             _get_embedding_for_node(currTaxo).add_exit_embeddings(treeIndex, pathPairPtr);
            currTaxo = currTaxo->get_parent();
            if (currTaxo == nullptr) {
                break;
            }
        }
    } else {
        ne.add_loop_embeddings(treeIndex, pathPairPtr);
    }
    return pathPairPtr;
}

void EmbeddedTree::embedTree(TreeMappedWithSplits & scaffold_tree,
                                      TreeMappedWithSplits & tree,
                                      std::size_t treeIndex,
                                      bool isScaffoldClone) {
    // do embedding
    std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
    std::set<NodePairingWithSplits *> tipPairings;
    for (auto nd : iter_post(tree)) {
        auto par = nd->get_parent();
        if (par == nullptr) {
            continue;
        }
        NodePairingWithSplits * ndPairPtr = nullptr;
        NodeWithSplits * taxoDes = nullptr;
        if (nd->is_tip()) {
             // TMP, Remove this next assert to save time?
            assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end());
            assert(nd->has_ott_id());
            auto ottId = nd->get_ott_id();
            taxoDes = scaffold_tree.get_data().get_node_by_ott_id(ottId);
            assert(taxoDes != nullptr);
            ndPairPtr = _add_node_mapping(taxoDes, nd, treeIndex);
            if (!isScaffoldClone) {
                for (auto former : tipPairings) {
                    if (are_linearly_related(taxoDes, former->scaffold_node)) {
                        std::string m = "Repeated or nested OTT ID in tip mapping of an input tree: \"";
                        m += nd->get_name();
                        m += "\" and \"";
                        m += former->phylo_node->get_name();
                        m += "\" found.";
                        throw OTCError(m);
                    }
                }
                tipPairings.insert(ndPairPtr);
            }
            currTreeNodePairings[nd] = ndPairPtr;
        } else {
            auto reuseNodePairingIt = currTreeNodePairings.find(nd);
            assert(reuseNodePairingIt != currTreeNodePairings.end());
            ndPairPtr = reuseNodePairingIt->second;
            taxoDes = ndPairPtr->scaffold_node;
            assert(taxoDes != nullptr);
        }
        NodePairingWithSplits * parPairPtr = nullptr;
        auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
        if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
            NodeWithSplits * taxoAnc = nullptr;
            if (isScaffoldClone) {
                // since it is a taxonomy, it will have internal node labels
                auto pottId = par->get_ott_id();
                taxoAnc = scaffold_tree.get_data().get_node_by_ott_id(pottId);
            
            } else {
                const auto & parDesIds = par->get_data().des_ids;
                taxoAnc = search_anc_for_mrca_of_des_ids(taxoDes, parDesIds);
            }
            assert(taxoAnc != nullptr);
            parPairPtr = _add_node_mapping(taxoAnc, par, treeIndex);
            currTreeNodePairings[par] = parPairPtr;
        } else {
            parPairPtr = prevAddedNodePairingIt->second;
        }
        _add_path_mapping(parPairPtr, ndPairPtr, treeIndex);
    }
}


void EmbeddedTree::write_dot_export(std::ostream & out,
                       const NodeEmbedding<NodeWithSplits, NodeWithSplits> & ,
                       const NodeWithSplits * nd,
                       const std::vector<TreeMappedWithSplits *> &t,
                       bool entireSubtree,
                       bool includeLastTree) const {
    write_dot_for_embedding(out, nd, t, scaffoldNdToNodeEmbedding, entireSubtree, includeLastTree);
}

} // namespace otc

