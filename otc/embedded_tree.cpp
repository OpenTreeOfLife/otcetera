#include <algorithm>
#include "otc/embedded_tree.h"
#include "otc/embedding.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"

namespace otc {

NodePairingWithSplits * EmbeddedTree::_addNodeMapping(NodeWithSplits *taxo, NodeWithSplits *nd, std::size_t treeIndex) {
    assert(taxo != nullptr);
    assert(nd != nullptr);
    nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
    auto ndPairPtr = &(*nodePairings.rbegin());
    auto & athreading = taxoToEmbedding[taxo];
    athreading.nodeEmbeddings[treeIndex].insert(ndPairPtr);
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
            taxoToEmbedding[currTaxo].edgeBelowEmbeddings[treeIndex].insert(pathPairPtr);
            currTaxo = currTaxo->getParent();
            if (currTaxo == nullptr) {
                break;
            }
        }
    } else {
        taxoToEmbedding[currTaxo].loopEmbeddings[treeIndex].insert(pathPairPtr);
    }
    return pathPairPtr;
}
void EmbeddedTree::embedNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
    // do threading
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

void EmbeddedTree::threadTaxonomyClone(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
    // do threading
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



void writeDOTContent(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> &,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree);
void writePathPairingToDOT(std::ostream & out, 
                           const PathPairingWithSplits & pp,
                           std::map<const NodeWithSplits *, std::string> & nd2name,
                           std::set<std::string> & names,
                           const std::string &style);


void writeNodeDOT(std::ostream & out,
                           const NodeWithSplits * pn,
                           std::map<const NodeWithSplits *, std::string> & nd2name,
                           std::set<std::string> & names,
                           const std::string &style) {
    if (!contains(nd2name, pn)) {
        const auto name = getDesignator(*pn);
        nd2name[pn] = name;
        if (!contains(names, name)) {
            out << "  \"" << nd2name[pn] << "\"";
            if (!style.empty()) {
                out << "[" << style <<']';
            }
            out << ";\n";
            names.insert(name);
        }
    }
}

const char * COLORS [] = {"crimson",
                    "blue",
                    "springgreen",
                    "magenta",
                    "darkorange",
                    "lightblue",
                    "goldenrod",
                    "brown",
                    "gray"};
constexpr unsigned LAST_COLOR_IND = 8;

void writePathPairingToDOT(std::ostream & out,
                           const PathPairingWithSplits & pp,
                           std::map<const NodeWithSplits *, std::string> & nd2name,
                           std::set<std::string> & names,
                           const std::string &style) {
    const auto * pn = pp.phyloParent;
    const auto * pd = pp.phyloChild;
    writeNodeDOT(out, pn, nd2name, names, style);
    writeNodeDOT(out, pd, nd2name, names, style);
    out << "    \"" << nd2name[pn] << "\" -> \"" << nd2name[pd] << '\"';
    if (!style.empty()) {
        out << "[" << style <<']';
    }
    out << ";\n";
}

void writeDOTContent(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> & tv,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree) {
    std::map<const NodeWithSplits *, std::string> nd2name;
    std::set<std::string> names;
    out << "digraph G{\n";
    std::string scafstyle = "fontcolor=\"black\" color=\"black\"";
    for (auto n : iter_pre_n_const(nd)) {
        writeNodeDOT(out, n, nd2name, names, scafstyle);
    }
    for (auto n : iter_pre_n_const(nd)) {
        auto p = n->getParent();
        if (p == nullptr) {
            continue;
        }
        out << "    " << nd2name[p] << " -> " << nd2name[n] << ";\n";
    }
    const NodeEmbeddingWithSplits & thr = eForNd.at(nd);
    const auto nt = tv.size();
    for (auto i = 0U; i < nt; ++i) {
        auto colorIndex = std::min(LAST_COLOR_IND, i);
        const char * color = COLORS[colorIndex];
        std::string style = "style=\"dashed\" fontcolor=\"";
        style.append(color);
        style.append("\" color=\"");
        style.append(color);
        style.append("\"");
        const auto eait = thr.edgeBelowEmbeddings.find(i);
        if (eait != thr.edgeBelowEmbeddings.end()) {
            for (const auto & pp : eait->second) {
                writePathPairingToDOT(out, *pp, nd2name, names, style);
            }
        }
        const auto lait = thr.loopEmbeddings.find(i);
        if (lait != thr.loopEmbeddings.end()) {
            for (const auto & pp : lait->second) {
                writePathPairingToDOT(out, *pp, nd2name, names, style);
            }
        }
    }

    out << "}\n";
}

void EmbeddedTree::writeDOTExport(std::ostream & out,
                       const NodeEmbedding<NodeWithSplits, NodeWithSplits> & thr,
                       const NodeWithSplits * nd,
                       const std::vector<TreeMappedWithSplits *> &t,
                       bool entireSubtree) const {
    writeDOTContent(out, nd, t, taxoToEmbedding, entireSubtree);
}






} // namespace otc

