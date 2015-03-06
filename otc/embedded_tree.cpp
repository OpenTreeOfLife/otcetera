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

typedef std::pair<const NodeWithSplits *, std::string> ToDotKey;
typedef std::map<ToDotKey, std::vector<std::string>> NodeToDotNames;

void writeDOTContent(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> &,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree);

void writeNodeDOT(std::ostream & out,
                           const NodeWithSplits * nd,
                           NodeToDotNames & nd2name,
                           std::set<std::string> & names,
                           const std::string &style, 
                           bool forceMRCANaming,
                           const std::string & prefix,
                           bool writeLabel,
                           bool pt) {
    const ToDotKey k{nd, prefix};
    std::string name;
    std::string unadorned;
    if (pt) {
        name = nd2name.at(k)[0];
        unadorned = name;
    } else {
        if (contains(nd2name, k)) {
            return;
        }
        if (nd->isTip() && !prefix.empty()) {
            return;
        } else {
            name.append(std::string(prefix));
            unadorned = forceMRCANaming ? getMRCADesignator(*nd) :  getDesignator(*nd);
            name.append(unadorned);
            nd2name[k].push_back(name);
        }
    }
    out << "  \"" << name << "\"[";
    if (pt) {
        out << "shape=point label=\"\" ";
    } else {
        if (writeLabel || nd->isTip()) {
            out << "label=\"" << unadorned << "\" ";
        } else {
            out << "label=\"\" ";
        }
    }
    if (!style.empty()) {
        out << style;
    }
    out << "];\n";
    names.insert(name);
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

void writeDOTEdge(std::ostream & out,
                  const ToDotKey anc,
                  const ToDotKey des,
                  NodeToDotNames & nd2name,
                  const std::string &style, 
                  bool toMidEdge) {
    out << "    \"" << nd2name[anc][0] << "\" -> \"" << nd2name[des][0] << '\"';
    if (!style.empty()) {
        out << "[" << style <<']';
    }
    out << ";\n";
}
void writeDOTCrossEdge(std::ostream & out,
                  const ToDotKey anc,
                  const ToDotKey des,
                  NodeToDotNames & nd2name,
                  const std::string &style) {
    out << "    \"" << nd2name[anc][0] << "\" -> \"" << nd2name[des][0] << '\"';
    if (!style.empty()) {
        out << "[arrowhead=none " << style <<']';
    }
    out << ";\n";
}

void writePathPairingToDOT(std::ostream & out,
                           const PathPairingWithSplits & pp,
                           NodeToDotNames & nd2name,
                           std::set<std::string> & names,
                           std::set<const PathPairingWithSplits *> & pathSet,
                           const char * color,
                           const char * prefix) {
    if (contains(pathSet, &pp)) {
        return;
    }
    pathSet.insert(&pp);
    std::string style = "fontcolor=\"";
    style.append(color);
    style.append("\" style=\"dashed\" color=\"");
    style.append(color);
    style.append("\"");
        
    const auto * pn = pp.phyloParent;
    const auto * pd = pp.phyloChild;
    writeNodeDOT(out, pn, nd2name, names, style, true, prefix, false, false);
    writeNodeDOT(out, pd, nd2name, names, style, false, prefix, false, false);
    ToDotKey pk{pd, "_path"};
    nd2name[pk].push_back("_" + std::to_string((long)(&pp)));
    writeNodeDOT(out, pd, nd2name, names, style, false, "_path", false, true);
    writeDOTEdge(out, ToDotKey{pn, prefix}, pk, nd2name, style, true);
    const auto * sn = pp.scaffoldAnc;
    const auto * sd = pp.scaffoldDes;
    if (pd->isTip()) {
        writeDOTEdge(out, pk, ToDotKey{sd, ""}, nd2name, style, false);
    } else {
        writeDOTEdge(out, pk, ToDotKey{pd, prefix}, nd2name, style, false);
    }
    ToDotKey sk{sd, "_path"};
    nd2name[sk].push_back("_" + std::to_string((long)(&pp)));
    writeNodeDOT(out, sn, nd2name, names, style, true, prefix, false, false);
    writeNodeDOT(out, sd, nd2name, names, style, false, prefix, false, false);
    writeNodeDOT(out, sd, nd2name, names, style, false, "_path", false, true);
    writeDOTEdge(out, ToDotKey{sn, ""}, sk, nd2name, style, true);
    writeDOTEdge(out, sk, ToDotKey{sd, ""}, nd2name, style, false);
    style.append(" style=\"dotted\"");
    writeDOTCrossEdge(out, pk, sk, nd2name, style);
}

void writeDOTEmbeddingForNode(std::ostream & out,
                              const NodeEmbeddingWithSplits & thr,
                              NodeToDotNames & nd2name,
                              std::set<std::string> & names,
                              std::set<const PathPairingWithSplits *> & pathSet,
                              const char * color,
                              const char * prefix,
                              std::size_t treeIndex
                              ) {
    const auto eait = thr.edgeBelowEmbeddings.find(treeIndex);
    if (eait != thr.edgeBelowEmbeddings.end()) {
        for (const auto & pp : eait->second) {
            writePathPairingToDOT(out, *pp, nd2name, names, pathSet, color, prefix);
        }
    }
    const auto lait = thr.loopEmbeddings.find(treeIndex);
    if (lait != thr.loopEmbeddings.end()) {
        for (const auto & pp : lait->second) {
            writePathPairingToDOT(out, *pp, nd2name, names, pathSet, color, prefix);
        }
    }

}
void writeDOTContent(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> & tv,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree) {
    NodeToDotNames nd2name;
    std::set<std::string> names;
    out << "digraph G{\n";
    for (auto n : iter_pre_n_const(nd)) {
        writeNodeDOT(out, n, nd2name, names, "", false, "", true, false);
    }
    for (auto n : iter_pre_n_const(nd)) {
        auto p = n->getParent();
        if (p == nullptr) {
            continue;
        }
        ToDotKey nk{n, ""};
        ToDotKey pk{p, ""};
        out << "    " << nd2name[pk][0] << " -> " << nd2name[nk][0] << ";\n";
    }
    std::set<const PathPairingWithSplits *> pathSet;
    const auto nt = tv.size();
    for (auto n : iter_pre_n_const(nd)) {
        const NodeEmbeddingWithSplits & thr = eForNd.at(n);
        for (auto i = 0U; i < nt; ++i) {
            const std::string tP = std::string("t") + std::to_string(i);
            auto colorIndex = std::min(LAST_COLOR_IND, i);
            const char * color = COLORS[colorIndex];
            writeDOTEmbeddingForNode(out, thr, nd2name, names, pathSet, color, tP.c_str(), i);
        }
        if (!entireSubtree) {
            break;
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

