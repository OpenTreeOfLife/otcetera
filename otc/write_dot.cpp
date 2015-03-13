#include "otc/write_dot.h"
#include <algorithm>
#include "otc/embedded_tree.h"
#include "otc/embedding.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/forest.h"

namespace otc {
typedef std::pair<const NodeWithSplits *, std::string> ToDotKey;
typedef std::pair<std::string, std::string> NamePair; // nodeName and nodeLabel
typedef std::map<ToDotKey, NamePair> NodeToDotNames;
void writeNodeDOT(std::ostream & out,
                           const ToDotKey & k,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           bool writeLabel,
                           bool pt);
void uncheckedWriteNodeDOT(std::ostream & out,
                           const ToDotKey & k,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           bool writeLabel,
                           bool pt);
void writeDOTEdge(std::ostream & out,
                  const ToDotKey anc,
                  const ToDotKey des,
                  NodeToDotNames & nd2name,
                  const std::string &style, 
                  bool toMidEdge);
void writeDOTCrossEdge(std::ostream & out,
                  const ToDotKey anc,
                  const ToDotKey des,
                  NodeToDotNames & nd2name,
                  const std::string &style);
void writePathPairingToDOT(std::ostream & out,
                           const NodeWithSplits *n, 
                           const PathPairingWithSplits & pp,
                           NodeToDotNames & nd2name,
                           std::set<const PathPairingWithSplits *> & pathSet,
                           const char * color,
                           const char * prefix) ;
void writeDOTEmbeddingForNode(std::ostream & out,
                              const NodeWithSplits *n, 
                              const NodeEmbeddingWithSplits & thr,
                              const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                              NodeToDotNames & nd2name,
                              std::set<const PathPairingWithSplits *> & pathSet,
                              const char * color,
                              const char * prefix,
                              std::size_t treeIndex
                              );
void writeOneSideOfPathPairingToDOT(std::ostream & out,
                                     const NodeWithSplits * pn,
                                     const NodeWithSplits * pd,
                                     const NodeWithSplits * sn,
                                     const NodeWithSplits * sd,
                                     const ToDotKey & midPointKey,
                                     const NamePair & namePair,
                                     NodeToDotNames & nd2name,
                                     const NodeWithSplits * focalNode,
                                     const std::string & style,
                                     const char * prefix);
void writeDOTForNodeWithoutEmbedding(std::ostream & out, const NodeWithSplits *nd , NodeToDotNames & nd2name);
std::string edgeStyleForColor(const char * color);
// IMPL
std::string edgeStyleForColor(const char * color) {
    std::string bstyle = "fontcolor=\"";
    bstyle.append(color);
    bstyle.append("\" color=\"");
    bstyle.append(color);
    bstyle.append("\"");
    return bstyle;
}
constexpr const char * COLORS [] = {"crimson",
                                    "blue",
                                    "springgreen",
                                    "magenta",
                                    "darkorange",
                                    "lightblue",
                                    "goldenrod",
                                    "brown",
                                    "gray"};
constexpr unsigned LAST_COLOR_IND = 8;

void writeNodeDOT(std::ostream & out,
                           const ToDotKey & k,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           bool writeLabel,
                           bool pt) {
    if (!contains(nd2name, k)) {
        uncheckedWriteNodeDOT(out, k, nd2name, style, forceMRCANaming, writeLabel, pt);
    }
}

void writeDOTForNodeWithoutEmbedding(std::ostream & , const NodeWithSplits * , NodeToDotNames & ) {
}

void uncheckedWriteNodeDOT(std::ostream & out,
                           const ToDotKey & k,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           bool writeLabel,
                           bool pt) {
    const NodeWithSplits * nd{k.first};
    const std::string & prefix{k.second};
    std::string name{prefix};
    const std::string unadorned = forceMRCANaming ? getMRCADesignator(*nd) :  getDesignator(*nd);
    name.append(unadorned);
    nd2name.emplace(k, NamePair{name, unadorned});
    assert(contains(nd2name, k));
    const NamePair & np = nd2name.at(k);
    out << "  \"" << np.first << "\"[";
    if (pt) {
        out << "shape=point ";
    }
    if ((!pt) && writeLabel) {
        out << "label=\"" << np.second << "\" ";
    } else {
        out << "label=\"\" ";
    }
    if (!style.empty()) {
        out << style;
    }
    out << "];\n";
}

void writeDOTEdge(std::ostream & out,
                  const ToDotKey anc,
                  const ToDotKey des,
                  NodeToDotNames & nd2name,
                  const std::string &style, 
                  bool /*toMidEdge*/) {
    assert(contains(nd2name, anc));
    assert(contains(nd2name, des));
    out << "    \"" << nd2name.at(anc).first << "\" -> \"" << nd2name.at(des).first << '\"';
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
    std::string s = "arrowhead=none ";
    s += style;
    writeDOTEdge(out, anc, des, nd2name, s, true);
}

void writeOneSideOfPathPairingToDOT(std::ostream & out,
                           const NodeWithSplits * sidePar,
                           const NodeWithSplits * sideDes,
                           const NodeWithSplits * scaffAnc,
                           const NodeWithSplits * scaffDes,
                           const ToDotKey & midPointKey,
                           const NamePair & namePair,
                           NodeToDotNames & nd2name,
                           const NodeWithSplits * focalNode,
                           const std::string & style,
                           const char * prefix) {
    const bool isScaffoldSide = (sideDes == scaffDes);
    constexpr bool mergeScPhAnc = false;
    const char * ancPref = (isScaffoldSide || (mergeScPhAnc && focalNode != scaffAnc) ? "" : prefix);
    const char * desPref = (isScaffoldSide || sideDes->isTip() ? "" : prefix);
    // If we have a tip or a node outside of this focal node, use the scaffold node for the 
    //  phylo side of the graph. Not accurate, but cuts down on the # of nodes.
    ToDotKey ancK{sidePar, ancPref};
    if (mergeScPhAnc && focalNode != scaffAnc) {
        ancK = ToDotKey{scaffAnc, ""};
    }
    ToDotKey desK{sideDes, desPref};
    if (scaffDes->isTip()) {
        desK = ToDotKey{scaffDes, ""};
    }
    
    writeNodeDOT(out, ancK, nd2name, style, true, false, false);
    writeNodeDOT(out, desK, nd2name, style, false, false, false);
    // always write the point node along the path
    nd2name[midPointKey] = namePair;
    uncheckedWriteNodeDOT(out, midPointKey, nd2name, style, false, false, true);
    // write the edges anck -> midpoint -> desk
    writeDOTEdge(out, ancK, midPointKey, nd2name, style, true);
    writeDOTEdge(out, midPointKey, desK, nd2name, style, false);
}

void writePathPairingToDOT(std::ostream & out,
                           const NodeWithSplits *n, 
                           const PathPairingWithSplits & pp,
                           NodeToDotNames & nd2name,
                           std::set<const PathPairingWithSplits *> & pathSet,
                           const char * color,
                           const char * prefix) {
    if (contains(pathSet, &pp)) {
        return;
    }
    const std::string emptyStr;
    pathSet.insert(&pp);
    std::string bstyle = edgeStyleForColor(color);
    std::string style = bstyle;
    style.append(" style=\"dashed\"");
    // The midpoint nodes are unlabeled dots with IDs that are _addr or __addr
    const std::string pname ="_" + std::to_string((long)(&pp));
    const std::string sname ="__" + std::to_string((long)(&pp));
    const NamePair pv{pname, emptyStr};
    const NamePair sv{sname, emptyStr};
    const auto * pn = pp.phyloParent;
    const auto * pd = pp.phyloChild;
    const auto * sd = pp.scaffoldDes;
    const auto * sn = pp.scaffoldAnc;
    const ToDotKey pk{pd, "_phpath"};
    writeOneSideOfPathPairingToDOT(out, pn, pd, sn, sd, pk, pv, nd2name, n, style, prefix);
    style = bstyle;
    style.append(" style=\"solid\"");
    const ToDotKey sk{sd, "_scpath"};
    writeOneSideOfPathPairingToDOT(out, sn, sd, sn, sd, sk, sv, nd2name, n, style, prefix);
    style = bstyle;
    style.append(" style=\"dotted\"");
    writeDOTCrossEdge(out, pk, sk, nd2name, style);
}

void writeDOTEmbeddingForNode(std::ostream & out,
                              const NodeWithSplits * nd,
                              const NodeEmbeddingWithSplits & thr,
                              const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                              NodeToDotNames & nd2name,
                              std::set<const PathPairingWithSplits *> & pathSet,
                              const char * color,
                              const char * prefix,
                              std::size_t treeIndex
                              ) {
    const auto eait = thr.edgeBelowEmbeddings.find(treeIndex);
    if (eait != thr.edgeBelowEmbeddings.end()) {
        for (const auto & pp : eait->second) {
            writePathPairingToDOT(out, nd, *pp, nd2name, pathSet, color, prefix);
        }
    }
    const auto incoming = thr.getAllIncomingPathPairs(eForNd, treeIndex);
    for (auto pp : incoming) {
        writePathPairingToDOT(out, nd, *pp, nd2name, pathSet, color, prefix);
    }
}

void writeDOTForEmbedding(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> & tv,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree,
                     bool includeLastTree) {
    NodeToDotNames nd2name;
    std::string emptyStr;
    out << "digraph G{\n";
    for (auto n : iter_pre_n_const(nd)) {
        const ToDotKey k{n, ""};
        writeNodeDOT(out, k, nd2name, "", false, true, false);
    }
    auto ndp = nd->getParent();
    if (ndp != nullptr) {
        const ToDotKey k{ndp, ""};
        writeNodeDOT(out, k, nd2name, "", false, true, false);
    }
    for (auto n : iter_pre_n_const(nd)) {
        auto p = n->getParent();
        if (p == nullptr) {
            continue;
        }
        const ToDotKey nk{n, ""};
        const ToDotKey pk{p, ""};
        writeDOTEdge(out, pk, nk, nd2name, emptyStr, false);
    }
    std::set<const PathPairingWithSplits *> pathSet;
    const auto nt = tv.size() - (includeLastTree ? 0U : 1U);
    for (auto n : iter_pre_n_const(nd)) {
        const auto emIt = eForNd.find(n);
        if (emIt == eForNd.end()) {
            writeDOTForNodeWithoutEmbedding(out, n, nd2name);
            continue;
        }
        const NodeEmbeddingWithSplits & thr = eForNd.at(n);
        for (auto i = 0U; i < nt; ++i) {
            const std::string tP = std::string("t") + std::to_string(i);
            auto colorIndex = std::min(LAST_COLOR_IND, i);
            const char * color = COLORS[colorIndex];
            writeDOTEmbeddingForNode(out, n, thr, eForNd, nd2name, pathSet, color, tP.c_str(), i);
        }
        if (!entireSubtree) {
            break;
        }
    }
    out << "}\n";
}


template<typename T>
using GCList = std::list<std::pair<T *, PhyloStatementSource> >;
template<typename T>
using NdToGCListMapping = std::map<T*, GCList<T> >;
template<typename T>
using NdToGCMapping = std::map<T*, std::pair<T *, PhyloStatementSource> >;
template<typename T>
void writeDOTGroupingConstraints(std::ostream & out,
                                 const NdToGCListMapping<T> & ndToGCList,
                                 const std::string & prefix,
                                 const std::string & connPrefix,
                                 NodeToDotNames & nd2name,
                                 const std::string & ndStyle,
                                 const std::string & style) {
    for (const auto & incNdGCPair : ndToGCList ) {
        const auto & ind = incNdGCPair.first;
        const auto & incList = incNdGCPair.second;
        if (incList.empty()) {
            continue;
        }
        const ToDotKey desK{incNdGCPair.first, prefix};
        writeNodeDOT(out, desK, nd2name, ndStyle, false, true, false);
        for (const auto & igc : incList) {
            const ToDotKey ancK{igc.first, connPrefix}; // use connPrefix for the prefix to connect to a connected node
            writeDOTEdge(out, ancK, desK, nd2name, style, false);
        }
    }
}
template<typename T>
void writeDOTGroupingConstraintSingle(std::ostream & out,
                                 const NdToGCMapping<T> & ndToGCList,
                                 const std::string & prefix,
                                 const std::string & connPrefix,
                                 NodeToDotNames & nd2name,
                                 const std::string & ndStyle,
                                 const std::string & style) {
    for (const auto & incNdGCPair : ndToGCList ) {
        const auto & ind = incNdGCPair.first;
        const ToDotKey desK{incNdGCPair.first, prefix};
        writeNodeDOT(out, desK, nd2name, ndStyle, false, true, false);
        const ToDotKey ancK{incNdGCPair.second.first, connPrefix}; // use connPrefix for the prefix to connect to a connected node
        writeDOTEdge(out, ancK, desK, nd2name, style, false);
    }
}

void writeDOTForFtree(std::ostream & out,
                              const FTree<RTSplits, MappedWithSplitsData> & tree, 
                              NodeToDotNames & nd2name,
                              const char * color,
                              std::size_t treeIndex
                              ) {
    std::string esty = edgeStyleForColor(color);
    const std::string nsty = esty;
    std::string tn = "t";
    tn += std::to_string(treeIndex);
    auto r = tree.getRoot();
    for (auto n : iter_pre_n_const(r)) {
        auto p = n->getParent();
        if (p == nullptr) {
            continue;
        }
        const ToDotKey ancK{p, tn};
        ToDotKey desK{n, tn};
        if (n->isTip()) {
            desK = ToDotKey{n, ""};
        }
        writeNodeDOT(out, ancK, nd2name, nsty, false, true, false);
        writeNodeDOT(out, desK, nd2name, nsty, false, true, false);
        writeDOTEdge(out, ancK, desK, nd2name, esty, false);
    }
    const std::string itn = std::string("inc") + tn;
    const auto & inc = tree.getIncluded2ConstraintMap();
    std::string incStyle = "shape=invtriangle color=\"";
    incStyle += color;
    incStyle += "\"";
    writeDOTGroupingConstraintSingle(out, inc, itn, tn, nd2name, incStyle, esty);
    const std::string etn = std::string("exc") + tn;
    const auto & exc = tree.getExcluded2ConstraintMap();
    std::string excStyle = "shape=octagon color=\"";
    excStyle += color;
    excStyle += "\"";
    writeDOTGroupingConstraints(out, exc, etn, tn, nd2name, excStyle, esty);
}

void writeDOTForest(std::ostream & out, const RootedForest<RTSplits, MappedWithSplitsData> &forest) {
    const auto & o2n = forest.getOttIdToNodeMapping();
    NodeToDotNames nd2name;
    std::string emptyStr;
    out << "digraph G{\n";
    for (auto & oidNodePair : o2n) {
        auto n = oidNodePair.second;
        bool writePlainNd = forest.isAttached(n->getOttId());
        if (!writePlainNd && (!forest.isMentionedInInclude(n)) && (!forest.isMentionedInExclude(n))) {
            writePlainNd = true;
        }
        if (writePlainNd) {
            const ToDotKey k{n, ""};
            writeNodeDOT(out, k, nd2name, emptyStr, false, true, false);
        }
    }
    const auto & trees = forest.getTrees();
    auto i = 0U;
    for (const auto & tiTrPair : trees) {
        const std::string tP = std::string("t") + std::to_string(i);
        auto colorIndex = std::min(LAST_COLOR_IND, i);
        const char * color = COLORS[colorIndex];
        const auto & tree = tiTrPair.second;
        writeDOTForFtree(out, tree, nd2name, color, i);
        ++i;
    }
    out << "}\n";
}

} // namespace otc

