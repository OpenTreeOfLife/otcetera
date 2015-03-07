#include "otc/write_dot.h"
#include <algorithm>
#include "otc/embedded_tree.h"
#include "otc/embedding.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"

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
                           const PathPairingWithSplits & pp,
                           NodeToDotNames & nd2name,
                           std::set<const PathPairingWithSplits *> & pathSet,
                           const char * color,
                           const char * prefix) ;
void writeDOTEmbeddingForNode(std::ostream & out,
                              const NodeEmbeddingWithSplits & thr,
                              NodeToDotNames & nd2name,
                              std::set<const PathPairingWithSplits *> & pathSet,
                              const char * color,
                              const char * prefix,
                              std::size_t treeIndex
                              );
void writeOneSideOfPathPairingToDOT(std::ostream & out,
                                     const NodeWithSplits * pn,
                                     const NodeWithSplits * pd,
                                     const NodeWithSplits * sd,
                                     const ToDotKey & midPointKey,
                                     const NamePair & namePair,
                                     NodeToDotNames & nd2name,
                                     const std::string & style,
                                     const char * prefix);
// IMPL
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
                           const NodeWithSplits * pn,
                           const NodeWithSplits * pd,
                           const NodeWithSplits * sd,
                           const ToDotKey & midPointKey,
                           const NamePair & namePair,
                           NodeToDotNames & nd2name,
                           const std::string & style,
                           const char * prefix) {
    const bool isScaffoldSide = (pd == sd);
    const char * ptu = (isScaffoldSide ? "" : prefix);
    const ToDotKey nk{pn, ptu};
    const ToDotKey dk{pd, ptu};
    writeNodeDOT(out, nk, nd2name, style, true, false, false);
    writeNodeDOT(out, dk, nd2name, style, false, false, false);
    nd2name[midPointKey] = namePair;
    uncheckedWriteNodeDOT(out, midPointKey, nd2name, style, false, false, true);
    writeDOTEdge(out, ToDotKey{pn, ptu}, midPointKey, nd2name, style, true);
    if (pd->isTip() && sd != pd) {
        writeDOTEdge(out, midPointKey, ToDotKey{sd, ""}, nd2name, style, false);
    } else {
        writeDOTEdge(out, midPointKey, dk, nd2name, style, false);
    }
    
}

void writePathPairingToDOT(std::ostream & out,
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
    std::string bstyle = "fontcolor=\"";
    bstyle.append(color);
    bstyle.append("\" color=\"");
    bstyle.append(color);
    bstyle.append("\"");
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
    const ToDotKey pk{pd, "_phpath"};
    writeOneSideOfPathPairingToDOT(out, pn, pd, sd, pk, pv, nd2name, style, prefix);
    const auto * sn = pp.scaffoldAnc;
    style = bstyle;
    style.append(" style=\"solid\"");
    const ToDotKey sk{sd, "_scpath"};
    writeOneSideOfPathPairingToDOT(out, sn, sd, sd, sk, sv, nd2name, style, prefix);
    style = bstyle;
    style.append(" style=\"dotted\"");
    writeDOTCrossEdge(out, pk, sk, nd2name, style);
}

void writeDOTEmbeddingForNode(std::ostream & out,
                              const NodeEmbeddingWithSplits & thr,
                              NodeToDotNames & nd2name,
                              std::set<const PathPairingWithSplits *> & pathSet,
                              const char * color,
                              const char * prefix,
                              std::size_t treeIndex
                              ) {
    const auto eait = thr.edgeBelowEmbeddings.find(treeIndex);
    if (eait != thr.edgeBelowEmbeddings.end()) {
        for (const auto & pp : eait->second) {
            writePathPairingToDOT(out, *pp, nd2name, pathSet, color, prefix);
        }
    }
    const auto lait = thr.loopEmbeddings.find(treeIndex);
    if (lait != thr.loopEmbeddings.end()) {
        for (const auto & pp : lait->second) {
            writePathPairingToDOT(out, *pp, nd2name, pathSet, color, prefix);
        }
    }
}

void writeDOTForEmbedding(std::ostream & out,
                     const NodeWithSplits * nd,
                     const std::vector<TreeMappedWithSplits *> & tv,
                     const std::map<const NodeWithSplits *, NodeEmbeddingWithSplits> & eForNd,
                     bool entireSubtree) {
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
    const auto nt = tv.size();
    for (auto n : iter_pre_n_const(nd)) {
        const NodeEmbeddingWithSplits & thr = eForNd.at(n);
        for (auto i = 0U; i < nt; ++i) {
            const std::string tP = std::string("t") + std::to_string(i);
            auto colorIndex = std::min(LAST_COLOR_IND, i);
            const char * color = COLORS[colorIndex];
            writeDOTEmbeddingForNode(out, thr, nd2name, pathSet, color, tP.c_str(), i);
        }
        if (!entireSubtree) {
            break;
        }
    }
    out << "}\n";
}

} // namespace otc

