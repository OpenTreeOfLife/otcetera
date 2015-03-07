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
                           const NodeWithSplits * nd,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           const std::string & prefix,
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
                           const NodeWithSplits * nd,
                           NodeToDotNames & nd2name,
                           const std::string &style, 
                           bool forceMRCANaming,
                           const std::string & prefix,
                           bool writeLabel,
                           bool pt) {
    const ToDotKey k{nd, prefix};
    if (!pt) {
        if (contains(nd2name, k)) {
            return;
        }
        if (nd->isTip() && !prefix.empty()) {
            return;
        }
        std::string name{prefix};
        const std::string unadorned = forceMRCANaming ? getMRCADesignator(*nd) :  getDesignator(*nd);
        name.append(unadorned);
        nd2name.emplace(k, NamePair{name, unadorned});
    }
    const NamePair & np = nd2name.at(k);
    out << "  \"" << np.first << "\"[";
    if (pt) {
        out << "shape=point ";
    }
    if ((!pt) && (writeLabel || nd->isTip())) {
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
                                     const char * prefix);

void writeOneSideOfPathPairingToDOT(std::ostream & out,
                           const NodeWithSplits * pn,
                           const NodeWithSplits * pd,
                           const NodeWithSplits * sd,
                           const ToDotKey & midPointKey,
                           const NamePair & namePair,
                           NodeToDotNames & nd2name,
                           const std::string & style,
                           const char * prefix) {
    writeNodeDOT(out, pn, nd2name, style, true, prefix, false, false);
    writeNodeDOT(out, pd, nd2name, style, false, prefix, false, false);
    nd2name[midPointKey] = namePair;
    writeNodeDOT(out, pd, nd2name, style, false, "_path", false, true);
    writeDOTEdge(out, ToDotKey{pn, prefix}, midPointKey, nd2name, style, true);
    if (pd->isTip() && sd != pd) {
        writeDOTEdge(out, midPointKey, ToDotKey{sd, ""}, nd2name, style, false);
    } else {
        writeDOTEdge(out, midPointKey, ToDotKey{pd, prefix}, nd2name, style, false);
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
    std::string style = "fontcolor=\"";
    style.append(color);
    style.append("\" style=\"dashed\" color=\"");
    style.append(color);
    style.append("\"");
    // The midpoint nodes are unlabeled dots with IDs that are _addr or __addr
    const std::string pname ="_" + std::to_string((long)(&pp));
    const std::string sname ="__" + std::to_string((long)(&pp));
    const NamePair pv{pname, emptyStr};
    const NamePair sv{pname, emptyStr};
    const auto * pn = pp.phyloParent;
    const auto * pd = pp.phyloChild;
    const auto * sd = pp.scaffoldDes;
    const ToDotKey pk{pd, "_phpath"};
    writeOneSideOfPathPairingToDOT(out, pn, pd, sd, pk, pv, nd2name, color, prefix);
    const auto * sn = pp.scaffoldAnc;
    const ToDotKey sk{sd, "_scpath"};
    writeOneSideOfPathPairingToDOT(out, sn, sd, sd, sk, sv, nd2name, color, prefix);
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
        writeNodeDOT(out, n, nd2name, "", false, "", true, false);
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

