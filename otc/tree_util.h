#ifndef OTCETERA_TREE_UTIL_H
#define OTCETERA_TREE_UTIL_H
// Very simple functions that operate on trees or treenode
//  without requiring iteration (contrast w/tree_operations.h)
// Depends on: tree.h 
// Depended on by: tree_iter.h tree_operation.h 

#include "otc/otc_base_includes.h"
#include "otc/tree.h"

namespace otc {

template<typename T>
bool isInternalNode(const T & nd);
template<typename T>
bool isLeaf(const T & nd);
template<typename T>
T findLeftmostInSubtree(T nd);
template<typename T>
T findRightmostInSubtree(T nd);

// alias for nd.isInternalNode
template<typename T>
inline bool isInternalNode(const T & nd) {
    return nd.isInternal();
}

// alias for nd.isTip
template<typename T>
inline bool isLeaf(const T & nd) {
    return nd.isTip();
}

/// Walks through the tree (using getFirstChild) to find the rightmost child
template<typename T>
inline T findLeftmostInSubtree(T nd) {
    if (nd == nullptr) {
        return nullptr;
    }
    auto next = nd->getFirstChild();
    while (next != nullptr) {
        nd = next;
        next = nd->getFirstChild();
    }
    return nd;
}

/// Walks through the tree (using getLastChild) to find the rightmost child
template<typename T>
inline T findRightmostInSubtree(T nd) {
    if (nd == nullptr) {
        return nullptr;
    }
    auto next = nd->getLastChild();
    while (next != nullptr) {
        nd = next;
        next = nd->getLastChild();
    }
    return nd;
}

/// Returns a pair of OTT Ids corresponding to:
//      the "leftmost" descendant of `nd` and the leftmost descendant of the rightmost child of `nd`
//      OR the Id of `nd` twice (if `nd` is a tip)
template<typename T>
inline std::pair<long, long> getMRCAOttIdPair(T nd) {
    assert(nd != nullptr);
    if (nd->isTip()) {
        return std::pair<long, long>(nd->get_ott_id(), nd->get_ott_id());
    }
    T f = findLeftmostInSubtree(nd);
    T r = nd->getLastChild();
    if (!r->isTip()) {
        r = findLeftmostInSubtree(r);
    }
    assert(f->isTip());
    assert(r->isTip());
    if (f->has_ott_id() && r->has_ott_id()) {
        return std::pair<long, long>(f->get_ott_id(), r->get_ott_id());
    }
    long fl = (f->has_ott_id() ? f->get_ott_id() : -1);
    long rl = (r->has_ott_id() ? r->get_ott_id() : -2);
    return std::pair<long, long>(fl, rl);
}

/// returns a string of "ott#" for any node with an OTT ID, "EMPTY_NODE" for a tip without an ID
//      or "MRCA(ott#,ott#) for an internal"
template<typename T>
inline std::string getDesignator(const T &nd) {
    if (nd.has_ott_id()) {
        std::string r = "ott";
        return r + std::to_string(nd.get_ott_id());
    }
    if (nd.isTip()) {
        return std::string("EMPTY_NODE");
    }
    auto p = getMRCAOttIdPair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}

/// returns a string of "ott#" for a tip or "MRCA(ott#,ott#) for an internal"
template<typename T>
inline std::string getMRCADesignator(const T &nd) {
    if (nd.isTip()) {
        std::string r = "ott";
        return r + std::to_string(nd.get_ott_id());
    }
    auto p = getMRCAOttIdPair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}

/// Writes a representation of the `extra` or `missing` IDs from the node `ndRef`. 
//  getDesignator is used to refer to the node in the output stream
template<typename T>
void emitConflictDetails(std::ostream & out, const T & ndRef, const std::set<long> & extras,  const std::set<long> & missing) {
    out << "    split: " << getDesignator(ndRef);
    out << ";    extras in phylo: ";
    for (auto o : extras) {
        out << o << ' ';
    }
    out << "    missing in phylo: ";
    for (auto o : missing) {
        out << o << ' ';
    }
    out << "\n";
}

} // namespace otc
#endif

