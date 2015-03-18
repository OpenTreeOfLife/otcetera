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
T findLeftmostInSubtreeM(T nd);
template<typename T>
T findRightmostInSubtreeM(T nd);

template<typename T>
inline bool isInternalNode(const T & nd) {
    return nd.isInternal();
}
template<typename T>
inline bool isLeaf(const T & nd) {
    return nd.isTip();
}

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

template<typename T>
inline std::pair<long, long> getMRCAOttIdPair(T nd) {
    assert(nd != nullptr);
    if (nd->isTip()) {
        return std::pair<long, long>(nd->getOttId(), nd->getOttId());
    }
    T f = findLeftmostInSubtree(nd);
    T r = nd->getLastChild();
    if (!r->isTip()) {
        r = findLeftmostInSubtree(r);
    }
    assert(f->isTip());
    assert(f->hasOttId());
    assert(r->isTip());
    assert(r->hasOttId());
    return std::pair<long, long>(f->getOttId(), r->getOttId());
}

template<typename T>
inline std::string getDesignator(const T &nd) {
    if (nd.hasOttId()) {
        std::string r = "ott";
        return r + std::to_string(nd.getOttId());
    }
    if (nd.isTip()) {
        return std::string("EMPTY_NODE");
    }
    auto p = getMRCAOttIdPair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}
template<typename T>
inline std::string getMRCADesignator(const T &nd) {
    if (nd.isTip()) {
        std::string r = "ott";
        return r + std::to_string(nd.getOttId());
    }
    auto p = getMRCAOttIdPair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}


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

