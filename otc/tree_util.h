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
bool is_internal_node(const T & nd);
template<typename T>        
bool is_leaf(const T & nd);
template<typename T>
T find_leftmost_in_subtree(T nd);
template<typename T>
T find_rightmost_in_subtree(T nd);

// alias for nd.is_internal_node
template<typename T>
inline bool is_internal_node(const T & nd) {
    return nd.is_internal();
}

// alias for nd.is_tip
template<typename T>
inline bool is_leaf(const T & nd) {
    return nd.is_tip();
}

/// Walks through the tree (using get_first_child) to find the rightmost child
template<typename T>
inline T find_leftmost_in_subtree(T nd) {
    if (nd == nullptr) {
        return nullptr;
    }
    auto next = nd->get_first_child();
    while (next != nullptr) {
        nd = next;
        next = nd->get_first_child();
    }
    return nd;
}

/// Walks through the tree (using get_last_child) to find the rightmost child
template<typename T>
inline T find_rightmost_in_subtree(T nd) {
    if (nd == nullptr) {
        return nullptr;
    }
    auto next = nd->get_last_child();
    while (next != nullptr) {
        nd = next;
        next = nd->get_last_child();
    }
    return nd;
}

/// Returns a pair of OTT Ids corresponding to:
//      the "leftmost" descendant of `nd` and the leftmost descendant of the rightmost child of `nd`
//      OR the Id of `nd` twice (if `nd` is a tip)
template<typename T>
inline std::pair<long, long> get_mrca_ott_id_pair(T nd) {
    assert(nd != nullptr);
    if (nd->is_tip()) {
        return std::pair<long, long>(nd->get_ott_id(), nd->get_ott_id());
    }
    T f = find_leftmost_in_subtree(nd);
    T r = nd->get_last_child();
    if (!r->is_tip()) {
        r = find_leftmost_in_subtree(r);
    }
    assert(f->is_tip());
    assert(r->is_tip());
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
inline std::string get_designator(const T &nd) {
    if (nd.has_ott_id()) {
        std::string r = "ott";
        return r + std::to_string(nd.get_ott_id());
    }
    if (nd.is_tip()) {
        return std::string("EMPTY_NODE");
    }
    auto p = get_mrca_ott_id_pair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}

/// returns a string of "ott#" for a tip or "MRCA(ott#,ott#) for an internal"
template<typename T>
inline std::string get_mrca_designator(const T &nd) {
    if (nd.is_tip()) {
        std::string r = "ott";
        return r + std::to_string(nd.get_ott_id());
    }
    auto p = get_mrca_ott_id_pair(&nd);
    std::string pf = std::to_string(p.first);
    std::string ps = std::to_string(p.second);
    return std::string("MRCA(ott") + pf + std::string{", ott"} + ps + std::string{")"};
}

/// Writes a representation of the `extra` or `missing` IDs from the node `ndRef`. 
//  get_designator is used to refer to the node in the output stream
template<typename T>
void emit_conflict_details(std::ostream & out, const T & ndRef, const std::set<long> & extras,  const std::set<long> & missing) {
    out << "    split: " << get_designator(ndRef);
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

