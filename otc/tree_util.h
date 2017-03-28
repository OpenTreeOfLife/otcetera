#ifndef OTCETERA_TREE_UTIL_H
#define OTCETERA_TREE_UTIL_H
// Very simple functions that operate on trees or treenode
//  without requiring iteration (contrast w/tree_operations.h)
// Depends on: tree.h 
// Depended on by: tree_iter.h tree_operation.h 
#include <deque>
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
template <typename N>
N * find_mrca_via_traversing(const std::set<N *> & tip_node_set, std::set<N*> * induced_nodes = nullptr);


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
inline std::pair<OttId, OttId> get_mrca_ott_id_pair(T nd) {
    assert(nd != nullptr);
    if (nd->is_tip()) {
        return std::pair<OttId, OttId>(nd->get_ott_id(), nd->get_ott_id());
    }
    T f = find_leftmost_in_subtree(nd);
    T r = nd->get_last_child();
    if (!r->is_tip()) {
        r = find_leftmost_in_subtree(r);
    }
    assert(f->is_tip());
    assert(r->is_tip());
    if (f->has_ott_id() && r->has_ott_id()) {
        return std::pair<OttId, OttId>(f->get_ott_id(), r->get_ott_id());
    }
    OttId fl = (f->has_ott_id() ? f->get_ott_id() : -1);
    OttId rl = (r->has_ott_id() ? r->get_ott_id() : -2);
    return std::pair<OttId, OttId>(fl, rl);
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
void emit_conflict_details(std::ostream & out, const T & ndRef, const OttIdSet & extras,  const OttIdSet & missing) {
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



// A slow method for finding the MRCA of all of the nodes pointed to by `tip_node_set`
//  If `induced_nodes` is provided, then on exit it will hold all of the nodes traversed
//    by walking from the tips to the returned ancestor.
// \returns nullptr if tip_node_set is empty
// Assumes that all nodes in tip_node_set are connected (in the same tree)
// Makes use of no info in the `data` field of the nodes
template <typename N>
inline N * find_mrca_via_traversing(const std::set<N *> & tip_node_set,
                                    std::set<N*> * induced_nodes) {
    if (tip_node_set.size() < 2)  {
        if (tip_node_set.empty()) {
            return nullptr;
        }
        N * r = *tip_node_set.begin();
        if (induced_nodes) {
            induced_nodes->insert(r);
        }
        return r;
    }
    // We walk rootward from each tip nodes until we intersect with a part
    //    of the tree we've visited (based on the node's presence in `seen_nodes`)
    auto nit = tip_node_set.begin();
    N * curr = *nit++;
    std::deque<N *> root_to_mrca;
    std::set<N *> local_induced_nodes;
    std::set<N *> * seen_nodes = (induced_nodes == nullptr ? &local_induced_nodes : induced_nodes);
    std::set<N *> dequed_nodes; // fast search for nodes in the deque
    // start by filling in the path from the first leaf to the root
    while (curr != nullptr) {
        root_to_mrca.push_front(curr);
        seen_nodes->insert(curr);
        dequed_nodes.insert(curr);
        curr = curr->get_parent();
    }
    // root_to_mrca has the root at the first pos and the leaf at the end
    // Now we walk through the other tips. Every time we hit a seen node, we 
    //    can stop retraversing.
    // If the seen node is in the deque, and it is not the last node in the deque
    //    then we've just walked from a tip that connects deeper, so we
    //    pop the shallower nodes off the end of the root_to_mrca deque.
    for (; nit != tip_node_set.end(); ++nit) {
        curr = *nit;
        while (curr != nullptr) {
            if (seen_nodes->count(curr) != 0) {
                if (dequed_nodes.count(curr) != 0) {
                    while (root_to_mrca.back() != curr) {
                        auto td = root_to_mrca.back();                
                        root_to_mrca.pop_back();
                        dequed_nodes.erase(td);
                        if (root_to_mrca.size() == 1) {
                            return root_to_mrca.back();
                        }
                    }
                }
                break;
            } else {
                seen_nodes->insert(curr);
            }
            curr = curr->get_parent();
        }
    }
    N * mrca = root_to_mrca.back();
    if (induced_nodes != nullptr) {
        // Now we need to remove any ancestors of the mrca from induced_nodes
        N * p = mrca->get_parent();
        while (p != nullptr) {
            if (induced_nodes->count(p) > 0) {
                induced_nodes->erase(p);
            } else {
                break;
            }
            p = p->get_parent();
        }
    }
    return mrca;
}

} // namespace otc
#endif

