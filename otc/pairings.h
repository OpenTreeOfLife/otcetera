#ifndef OTCETERA_PAIRINGS_H
#define OTCETERA_PAIRINGS_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
namespace otc {
template<typename T, typename U> class NodeEmbedding;

template<typename T, typename U>
void update_ancestral_path_ott_id_set(T * nd,
                                 const OttIdSet & oldEls,
                                 const OttIdSet & newEls,
                                 std::map<const T *, NodeEmbedding<T, U> > & m);

/* a pair of aligned nodes from an embedding of a phylogeny onto a scaffold
   In NodeEmbedding objects two forms of these pairings are created:
        1. tips of the tree are mapped to the scaffold_node assigned the same
            OTT Id,
        2. internal nodes in phylo tree are mapped to the least inclusive node
            in the scaffold tree that has all of the descendant OTT Ids from
            the phylo tree. (so the phylo_node->get_data().des_ids will be a subset
            of scaffold_node->get_data().des_ids)
*/
template<typename T, typename U>
class NodePairing {
    public:
    T * scaffold_node;
    U * phylo_node;
    NodePairing(T *taxo, U *phylo)
        :scaffold_node(taxo),
        phylo_node(phylo) {
        assert(taxo != nullptr);
        assert(phylo != nullptr);
    }
};

/* Represents the mapping of an edge from phylo_parent -> phylo_child onto
    a scaffold tree. The endpoints will be pairs of nodes that were aligned.
    Note that the phylo_child is a child of phylo_parent, but scaffold_des can
    be the same node as scaffold_anc or it can be any descendant of that node.
    Thus an directed edge in the phylo tree pairs with a directed path in the
    scaffold.
*/
template<typename T, typename U>
class PathPairing {
    public:
    T * scaffold_des;
    T * scaffold_anc;
    U * const phylo_child;
    U * phylo_parent; //@TMP used to be const, but the resolving step edits
                     // rather than creating a new path pairing... 
                     // that is probably the best way to do this (const was
                     // probably to retrictive). But more careful consideration may be needed
    OttIdSet curr_child_ott_id_set;
    bool path_is_now_trivial() {
        return curr_child_ott_id_set.size() == 1;
    }
    void set_ott_id_set(long oid,
                     std::map<const U *, NodeEmbedding<T, U> > & m) {
        if (curr_child_ott_id_set.size() == 1 && *curr_child_ott_id_set.begin() == oid) {
            return;
        }
        LOG(DEBUG) << "set_ott_id_set to " << oid << "  for path " << reinterpret_cast<long>(this);
        db_write_ott_id_set("prev curr_child_ott_id_set = ", curr_child_ott_id_set);
        OttIdSet n;
        OttIdSet oldIds;
        std::swap(oldIds, curr_child_ott_id_set);
        if (contains(oldIds, oid)) {
            oldIds.erase(oid);
        }
        n.insert(oid);
        update_des_ids_for_self_and_anc(oldIds, n, m);
    }
    void update_des_ids_for_self_and_anc(const OttIdSet & oldIds,
                                   const OttIdSet & newIds,
                                   std::map<const U *, NodeEmbedding<T, U> > & m) {
        update_ancestral_path_ott_id_set(scaffold_des, oldIds, newIds, m);
        curr_child_ott_id_set = newIds;
        db_write_ott_id_set(" update_des_ids_for_self_and_anc onExit curr_child_ott_id_set = ", curr_child_ott_id_set);
    }
    bool update_ott_id_set_no_traversal(const OttIdSet & oldEls, const OttIdSet & newEls);
    PathPairing(const NodePairing<T, U> & parent,
                const NodePairing<T, U> & child)
        :scaffold_des(child.scaffold_node),
        scaffold_anc(parent.scaffold_node),
        phylo_child(child.phylo_node),
        phylo_parent(parent.phylo_node),
        curr_child_ott_id_set(child.phylo_node->get_data().des_ids) {
        assert(phylo_child->get_parent() == phylo_parent);
        assert(scaffold_anc == scaffold_des || is_ancestor_des_no_iter(scaffold_anc, scaffold_des));
    }
    PathPairing(T * scafPar,
                U * phyPar,
                const NodePairing<T, U> & child)
        :scaffold_des(child.scaffold_node),
        scaffold_anc(scafPar),
        phylo_child(child.phylo_node),
        phylo_parent(phyPar),
        curr_child_ott_id_set(child.phylo_node->get_data().des_ids) {
        assert(phylo_child->get_parent() == phylo_parent);
        assert(scaffold_anc == scaffold_des || is_ancestor_des_no_iter(scaffold_anc, scaffold_des));
    }
    // as Paths get paired back deeper in the tree, the ID may be mapped to a higher
    // taxon. The curr_child_ott_id_set starts out identical to the phylogenetic node's 
    // descendant Ott Id set. But may change to reflect this remapping to the effective
    // set of IDs that include the tip.
    const OttIdSet & get_ott_id_set() const {
        return curr_child_ott_id_set;
    }
    const OttIdSet & get_phylo_child_des_id() const {
        return phylo_child->get_data().des_ids;
    }
};


} // namespace
#endif
