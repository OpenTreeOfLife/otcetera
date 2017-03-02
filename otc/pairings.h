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
        1. tips of the tree are mapped to the scaffoldNode assigned the same
            OTT Id,
        2. internal nodes in phylo tree are mapped to the least inclusive node
            in the scaffold tree that has all of the descendant OTT Ids from
            the phylo tree. (so the phyloNode->get_data().desIds will be a subset
            of scaffoldNode->get_data().desIds)
*/
template<typename T, typename U>
class NodePairing {
    public:
    T * scaffoldNode;
    U * phyloNode;
    NodePairing(T *taxo, U *phylo)
        :scaffoldNode(taxo),
        phyloNode(phylo) {
        assert(taxo != nullptr);
        assert(phylo != nullptr);
    }
};

/* Represents the mapping of an edge from phyloParent -> phyloChild onto
    a scaffold tree. The endpoints will be pairs of nodes that were aligned.
    Note that the phyloChild is a child of phyloParent, but scaffoldDes can
    be the same node as scaffoldAnc or it can be any descendant of that node.
    Thus an directed edge in the phylo tree pairs with a directed path in the
    scaffold.
*/
template<typename T, typename U>
class PathPairing {
    public:
    T * scaffoldDes;
    T * scaffoldAnc;
    U * const phyloChild;
    U * phyloParent; //@TMP used to be const, but the resolving step edits
                     // rather than creating a new path pairing... 
                     // that is probably the best way to do this (const was
                     // probably to retrictive). But more careful consideration may be needed
    OttIdSet currChildOttIdSet;
    bool path_is_now_trivial() {
        return currChildOttIdSet.size() == 1;
    }
    void set_ott_id_set(long oid,
                     std::map<const U *, NodeEmbedding<T, U> > & m) {
        if (currChildOttIdSet.size() == 1 && *currChildOttIdSet.begin() == oid) {
            return;
        }
        LOG(DEBUG) << "set_ott_id_set to " << oid << "  for path " << reinterpret_cast<long>(this);
        dbWriteOttSet("prev currChildOttIdSet = ", currChildOttIdSet);
        OttIdSet n;
        OttIdSet oldIds;
        std::swap(oldIds, currChildOttIdSet);
        if (contains(oldIds, oid)) {
            oldIds.erase(oid);
        }
        n.insert(oid);
        update_des_ids_for_self_and_anc(oldIds, n, m);
    }
    void update_des_ids_for_self_and_anc(const OttIdSet & oldIds,
                                   const OttIdSet & newIds,
                                   std::map<const U *, NodeEmbedding<T, U> > & m) {
        update_ancestral_path_ott_id_set(scaffoldDes, oldIds, newIds, m);
        currChildOttIdSet = newIds;
        dbWriteOttSet(" update_des_ids_for_self_and_anc onExit currChildOttIdSet = ", currChildOttIdSet);
    }
    bool update_ott_id_set_no_traversal(const OttIdSet & oldEls, const OttIdSet & newEls);
    PathPairing(const NodePairing<T, U> & parent,
                const NodePairing<T, U> & child)
        :scaffoldDes(child.scaffoldNode),
        scaffoldAnc(parent.scaffoldNode),
        phyloChild(child.phyloNode),
        phyloParent(parent.phyloNode),
        currChildOttIdSet(child.phyloNode->get_data().desIds) {
        assert(phyloChild->getParent() == phyloParent);
        assert(scaffoldAnc == scaffoldDes || isAncestorDesNoIter(scaffoldAnc, scaffoldDes));
    }
    PathPairing(T * scafPar,
                U * phyPar,
                const NodePairing<T, U> & child)
        :scaffoldDes(child.scaffoldNode),
        scaffoldAnc(scafPar),
        phyloChild(child.phyloNode),
        phyloParent(phyPar),
        currChildOttIdSet(child.phyloNode->get_data().desIds) {
        assert(phyloChild->getParent() == phyloParent);
        assert(scaffoldAnc == scaffoldDes || isAncestorDesNoIter(scaffoldAnc, scaffoldDes));
    }
    // as Paths get paired back deeper in the tree, the ID may be mapped to a higher
    // taxon. The currChildOttIdSet starts out identical to the phylogenetic node's 
    // descendant Ott Id set. But may change to reflect this remapping to the effective
    // set of IDs that include the tip.
    const OttIdSet & get_ott_id_set() const {
        return currChildOttIdSet;
    }
    const OttIdSet & get_phylo_child_des_id() const {
        return phyloChild->get_data().desIds;
    }
};


} // namespace
#endif
