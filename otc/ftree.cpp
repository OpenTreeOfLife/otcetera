#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/debug.h"

#include "otc/write_dot.h"
namespace otc {

void PhyloStatement::write_as_newick(std::ofstream & out) const {
    if (!exclude_group.empty()) {
        out << '(';
    }
    out << "(";
    bool first = true;
    for (const auto & oid : include_group) {
        if (!first) {
            out << ',';
        }
        out << "ott" << oid;
        first = false;
    }
    out << ')';
    if (!exclude_group.empty()) {
        for (const auto & oid : exclude_group) {
            out << ",ott" << oid;
        }
        out << ")";
    }
    out << ';';
}

template<typename T>
bool ExcludeConstraints<T>::is_excluded_from(const node_type * ndToCheck,
                                           const node_type * potentialAttachment,
                                           const std::map<long, node_type*> * o2n) const {
    const auto & ndi = ndToCheck->get_data().desIds;
    if (ndi.size() == 1) {
        auto nit = by_exclude_node.find(ndToCheck);
        if (nit == by_exclude_node.end()) {
            return false;
        }
        return contains(nit->second, potentialAttachment);
    }
    assert(o2n != nullptr);
    for (auto oid : ndi) {
        auto tn = o2n->at(oid);
        auto nit = by_exclude_node.find(tn);
        if (nit != by_exclude_node.end() && contains(nit->second, potentialAttachment)) {
            return true;
        }
    }
    return false;
}

template<typename T>
bool ExcludeConstraints<T>::add_exclude_statement(const node_type * nd2Exclude,
                                                const node_type * forbiddenAttach) {
    if (is_excluded_from(nd2Exclude, forbiddenAttach, nullptr)) {
        return false; // already excluded from this subtree
    }
    // If any of the descendants exclude this node, we can remove those exclude statements,
    //  because they'll be "dominated by this one"
    std::list<node_pair> toRemove;
    auto pIt = by_exclude_node.find(nd2Exclude);
    if (pIt != by_exclude_node.end()) {
        for (const auto & cp : pIt->second) {
            if (isAncestorDesNoIter(forbiddenAttach, cp)) {
                toRemove.push_back(node_pair(nd2Exclude, cp));
            }
        }
    }
    for (auto & np : toRemove) {
        purge_exclude_raw(np.first, np.second);
    }
    ingest_exclude_raw(nd2Exclude, forbiddenAttach);
    return true;
}

template<typename T>
void ExcludeConstraints<T>::purge_exclude_raw(const node_type * nd2Exclude,
                                            const node_type * forbiddenAttach)  {
    cnode_set & bev = by_exclude_node.at(nd2Exclude);
    assert(contains(bev, forbiddenAttach));
    bev.erase(forbiddenAttach);
    if (bev.empty()) {
        by_exclude_node.erase(nd2Exclude);
    }
    cnode_set & bcv = by_node_with_constraints.at(forbiddenAttach);
    assert(contains(bcv, nd2Exclude));
    bcv.erase(nd2Exclude);
    if (bcv.empty()) {
        by_node_with_constraints.erase(forbiddenAttach);
    }
}

template<typename T>
void ExcludeConstraints<T>::ingest_exclude_raw(const node_type * nd2Exclude,
                                             const node_type * forbiddenAttach) {
    by_exclude_node[nd2Exclude].insert(forbiddenAttach);
    by_node_with_constraints[forbiddenAttach].insert(nd2Exclude);
}

template<typename T>
void InterTreeBandBookkeeping<T>::reassign_attachment_node(InterTreeBand<T> * b,
                                         RootedTreeNode<T> * oldAnc,
                                         RootedTreeNode<T> * newAnc,
                                         const PhyloStatement & ) {
    assert(b != nullptr);
    assert(oldAnc != nullptr);
    assert(newAnc != nullptr);
    band_to_node[b] = newAnc;
    node_to_band.at(oldAnc).erase(b);
    node_to_band[newAnc].insert(b);
    b->reassign_attachment_node(oldAnc, newAnc);
}

template<typename T, typename U>
void FTree<T, U>::update_to_reflect_resolution(node_type *oldAnc,
                                            node_type * newAnc,
                                            const std::set<node_type *> & movedNodes,
                                                            const PhyloStatement & ps) {
    std::set<node_type *> bn;
    for (auto np : movedNodes) {
        if (np->isTip()
            && np->has_ott_id()
            && contains(ps.include_group, np->get_ott_id())
            && (!ott_id_is_connected(np->get_ott_id()))) {
            bn.insert(np);
        }
    }
    if (bn.empty()) {
        return;
    }
    auto relevantBands = bands.get_bands_for_node(oldAnc);
    for (auto & b : relevantBands) {
        if (b->is_the_set_of_phantom_nodes(oldAnc, bn)) {
            bands.reassign_attachment_node(b, oldAnc, newAnc, ps);
        } else {
            auto nitbp = forest._createNewBand(*this, *newAnc, ps);
            assert(nitbp != nullptr);
            bands._add_ref_to_band(nitbp, newAnc);
            nitbp->insert_set(newAnc, bn);
        }
    }
}

bool PhyloStatement::debug_check() const {
#ifdef DEBUGGING_PHYLO_STATEMENTS
    const OttIdSet ie = set_union_as_set(include_group, exclude_group);
    if (ie != leaf_set) {
        dbWriteOttSet(" include_group ",include_group);
        dbWriteOttSet(" exclude_group ", exclude_group);
        dbWriteOttSet(" leaf_set ", leaf_set);
        assert(false);
    }
#endif
    return true;
}

template<typename T, typename U>
void FTree<T, U>::create_deeper_root() {
    auto nr = forest.create_node(nullptr, this);
    forest.add_and_update_child(nr, root, *this);
    root = nr;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T, U>::create_deeper_node(RootedTreeNode<T> *nd) {
    if (nd == root) {
        create_deeper_root();
        return root;
    }
    auto p = nd->getParent();
    assert(p != nullptr);
    // manually place the new node nn in the same spot as nd occupies
    auto nn = forest.create_node(nullptr, this);
    nd->addSibOnRight(nd);
    nd->detachThisNode();
    // we can add nd to nn in the normal way
    forest.add_and_update_child(nd, nn, *this);
    // but we placed nn in the exact spot that nd used to occupy, so
    //  we don't want to call an addChild... method.
    forest._and_update_child(nn, p, *this);
    return nn;
}

// profiling indicates that this is a big problem spot. We should cache this...
template<typename T, typename U>
const OttIdSet FTree<T, U>::get_connected_ott_ids() const {
    OttIdSet r;
    // root can be a tip but not a named node, in the process of stealing
    //  children from one tree in the merging of the forests
    if (root == nullptr || (root->isTip() && !root->has_ott_id())) {
        return r;
    }
    for (auto t : iter_leaf_n_const(*root)) {
        assert(isAncestorDesNoIter(root, t));
        if (t->has_ott_id()) {
            r.insert(t->get_ott_id());
        }
    }
    return r;
}



template<typename T, typename U>
bool FTree<T, U>::any_excluded_at_node(const node_type * nd, const OttIdSet &ott_id_set) const {
    for (auto oid : ott_id_set) {
        if (exclude.is_excluded_from(ott_id_to_node_map.at(oid), nd, nullptr)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool FTree<T, U>::any_included_at_node(const node_type * nd, const OttIdSet &ott_id_set) const {
    if (!areDisjoint(nd->get_data().desIds, ott_id_set)) {
        return true;
    }
    auto c = add_phantom_nodes_at_node(nd, ott_id_set);
    assert(c == false);// if we are correctly updating desIds we don't need this branch.... TMP
    return c;
}

template<typename T, typename U>
bool FTree<T, U>::add_phantom_nodes_at_node(const node_type * nd, const OttIdSet &ott_id_set) const {
    const auto & b = bands.get_bands_for_node(nd);
    if (!b.empty()) {
       const auto ns = ott_id_set_to_node_set(ott_id_set);
        for (auto ob : b) {
            if (ob->band_point_has_any(nd, ns)) {
                return true;
            }
        }
    }
    return false;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T, U>::add_leaf_no_des_update(RootedTreeNode<T> * par, long ottId) {
    //connectedIds.insert(ottId);
    return forest.create_leaf(par, ottId, this);
}


template<typename T, typename U>
RootedTreeNode<T> * FTree<T, U>::resolve_to_create_clade_of_included(RootedTreeNode<T> * par,
                                                               const PhyloStatement & ps) {
    const OttIdSet & oids = ps.include_group;
    dbWriteOttSet("  ott_id_is_connected oids = ", oids);
    dbWriteOttSet("                                 nd->get_data().desIds = ", par->get_data().desIds);
    std::set<RootedTreeNode<T> *> cToMove;
    std::list<RootedTreeNode<T> *> orderedToMove;
    for (auto oid : oids) {
        auto n = ott_id_to_node_map.at(oid);
        bool connectionFound = false;
        if (n->getParent() == par) {
            cToMove.insert(n);
            orderedToMove.push_back(n);
            connectionFound = true;
        } else {
            for (auto anc : iter_anc(*n)) {
                if (anc->getParent() == par) {
                    if (!contains(cToMove, anc)) {
                        cToMove.insert(anc);
                        orderedToMove.push_back(anc);
                        connectionFound = true;
                        break;
                    }
                }
            }
        }
        if (connectionFound) {
            continue;
        }
    }
    auto newNode = forest.create_node(par, this); // parent of include_group
    for (auto c : orderedToMove) {
        c->detachThisNode();
        forest.add_and_update_child(newNode, c, *this);
    }
    update_to_reflect_resolution(par, newNode, cToMove, ps);
    LOG(DEBUG) << "resolve_to_create_clade_of_included postcondition check";
    forest.debug_invariants_check();
    LOG(DEBUG) << "resolve_to_create_clade_of_included exiting";
    return newNode;
}

template<typename T, typename U>
bool FTree<T, U>::insert_into_band_no_des_update(InterTreeBand<T> * itbp,
                                    RootedTreeNode<T> * connectedNode,
                                    long phantomID) {
    assert(connectedNode != nullptr);
    assert(itbp != nullptr);
    auto noid = ott_id_to_node_map.at(phantomID);
    itbp->insert(connectedNode, noid);
    return true;
}

template<typename T>
T * rootToTipSearchByDesIds(T * nd, const OttIdSet &oids) {
    assert(nd);
    T * curr = nd;
    for (;;) {
        const auto & di = curr->get_data().desIds;
        assert(isSubset(oids, di));
        if (curr->isTip()) {
            return curr;
        }
        T * nextNd = nullptr;
        for (auto c : iter_child(*curr)) {
            const auto & chdi = c->get_data().desIds;
            if (!areDisjoint(oids, chdi)) {
                if (nextNd == nullptr) {
                    nextNd = &(*c);
                } else {
                    return curr;
                }
            }
        }
        if (nextNd == nullptr) {
            return curr;
        }
        curr = nextNd;
    }
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T, U>::get_mrca(const OttIdSet &ott_id_set) {
    if (ott_id_set.empty()) {
        assert(false);
        throw OTCError("empty MRCA");
    }
    check_all_node_pointers_iter(*root);
    const auto con = get_connected_ott_ids();
    const auto rel = set_intersection_as_set(ott_id_set, con);
    dbWriteOttSet(" get_mrca ingroup", ott_id_set);
    dbWriteOttSet(" get_mrca connected ", con);
    dbWriteOttSet(" get_mrca connected ingroup", rel);
    const auto & relCheck = root->get_data().desIds;
    dbWriteOttSet(" get_mrca relCheck", relCheck);
    assert(isSubset(rel, relCheck));
    for (auto nextOttId : rel) {
        auto x = ott_id_to_node_map.find(nextOttId);
        assert(x != ott_id_to_node_map.end());
        node_type * aTip = x->second;
        assert(forest.get_tree_for_node(aTip) == this);
        if (!isAncestorDesNoIter(root, aTip)) {
            LOG(ERROR) << "aTip->get_ott_id() = " << aTip->get_ott_id();
            for (auto a : iter_anc(*aTip)) {
                LOG(ERROR) << " anc address =  " << long(a);
            }
            LOG(ERROR) << " root address = " << long(root);
            assert(false);
        }
        assert(aTip != nullptr);
        if (ott_id_set.size() == 1) {
            return aTip;
        }
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    const auto referredTo = set_intersection_as_set(ott_id_set, root->get_data().desIds);
    assert(!referredTo.empty());
    return rootToTipSearchByDesIds(root, referredTo);
}

template<typename T, typename U>
void FTree<T, U>::mirror_phylo_statement(const PhyloStatement &ps) {
    assert(root == nullptr);
    root = forest.create_node(nullptr, this);
    add_phylo_statement_as_child_of_root(ps);
}

template<typename T, typename U>
void FTree<T, U>::steal_exclusion_statements(node_type * newPar,
                                           node_type * srcNode,
                                           FTree<T, U>  & donorTree) {
    auto e = donorTree.exclude.steal_exclusions(srcNode);
    for (auto nd : e) {
        this->exclude.add_exclude_statement(nd, newPar);
    }
}

template<typename T, typename U>
OttIdSet FTree<T, U>::steal_inclusion_statements(node_type * newPar,
                                               node_type * srcNode,
                                               FTree<T, U>  & donorTree,
                                               InterTreeBand<T> * bandToSkip) {
    assert(newPar != nullptr);
    auto bandSet = donorTree.bands.steal_bands(srcNode);
    OttIdSet r;
    for (auto bandPtr : bandSet) {
        if (bandToSkip == bandPtr) {
            continue;
        }
        const OttIdSet pis = bandPtr->get_phantom_ids(srcNode);
        r.insert(begin(pis), end(pis));
        bandPtr->reassign_attachment_node(srcNode, newPar);
        this->bands._add_ref_to_band(bandPtr, newPar);
    }
    return r;
}

// moves the exclusion statements to a different FTree (but with the same nodes)
template<typename T, typename U>
void FTree<T, U>::register_exclusion_statement_for_transferring_node(node_type * srcNode,
                                                                FTree<T, U>  & donorTree) {
    auto e = donorTree.exclude.steal_exclusions(srcNode);
    for (auto nd : e) {
        this->exclude.add_exclude_statement(nd, srcNode);
    }
}
// moves the inclusion statements to a different FTree (but with the same nodes)
template<typename T, typename U>
void FTree<T, U>::register_inclusion_statement_for_transferring_node(node_type * srcNode,
                                                                FTree<T, U>  & donorTree) {
    auto bandSet = donorTree.bands.steal_bands(srcNode);
    for (auto bandPtr : bandSet) {
        this->bands._add_ref_to_band(bandPtr, srcNode);
    }
}

template<typename T, typename U>
void FTree<T, U>::add_phylo_statement_as_child_of_root(const PhyloStatement &ps) {
    dbWriteOttSet(" add_phylo_statement_as_child_of_root", ps.include_group);
    assert(root != nullptr);
    if (!root->isTip()) {
        forest.debug_invariants_check();
    }
    if (any_excluded_at_node(root, ps.include_group)) {
        create_deeper_root();
        assert(!root->isTip());
    }
    assert(root != nullptr);
    auto parOfIncGroup = forest.create_node(root, this); // parent of include_group
    assert(ps.exclude_group.size() > 0);
    for (auto i : ps.exclude_group) {
        if (!forest.is_attached(i)) { // greedy
            add_leaf_no_des_update(root, i);
            root->get_data().desIds.insert(i);
        } else {
            add_exclude_statement(i, parOfIncGroup, ps.provenance);
        }
    }
    for (auto i : ps.include_group) {
        assert(!forest.is_attached(i));
        add_leaf_no_des_update(parOfIncGroup, i);
    }
    root->get_data().desIds.insert(begin(ps.include_group), end(ps.include_group));
    parOfIncGroup->get_data().desIds = ps.include_group;
    forest.debug_invariants_check();
    LOG(DEBUG) << "Leaving add_phylo_statement_as_child_of_root";
}

template<typename T>
std::set<T *> getAncSet(T *nd) {
    std::set<T *> r;
    T * p = nd->getParent();
    while (p != nullptr) {
        r.insert(p);
        p = p->getParent();
    }
    return r;
}


template<typename T, typename U>
OttIdSet FTree<T, U>::add_phylo_statement_at_node(const PhyloStatement & ps, 
                                             RootedTreeNode<T> * includeGroupA,
                                             const OttIdSet & attachedElsewhere,
                                             InterTreeBand<T> * itbp) {
    dbWriteOttSet(" FTree<T, U>::add_phylo_statement_at_node inc", ps.include_group);
    LOG(DEBUG) << "includeGroupA = " << (long)(includeGroupA)  << " " << std::hex << (long)(includeGroupA) << std::dec;
    LOG(DEBUG) << "itbp = " << (long)(itbp) << " " << std::hex << (long)(itbp) << std::dec;
    LOG(DEBUG) << "FTree = " << (long)(this) << " " << std::hex << (long)(this) << std::dec;
    dbWriteOttSet("    includeGroupA->get_data().desIds", includeGroupA->get_data().desIds);
    assert(forest.get_tree_for_node(includeGroupA) == this);
    OttIdSet r;
    if (itbp != nullptr) {
        bands._add_ref_to_band(itbp, includeGroupA);
    }
    for (auto oid : ps.include_group) {
        if (!ott_id_is_connected(oid)) {
            LOG(DEBUG) << " not connected " << oid;
            if (contains(attachedElsewhere, oid)) {
                assert(itbp != nullptr);
                LOG(DEBUG) << " attachedElsewhere " << oid;
                insert_into_band_no_des_update(itbp, includeGroupA, oid);
            } else {
                LOG(DEBUG) << " adding leaf " << oid;
                add_leaf_no_des_update(includeGroupA, oid);
                r.insert(oid);
            }
        } else {
            LOG(DEBUG) << " connected " << oid;
            assert(contains(includeGroupA->get_data().desIds, oid));
        }
    }
    add_des_ids_to_node_and_anc(includeGroupA, ps.include_group);
    dbWriteOttSet("    later includeGroupA->get_data().desIds", includeGroupA->get_data().desIds);
    for (auto oid : ps.exclude_group) {
        if (!ott_id_is_connected(oid)) {
            add_exclude_statement(oid, includeGroupA, ps.provenance);
        }
    }
    LOG(DEBUG) << "before add_phylo_statement_at_node exit";
    debug_invariants_check_ft();
    LOG(DEBUG) << "foreset level check";
    forest.debug_invariants_check();
    LOG(DEBUG) << "bout to add_phylo_statement_at_node exit";
    return r;
}
#if defined(DO_DEBUG_CHECKS)
template<typename T, typename U>
void FTree<T, U>::debug_verify_des_ids_assuming_des(const OttIdSet &s, const RootedTreeNode<T> *nd) const{
    OttIdSet ois;
    if (nd->isTip()) {
        if (nd->has_ott_id()) {
            ois.insert(nd->get_ott_id());
        }
    } else {
        for (auto c : iter_child_const(*nd)) {
            const auto & coids = c->get_data().desIds;
            assert(coids.find(LONG_MAX) == coids.end());
            ois.insert(begin(coids), end(coids));
        }
    }
    const auto pids = bands.get_phantom_ids(nd);
    assert(pids.find(LONG_MAX) == pids.end());
    ois.insert(pids.begin(), pids.end());
    if(s != ois) {
        LOG(DEBUG) << "FTree = " << (long)(this) << " " << std::hex << (long)(this) << std::dec;
        LOG(DEBUG) << "nd = " << (long)(nd) << " " << std::hex << (long)(nd) << std::dec;
        dbWriteNewick(nd);
        dbWriteOttSet("debug_verify_des_ids_assuming_des incoming", s);
        dbWriteOttSet("calculated:", ois);
        const auto extras = set_difference_as_set(s, ois);
        dbWriteOttSet("inc - calc:", extras);
        const auto missing = set_difference_as_set(ois, s);
        dbWriteOttSet("calc - inc:", missing);
        for (auto m : missing) {
            if (contains(pids, m)) {
                LOG(DEBUG) << m << " is from the phantomIDs";
            } else {
                LOG(DEBUG) << m << " is from a descendant ott_id_is_connected returns " << ott_id_is_connected(m);
                for (auto c : iter_child_const(*nd)) {
                    LOG(DEBUG) << "child " << getDesignator(*c);
                    dbWriteOttSet("   a child desIds", c->get_data().desIds);
                }
            }
        }
        assert(s == ois);
    }
}
template<typename T, typename U>
void FTree<T, U>::debug_invariants_check_ft() const {
    check_all_node_pointers_iter(*root);
    for (auto n : iter_post_n_const(*root)) {
         auto & b = bands.get_bands_for_node(n);
         if (!b.empty()) {
            for (auto ob : b) {
                assert(contains(ob->get_banded_nodes(), n));
            }
        }
        if(forest.get_tree_for_node(n) != this) {
            long x = (long)(forest.get_tree_for_node(n));
            long p = (long)(n->getParent());
            LOG(DEBUG) << getDesignator(*n) << " in the wrong tree reporting " << x << " " << std::hex << x << std::dec << " p = " << p << " " << std::hex << p << std::dec;
            if (n->getParent() != nullptr) {
                std::cerr << "the parent newick\n";
                dbWriteNewick(n->getParent());
                std::cerr << std::endl;
            }
            assert(false);
        }
        OttIdSet noids;
        if (n->isTip()) {
            if (n->has_ott_id()) {
                const auto o = n->get_ott_id();
                assert(ott_id_to_node_map.at(o) == n);
            }
            // Make sure that our ancestors do not exclude us.
            const std::set<const node_type *> ancSet = getAncSet(n);
            for (auto a : ancSet) {
                assert(!exclude.is_excluded_from(n, a, &ott_id_to_node_map));
            }
        } else {
            assert(!n->has_ott_id());
        }
        if (n != root) {
            assert(isAncestorDesNoIter(root, n));
        }
        debug_verify_des_ids_assuming_des(n->get_data().desIds, n);
    }
}
#endif

template class ExcludeConstraints<RTSplits>; // force explicit instantiaion of this template.
template class InterTreeBand<RTSplits>; // force explicit instantiaion of this template.
template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
