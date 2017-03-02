#include "otc/greedy_forest.h"
#include "otc/node_embedding.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/write_dot.h"
#include "debug.h"
namespace otc {

template<typename T, typename U>
bool GreedyBandedForest<T, U>::attempt_to_add_grouping(const OttIdSet & incGroup,
                                                        const OttIdSet & leaf_set,
                                                        int treeIndex,
                                                        long groupIndex,
                                                        SupertreeContextWithSplits *sc) {
    if (incGroup.size() == 1) {
        add_leaf(incGroup, leaf_set, treeIndex, groupIndex, sc);
        return true;
    }
    return create_and_add_phylo_statement(incGroup, leaf_set, treeIndex, groupIndex);
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::add_leaf(
                      const OttIdSet & incGroup,
                      const OttIdSet & ,
                      int ,
                      long ,
                      SupertreeContextWithSplits *) {
    assert(incGroup.size() == 1);
    const auto ottId = *incGroup.begin();
    register_leaf(ottId);
    return true;
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::create_and_add_phylo_statement(
                const OttIdSet & incGroup,
                const OttIdSet & leaf_set,
                int treeIndex,
                long groupIndex) {
    const auto & i = *(encountered.insert(incGroup).first);
    const auto & ls = (leaf_set.empty() ? i : *(encountered.insert(leaf_set).first));
    const OttIdSet exc = set_difference_as_set(ls, i);
    const auto & e = *(encountered.insert(exc).first);
    PhyloStatement ps(i, e, ls, PhyloStatementSource(treeIndex, groupIndex));
    LOG(DEBUG) << " GPF calling add_phylo_statement for tree=" << treeIndex << " group=" << groupIndex;
    return add_phylo_statement(ps);
}


// 1. finalize_tree
// 2. copy the structure into the scaffold tree
// 3. update all of the outgoing paths so that they map
//      to this taxon
template<typename T, typename U>
void GreedyBandedForest<T, U>::finish_resolution_of_embedded_clade(U & scaffoldNode,
                                                                    NodeEmbedding<T, U> * embedding,
                                                                    SupertreeContextWithSplits * sc) {
    assert(sc != nullptr);
    const auto snoid = scaffoldNode.get_ott_id();
    LOG(DEBUG) << "finish_resolution_of_embedded_clade for " << snoid;
    debug_invariants_check();
    finalize_tree(sc);
    debug_invariants_check();
    assert(trees.size() == 1);
    auto & resolvedTree = begin(trees)->second;
    const auto beforePar = scaffoldNode.getParent();
    check_all_node_pointers_iter(scaffoldNode);
    copyStructureToResolvePolytomy(resolvedTree.get_root(), sc->scaffoldTree, &scaffoldNode, sc);
    check_all_node_pointers_iter(scaffoldNode);
    assert(beforePar == scaffoldNode.getParent());
    // merge multiple exit paths (if needed) and remap all path pairings out of this node...
    embedding->merge_exit_embeddings_if_multiple();
    embedding->set_ott_id_for_exit_embeddings(&scaffoldNode, snoid, sc->scaffold2NodeEmbedding);
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::finalize_tree(SupertreeContextWithSplits *sc) {
    LOG(DEBUG) << "finalize_tree for a forest with " << trees.size() << " roots:";
    debug_invariants_check();
    auto roots = get_roots();
    for (auto r : roots) {
        LOG(DEBUG) << " tree-in-forest = "; dbWriteNewick(r);
    }
    if (trees.size() > 1) {
        LOG(WARNING) << "finalize_tree is not well thought out. merging of multiple trees is questionable.";
        const char * dbfn = "real-forest-in-finalize_tree.dot";
        LOG(WARNING) << "  should be writing DOT to " << dbfn;
        //std::ofstream outf(dbfn);
        //writeDOTForest(outf, *this);
        //outf.close();
        merge_forest(sc);
        LOG(WARNING) << "finished questionable merge_forest.";
    }
    debug_invariants_check();
    roots = get_roots();
    assert(roots.size() < 2);
    attach_all_detached_tips();
    debug_invariants_check();
    roots = get_roots();
    assert(roots.size() == 1);
    auto onlyRoot = *roots.begin();
    onlyRoot->get_data().desIds = ott_id_set;
    LOG(DEBUG)<< " finalized-tree-from-forest = "; dbWriteNewick(onlyRoot);
}


template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly,
                                    SupertreeContextWithSplits * sc) {
    assert(sc != nullptr);
    std::map<const T *, typename U::node_type *> gpf2scaff;
    std::map<long, typename U::node_type *> & dOttIdToNode = destTree.get_data().ottIdToNode;
    LOG(DEBUG) << " adding " << srcPoly;
    LOG(DEBUG) << " copying structure to resolve " << destPoly->get_ott_id();
    gpf2scaff[srcPoly] = destPoly;
    for (auto sn : iter_pre_n_const(srcPoly)) {
        if (sn == srcPoly) {
            continue;
        }
        auto sp = sn->getParent();
        auto dp = gpf2scaff.at(sp);
        typename U::node_type * dn;
        if (sn->has_ott_id() && sn->get_ott_id() > 0) { // might assign negative number to nodes created in synth...
        auto nid = sn->get_ott_id();
            //LOG(DEBUG) << " node in src has ID " << nid;
            dn = dOttIdToNode.at(nid);
            //LOG(DEBUG) << " in dest node, that ID maps to a node with id:  " << dn->get_ott_id();
            assert(dn != destPoly);
            if (contains(sc->detachedScaffoldNodes, dn)) {
                assert(not dn->getFirstChild());
                assert(not dn->getNextSib());
                sc->detachedScaffoldNodes.erase(dn);
                sc->scaffoldTree.markAsAttached(dn);
            }
            if (dn->getParent() != dp) {
                if (dn->getParent() != nullptr) {
                    dn->detachThisNode();
                    sc->scaffoldTree.markAsDetached(dn);
                }
                assert(dn->getNextSib() == nullptr);
                dp->addChild(dn);
            }
        } else {
            dn = destTree.create_child(dp);
        if (sn->has_ott_id())
          dn->set_ott_id(sn->get_ott_id());
        }
        //LOG(DEBUG) << " adding " << sn;
        gpf2scaff[sn] = dn;
    }
}

template<typename T, typename U>
std::vector<T *> GreedyBandedForest<T, U>::get_roots(){
    std::vector<T *> r;
    r.reserve(trees.size());
    for (auto & t : trees) {
        r.push_back(t.second.get_root());
    }
    return r;
}

template<typename T, typename U>
std::set<InterTreeBand<typename T::data_type> * > collectBandsForSubtree(U & tree, T * node) {
    std::set<InterTreeBand<typename T::data_type> *> r;
    for (auto nd : iter_child(*node)) {
        const auto & bs = tree.get_bands_for_node(&(*nd));
        for (auto & b : bs) {
            if (!b->is_single_tree_band()) {
                r.insert(b);
            }
        }
    }
    return r;
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::merge_forest(SupertreeContextWithSplits *sc) {
    if (trees.size() == 1) {
        return;
    }
    using FTreeType = FTree<RTSplits, MappedWithSplitsData>;
    std::vector<FTreeType *> sortedTrees;
    sortedTrees.reserve(trees.size());
    for (auto & t : trees) {
        FTreeType & tre = t.second;
        sortedTrees.push_back(&tre);
    }
    const std::size_t nTrees = sortedTrees.size();
    bool hasLeafInit = true;
    std::vector<bool> stillHasLeaves (nTrees, hasLeafInit);
    for (auto treeInd = sortedTrees.size() - 1 ; treeInd > 0; --treeInd) {
        FTreeType & toDie = *sortedTrees.at(treeInd);
        std::set<InterTreeBand<typename T::data_type> *> itbSet = collectBandsForSubtree(toDie, toDie.get_root());
        if (!itbSet.empty()) {
            if (itbSet.size() == 1) {
                stillHasLeaves[treeInd] = perform_single_band_merge(treeInd, *itbSet.begin(), sortedTrees, sc);
            } else {
                stillHasLeaves[treeInd] = perform_set_of_single_band_merges(treeInd, itbSet, sortedTrees, sc);
            }
        }
    }
    // clean up all trees with no leaves...
    assert(stillHasLeaves[0]);
    
    auto trIt = begin(trees);
    for (unsigned i = 0; trIt != end(trees); ++i) {
        if (stillHasLeaves.at(i)) {
            ++trIt;
        } else {
            trIt = trees.erase(trIt);
        }
    }
    merge_trees_to_first_post_band_handling(sc);
    LOG(DEBUG) << "exiting merge_forest";
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::perform_single_band_merge(
            std::size_t treeInd,
            InterTreeBand<RTSplits> * itb,
            const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
            SupertreeContextWithSplits *sc) {
    assert(itb != nullptr);
    auto btm = get_tree_to_node_map_for_band(*itb);
    FTree<RTSplits, MappedWithSplitsData> * toMergeTo = nullptr;
    for (auto j = 0U; j < treeInd; ++j) {
        auto btp = sortedTrees.at(j);
        if (contains(btm, btp)) {
            toMergeTo = btp;
            break;
        }
    }
    assert(toMergeTo != nullptr);
    FTree<RTSplits, MappedWithSplitsData> & toDie = *sortedTrees.at(treeInd);
    return merge_single_banded_tree(toDie, itb, *toMergeTo, sc);
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::perform_set_of_single_band_merges(
            std::size_t treeInd,
            std::set<InterTreeBand<typename T::data_type> *> & itbSet,
            const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
            SupertreeContextWithSplits *sc) {
    LOG(WARNING) << itbSet.size() <<  " bands for a tree. It is not a great idea to merge these one at at time...";
    LOG(DEBUG) << itbSet.size() <<  " bands for a tree. It is not a great idea to merge these one at at time...";
    //            write_forest_dot_to_fn("writingForestMerge.dot");
    //            assert("not implemented"[0] == 'f');;
    std::vector<InterTreeBand<typename T::data_type> * > postOrd;
    postOrd.reserve(itbSet.size());
    FTree<RTSplits, MappedWithSplitsData> & toDie = *sortedTrees.at(treeInd);
    std::set<InterTreeBand<typename T::data_type> *> toOrganize = itbSet;
    for (auto nd : iter_post(toDie)) {
        if (toOrganize.empty()) {
            break;
        }
        auto itbIt = begin(toOrganize);
        for (; itbIt != end(toOrganize);) {
            InterTreeBand<typename T::data_type> * ip = *itbIt;
            if (ip->is_a_banded_node_in_this(nd)) {
                postOrd.push_back(ip);
                itbIt = toOrganize.erase(itbIt);
            } else {
                ++itbIt;
            }
        }
    }
    assert(toOrganize.empty());
    for (auto sit : postOrd) {
        if (!perform_single_band_merge(treeInd, sit, sortedTrees, sc)) {
            return false;
        }
    }
    return true;
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::merge_trees_to_first_post_band_handling(SupertreeContextWithSplits *) {
    if (trees.size() == 1) {
        return;
    }
    assert(trees.size() > 0);
    auto trIt = begin(trees);
    auto & firstTree = trIt->second;
    auto firstTreeRoot = firstTree.get_root();
    OttIdSet idsIncluded = firstTreeRoot->get_data().desIds;
    for (++trIt; trIt != end(trees);) {
        debug_invariants_check();
        auto & currTree = trIt->second;
        auto currTreeRoot = currTree.get_root();
        assert(currTreeRoot->getParent() == nullptr);
        assert(!currTreeRoot->isTip());
        auto p = move_all_children(currTreeRoot, currTree, firstTreeRoot, firstTree, nullptr);
        firstTreeRoot = p.first;
        currTreeRoot = currTree.get_root();
        for (auto j : node_to_tree) {
            assert((j.first == currTreeRoot || j.second != &currTree));
        }
        LOG(DEBUG) << "Before erase";
        debug_invariants_check();
        trIt = trees.erase(trIt);
        register_tree_for_node(currTreeRoot, nullptr);
        LOG(DEBUG) << "After erase";
        debug_invariants_check();
        LOG(DEBUG) << "After erase check";
    }
    LOG(DEBUG) << "check before exit of merge_trees_to_first_post_band_handling";
    debug_invariants_check();
    LOG(DEBUG) << "exiting merge_trees_to_first_post_band_handling";
}


template<typename T, typename U>
bool GreedyBandedForest<T, U>::merge_single_banded_tree(
            FTree<RTSplits, MappedWithSplitsData> &donorTree,
            InterTreeBand<RTSplits> * band,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree,
            SupertreeContextWithSplits *sc) {
    LOG(DEBUG) << "pre merge_single_banded_tree check";
    debug_invariants_check();
    const auto t2n = get_tree_to_node_map_for_band(*band);
    auto dn = const_cast<NodeWithSplits *>(t2n.at(&donorTree));
    auto dp = dn->getParent();
    auto rn = const_cast<NodeWithSplits *>(t2n.at(&recipientTree));
    const auto & beforePhantomsRN = band->get_phantom_nodes(rn);
    const auto & beforePhantomsDN = band->get_phantom_nodes(dn);
    // any phantom nodes in common between dn and rn must be attached
    //  to a different tree. So they will remain "phantom" wrt rn
    auto movedNodeSet = set_difference_as_set(beforePhantomsRN, beforePhantomsDN);
    auto p = move_all_children(dn, donorTree, rn, recipientTree, band);
    OttIdSet movedIDSet;
    for (auto mel : movedNodeSet) {
        movedIDSet.insert(mel->get_ott_id());
    }
    for (auto mel : beforePhantomsDN) {
        movedIDSet.insert(mel->get_ott_id());
    }
    removeDesIdsToNdAndAnc(dn, movedIDSet);
    dbWriteOttSet(" movedIDSet =", movedIDSet);
    rn = p.first;
    band->remove_node(dn);
    band->remove_from_set(rn, movedNodeSet);
    if (p.second != rn) {
        assert("not implemented"[0] == 'f');;
    }
    auto anc = find_grandparent_that_is_root_or_band_sharing(donorTree, dn, recipientTree);
    bool r = true;
    if (anc.first == nullptr) {
        // dn is a child of the root
        register_tree_for_node(dn, nullptr);
        if (dp == nullptr) {
            donorTree._set_root(nullptr);
            r = false;
        } else {
            dn->detachThisNode();
        }
    } else {
        r = zip_paths_from_barren_node(donorTree,
                                   dn,
                                   anc.first,
                                   recipientTree,
                                   p.first,
                                   anc.second,
                                   sc);
    }
    LOG(DEBUG) << "post merge_single_banded_tree check";
    debug_invariants_check();
    LOG(DEBUG) << "  exiting merge_single_banded_tree";
    return r;
}

template<typename T, typename U>
std::pair<NodeWithSplits *, NodeWithSplits*> GreedyBandedForest<T, U>::move_all_children(NodeWithSplits * donorParent,
                                                     FTree<RTSplits, MappedWithSplitsData> &donorTree,
                                                     NodeWithSplits * recipientNode,
                                                     FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                                                     InterTreeBand<RTSplits> * bandBeingMerged) {
    if (recipientTree.is_excluded_from(donorParent, recipientNode)) {
        if (recipientNode == recipientTree.get_root()) {
            recipientTree.create_deeper_root();
            recipientNode = recipientTree.get_root();
        } else {
            recipientNode = recipientTree.create_deeper_node(recipientNode);
        }
        debug_invariants_check();
    }
    auto attachmentPoint = recipientNode;
    if (donorTree.is_excluded_from(recipientNode, donorParent)) {
        attachmentPoint = create_node(recipientNode, &recipientTree);
        debug_invariants_check();
    }
    std::list<node_type *> rc;
    for (auto currChild : iter_child(*donorParent)) {
        rc.push_back(currChild);
    }
    for (auto currChild : rc) {
        if (bandBeingMerged == nullptr) {
            debug_invariants_check();
        }
        transfer_subtree_in_forest(currChild, donorTree, attachmentPoint, recipientTree, &donorTree, bandBeingMerged);
        LOG(DEBUG) << "Back from transfer_subtree_in_forest";
    }
    assert(not donorParent->getFirstChild());
    return std::pair<NodeWithSplits *, NodeWithSplits*>{recipientNode, attachmentPoint};
}

// If the parent of a banded node that was just removed is also banded to
//  the same tree, or is the root, then the subsequent handling of nodes attached
//  to a band or the root should work.
// But if there are intervening unbannded nodes, we want to zip the path from
//  nd to the next relevant ancestor onto the existing path...
template<typename T, typename U>
std::pair<NodeWithSplits *, NodeWithSplits *>
GreedyBandedForest<T, U>::find_grandparent_that_is_root_or_band_sharing(
            FTree<RTSplits, MappedWithSplitsData> & donorTree,
            NodeWithSplits * nd,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree) {
    assert(nd != nullptr);
    std::pair<NodeWithSplits *, NodeWithSplits *> r{nullptr, nullptr};
    auto anc = nd->getParent();
    if (anc == nullptr || anc->getParent() == nullptr) {
        return r;
    }
    for (auto a : iter_anc(*nd)) {
        if (a->getParent() == nullptr) {
            r.first = a;
            r.second = recipientTree.get_root();
            return r;
        }
        const auto & bfn = donorTree.get_bands_for_node(a);
        for (auto & b : bfn) {
            const auto m = get_tree_to_node_map_for_band(*b);
            if (contains(m, &recipientTree)) {
                r.first = a;
                r.second = const_cast<NodeWithSplits *>(m.at(&recipientTree));
                return r;
            }
        }
    }
    OTC_UNREACHABLE;
}


template<typename T, typename U>
void GreedyBandedForest<T, U>::transfer_subtree_in_forest(
                NodeWithSplits * des,
                FTree<RTSplits, MappedWithSplitsData> & possDonor,
                NodeWithSplits * newPar,
                FTree<RTSplits, MappedWithSplitsData> & recipientTree, 
                FTree<RTSplits, MappedWithSplitsData> *donorTree,
                InterTreeBand<RTSplits> * bandBeingMerged) {
    assert(des != nullptr);
    auto oldPar = des->getParent();
    //LOG(DEBUG) << "top transfer_subtree_in_forest pre des == " << getDesignator(*des) << " " << (long) des << " " << std::hex << (long) des << std::dec;
    //LOG(DEBUG) << "                            newPar " << (long) newPar << " " << std::hex << (long) newPar << std::dec;
    //LOG(DEBUG) << "                            oldPar " << (long) oldPar << " " << std::hex << (long) oldPar << std::dec;
    //LOG(DEBUG) << "    newick of newPar before actions of transfer_subtree_in_forest";
    //dbWriteNewick(newPar);
    //LOG(DEBUG) << "    newick of oldPar before actions of transfer_subtree_in_forest";
    //dbWriteNewick(oldPar);
    //dbWriteOttSet("    on entry des->desIds", des->get_data().desIds);
    //dbWriteOttSet("    on entry newPar->desIds", newPar->get_data().desIds);
    //dbWriteOttSet("    on entry oldPar->desIds", oldPar->get_data().desIds);
    if (bandBeingMerged == nullptr) {
        debug_invariants_check();
    }
    assert(des != nullptr);
    assert(newPar != nullptr);
    if (donorTree == nullptr) {
        if (isAncestorDesNoIter(possDonor.get_root(), des)) {
            donorTree = &possDonor;
        } else {
            assert(isAncestorDesNoIter(recipientTree.get_root(), des));
            donorTree = &recipientTree;
        }
    }
    assert(!recipientTree.is_excluded_from(des, newPar));
    assert(get_tree_for_node(des) == donorTree);
    des->detachThisNode();
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" des pre add_and_update_child", des->get_data().desIds);
        dbWriteOttSet(" newPar pre add_and_update_child", newPar->get_data().desIds);
    }
    add_and_update_child(newPar, des, recipientTree);
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" des pre loop", des->get_data().desIds);
        dbWriteOttSet(" newPar pre loop", newPar->get_data().desIds);
    }
    if (oldPar != nullptr) {
        removeDesIdsToNdAndAnc(oldPar, des->get_data().desIds);
    }
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" oldPar post removeDesIdsToNdAndAnc", oldPar->get_data().desIds);
    }
    if (donorTree != &recipientTree) {
        recipientTree.steal_exclusion_statements(newPar, oldPar, *donorTree);
        const auto oids = recipientTree.steal_inclusion_statements(newPar, oldPar, *donorTree, bandBeingMerged);
        if (oldPar != nullptr) {
            removeDesIdsToNdAndAnc(oldPar, oids);
        }
        newPar->get_data().desIds.insert(begin(oids), end(oids));
        for (auto nd : iter_pre_n(des)) {
            register_tree_for_node(nd, &recipientTree);
            recipientTree.register_exclusion_statement_for_transferring_node(nd, *donorTree);
            recipientTree.register_inclusion_statement_for_transferring_node(nd, *donorTree);
        }
    }
    if (newPar->getParent()) {
        addDesIdsToNdAndAnc(newPar->getParent(), newPar->get_data().desIds);
    }
    if (bandBeingMerged == nullptr) {
        LOG(DEBUG) << "after actions of transfer_subtree_in_forest";
        dbWriteNewick(newPar);
        LOG(DEBUG) << "transfer_subtree_in_forest post";
        debug_invariants_check();
        LOG(DEBUG) << "transfer_subtree_in_forest exiting";
    }
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::zip_paths_from_barren_node(
            FTree<RTSplits, MappedWithSplitsData> &donorTree,
            NodeWithSplits * donorDes,
            NodeWithSplits * donorAnc,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree,
            NodeWithSplits * recipientDes,
            NodeWithSplits * recipientAnc,
            SupertreeContextWithSplits *sc) {
    auto currDoomedChild = donorDes;
    auto currAttachPoint = recipientDes;
    bool hitDeepest = false;
    while (currDoomedChild != donorAnc) {
        auto nextDoomedChild = currDoomedChild->getParent();
        if (hitDeepest || currAttachPoint == recipientAnc) {
            hitDeepest = true;
            currAttachPoint = recipientTree.create_deeper_node(currAttachPoint);
        } else {
            currAttachPoint = currAttachPoint->getParent();
            assert(currAttachPoint != nullptr);
        }
        currAttachPoint = move_all_sibs(currDoomedChild, donorTree, currAttachPoint, recipientTree, sc);
        if (nextDoomedChild == nullptr) {
            assert(false);
        }
        currDoomedChild = nextDoomedChild;
    }
    if (hitDeepest || currAttachPoint == recipientAnc) {
        hitDeepest = true;
        currAttachPoint = recipientTree.create_deeper_node(currAttachPoint);
    } else {
        currAttachPoint = currAttachPoint->getParent();
        assert(currAttachPoint != nullptr);
    }
    auto dp = donorAnc->getParent();
    OttIdSet dpoids;
    for (auto dac : iter_anc(*donorAnc)) {
        auto & dacdi = dac->get_data().desIds;
        dpoids.insert(begin(dacdi), end(dacdi));
    }
    move_all_children(donorAnc, donorTree, currAttachPoint, recipientTree, nullptr);
    dbWriteOttSet("   donorAnc = ", dpoids);
    removeDesIdsToNdAndAnc(donorAnc, dpoids);
    if (dp == nullptr) {
        donorTree._set_root(nullptr);
        register_tree_for_node(dp, nullptr);
        return false;
    }
    return true;
}

template<typename T, typename U>
NodeWithSplits * GreedyBandedForest<T, U>::move_all_sibs(
            NodeWithSplits * donorC,
            FTree<RTSplits, MappedWithSplitsData> &donorTree,
            NodeWithSplits * attachPoint,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree,
            SupertreeContextWithSplits *) {
    auto dp = donorC->getParent();
    assert(dp != nullptr);
    OttIdSet dpoids;
    for (auto c :iter_child(*dp)) {
        dpoids.insert(begin(c->get_data().desIds), end(c->get_data().desIds));
    }
    donorC->detachThisNode();
    auto p = move_all_children(dp, donorTree, attachPoint, recipientTree, nullptr);
    assert(not dp->getFirstChild());
    dbWriteOttSet("   dpoids =", dpoids);
    removeDesIdsToNdAndAnc(dp, dpoids);
    register_tree_for_node(donorC, nullptr);
    if (dp == nullptr) {
        donorTree._set_root(nullptr);
        register_tree_for_node(dp, nullptr);
    }
    return p.first;
}



template class GreedyBandedForest<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

} // namespace otc

