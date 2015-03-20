#include "otc/greedy_forest.h"
#include "otc/embedding.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/write_dot.h"
#include "debug.h"
namespace otc {

template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly,
                                    SupertreeContextWithSplits * sc) {
    assert(sc != nullptr);
    std::map<const T *, typename U::node_type *> gpf2scaff;
    std::map<long, typename U::node_type *> & dOttIdToNode = destTree.getData().ottIdToNode;
    LOG(DEBUG) << " adding " << srcPoly;
    LOG(DEBUG) << " copying structure to resolve " << destPoly->getOttId();
    gpf2scaff[srcPoly] = destPoly;
    for (auto sn : iter_pre_n_const(srcPoly)) {
        if (sn == srcPoly) {
            continue;
        }
        auto sp = sn->getParent();
        auto dp = gpf2scaff.at(sp);
        typename U::node_type * dn;
        auto nid = sn->getOttId();
        if (sn->hasOttId() && nid > 0) { // might assign negative number to nodes created in synth...
            //LOG(DEBUG) << " node in src has ID " << nid;
            dn = dOttIdToNode.at(nid);
            //LOG(DEBUG) << " in dest node, that ID maps to a node with id:  " << dn->getOttId();
            assert(dn != destPoly);
            if (contains(sc->detachedScaffoldNodes, dn)) {
                dn->_setFirstChild(nullptr);
                dn->_setNextSib(nullptr);
                sc->detachedScaffoldNodes.erase(dn);
                sc->scaffoldTree.markAsAttached(dn);
            }
            if (dn->getParent() != dp) {
                if (dn->getParent() != nullptr) {
                    dn->_detachThisNode();
                    sc->scaffoldTree.markAsDetached(dn);
                }
                assert(dn->getNextSib() == nullptr);
                dp->addChild(dn);
            }
        } else {
            dn = destTree.createChild(dp);
            dn->setOttId(sn->getOttId());
        }
        //LOG(DEBUG) << " adding " << sn;
        gpf2scaff[sn] = dn;
    }
}

template<typename T>
bool canBeResolvedToDisplayExcGroup(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup) {
    for (auto c : iter_child_const(*nd)) {
        if (haveIntersection(incGroup, c->getData().desIds) && haveIntersection(excGroup, c->getData().desIds)) {
            return false;
        }
    }
    return true;
}

template<typename T, typename U>
void GreedyBandedForest<T,U>::writeFirstTree(std::ostream & treeFileStream) {
    writeNewick(treeFileStream, getRoots().at(0));
}

// 1. finalizeTree
// 2. copy the structure into the scaffold tree
// 3. update all of the outgoing paths so that they map
//      to this taxon
template<typename T, typename U>
void GreedyBandedForest<T,U>::finishResolutionOfEmbeddedClade(U & scaffoldNode,
                                                                    NodeEmbedding<T, U> * embedding,
                                                                    SupertreeContextWithSplits * sc) {
    assert(sc != nullptr);
    const auto snoid = scaffoldNode.getOttId();
    LOG(DEBUG) << "finishResolutionOfEmbeddedClade for " << snoid;
    debugInvariantsCheck();
    finalizeTree(sc);
    debugInvariantsCheck();
    assert(trees.size() == 1);
    auto & resolvedTree = begin(trees)->second;
    const auto beforePar = scaffoldNode.getParent();
    checkAllNodePointersIter(scaffoldNode);
    copyStructureToResolvePolytomy(resolvedTree.getRoot(), sc->scaffoldTree, &scaffoldNode, sc);
    checkAllNodePointersIter(scaffoldNode);
    assert(beforePar == scaffoldNode.getParent());
    // remap all path pairings out of this node...
    for (auto treeInd2eout : embedding->edgeBelowEmbeddings) {
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " setOttId(" << snoid<< ')';
            eout->scaffoldDes = &scaffoldNode;
            eout->setOttIdSet(snoid, sc->scaffold2NodeEmbedding);
        }
    }
}

template<typename T, typename U>
bool GreedyBandedForest<T,U>::createAndAddPhyloStatement(
                const OttIdSet & incGroup,
                const OttIdSet & leafSet,
                int treeIndex,
                long groupIndex) {
    const auto & i = *(encountered.insert(incGroup).first);
    const auto & ls = (leafSet.empty() ? i : *(encountered.insert(leafSet).first));
    const OttIdSet exc = set_difference_as_set(ls, i);
    const auto & e = *(encountered.insert(exc).first);
    PhyloStatement ps(i, e, ls, PhyloStatementSource(treeIndex, groupIndex));
    LOG(DEBUG) << " GPF calling addPhyloStatement for tree=" << treeIndex << " group=" << groupIndex;
    return addPhyloStatement(ps);
}

template<typename T, typename U>
bool GreedyBandedForest<T,U>::addLeaf(
                      const OttIdSet & incGroup,
                      const OttIdSet & ,
                      int ,
                      long ,
                      SupertreeContextWithSplits *) {
    assert(incGroup.size() == 1);
    const auto ottId = *incGroup.begin();
    registerLeaf(ottId);
    return true;
}

template<typename T, typename U>
bool GreedyBandedForest<T,U>::attemptToAddGrouping(const OttIdSet & incGroup,
                                                        const OttIdSet & leafSet,
                                                        int treeIndex,
                                                        long groupIndex,
                                                        SupertreeContextWithSplits *sc) {
    if (incGroup.size() == 1) {
        addLeaf(incGroup, leafSet, treeIndex, groupIndex, sc);
        return true;
    }
    return createAndAddPhyloStatement(incGroup, leafSet, treeIndex, groupIndex);
}


// returns:
//      false, nullptr, nullptr if incGroup/leafset can't be added. 
//      true, nullptr, nullptr if there is no intersection with the leafset
//      true, nullptr, excGroup-MRCA if there is no intersection with the leafset, but not the incGroup
//      true, incGroup-MRCA, excGroup-MRCA if there is an intersection the leafset and the incGroup and the incGroup could be added
//          to this tree of the forest
template<typename T, typename U>
CouldAddResult GreedyBandedForest<T,U>::couldAddToTree(NodeWithSplits *root, const OttIdSet & incGroup, const OttIdSet & leafSet) {
    if (areDisjoint(root->getData().desIds, leafSet)) {
        return std::make_tuple(true, nullptr, nullptr);
    }
    const OttIdSet ointers = set_intersection_as_set(root->getData().desIds, leafSet);
    if (ointers.empty()) {
        return std::make_tuple(true, nullptr, nullptr);
    }
    const OttIdSet inters = set_intersection_as_set(root->getData().desIds, incGroup);
    if (inters.empty()) {
        auto aoLOttId = *(ointers.begin());
        auto aoLeaf = nodeSrc.getData().ottIdToNode[aoLOttId];
        assert(aoLeaf != nullptr);
        auto oNd = searchAncForMRCAOfDesIds(aoLeaf, ointers);
        return std::make_tuple(true, nullptr, oNd);
    }
    auto aLOttId = *(inters.begin());
    auto aLeaf = nodeSrc.getData().ottIdToNode[aLOttId];
    assert(aLeaf != nullptr);
    auto iNd = searchAncForMRCAOfDesIds(aLeaf, inters);
    if (ointers.size() == inters.size()) {
        return std::make_tuple(true, iNd, root);
    }
    auto oNd = searchAncForMRCAOfDesIds(iNd, ointers);
    if (iNd == oNd && !canBeResolvedToDisplay(iNd, inters, ointers)) {
        return std::make_tuple(false, iNd, iNd);
    }
    return std::make_tuple(true, iNd, oNd);
}

template<typename T, typename U>
std::vector<T *> GreedyBandedForest<T,U>::getRoots(){
    std::vector<T *> r;
    r.reserve(trees.size());
    for (auto & t : trees) {
        r.push_back(t.second.getRoot());
    }
    return r;
}

template<typename T, typename U>
std::map<T, std::set<T> > invertGCMap(const std::map<T, std::pair<T, U> > & k2pl) {
    std::map<T, std::set<T> > r;
    for (const auto & kIt : k2pl) {
        const T & k{kIt.first};
        r[kIt.second.first].insert(k);
    }
    return r;
}

template<typename T, typename U>
std::map<T, std::set<T> > invertGCListMap(const std::map<T, std::list<std::pair<T, U> > > & k2pl) {
    std::map<T, std::set<T> > r;
    for (const auto & kIt : k2pl) {
        const T & k{kIt.first};
        for (const auto & v : kIt.second) {
            r[v.first].insert(k);
        }
    }
    return r;
}

template<typename T, typename U>
std::set<InterTreeBand<typename T::data_type> * > collectBandsForSubtree(U & tree, T * node) {
    std::set<InterTreeBand<typename T::data_type> *> r;
    for (auto nd : iter_child(*node)) {
        const auto & bs = tree.getBandsForNode(&(*nd));
        for (auto & b : bs) {
            if (!b->isSingleTreeBand()) {
                r.insert(b);
            }
        }
    }
    return r;
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::transferSubtreeInForest(
                NodeWithSplits * des,
                FTree<RTSplits, MappedWithSplitsData> & possDonor,
                NodeWithSplits * newPar,
                FTree<RTSplits, MappedWithSplitsData> & recipientTree, 
                FTree<RTSplits, MappedWithSplitsData> *donorTree,
                InterTreeBand<RTSplits> * bandBeingMerged) {
    assert(des != nullptr);
    auto oldPar = des->getParent();
    LOG(DEBUG) << "top transferSubtreeInForest pre des == " << getDesignator(*des) << " " << (long) des << " " << std::hex << (long) des << std::dec;
    LOG(DEBUG) << "                            newPar " << (long) newPar << " " << std::hex << (long) newPar << std::dec;
    LOG(DEBUG) << "                            oldPar " << (long) oldPar << " " << std::hex << (long) oldPar << std::dec;
    LOG(DEBUG) << "    newick of newPar before actions of transferSubtreeInForest";
    dbWriteNewick(newPar);
    LOG(DEBUG) << "    newick of oldPar before actions of transferSubtreeInForest";
    dbWriteNewick(oldPar);
    dbWriteOttSet("    on entry des->desIds", des->getData().desIds);
    dbWriteOttSet("    on entry newPar->desIds", newPar->getData().desIds);
    dbWriteOttSet("    on entry oldPar->desIds", oldPar->getData().desIds);
    if (bandBeingMerged == nullptr) {
        debugInvariantsCheck();
    }
    assert(des != nullptr);
    assert(newPar != nullptr);
    if (donorTree == nullptr) {
        if (isAncestorDesNoIter(possDonor.getRoot(), des)) {
            donorTree = &possDonor;
        } else {
            assert(isAncestorDesNoIter(recipientTree.getRoot(), des));
            donorTree = &recipientTree;
        }
    }
    assert(!recipientTree.isExcludedFrom(des, newPar));
    assert(getTreeForNode(des) == donorTree);
    des->_detachThisNode();
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" des pre addAndUpdateChild", des->getData().desIds);
        dbWriteOttSet(" newPar pre addAndUpdateChild", newPar->getData().desIds);
    }
    addAndUpdateChild(newPar, des, recipientTree);
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" des pre loop", des->getData().desIds);
        dbWriteOttSet(" newPar pre loop", newPar->getData().desIds);
    }
    if (oldPar != nullptr) {
        removeDesIdsToNdAndAnc(oldPar, des->getData().desIds);
    }
    if (bandBeingMerged == nullptr) {
        dbWriteOttSet(" oldPar post removeDesIdsToNdAndAnc", oldPar->getData().desIds);
    }
    if (donorTree != &recipientTree) {
        recipientTree.stealExclusionStatements(newPar, oldPar, *donorTree);
        const auto oids = recipientTree.stealInclusionStatements(newPar, oldPar, *donorTree, bandBeingMerged);
        if (oldPar != nullptr) {
            removeDesIdsToNdAndAnc(oldPar, oids);
        }
        newPar->getData().desIds.insert(begin(oids), end(oids));
        for (auto nd : iter_pre_n(des)) {
            registerTreeForNode(nd, &recipientTree);
            recipientTree.registerExclusionStatementForTransferringNode(nd, *donorTree);
            recipientTree.registerInclusionStatementForTransferringNode(nd, *donorTree);
        }
    }
    if (newPar->getParent()) {
        addDesIdsToNdAndAnc(newPar->getParent(), newPar->getData().desIds);
    }
    if (bandBeingMerged == nullptr) {
        LOG(DEBUG) << "after actions of transferSubtreeInForest";
        dbWriteNewick(newPar);
        LOG(DEBUG) << "transferSubtreeInForest post";
        debugInvariantsCheck();
        LOG(DEBUG) << "transferSubtreeInForest exiting";
    }
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::performSingleBandMerge(
            std::size_t treeInd,
            InterTreeBand<RTSplits> * itb,
            const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
            SupertreeContextWithSplits *sc) {
    assert(itb != nullptr);
    auto btm = getTreeToNodeMapForBand(*itb);
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
    return mergeSingleBandedTree(toDie, itb, *toMergeTo, sc);
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::performSetOfSingleBandMerges(
            std::size_t treeInd,
            std::set<InterTreeBand<typename T::data_type> *> & itbSet,
            const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
            SupertreeContextWithSplits *sc) {
    LOG(WARNING) << itbSet.size() <<  " bands for a tree. It is not a great idea to merge these one at at time...";
    LOG(DEBUG) << itbSet.size() <<  " bands for a tree. It is not a great idea to merge these one at at time...";
    //            writeForestDOTToFN("writingForestMerge.dot");
    //            NOT_IMPLEMENTED;
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
            if (ip->isABandedNodeInThis(nd)) {
                postOrd.push_back(ip);
                itbIt = toOrganize.erase(itbIt);
            } else {
                ++itbIt;
            }
        }
    }
    assert(toOrganize.empty());
    for (auto sit : postOrd) {
        if (!performSingleBandMerge(treeInd, sit, sortedTrees, sc)) {
            return false;
        }
    }
    return true;
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::mergeForest(SupertreeContextWithSplits *sc) {
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
        std::set<InterTreeBand<typename T::data_type> *> itbSet = collectBandsForSubtree(toDie, toDie.getRoot());
        if (!itbSet.empty()) {
            if (itbSet.size() == 1) {
                stillHasLeaves[treeInd] = performSingleBandMerge(treeInd, *itbSet.begin(), sortedTrees, sc);
            } else {
                stillHasLeaves[treeInd] = performSetOfSingleBandMerges(treeInd, itbSet, sortedTrees, sc);
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
    mergeTreesToFirstPostBandHandling(sc);
    LOG(DEBUG) << "exiting mergeForest";
}

template<typename T, typename U>
NodeWithSplits * GreedyBandedForest<T, U>::moveAllSibs(
            NodeWithSplits * donorC,
            FTree<RTSplits, MappedWithSplitsData> &donorTree,
            NodeWithSplits * attachPoint,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree,
            SupertreeContextWithSplits *) {
    auto dp = donorC->getParent();
    assert(dp != nullptr);
    OttIdSet dpoids;
    for (auto c :iter_child(*dp)) {
        dpoids.insert(begin(c->getData().desIds), end(c->getData().desIds));
    }
    donorC->_detachThisNode();
    auto p = moveAllChildren(dp, donorTree, attachPoint, recipientTree, nullptr);
    dp->_setFirstChild(nullptr);
    dbWriteOttSet("   dpoids =", dpoids);
    removeDesIdsToNdAndAnc(dp, dpoids);
    registerTreeForNode(donorC, nullptr);
    if (dp == nullptr) {
        donorTree._setRoot(nullptr);
        registerTreeForNode(dp, nullptr);
    }
    return p.first;
}


template<typename T, typename U>
bool GreedyBandedForest<T, U>::zipPathsFromBarrenNode(
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
            currAttachPoint = recipientTree.createDeeperNode(currAttachPoint);
        } else {
            currAttachPoint = currAttachPoint->getParent();
            assert(currAttachPoint != nullptr);
        }
        currAttachPoint = moveAllSibs(currDoomedChild, donorTree, currAttachPoint, recipientTree, sc);
        if (nextDoomedChild == nullptr) {
            assert(false);
        }
        currDoomedChild = nextDoomedChild;
    }
    if (hitDeepest || currAttachPoint == recipientAnc) {
        hitDeepest = true;
        currAttachPoint = recipientTree.createDeeperNode(currAttachPoint);
    } else {
        currAttachPoint = currAttachPoint->getParent();
        assert(currAttachPoint != nullptr);
    }
    auto dp = donorAnc->getParent();
    OttIdSet dpoids;
    for (auto dac : iter_anc(*donorAnc)) {
        auto & dacdi = dac->getData().desIds;
        dpoids.insert(begin(dacdi), end(dacdi));
    }
    moveAllChildren(donorAnc, donorTree, currAttachPoint, recipientTree, nullptr);
    dbWriteOttSet("   donorAnc = ", dpoids);
    removeDesIdsToNdAndAnc(donorAnc, dpoids);
    if (dp == nullptr) {
        donorTree._setRoot(nullptr);
        registerTreeForNode(dp, nullptr);
        return false;
    }
    return true;
}


// If the parent of a banded node that was just removed is also banded to
//  the same tree, or is the root, then the subsequent handling of nodes attached
//  to a band or the root should work.
// But if there are intervening unbannded nodes, we want to zip the path from
//  nd to the next relevant ancestor onto the existing path...
template<typename T, typename U>
std::pair<NodeWithSplits *, NodeWithSplits *>
GreedyBandedForest<T, U>::findGrandparentThatIsRootOrBandSharing(
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
            r.second = recipientTree.getRoot();
            return r;
        }
        const auto & bfn = donorTree.getBandsForNode(a);
        for (auto & b : bfn) {
            const auto m = getTreeToNodeMapForBand(*b);
            if (contains(m, &recipientTree)) {
                r.first = a;
                r.second = const_cast<NodeWithSplits *>(m.at(&recipientTree));
                return r;
            }
        }
    }
    UNREACHABLE;
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::mergeSingleBandedTree(
            FTree<RTSplits, MappedWithSplitsData> &donorTree,
            InterTreeBand<RTSplits> * band,
            FTree<RTSplits, MappedWithSplitsData> &recipientTree,
            SupertreeContextWithSplits *sc) {
    LOG(DEBUG) << "pre mergeSingleBandedTree check";
    debugInvariantsCheck();
    const auto t2n = getTreeToNodeMapForBand(*band);
    auto dn = const_cast<NodeWithSplits *>(t2n.at(&donorTree));
    auto dp = dn->getParent();
    auto rn = const_cast<NodeWithSplits *>(t2n.at(&recipientTree));
    const auto & beforePhantomsRN = band->getPhantomNodes(rn);
    const auto & beforePhantomsDN = band->getPhantomNodes(dn);
    // any phantom nodes in common between dn and rn must be attached
    //  to a different tree. So they will remain "phantom" wrt rn
    auto movedNodeSet = set_difference_as_set(beforePhantomsRN, beforePhantomsDN);
    auto p = moveAllChildren(dn, donorTree, rn, recipientTree, band);
    OttIdSet movedIDSet;
    for (auto mel : movedNodeSet) {
        movedIDSet.insert(mel->getOttId());
    }
    for (auto mel : beforePhantomsDN) {
        movedIDSet.insert(mel->getOttId());
    }
    removeDesIdsToNdAndAnc(dn, movedIDSet);
    dbWriteOttSet(" movedIDSet =", movedIDSet);
    rn = p.first;
    band->removeNode(dn);
    band->removeFromSet(rn, movedNodeSet);
    if (p.second != rn) {
        NOT_IMPLEMENTED;
    }
    auto anc = findGrandparentThatIsRootOrBandSharing(donorTree, dn, recipientTree);
    bool r = true;
    if (anc.first == nullptr) {
        // dn is a child of the root
        registerTreeForNode(dn, nullptr);
        if (dp == nullptr) {
            donorTree._setRoot(nullptr);
            r = false;
        } else {
            dn->_detachThisNode();
        }
    } else {
        r = zipPathsFromBarrenNode(donorTree,
                                   dn,
                                   anc.first,
                                   recipientTree,
                                   p.first,
                                   anc.second,
                                   sc);
    }
    LOG(DEBUG) << "post mergeSingleBandedTree check";
    debugInvariantsCheck();
    LOG(DEBUG) << "  exiting mergeSingleBandedTree";
    return r;
}

template<typename T, typename U>
std::pair<NodeWithSplits *, NodeWithSplits*> GreedyBandedForest<T, U>::moveAllChildren(NodeWithSplits * donorParent,
                                                     FTree<RTSplits, MappedWithSplitsData> &donorTree,
                                                     NodeWithSplits * recipientNode,
                                                     FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                                                     InterTreeBand<RTSplits> * bandBeingMerged) {
    if (recipientTree.isExcludedFrom(donorParent, recipientNode)) {
        if (recipientNode == recipientTree.getRoot()) {
            recipientTree.createDeeperRoot();
            recipientNode = recipientTree.getRoot();
        } else {
            recipientNode = recipientTree.createDeeperNode(recipientNode);
        }
        debugInvariantsCheck();
    }
    auto attachmentPoint = recipientNode;
    if (donorTree.isExcludedFrom(recipientNode, donorParent)) {
        attachmentPoint = createNode(recipientNode, &recipientTree);
        debugInvariantsCheck();
    }
    std::list<node_type *> rc;
    for (auto currChild : iter_child(*donorParent)) {
        rc.push_back(currChild);
    }
    for (auto currChild : rc) {
        if (bandBeingMerged == nullptr) {
            debugInvariantsCheck();
        }
        transferSubtreeInForest(currChild, donorTree, attachmentPoint, recipientTree, &donorTree, bandBeingMerged);
        LOG(DEBUG) << "Back from transferSubtreeInForest";
    }
    donorParent->_setFirstChild(nullptr);
    return std::pair<NodeWithSplits *, NodeWithSplits*>{recipientNode, attachmentPoint};
}

template<typename T, typename U>
void GreedyBandedForest<T, U>::mergeTreesToFirstPostBandHandling(SupertreeContextWithSplits *) {
    if (trees.size() == 1) {
        return;
    }
    assert(trees.size() > 0);
    auto trIt = begin(trees);
    auto & firstTree = trIt->second;
    auto firstTreeRoot = firstTree.getRoot();
    OttIdSet idsIncluded = firstTreeRoot->getData().desIds;
    for (++trIt; trIt != end(trees);) {
        debugInvariantsCheck();
        auto & currTree = trIt->second;
        auto currTreeRoot = currTree.getRoot();
        assert(currTreeRoot->getParent() == nullptr);
        assert(!currTreeRoot->isTip());
        auto p = moveAllChildren(currTreeRoot, currTree, firstTreeRoot, firstTree, nullptr);
        firstTreeRoot = p.first;
        currTreeRoot = currTree.getRoot();
        for (auto j : nd2Tree) {
            assert((j.first == currTreeRoot || j.second != &currTree));
        }
        LOG(DEBUG) << "Before erase";
        debugInvariantsCheck();
        trIt = trees.erase(trIt);
        registerTreeForNode(currTreeRoot, nullptr);
        LOG(DEBUG) << "After erase";
        debugInvariantsCheck();
        LOG(DEBUG) << "After erase check";
    }
    LOG(DEBUG) << "check before exit of mergeTreesToFirstPostBandHandling";
    debugInvariantsCheck();
    LOG(DEBUG) << "exiting mergeTreesToFirstPostBandHandling";
}

template<typename T, typename U>
void GreedyBandedForest<T,U>::finalizeTree(SupertreeContextWithSplits *sc) {
    LOG(DEBUG) << "finalizeTree for a forest with " << trees.size() << " roots:";
    debugInvariantsCheck();
    auto roots = getRoots();
    for (auto r : roots) {
        LOG(DEBUG) << " tree-in-forest = "; dbWriteNewick(r);
    }
    if (trees.size() > 1) {
        LOG(WARNING) << "finalizeTree is not well thought out. merging of multiple trees is questionable.";
        const char * dbfn = "real-forest-in-finalizeTree.dot";
        LOG(WARNING) << "  should be writing DOT to " << dbfn;
        //std::ofstream outf(dbfn);
        //writeDOTForest(outf, *this);
        //outf.close();
        mergeForest(sc);
        LOG(WARNING) << "finished questionable mergeForest.";
    }
    debugInvariantsCheck();
    roots = getRoots();
    assert(roots.size() < 2);
    attachAllDetachedTips();
    debugInvariantsCheck();
    roots = getRoots();
    assert(roots.size() == 1);
    auto onlyRoot = *roots.begin();
    onlyRoot->getData().desIds = ottIdSet;
    LOG(DEBUG)<< " finalized-tree-from-forest = "; dbWriteNewick(onlyRoot);
}

template class GreedyBandedForest<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

} // namespace otc

