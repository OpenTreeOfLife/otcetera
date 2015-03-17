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
                                    SupertreeContextWithSplits & sc) {
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
            if (contains(sc.detachedScaffoldNodes, dn)) {
                dn->_setLChild(nullptr);
                dn->_setNextSib(nullptr);
                sc.detachedScaffoldNodes.erase(dn);
                sc.scaffoldTree.markAsAttached(dn);
            }
            if (dn->getParent() != dp) {
                if (dn->getParent() != nullptr) {
                    dn->_detachThisNode();
                    sc.scaffoldTree.markAsDetached(dn);
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


// 1. finalizeTree
// 2. copy the structure into the scaffold tree
// 3. update all of the outgoing paths so that they map
//      to this taxon
template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::finishResolutionOfEmbeddedClade(U & scaffoldNode,
                                                                    NodeEmbedding<T, U> * embedding,
                                                                    SupertreeContextWithSplits & sc) {
    const auto snoid = scaffoldNode.getOttId();
    LOG(DEBUG) << "finishResolutionOfEmbeddedClade for " << snoid;
    debugInvariantsCheck();
    finalizeTree(sc);
    debugInvariantsCheck();
    assert(trees.size() == 1);
    auto & resolvedTree = begin(trees)->second;
    const auto beforePar = scaffoldNode.getParent();
    checkAllNodePointersIter(scaffoldNode);
    copyStructureToResolvePolytomy(resolvedTree.getRoot(), sc.scaffoldTree, &scaffoldNode, sc);
    checkAllNodePointersIter(scaffoldNode);
    assert(beforePar == scaffoldNode.getParent());
    // remap all path pairings out of this node...
    for (auto treeInd2eout : embedding->edgeBelowEmbeddings) {
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " setOttId(" << snoid<< ')';
            eout->scaffoldDes = &scaffoldNode;
            eout->setOttIdSet(snoid, sc.scaffold2NodeEmbedding);
        }
    }
}

template<typename T, typename U>
bool GreedyPhylogeneticForest<T,U>::addGroupToNewTree(const OttIdSet & incGroup,
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
bool GreedyPhylogeneticForest<T,U>::addLeaf(PathPairing<T, U> * ,
                      const OttIdSet & incGroup,
                      const OttIdSet & ,
                      int ,
                      long ,
                      SupertreeContextWithSplits &) {
    assert(incGroup.size() == 1);
    const auto ottId = *incGroup.begin();
    registerLeaf(ottId);
    return true;
}

template<typename T, typename U>
bool GreedyPhylogeneticForest<T,U>::attemptToAddGrouping(PathPairing<T, U> * ppptr,
                                                        const OttIdSet & incGroup,
                                                        const OttIdSet & leafSet,
                                                        int treeIndex,
                                                        long groupIndex,
                                                        SupertreeContextWithSplits &sc) {
    if (incGroup.size() == 1) {
        addLeaf(ppptr, incGroup, leafSet, treeIndex, groupIndex, sc);
        return true;
    }
    return addGroupToNewTree(incGroup, leafSet, treeIndex, groupIndex);
    /*
    if (this->empty() || areDisjoint(ottIdSet, leafSet)) { // first grouping, always add...
        sc.log(CLADE_CREATES_TREE, ppptr->phyloChild);
        addGroupToNewTree(incGroup, leafSet, treeIndex, groupIndex);
        return true;
    }
    std::list<std::tuple<NodeWithSplits *, NodeWithSplits *, NodeWithSplits *> > rootIngroupPairs;
    for (auto & itIt : trees) {
        auto & t = itIt.second;
        auto r = t.getRoot();
        auto srca = couldAddToTree(r, incGroup, leafSet);
        if (!std::get<0>(srca)) {
            sc.log(CLADE_REJECTED, ppptr->phyloChild);
            return false;
        }
        auto incGroupMRCA = std::get<1>(srca);
        if (incGroupMRCA != nullptr) {
            rootIngroupPairs.push_back(std::make_tuple(r, incGroupMRCA, std::get<2>(srca)));
        }
    }
    assert(!rootIngroupPairs.empty());
    auto rit = rootIngroupPairs.begin();
    const auto & rip = *rit;
    auto retainedRoot = std::get<0>(rip);
    auto ing = std::get<1>(rip);
    if (rootIngroupPairs.size() == 1) {
        sc.log(CLADE_ADDED_TO_TREE, ppptr->phyloChild);
        auto outg = std::get<1>(rip);
        addIngroupAtNode(retainedRoot, ing, outg, incGroup, leafSet);
    } else {
        for (++rit; rit != rootIngroupPairs.end(); ++rit) {
            auto dr = std::get<0>(*rit);
            auto di = std::get<1>(*rit);
            auto dout = std::get<2>(*rit);
            graftTreesTogether(retainedRoot, ing, dr, di, dout, incGroup, leafSet);
        }
    }
    return true; */
}


// returns:
//      false, nullptr, nullptr if incGroup/leafset can't be added. 
//      true, nullptr, nullptr if there is no intersection with the leafset
//      true, nullptr, excGroup-MRCA if there is no intersection with the leafset, but not the incGroup
//      true, incGroup-MRCA, excGroup-MRCA if there is an intersection the leafset and the incGroup and the incGroup could be added
//          to this tree of the forest
template<typename T, typename U>
CouldAddResult GreedyPhylogeneticForest<T,U>::couldAddToTree(NodeWithSplits *root, const OttIdSet & incGroup, const OttIdSet & leafSet) {
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
std::vector<T *> GreedyPhylogeneticForest<T,U>::getRoots(){
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
MergeStartInfo calcTreePairMergeStartInfo(T & tdt, U* toDie, const std::set<InterTreeBand<T> *> &, T & tt, U* target);
std::set<NodeWithSplits *> detectExclusionsToMoveDown(FTree<RTSplits, MappedWithSplitsData> & donor,
                                                    NodeWithSplits *dr,
                                                    FTree<RTSplits, MappedWithSplitsData> & recipient,
                                                    NodeWithSplits *rr);

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
MergeStartInfo calcTreePairMergeStartInfo(T & tdt, U* toDie, const std::set<InterTreeBand<T> *> & bands, T & tt, U* target) {
    const std::set<NodeWithSplits *> emptySet;
    std::list<MRCABandedPaths> emptyList;
    MergeStartInfo r{false, emptyList, emptySet};
    std::list<MRCABandedPaths> & desBands = std::get<1>(r);
    std::set<NodeWithSplits *> & unbandedChildren = std::get<2>(r);
    std::map<NodeWithSplits *, InterTreeBand<T> * > targetNd2ITB;
    for (auto b : bands) {
        auto n = b->nodeForTree(tt);
        if (n != nullptr) {
            targetNd2ITB[n] = b;
        }
    }
    std::map<NodeWithSplits *, std::set<InterTreeBand<T> *> > childNd2DesITB;
    for (auto nd : iter_child(target)) {
        childNd2DesITB = std::set<InterTreeBand<T> *>{};
    }
    for (auto tp : targetNd2ITB) {
        for (auto & cp : childNd2DesITB) {
            if (isAncestorDesNoIter(cp.first, tp.first)) {
                cp.second.insert(tp.second);
            }
        }
    }
    NOT_IMPLEMENTED
    /*for (auto & cp :childNd2DesITB) {
        if (cp.second.empty()) {
            unbandedChildren.insert(cp.first);
        } else {
            MRCABandedPaths mbp;
            std::set<NodeWithSplits *> tdns;
            std::set<NodeWithSplits *> ttns;
            for (auto bpp : cp.second) {
                auto tdn = bpp->nodeForTree(tdt);
                auto n = bpp->nodeForTree(tt);
                tdns.insert(tdn);
                ttns.insert(n);
                mbp.bandedPairs.insert(std::pair<NodeWithSplits *, NodeWithSplits *>(tdn, n));
            }
            mbp.mrcaF = getMRCA(tdns);
            mbp.mrcaS = getMRCA(ttns);
            desBands.push_back(mbp);
        }
    }*/
    return r;
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T, U>::mergePathToNextBand(
            FTree<RTSplits, MappedWithSplitsData> & donor,
            NodeWithSplits * spikeDes,
            FTree<RTSplits, MappedWithSplitsData> & recipient, 
            SupertreeContextWithSplits &sc) {
    NOT_IMPLEMENTED;
}

std::set<NodeWithSplits *> detectExclusionsToMoveDown(FTree<RTSplits, MappedWithSplitsData> & donor,
                                                    NodeWithSplits *dr,
                                                    FTree<RTSplits, MappedWithSplitsData> & recipient,
                                                    NodeWithSplits *rr) {
    assert(dr != nullptr);
    assert(rr != nullptr);
    std::set<NodeWithSplits *> r;
    for (auto c : iter_child(*rr)) {
        if (donor.anyExcludedAtNode(dr, c->getData().desIds)) {
            r.insert(c);
        }
    }for (auto c : iter_child(*dr)) {
        if (recipient.anyExcludedAtNode(rr, c->getData().desIds)) {
            r.insert(c);
        }
    }
    return r;
}


template<typename T, typename U>
void GreedyPhylogeneticForest<T, U>::transferSubtreeInForest(NodeWithSplits * des,
                                                    FTree<RTSplits, MappedWithSplitsData> & possDonor,
                                                    NodeWithSplits * newPar,
                                                 FTree<RTSplits, MappedWithSplitsData> & recipientTree, 
                                                 FTree<RTSplits, MappedWithSplitsData> *donorTree) {
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
    if (donorTree != &recipientTree) {
        auto p = des->getParent();
        recipientTree.stealExclusionStatements(newPar, p, *donorTree);
        for (auto nd : iter_pre_n(des)) {
            recipientTree.registerExclusionStatementForTransferringNode(nd, *donorTree);
        }
    }
    des->_detachThisNode();
    newPar->addChild(des);
    addDesIdsToNdAndAnc(newPar, des->getData().desIds);
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T, U>::mergeBandedTrees(FTree<RTSplits, MappedWithSplitsData> & donor,
                                                 FTree<RTSplits, MappedWithSplitsData> & recipient, 
                                                 const MergeStartInfo & mi,
                                                 SupertreeContextWithSplits &sc) {
    NOT_IMPLEMENTED;
    /*
    const std::set<NodeWithSplits *> & desBands = std::get<1>(mi);
    const std::set<NodeWithSplits *> & unbandedChildren = std::get<2>(mi);
    if (std::get<0>(mi)) {
        //merge root of donor recipient
        auto dr = donor.getRoot();
        auto rr = recipient.getRoot();
        auto & rdi = rr->getData().desIds;
        const auto toPushDown = detectExclusionsToMoveDown(donor, dr, recipient, rr);
        NodeWithSplits * rrr = nullptr;
        if (!toPushDown.empty()) {
            recipient.createDeeperRoot();
            rrr = recipient.getRoot();
            assert(rr->getParent() == rrr);
            rrr->getData().desIds = rdi;
            for (auto nd : toPushDown) {
                transferSubtreeInForest(nd, donor, rrr, recipient, nullptr);
            }
        }
        for (auto spikeDes : desBands) {
            mergePathToNextBand(donor, spikeDes, recipient, sc);
        }
        for (auto c : iter_child(*dr)) {
            if (contains(toPushDown, c)) {
                continue;
            }
            transferSubtreeInForest(c, donor, rrr, recipient, &donor);
        }
        if (rrr != nullptr) {
            rrr->getData().desIds.insert(begin(rdi), end(rdi));
        }
        return;
    }*/
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T, U>::mergeForest(SupertreeContextWithSplits &sc) {
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
    for (auto i = sortedTrees.size() -1 ; i > 0; --i) {
        FTreeType & toDie = *sortedTrees[i];
        std::set<InterTreeBand<typename T::data_type> *> itbSet = collectBandsForSubtree(toDie, toDie.getRoot());
        if (!itbSet.empty()) {
            dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
            NOT_IMPLEMENTED;
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
    /*
    const std::set<NodeWithSplits *> emptySet;
    const std::set<NodeWithSplits *> emptyList;
    for (auto i = sortedTrees.size() -1 ; i > 0; --i) {
        FTreeType & toDie = *sortedTrees[i];

        std::set<InterTreeBand<T> *> itbSet = collectBandsForSubtree(toDie, toDie.getRoot());
        MergeStartInfo bmi{false, emptySet, emptySet};
        FTreeType * bestTarget  = nullptr;
        for (auto j = 0U; j < i; ++j) {
            FTreeType & target = *sortedTrees[i];
            const auto mi = calcTreePairMergeStartInfo(toDie, toDie.getRoot(), itbSet, target, target.getRoot());
            if (j == 0 || std::get<0>(mi) || std::get<1>(mi).size() >  std::get<1>(bmi).size()) {
                bestTarget = &target;
                bmi = mi;
                if (std::get<0>(mi)) {
                    break;
                }
            }
        }
        assert(bestTarget != nullptr);
        mergeBandedTrees(toDie, *bestTarget, bmi, sc);
    }
    auto tmIt = trees.begin();
    ++tmIt; // keep the first tree this is the final target
    while (tmIt != trees.end()) {
        tmIt = trees.erase(tmIt);
    }
    */
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T, U>::mergeTreesToFirstPostBandHandling(SupertreeContextWithSplits &sc) {
    if (trees.size() == 1) {
        return;
    }
    assert(trees.size() > 0);
    auto trIt = begin(trees);
    auto & firstTree = trIt->second;
    auto firstRoot = firstTree.getRoot();
    OttIdSet idsIncluded = firstRoot->getData().desIds;
    for (++trIt; trIt != end(trees); trIt = trees.erase(trIt)) {
        auto & currTree = trIt->second;
        auto currRoot = currTree.getRoot();
        assert(currRoot->getParent() == nullptr);
        assert(!currRoot->isTip());
        std::list<node_type *> rc;
        for (auto currChild : iter_child(*currRoot)) {
            rc.push_back(currChild);
        }
        for (auto currChild : rc) {
            transferSubtreeInForest(currChild, currTree, firstRoot, firstTree, &currTree);
        }
        currRoot->_setLChild(nullptr);
    }
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::finalizeTree(SupertreeContextWithSplits &sc) {
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

// addIngroupAtNode adds leaves for all "new" ottIds to:
//      the root if they are in the excGroup only, OR
//      ing if they are in the incGroup.
// The caller has guaranteed that:
//   outg is a part of the only tree in the forest that intersects with leafSet, AND
//   outg is the MRCA of all of the taxa in leafset in this tree.
//   ing is the MRCA of the incGroup in this tree
//
template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::addIngroupAtNode(NodeWithSplits *, //delete root param?
                                                     NodeWithSplits *ing,
                                                     NodeWithSplits *outg,
                                                     const OttIdSet & incGroup,
                                                     const OttIdSet & leafSet) {
    assert(outg != nullptr);
    if (ing == nullptr) {
        // there was no intersection with the incGroup...
        // So: create a new node at the MRCA of the excGroup, and add the incGroup to that node.
        assert(outg != nullptr);
        ing = nodeSrc.createChild(outg);
        ing->getData().desIds = incGroup;
        for (auto io : incGroup) {
            addChildForOttId(*ing, io, nodeSrc);
        }
    } else {
        const auto newIngLeaves = set_difference_as_set(incGroup, ottIdSet);
        if (!newIngLeaves.empty()) {
            ing->getData().desIds.insert(begin(newIngLeaves), end(newIngLeaves));
            for (auto io : newIngLeaves) {
                addChildForOttId(*ing, io, nodeSrc);
            }
            for (auto ianc : iter_anc(*ing)) {
                if (ianc == outg) {
                    break; // we'll update outg and below in the code below
                }
                ianc->getData().desIds.insert(begin(newIngLeaves), end(newIngLeaves));
            }
        }
    }
    // attach any new excGroup leaves as children of outg, creating/expanding a polytomy.
    const auto newLeaves = set_difference_as_set(leafSet, ottIdSet);
    ottIdSet.insert(begin(newLeaves), end(newLeaves));
    const auto newOutLeaves = set_difference_as_set(newLeaves, incGroup);
    for (auto oo : newOutLeaves) {
        addChildForOttId(*outg, oo, nodeSrc);
    }
    // add new ottids to desIDs in outg and its ancestors
    auto f = outg->getFirstChild();
    assert(f != nullptr);
    for (auto oanc : iter_anc(*f)) {
        oanc->getData().desIds.insert(begin(newLeaves), end(newLeaves));
    }
}

template<typename T, typename U>
void GreedyPhylogeneticForest<T,U>::graftTreesTogether(NodeWithSplits *, //rr,
                                                       NodeWithSplits *, //ri,
                                                       NodeWithSplits *, //delr,
                                                       NodeWithSplits *, //deli,
                                                       NodeWithSplits *, //delo,
                                                       const OttIdSet & , //incGroup,
                                                       const OttIdSet & ){  //leafSet
    NOT_IMPLEMENTED; // refactoring to banded
}

template class GreedyPhylogeneticForest<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

} // namespace otc

