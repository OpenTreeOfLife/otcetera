Code that can probably be thrown away
This file is NOT part of the build!
It was helpful to have it around for reference

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
void GreedyBandedForest<T, U>::mergePathToNextBand(
            FTree<RTSplits, MappedWithSplitsData> & donor,
            NodeWithSplits * spikeDes,
            FTree<RTSplits, MappedWithSplitsData> & recipient, 
            SupertreeContextWithSplits *sc) {
    NOT_IMPLEMENTED;
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
void GreedyBandedForest<T, U>::addIngroupAtNode(NodeWithSplits *, //delete root param?
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
            auto nn = addChildForOttId(*ing, io, nodeSrc);
            registerTreeForNode(nn, )

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
void GreedyBandedForest<T, U>::graftTreesTogether(NodeWithSplits *, //rr,
                                                       NodeWithSplits *, //ri,
                                                       NodeWithSplits *, //delr,
                                                       NodeWithSplits *, //deli,
                                                       NodeWithSplits *, //delo,
                                                       const OttIdSet & , //incGroup,
                                                       const OttIdSet & ){  //leafSet
    NOT_IMPLEMENTED; // refactoring to banded
}

template<typename T, typename U>
bool GreedyBandedForest<T, U>::attemptToAddGrouping(const OttIdSet & incGroup,
                                                        const OttIdSet & leafSet,
                                                        int treeIndex,
                                                        long groupIndex,
                                                        SupertreeContextWithSplits *sc) {
    if (incGroup.size() == 1) {
        addLeaf(incGroup, leafSet, treeIndex, groupIndex, sc);
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


template<typename T, typename U>
void FTree<T, U>::addSubtree(RootedTreeNode<T> * subtreeRoot,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedInc,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedExc,
                const std::map<node_type *, std::pair<node_type*, PhyloStatementSource> > & otherIC,
                const std::map<node_type *, std::list<std::pair<node_type*, PhyloStatementSource> > > & otherEC) {
    assert(subtreeRoot);
    auto par = subtreeRoot->getParent();
    assert(par);
    const OttIdSet & sroids = subtreeRoot->getData().desIds;
    OttIdSet overlapIds = set_intersection_as_set(root->getData().desIds, sroids);
    if (overlapIds.empty()) {
        //the new tree does not have any leaves with include statements in the current tree...
        bool needDeeperRoot = false;
        for (auto oid : sroids) {
            if (ottIdIsExcludedFromRoot(oid)) {
                needDeeperRoot = true;
                break;
            }
        }
        if (needDeeperRoot) {
            createDeeperRoot();
        }
        addAndUpdateChild(root, subtreeRoot, *this) UNTESTED
        for (auto c : iter_leaf_n_const(*subtreeRoot)) {
            //connectedIds.insert(c->getOttId());
        }
        root->getData().desIds.insert(begin(subtreeRoot->getData().desIds), end(subtreeRoot->getData().desIds));
        const PhyloStatementSource bogusPSS{-1, -1};
        for (auto sn : iter_pre_n(subtreeRoot)) { //TMP need iter_node_n
            const auto esn = otherInvertedExc.find(sn);
            if (esn != otherInvertedExc.end()) {
                for (auto eValNd : esn->second) {
                    excludesConstraints[eValNd].push_back(GroupingConstraint{sn, bogusPSS});
                }
            }
            const auto isn = otherInvertedInc.find(sn);
            if (isn != otherInvertedInc.end()) {
                for (auto eValNd : isn->second) {
                    auto pinP = includesConstraints.find(eValNd);
                    if (pinP == includesConstraints.end() || isAncestorDesNoIter(pinP->first, eValNd)) {
                        includesConstraints.emplace(eValNd, GroupingConstraint{sn, bogusPSS});
                    }
                }
            }
        }
        return;
    }
    assert(false);
}

template<typename T, typename U>
MergeStartInfo calcTreePairMergeStartInfo(T & tdt, U* toDie, const std::set<InterTreeBand<T> *> &, T & tt, U* target);
std::set<NodeWithSplits *> detectExclusionsToMoveDown(FTree<RTSplits, MappedWithSplitsData> & donor,
                                                    NodeWithSplits *dr,
                                                    FTree<RTSplits, MappedWithSplitsData> & recipient,
                                                    NodeWithSplits *rr);

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
    for (auto i = sortedTrees.size() - 1 ; i > 0; --i) {
        FTreeType & toDie = *sortedTrees.at(i);
        std::set<InterTreeBand<typename T::data_type> *> itbSet = collectBandsForSubtree(toDie, toDie.getRoot());
        if (!itbSet.empty()) {
            if (itbSet.size() == 1) {
                InterTreeBand<typename T::data_type> * itb = *itbSet.begin();
                FTreeType * toMergeTo = nullptr;
                for (auto j = 0U; j < i; ++j) {
                    auto btp = sortedTrees.at(j);
                    if (bts.find(ntp) != bts.end()) {
                        toMergeTo = btp;
                        break;
                    }
                }
                assert(toMergeTo != nullptr);
                stillHasLeaves[i] = mergeSingleBandedTree(toDie, itb, *toMergeTo, sc);
            } else {
                LOG(DEBUG) <<  << "bands for a tree";
                writeForestDOTToFN("writingForestMerge.dot");
                NOT_IMPLEMENTED;
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
void GreedyBandedForest<T, U>::mergeBandedTrees(FTree<RTSplits, MappedWithSplitsData> & donor,
                                                 FTree<RTSplits, MappedWithSplitsData> & recipient, 
                                                 const MergeStartInfo & mi,
                                                 SupertreeContextWithSplits *sc) {
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

