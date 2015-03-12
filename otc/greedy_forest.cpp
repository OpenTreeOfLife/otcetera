#include "otc/greedy_forest.h"
#include "otc/embedding.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
namespace otc {

template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly) {
    std::map<const T *, typename U::node_type *> gpf2scaff;
    std::map<long, typename U::node_type *> & dOttIdToNode = destTree.getData().ottIdToNode;
    //LOG(DEBUG) << " adding " << srcPoly;
    gpf2scaff[srcPoly] = destPoly;
    for (auto sn : iter_pre_n_const(srcPoly)) {
        if (sn == srcPoly) {
            continue;
        }
        auto sp = sn->getParent();
        //LOG(DEBUG) << " looking for " << sp << " the parent of " << sn;
        auto dp = gpf2scaff.at(sp);
        typename U::node_type * dn;
        auto nid = sn->getOttId();
        if (sn->hasOttId() && nid > 0) { // might assign negative number to nodes created in synth...
            auto oid = sn->getOttId();
            dn = dOttIdToNode.at(oid);
            if (dn->getParent() != dp) {
                dn->_detachThisNode();
                dn->_setNextSib(nullptr);
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
    finalizeTree(sc);
    assert(trees.size() == 1);
    auto & resolvedTree = begin(trees)->second;
    copyStructureToResolvePolytomy(resolvedTree.getRoot(), sc.scaffoldTree, &scaffoldNode);
    // remap all path pairings out of this node...
    for (auto treeInd2eout : embedding->edgeBelowEmbeddings) {
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " setOttId(" << snoid<< ')';
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
void GreedyPhylogeneticForest<T,U>::finalizeTree(SupertreeContextWithSplits &) {
    LOG(DEBUG) << "finalizeTree for a forest with " << trees.size() << " roots:";
    auto roots = getRoots();
    if (roots.size() == 0) {
        attachAllKnownTipsAsNewTree();
        return;
    } else {
        for (auto r : roots) {
            std::cerr << " tree-in-forest = "; writeNewick(std::cerr, r); std::cerr << '\n';
        }
    }
    if (trees.size() < 2) {
        return;
    }
    LOG(WARNING) << "finalizeTree is not correctly merging";
    auto rit = begin(roots);
    auto firstRoot = *rit;
    for (++rit; rit != end(roots); rit = roots.erase(rit)) {
        if ((*rit)->isTip()) {
            firstRoot->addChild(*rit);
        } else {
            assert((*rit)->getParent() == nullptr);
            for (auto c : iter_child(**rit)) {
                assert(c->isTip());
                firstRoot->addChild(c);
            }
        }
    }
    firstRoot->getData().desIds = ottIdSet;
    std::cerr << " finalized-tree-from-forest = "; writeNewick(std::cerr, firstRoot); std::cerr << '\n';
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
    assert(false);
}

template class GreedyPhylogeneticForest<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

} // namespace otc

