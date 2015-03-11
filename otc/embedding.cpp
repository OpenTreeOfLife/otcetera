#include <queue>
#include "otc/greedy_forest.h"
#include "otc/embedding.h"
#include "otc/util.h"
#include "otc/supertree_context.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/error.h"
namespace otc {
constexpr bool COLLAPSE_IF_CONFLICT = true;


bool culledAndCompleteIncompatWRTLeafSet(const OttIdSet & culled,
                                                const OttIdSet & complete,
                                                const OttIdSet & leafSet) {
    //TMP this could be more efficient. See areCompatibleDesIdSets
    const OttIdSet inter = set_intersection_as_set(culled, complete);
    if (inter.empty()) {
        return false;
    }
    if (inter == culled) {
        return false;
    }
    const OttIdSet compCulled = set_intersection_as_set(complete, leafSet);
    return (inter != compCulled);
}


// Returns all loop paths for nd and all edgeBelowEmbeddings of its children
template<typename T, typename U>
std::vector<const PathPairing<T, U> *>
NodeEmbedding<T, U>::getAllIncomingPathPairs(const T *nd,
                                                 const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                 std::size_t treeIndex) const {
    std::vector<const PathPairingWithSplits *> r;
    const auto lait = loopEmbeddings.find(treeIndex);
    if (lait != loopEmbeddings.end()) {
        for (const auto pp : lait->second) {
            r.push_back(pp);
        }
    }
    for (auto c : iter_child_const(*nd)) {
        //LOG(DEBUG) << "    getAllIncomingPathPairs c = " << getDesignator(*c);
        const auto cembed = eForNd.find(c);
        if (cembed == eForNd.end()) {
            //LOG(DEBUG) << "     No embedding found";
            continue;
        }
        const auto & emb = cembed->second;
        const auto ceait = emb.edgeBelowEmbeddings.find(treeIndex);
        if (ceait != emb.edgeBelowEmbeddings.end()) {
            for (const auto & e : ceait->second) {
                //std::cerr << "     Found some paths e.currChildOttIdSet ="; writeOttSet(std::cerr, "      ", e->currChildOttIdSet, " "); std::cerr << '\n';
    
                r.push_back(e);
            }
        } else {
                //LOG(DEBUG) << "     No embedding found for this tree";
        }
    }
    return r;
}



template<typename T, typename U>
bool PathPairing<T,U>::updateOttIdSetNoTraversal(const OttIdSet & oldEls, const OttIdSet & newEls) {
    if (debuggingOutputEnabled) {
        std::cerr << "  updateOttIdSetNoTraversal for " << (long)this << " in currChildOttIdSet"; writeOttSet(std::cerr, " ", currChildOttIdSet, " "); std::cerr << '\n';
        std::cerr << "  updateOttIdSetNoTraversal for " << (long)this << " in oldEls"; writeOttSet(std::cerr, " ", oldEls, " "); std::cerr << '\n';
        std::cerr << "  updateOttIdSetNoTraversal for " << (long)this << " in newEls"; writeOttSet(std::cerr, " ", newEls, " "); std::cerr << '\n';
    }
    auto i = set_intersection_as_set(oldEls, currChildOttIdSet);
    if (i.size() < oldEls.size()) {
        return false;
    }
    if (!i.empty()) {
        for (auto o : i) {
            currChildOttIdSet.erase(o);
        }
    }
    currChildOttIdSet.insert(begin(newEls), end(newEls));
    if (debuggingOutputEnabled) {
        std::cerr << "  updateOttIdSetNoTraversal for " << (long)this << " on updateOttIdSetNoTraversal exit "; writeOttSet(std::cerr, " ", currChildOttIdSet, " "); std::cerr << '\n';
    }
    return true;
}

template<typename T, typename U>
inline const OttIdSet & NodeEmbedding<T, U>::getRelevantDesIdsFromPath(const PathPairing<T, U> & path) {
    return path.getOttIdSet();
}
template<typename T, typename U>
OttIdSet NodeEmbedding<T, U>::getRelevantDesIdsFromPathPairSet(const PathPairSet & pps) {
    OttIdSet relevantIds;
    for (auto path : pps) {
        const auto & cdi = getRelevantDesIdsFromPath(*path);
        relevantIds.insert(begin(cdi), end(cdi));
    }
    //std::cerr << " getRelevantDesIdsFromPathPairSet returning"; writeOttSet(std::cerr, " ", relevantIds, " "); std::cerr << '\n';
    return relevantIds;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseSourceEdge(const T * , //phyloParent,
                                                 PathPairing<T, U> * ) { //path
    assert(false);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseSourceEdgesToForceOneEntry(U & ,
                                                             PathPairSet & pps,
                                                             std::size_t treeIndex,
                                                             SupertreeContextWithSplits & sc) {
    if (pps.size() < 2) {
        return;
    }
    auto relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeIndex);
    PathPairing<T, U> * firstPairing = *pps.begin();
    const T * onePhyloPar = firstPairing->phyloParent;
    const T * phyloMrca = searchAncForMRCAOfDesIds(onePhyloPar, relevantIds);
    std::set<const T *> prevCollapsed; 
    prevCollapsed.insert(phyloMrca); // we don't actually collapse this edge, we just add it to the set so we don't collapse it below....
    for (auto path : pps) {
        const auto pp = path->phyloParent;
        if (!contains(prevCollapsed, pp)) {
            collapseSourceEdge(pp, path);
            prevCollapsed.insert(pp);
        }
    }
}
template<typename T, typename U>
void NodeEmbedding<T, U>::resolveGivenContestedMonophyly(U & scaffoldNode,
                                                         SupertreeContextWithSplits & sc) {
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto ebaIt = edgeBelowEmbeddings.find(treeInd);
        if (ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        PathPairSet & pps = ebaIt->second;
        collapseSourceEdgesToForceOneEntry(scaffoldNode, pps, treeInd, sc);
    }
    resolveGivenUncontestedMonophyly(scaffoldNode, sc);
}
template<typename T, typename U>
std::set<PathPairing<T, U> *> NodeEmbedding<T, U>::getAllChildExitPaths(U & scaffoldNode, SupertreeContextWithSplits & sc) {
    std::set<PathPairing<T, U> *> r;
    for (auto c : iter_child(scaffoldNode)) {
        const auto & thr = sc.scaffold2NodeEmbedding.at(c);
        for (auto te : thr.edgeBelowEmbeddings) {
            r.insert(begin(te.second), end(te.second));
        }
    }
    return r;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::resolveGivenUncontestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc) {
    const OttIdSet EMPTY_SET;
    LOG(DEBUG) << "resolveGivenUncontestedMonophyly for " << scaffoldNode.getOttId();
    GreedyPhylogeneticForest<T,U> gpf;
    std::set<PathPairing<T, U> *> considered;
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto laIt = loopEmbeddings.find(treeInd);
        if (laIt == loopEmbeddings.end()) {
            continue;
        }
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeInd);
        PathPairSet & pps = laIt->second;
        // leaf set of this tree for this subtree
        // for repeatability, we'll try to add groupings in reverse order of desIds sets (deeper first)
        std::map<OttIdSet, PathPairing<T,U> *> mapToProvideOrder;
        for (auto pp : pps) {
            mapToProvideOrder[pp->getOttIdSet()] = pp;
        }
        long bogusGroupIndex = 0; // should get this from the node!
        typedef std::pair<const OttIdSet *, PathPairing<T,U> *>  q_t;
        std::queue<q_t> trivialQ;
        for (auto mpoIt = mapToProvideOrder.rbegin(); mpoIt != mapToProvideOrder.rend(); ++mpoIt) {
            auto ppptr = mpoIt->second;
            if (ppptr->pathIsNowTrivial()) {
                LOG(DEBUG) << "pathIsNowTrivial" ; writeOttSet(std::cerr, " ", mpoIt->first, " "); std::cerr << std::endl;
                const q_t toQ{&(mpoIt->first), ppptr};
                trivialQ.push(toQ);
            } else {
                const auto & d = mpoIt->first;
                gpf.attemptToAddGrouping(ppptr, d, relevantIds, static_cast<int>(treeInd), bogusGroupIndex++, sc);
                considered.insert(ppptr);
            }
        }
        while (!trivialQ.empty()) {
            const q_t triv = trivialQ.front();
            const OttIdSet * inc = triv.first;
            const auto ppptr = triv.second;
            gpf.addLeaf(ppptr, *inc, relevantIds, static_cast<int>(treeInd), bogusGroupIndex++, sc);
            considered.insert(ppptr);
            trivialQ.pop();
        }
    }
    // we might have missed some descendants  - any child that is has
    //  "scaffoldNode" as its embedded parent, but which is not involved
    //  any loop or exiting edges.
    //  This means that we have no info on the placement of such nodes.
    //      so we'll just attach them here.
    //  First step: get the list of paths for the children.
    int bogusTreeIndex = 123456; // should get this from the node!
    long bogusGroupIndex = 100000; // should get this from the node!
    auto childExitPaths = getAllChildExitPaths(scaffoldNode, sc);
    for (auto pathPtr : childExitPaths) {
        if (!contains(considered, pathPtr)) {
            gpf.attemptToAddGrouping(pathPtr, pathPtr->getOttIdSet(), EMPTY_SET, bogusTreeIndex, bogusGroupIndex++, sc);
            considered.insert(pathPtr); // @TMP not needed
        }
    }
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        for (auto snc : iter_child(scaffoldNode)) {
            assert(snc != nullptr);
        }
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, sc);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseGroup(U & scaffoldNode, SupertreeContext<T,U> & sc) {
    sc.log(COLLAPSE_TAXON, scaffoldNode);
    U * p = scaffoldNode.getParent();
    assert(p != nullptr); // can't disagree with the root !
    // remap all nodes in NodePairing to parent
    for (auto nai : nodeEmbeddings) {
        for (auto np : nai.second) {
            np->scaffoldNode = p;
        }
    }
    NodeEmbedding<T, U>& parEmbedding = sc.scaffold2NodeEmbedding.at(p);
    // every loop for this node becomes a loop for its parent
    for (auto lai : loopEmbeddings) {
        for (auto lp : lai.second) {
            assert(lp->scaffoldDes == &scaffoldNode);
            assert(lp->scaffoldAnc == &scaffoldNode);
            lp->scaffoldDes = p;
            lp->scaffoldAnc = p;
            parEmbedding.loopEmbeddings[lai.first].insert(lp);
        }
    }
    // every exit edge for this node becomes a loop for its parent if it is not trivial
    for (auto ebai : edgeBelowEmbeddings) {
        for (auto lp : ebai.second) {
            if (lp->scaffoldAnc == p) {
                if (lp->scaffoldDes == &scaffoldNode) {
                    // this only happens if a terminal was mapped to this higher level taxon
                    // we don't know how to interpret this label any more, so we'll drop that 
                    // leaf. The taxa will be included by other relationships (the taxonomy as
                    // a last resort), so we don't need to worry about losing leaves by skipping this...
                    if (lp->phyloChild->hasOttId()) {
                        LOG(DEBUG) << "scaff = " << scaffoldNode.getOttId() << " == phylo " << lp->phyloChild->getOttId();
                        assert(scaffoldNode.getOttId() == lp->phyloChild->getOttId());
                    }
                    sc.log(IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON, *lp->phyloChild);
                } else {
                    parEmbedding.loopEmbeddings[ebai.first].insert(lp);
                }
            } else {
                // if the anc isn't the parent, then it must pass through scaffoldNode's par
                assert(contains(parEmbedding.edgeBelowEmbeddings[ebai.first], lp));
            }
        }
    }
    for (auto child : iter_child(scaffoldNode)) {
        auto cit = sc.scaffold2NodeEmbedding.find(child);
        if (cit == sc.scaffold2NodeEmbedding.end()) {
            continue;
        }
        NodeEmbedding<T, U>& childEmbedding = cit->second;
        for (auto ceabi : childEmbedding.edgeBelowEmbeddings) {
            for (auto clp : ceabi.second) {
                if (clp->scaffoldAnc == &scaffoldNode) {
                    clp->scaffoldAnc = p;
                }
            }
        }
    }
    pruneCollapsedNode(scaffoldNode, sc);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::pruneCollapsedNode(U & scaffoldNode, SupertreeContextWithSplits & sc) {
    LOG(DEBUG) << "collapsed paths from ott" << scaffoldNode.getOttId() << ", adding child to parent";
    // NOTE: it is important that we add the children of scaffoldNode the left of its location
    //  in the tree so that the postorder traversal will not iterate over them.
    auto lc = scaffoldNode.getPrevSib();
    auto rc = scaffoldNode.getNextSib();
    scaffoldNode._detachThisNode();
    sc.detachedScaffoldNodes.insert(&scaffoldNode);
    auto p = scaffoldNode.getParent();
    const auto cv = scaffoldNode.getChildren();
    if (cv.empty()) {
        assert(false);
        return;
    }
    for (auto c : cv) {
        c->_setParent(p);
    }
    if (lc == nullptr) {
        p->_setLChild(cv[0]);
    } else {
        lc->_setNextSib(cv[0]);
    }
    (*cv.rbegin())->_setNextSib(rc);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::constructPhyloGraphAndCollapseIfNecessary(U & scaffoldNode, SupertreeContextWithSplits & sc) {
    LOG(DEBUG) << "constructPhyloGraphAndCollapseIfNecessary for " << scaffoldNode.getOttId();
    LOG(DEBUG) << "TEMP collapsing if conflict..." ;
    if (COLLAPSE_IF_CONFLICT) {
        collapseGroup(scaffoldNode, sc);
        return;
    }
    GreedyPhylogeneticForest<T,U> gpf;
    gpf.setPossibleMonophyletic(scaffoldNode);
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto laIt = loopEmbeddings.find(treeInd);
        const auto ebaIt = edgeBelowEmbeddings.find(treeInd);
        if (laIt == loopEmbeddings.end() && ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        /* order the groupings */
        std::map<OttIdSet, PathPairing<T,U> *> mapToProvideOrder;
        for (auto pp : laIt->second) {
            mapToProvideOrder[pp->getOttIdSet()] = pp;
        }
        for (auto pp : ebaIt->second) {
            mapToProvideOrder[pp->getOttIdSet()] = pp;
        }
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeInd);
        /* try to add groups bail out when we know that the possible group is not monophyletic */
        long bogusGroupIndex = 200000; // should get this from the node!
        for (auto mpoIt : mapToProvideOrder) {
            const auto & d = mpoIt.first;
            auto ppptr = mpoIt.second;
            gpf.attemptToAddGrouping(ppptr, d, relevantIds, static_cast<int>(treeInd), bogusGroupIndex++, sc);
            if (!gpf.possibleMonophyleticGroupStillViable()) {
                collapseGroup(scaffoldNode, sc);
                return;
            }
        }
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, sc);
}

template<typename T, typename U>
OttIdSet NodeEmbedding<T, U>::getRelevantDesIds(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                std::size_t treeIndex) {
    /* find MRCA of the phylo nodes */
    auto ippV = getAllIncomingPathPairs(nodeWithEmbedding, eForNd, treeIndex);
    OttIdSet relevantIds;
    for (auto pIt : ippV) {
        const OttIdSet otherRelevantIds = getRelevantDesIdsFromPath(*pIt);
        //std::cerr << "    for tree " << treeIndex << " otherRelevantIds returning"; writeOttSet(std::cerr, " ", otherRelevantIds, " "); std::cerr << '\n';
        relevantIds.insert(otherRelevantIds.begin(), otherRelevantIds.end());
    }
    std::cerr << " for tree " << treeIndex << " getRelevantDesIds returning"; writeOttSet(std::cerr, " ", relevantIds, " "); std::cerr << '\n';
    return relevantIds;
}

template<typename T, typename U>
bool NodeEmbedding<T, U>::reportIfContested(std::ostream & out,
                       const U * nd,
                       const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                       const std::vector<NodeWithSplits *> & aliasedBy,
                       bool verbose) const {
    if (isContested()) {
        auto c = getContestingTrees();
        for (auto cti : c) {
            auto ctree = treePtrByIndex.at(cti);
            const OttIdSet ls = getOttIdSetForLeaves(*ctree);
            const auto & edges = getEdgesExiting(cti);
            const std::string prefix = getContestedPreamble(*nd, *ctree);
            if (verbose) {
                reportOnConflicting(out, prefix, nd, edges, ls);
            } else {
                out << prefix << '\n';
            }
            for (auto na : aliasedBy) {
                const std::string p2 = getContestedPreamble(*na, *ctree);
                if (verbose) {
                    reportOnConflicting(out, p2, na, edges, ls);
                } else {
                    out << prefix << '\n';
                }
            }
        }
        return true;
    }
    return false;
}

template<typename T, typename U>
void reportOnConflicting(std::ostream & out,
                        const std::string & prefix,
                        const T * scaffold,
                        const std::set<PathPairing<T, U> *> & exitPaths,
                        const OttIdSet & phyloLeafSet) {
    if (exitPaths.size() < 2) {
        assert(false);
        return;
    }
    const auto scaffoldDes = set_intersection_as_set(scaffold->getData().desIds, phyloLeafSet);
    auto epIt = begin(exitPaths);
    const PathPairing<T, U> * ep = *epIt;
    const U * phyloPar = ep->phyloParent;
    const U * deepestPhylo = nullptr;
    std::map<OttIdSet, const U *> desIdSet2NdConflicting;
    if (isProperSubset(scaffoldDes, phyloPar->getData().desIds)) {
        deepestPhylo = phyloPar;
    } else {
        desIdSet2NdConflicting[phyloPar->getData().desIds] = phyloPar;
        for (auto anc : iter_anc_const(*phyloPar)) {
            if (isProperSubset(scaffoldDes, anc->getData().desIds)) {
                deepestPhylo = anc;
                break;
            }
            desIdSet2NdConflicting[anc->getData().desIds] = anc;
        }
        assert(deepestPhylo != nullptr);
    }
    for (++epIt; epIt != end(exitPaths); ++epIt) {
        const U * phyloNd  = (*epIt)->phyloChild;
        assert(phyloNd != nullptr);
        for (auto anc : iter_anc_const(*phyloNd)) {
            if (anc == deepestPhylo) {
                break;
            }
            desIdSet2NdConflicting[anc->getData().desIds] = anc;
        }
    }
    if (desIdSet2NdConflicting.empty()) {
        out << "VERY ODD " << prefix << ", desIdSet2NdConflicting but desIdSet2NdConflicting is empty()!\n";
        //assert(false); // not reachable if we are calling isContested first as a test.
        return;
    }
    for (const auto & mIt : desIdSet2NdConflicting) {
        const auto & di = mIt.first;
        auto nd = mIt.second;
        const OttIdSet e = set_difference_as_set(di, scaffoldDes);
        const OttIdSet m = set_difference_as_set(scaffoldDes, di);
        out << prefix;
        emitConflictDetails(out, *nd, e, m);
    }
}

template class NodePairing<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.
template class PathPairing<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.
template class NodeEmbedding<NodeWithSplits, NodeWithSplits>; // force explicit instantiaion of this template.

}// namespace
