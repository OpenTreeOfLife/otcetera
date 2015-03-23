#include <queue>
#include "otc/greedy_forest.h"
#include "otc/node_embedding.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/debug.h"
#include "otc/supertree_context.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/error.h"
namespace otc {
constexpr bool COLLAPSE_IF_CONFLICT = true;

template<typename T, typename U>
bool NodeEmbedding<T, U>::debugNodeEmbedding(bool isContested,
                                            const std::map<const T *, NodeEmbedding<T, U> > &sn2ne) const {
    for (const auto & t2l : loopEmbeddings) {
        auto treeIndex = t2l.first;
        const auto childExitForThisTree = getAllChildExitPathsForTree(embeddedNode,
                                                                      treeIndex,
                                                                      sn2ne);
        auto lnd2par = getLoopedPhyloNd2Par(treeIndex);
        auto end2par = getExitPhyloNd2Par(treeIndex);
        // all of the loop phyloParents should be in a child of another
        //  loop or a child of an exit path.
        // If the node is uncontested, then 
        U * root = nullptr;
        std::set<U *> parentsOfExits;
        for (const auto & c2p : lnd2par) {
            auto p = c2p.second;
            if (p != root && !contains(lnd2par, p)) {
                if (contains(end2par, p)) {
                    parentsOfExits.insert(end2par.at(p));
                    assert(root == nullptr || !parentsOfExits.empty());
                    assert(isContested || parentsOfExits.size() < 2);
                } else {
                    LOG(DEBUG)  << " parentless " << getDesignator(*p) << ' ' << (long) p;
                    assert(p->getParent() == nullptr);
                    assert(parentsOfExits.empty());
                }
                root = p;
            }
        }
        if (parentsOfExits.empty()) {
            // no loops, the tree is just polytomy for this subproblem
            for (auto pathPtr : childExitForThisTree) {
                const auto & rids = pathPtr->getOttIdSet();
                assert(rids.size() == 1);
                assert(*rids.begin() != LONG_MAX);
            }
        } else {
            for (auto pp : childExitForThisTree) {
                if (pp->phyloParent != root) {
                    parentsOfExits.insert(pp->phyloParent);
                    assert(isContested || parentsOfExits.size() < 2);
                }
                assert(!contains(lnd2par, pp->phyloChild));
                if (contains(end2par, pp->phyloChild)) {
                    assert(!contains(lnd2par, pp->phyloParent));
                } else {
                    assert(contains(lnd2par, pp->phyloParent));
                }
                const auto & rids = pp->getOttIdSet();
                assert(rids.size() == 1);
                assert(rids.size() == 1);
                assert(*rids.begin() != LONG_MAX);
            }
        }
    }
    return true;
}



template<typename T, typename U>
void NodeEmbedding<T, U>::setOttIdForExitEmbeddings(
                    T * newScaffDes,
                    long ottId,
                    std::map<const T *, NodeEmbedding<T, U> > & n2ne) {
    for (auto treeInd2eout : edgeBelowEmbeddings) {
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " setOttId(" << ottId<< ')';
            eout->scaffoldDes = newScaffDes;
            eout->setOttIdSet(ottId, n2ne);
        }
    }
}

// Returns all loop paths for nd and all edgeBelowEmbeddings of its children
template<typename T, typename U>
std::vector<const PathPairing<T, U> *>
NodeEmbedding<T, U>::getAllIncomingPathPairs(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                 std::size_t treeIndex) const {
    const T *nd = embeddedNode;
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
                r.push_back(e);
            }
        } else {
                //LOG(DEBUG) << "     No embedding found for this tree";
        }
    }
    return r;
}

template<typename T, typename U>
bool PathPairing<T, U>::updateOttIdSetNoTraversal(const OttIdSet & oldEls, const OttIdSet & newEls) {
    if (false && debuggingOutputEnabled) {
        LOG(DEBUG) << "  updateOttIdSetNoTraversal for " << (long)this << " in ";
        dbWriteOttSet("currChildOttIdSet", currChildOttIdSet);
        dbWriteOttSet("oldEls", oldEls);
        dbWriteOttSet("newEls", newEls);
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
    if (false && debuggingOutputEnabled) {
        LOG(DEBUG) << "  updateOttIdSetNoTraversal for " << (long)this;
        dbWriteOttSet("updateOttIdSetNoTraversal exit ", currChildOttIdSet);
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
    return relevantIds;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseSourceEdge(const T * , //phyloParent,
                                                 PathPairing<T, U> * ) { //path
    NOT_IMPLEMENTED; // until we check for "high ranks preserve contested monophyly optimization"
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseSourceEdgesToForceOneEntry(T & ,
                                                             PathPairSet & pps,
                                                             std::size_t treeIndex,
                                                             SupertreeContextWithSplits & sc) {
    NOT_IMPLEMENTED; // until we check for "high ranks preserve contested monophyly optimization"
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
void NodeEmbedding<T, U>::resolveGivenContestedMonophyly(T & scaffoldNode,
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
std::set<PathPairing<T, U> *> NodeEmbedding<T, U>::getAllChildExitPaths(
                const T & scaffoldNode,
                const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    std::set<PathPairing<T, U> *> r;
    for (auto c : iter_child_const(scaffoldNode)) {
        const auto & thr = sn2ne.at(c);
        for (auto te : thr.edgeBelowEmbeddings) {
            r.insert(begin(te.second), end(te.second));
        }
    }
    return r;
}

enum DOTFileStep {
    INFORMATIVE_SPLIT,
    TRIVIAL_SPLIT,
    NONEMBEDDED_SPLIT
};
const std::string getForestDOTFilename(const std::string & prefix,
                                       const DOTFileStep step,
                                       std::size_t treeInd,
                                       long groupIndex);
const std::string getForestDOTFilename(const std::string & prefix,
                                       const DOTFileStep step,
                                       std::size_t treeInd,
                                       long groupIndex) {
    std::string base = prefix;
    if (step == INFORMATIVE_SPLIT) {
        base += "AfterInfTree";
    } else if (step == TRIVIAL_SPLIT) {
        base += "AfterTrivTree";
    } else if (step == NONEMBEDDED_SPLIT) {
        base += "AfterXChildSplit";
    }
    if (step == INFORMATIVE_SPLIT || step == TRIVIAL_SPLIT || step == NONEMBEDDED_SPLIT) {
        base += std::to_string(treeInd);
        base += "Split";
        base += std::to_string(groupIndex);
    }
    base += ".dot";
    return base;
}

template<typename T, typename U>
std::set<PathPairing<T, U> *> NodeEmbedding<T, U>::getAllChildExitPathsForTree(
                const T & scaffoldNode,
                std::size_t treeIndex,
                const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    std::set<PathPairing<T, U> *> r;
    for (auto c : iter_child_const(scaffoldNode)) {
        const auto & thr = sn2ne.at(c);
        const auto & tebeIt = thr.edgeBelowEmbeddings.find(treeIndex);
        if (tebeIt != thr.edgeBelowEmbeddings.end()) {
            r.insert(begin(tebeIt->second), end(tebeIt->second));
        }
    }
    return r;
}
// By "fake" the resolution here we mean that on exit, this
//  subproblem's outgoing edges should be labeled with the scaffoldNode's OTT ID, just
//  like the outgoing edges would look upon successful completion of the resolution of
//  the clade. This allows export of deeper subproblems to be performed.
template<typename T, typename U>
void NodeEmbedding<T, U>::exportSubproblemAndFakeResolution(
                            T & scaffoldNode,
                            const std::string & exportDir,
                            std::ostream * exportStream,
                            SupertreeContextWithSplits & sc) {
    debugNodeEmbedding(false, sc.scaffold2NodeEmbedding);
    const OttIdSet EMPTY_SET;
    LOG(DEBUG) << "exportSubproblemAndFakeResolution for " << scaffoldNode.getOttId();
    const auto scaffOTTId = scaffoldNode.getOttId();
    std::ofstream treeFileStream;
    std::ofstream provFileStream;
    std::ostream * treeExpStream = exportStream;
    std::ostream * provExpStream = exportStream;
    if (exportStream == nullptr) {
        std::string outFilename = exportDir;
        outFilename.append("/ott");
        outFilename += std::to_string(scaffOTTId);
        outFilename += ".tre";
        std::string provOutFilename = exportDir;
        provOutFilename.append("/ott");
        provOutFilename += std::to_string(scaffOTTId);
        provOutFilename += "-tree-names.txt";
        treeFileStream.open(outFilename);
        treeExpStream = &treeFileStream;
        provFileStream.open(provOutFilename);
        provExpStream = &provFileStream;
    }
    
    //TMP this could be done a lot more efficiently. Writing the trees through the GBF should
    // exercise some of the code for the scaffolded supertree operation. So this should help
    //  us find bugs in that code...
    
    OttIdSet totalLeafSet;
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto * treePtr = sc.treesByIndex.at(treeInd);
        assert(treePtr != nullptr);
        const auto childExitForThisTree = getAllChildExitPathsForTree(scaffoldNode,
                                                                      treeInd,
                                                                      sc.scaffold2NodeEmbedding);
        auto laIt = loopEmbeddings.find(treeInd);
        if (laIt == loopEmbeddings.end() && childExitForThisTree.empty()) {
            continue;
        }
        const auto firstLaIt = laIt;
        LOG(INFO) << "      treeInd = " << treeInd;
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeInd);
        totalLeafSet.insert(begin(relevantIds), end(relevantIds));
        auto nd2par = getLoopedPhyloNd2Par(treeInd);
        unsigned numParentsWOParentsInMap = 0;
        NodeWithSplits * root = nullptr;
        for (auto c2p : nd2par) {
            auto p = c2p.second;
            if (p != root && !contains(nd2par, p)) {
                ++numParentsWOParentsInMap;
                LOG(DEBUG)  << " parentless " << getDesignator(*p) << ' ' << (long) p;
                root = p;
            }
        }
        assert(numParentsWOParentsInMap < 2);
        if (numParentsWOParentsInMap == 0) {
            assert(firstLaIt == loopEmbeddings.end());
            // no loops, the tree is just polytomy for this subproblem
            *treeExpStream << '(';
            bool first = true;
            for (auto pathPtr : childExitForThisTree) {
                auto rids = pathPtr->getOttIdSet();
                assert(rids.size() == 1);
                long ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                totalLeafSet.insert(ottId);
                if (!first) {
                    *treeExpStream << ',';
                }
                *treeExpStream << "ott" << ottId;
                first = false;
            }
            *treeExpStream << ");\n";
        } else {
            NodeWithSplits * deeperNd = nullptr; // nodes passing through this taxon must all have the same parent (or the taxon would be contested)
            std::map<NodeWithSplits *, long> nd2id;
            LOG(DEBUG) << "Before adding child exits...";
            for (auto np : nd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << getDesignator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << getDesignator(*np.second);
            }
            for (auto pp : childExitForThisTree) {
                assert(!contains(nd2par, pp->phyloChild));
                if (pp->phyloParent == root) {
                    nd2par[pp->phyloChild] = root;
                } else {
                    if (!contains(nd2par, pp->phyloParent)) {
                        assert(deeperNd == nullptr || deeperNd == pp->phyloParent);
                        deeperNd = pp->phyloParent;
                        nd2par[pp->phyloChild] = root;
                    } else {
                        nd2par[pp->phyloChild] = pp->phyloParent;
                    }
                }
                auto rids = pp->getOttIdSet();
                assert(rids.size() == 1);
                long ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                totalLeafSet.insert(ottId);
                nd2id[pp->phyloChild] = ottId;
            }
            RootedTreeTopologyNoData toWrite;
            LOG(DEBUG) << "After adding child exits...";
            for (auto np : nd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << getDesignator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << getDesignator(*np.second);
            }
            copyTreeStructure(nd2par, nd2id, toWrite);
            writeTreeAsNewick(*treeExpStream, toWrite);
            *treeExpStream << "\n";
        }
        *provExpStream << treePtr->getName() << "\n";
    }
    GreedyBandedForest<T, U> gpf{scaffoldNode.getOttId()};
    gpf.attemptToAddGrouping(totalLeafSet, EMPTY_SET, 0, 1, &sc);
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        for (auto snc : iter_child(scaffoldNode)) {
            assert(snc != nullptr);
        }
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, &sc);
    if (exportStream == nullptr) {
        provFileStream.close();
        treeFileStream.close();
    }
}


template<typename T, typename U>
void NodeEmbedding<T, U>::resolveGivenUncontestedMonophyly(T & scaffoldNode,
                                                           SupertreeContextWithSplits & sc) {
    const OttIdSet EMPTY_SET;
    LOG(DEBUG) << "resolveGivenUncontestedMonophyly for " << scaffoldNode.getOttId();
    GreedyBandedForest<T, U> gpf{scaffoldNode.getOttId()};
    std::set<PathPairing<T, U> *> considered;
    const auto scaffOTTId = scaffoldNode.getOttId();
    std::string forestDOTfile = "forestForOTT";
    forestDOTfile += std::to_string(scaffOTTId);
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto laIt = loopEmbeddings.find(treeInd);
        if (laIt == loopEmbeddings.end()) {
            continue;
        }
        LOG(INFO) << "      treeInd = " << treeInd;
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeInd);
        PathPairSet & pps = laIt->second;
        // leaf set of this tree for this subtree
        // for repeatability, we'll try to add groupings in reverse order of desIds sets (deeper first)
        std::map<OttIdSet, PathPairing<T, U> *> mapToProvideOrder;
        for (auto pp : pps) {
            mapToProvideOrder[pp->getOttIdSet()] = pp;
        }
        long bogusGroupIndex = 0; // should get this from the node!
        typedef std::pair<const OttIdSet *, PathPairing<T, U> *>  q_t;
        std::queue<q_t> trivialQ;
        for (auto mpoIt = mapToProvideOrder.rbegin(); mpoIt != mapToProvideOrder.rend(); ++mpoIt) {
            auto ppptr = mpoIt->second;
            if (ppptr->pathIsNowTrivial()) {
                dbWriteOttSet("pathIsNowTrivial", mpoIt->first);
                const q_t toQ{&(mpoIt->first), ppptr};
                trivialQ.push(toQ);
            } else {
                const auto & d = mpoIt->first;
                gpf.debugInvariantsCheck();
                if (scaffOTTId == ottIDBeingDebugged) {
                    appendIncludeLeafSetAsNewick("phyloStatementAttempt", d, relevantIds);
                }
                LOG(INFO) << "        bogusGroupIndex = " << bogusGroupIndex << " out of " << mapToProvideOrder.size() << " (some of which may be skipped as trivial)";
                gpf.attemptToAddGrouping(d, relevantIds, static_cast<int>(treeInd), bogusGroupIndex, &sc);
                gpf.debugInvariantsCheck();
                if (scaffOTTId == ottIDBeingDebugged) {
                    gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                    gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, INFORMATIVE_SPLIT, treeInd, bogusGroupIndex -1));
                }
                considered.insert(ppptr);
            }
            bogusGroupIndex++;
        }
        while (!trivialQ.empty()) {
            const q_t triv = trivialQ.front();
            const OttIdSet * inc = triv.first;
            const auto ppptr = triv.second;
            if (scaffOTTId == ottIDBeingDebugged) {
                appendIncludeLeafSetAsNewick("phyloStatementAttempt", *inc, relevantIds);
            }
            gpf.addLeaf(*inc, relevantIds, static_cast<int>(treeInd), bogusGroupIndex++, &sc);
            if (scaffOTTId == ottIDBeingDebugged) {
                gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, TRIVIAL_SPLIT, treeInd, bogusGroupIndex -1));
            }
            gpf.debugInvariantsCheck();
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
    std::size_t bogusTreeIndex = 123456; // should get this from the node!
    long bogusGroupIndex = 100000; // should get this from the node!
    auto childExitPaths = getAllChildExitPaths(scaffoldNode, sc.scaffold2NodeEmbedding);
    for (auto pathPtr : childExitPaths) {
        if (!contains(considered, pathPtr)) {
            if (scaffOTTId == ottIDBeingDebugged) {
                appendIncludeLeafSetAsNewick("phyloStatementAttempt", pathPtr->getOttIdSet(), pathPtr->getOttIdSet());
            }
            gpf.attemptToAddGrouping(pathPtr->getOttIdSet(), EMPTY_SET, (int)bogusTreeIndex, bogusGroupIndex++, &sc);
            if (scaffOTTId == ottIDBeingDebugged) {
                gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, NONEMBEDDED_SPLIT, bogusTreeIndex, bogusGroupIndex -1));
            }
            gpf.debugInvariantsCheck();
            considered.insert(pathPtr); // @TMP not needed
        }
    }
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        for (auto snc : iter_child(scaffoldNode)) {
            assert(snc != nullptr);
        }
    }
    if (scaffOTTId == ottIDBeingDebugged) {
        std::string fn = forestDOTfile;
        fn += "BeforeFinalize.dot";
        gpf.writeForestDOTToFN(fn);
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, &sc);
    if (scaffOTTId == ottIDBeingDebugged) {
        std::string fn = forestDOTfile;
        fn += "AfterFinalize.dot";
        gpf.writeForestDOTToFN(fn);
    }
}

template<typename T, typename U>
std::map<std::size_t, std::set<PathPairing<T, U> *> > copyAllLoopPathPairing(const T *nd, const std::map<const T *, NodeEmbedding<T, U> > & eForNd) {
    const NodeEmbedding<T, U> & ne = eForNd.at(nd);
    return ne.loopEmbeddings;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseGroup(T & scaffoldNode, SupertreeContext<T, U> & sc) {
    sc.log(COLLAPSE_TAXON, scaffoldNode);
    assert(!scaffoldNode.isTip());
    U * p = scaffoldNode.getParent();
    assert(p != nullptr); // can't disagree with the root !
    // remap all nodes in NodePairing to parent
    for (auto nai : nodeEmbeddings) {
        for (auto np : nai.second) {
            np->scaffoldNode = p;
        }
    }
    NodeEmbedding<T, U>& parEmbedding = sc.scaffold2NodeEmbedding.at(p);
    //const auto beforePL = copyAllLoopPathPairing(p, sc.scaffold2NodeEmbedding);
    // every loop for this node becomes a loop for its parent
    for (auto lai : loopEmbeddings) {
        const auto & treeIndex = lai.first;
        for (auto lp : lai.second) {
            assert(lp->scaffoldDes == &scaffoldNode);
            assert(lp->scaffoldAnc == &scaffoldNode);
            lp->scaffoldDes = p;
            lp->scaffoldAnc = p;
            parEmbedding.loopEmbeddings[treeIndex].insert(lp);
        }
    }
    // every exit edge for this node becomes a loop for its parent if it is not trivial
    for (auto ebai : edgeBelowEmbeddings) {
        const auto & treeIndex = ebai.first;
        for (auto lp : ebai.second) {
            if (lp->scaffoldAnc == p) {
                if (lp->scaffoldDes == &scaffoldNode) {
                    if(lp->phyloChild->isTip()) {
                        // this only happens if a terminal was mapped to this higher level taxon
                        // we don't know how to interpret this label any more, so we'll drop that 
                        // leaf. The taxa will be included by other relationships (the taxonomy as
                        // a last resort), so we don't need to worry about losing leaves by skipping this...
                        LOG(DEBUG) << "IGNORING scaff = " << scaffoldNode.getOttId() << " == phylo " << lp->phyloChild->getOttId();
                        assert(scaffoldNode.getOttId() == lp->phyloChild->getOttId());
                        sc.log(IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON, *lp->phyloChild);
                        OttIdSet innerOTTId;
                        innerOTTId.insert(lp->phyloChild->getOttId());
                        OttIdSet n = scaffoldNode.getData().desIds; // expand the internal name to it taxonomic content
                        n.erase(lp->phyloChild->getOttId());
                        lp->updateDesIdsForSelfAndAnc(innerOTTId, n, sc.scaffold2NodeEmbedding);
                        // we'll let the path pairing be lost (since it is attached to a node that will be detached...)
                    } else {
                        lp->scaffoldDes = p;
                        parEmbedding.loopEmbeddings[treeIndex].insert(lp);
                    }
                }
            } else {
                // if the anc isn't the parent, then it must pass through scaffoldNode's par
                assert(contains(parEmbedding.edgeBelowEmbeddings[treeIndex], lp));
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
    /*
    const auto afterPL = copyAllLoopPathPairing(p, sc.scaffold2NodeEmbedding);

    for (auto bpl : beforePL) {
        const auto & afterVal = afterPL.at(bpl.first);
        const auto tv = parEmbedding.getAllIncomingPathPairs(sc.scaffold2NodeEmbedding, bpl.first);
        assert(isSubset(bpl.second, afterVal));
        for (auto pp : bpl.second) {
            const PathPairing<T, U> * cpp = pp;
            LOG(DEBUG) << "Checking tree" << bpl.first << " scaff" << p->getOttId() << " " << cpp->scaffoldAnc->getOttId() << " -> " <<  cpp->scaffoldDes->getOttId();
            assert(vcontains(tv, cpp));
        }
    }
    */
}

template<typename T, typename U>
void NodeEmbedding<T, U>::pruneCollapsedNode(T & scaffoldNode, SupertreeContextWithSplits & sc) {
    checkAllNodePointersIter(scaffoldNode);
    LOG(DEBUG) << "collapsed paths from ott" << scaffoldNode.getOttId() << ", adding child to parent";
    // NOTE: it is important that we add the children of scaffoldNode the left of its location
    //  in the tree so that the postorder traversal will not iterate over them.
    assert(!scaffoldNode.isTip());
    auto entrySib = scaffoldNode.getPrevSib();
    assert(entrySib == nullptr || entrySib->getNextSib() == &scaffoldNode);
    auto exitSib = scaffoldNode.getNextSib();
    auto p = scaffoldNode.getParent();
    const auto cv = scaffoldNode.getChildren();
    scaffoldNode._detachThisNode();
    sc.scaffoldTree.markAsDetached(&scaffoldNode);
    sc.detachedScaffoldNodes.insert(&scaffoldNode);
    if (cv.empty()) {
        assert(false);
        throw OTCError("disabled assert is false");
    }
    for (auto c : cv) {
        c->_setParent(p);
    }
    auto firstMovingChild = cv[0];
    auto lastMovingChild = *cv.rbegin();
    if (entrySib == nullptr) {
        p->_setFirstChild(firstMovingChild);
    } else {
        assert(entrySib->getNextSib() == exitSib); // _detachThisNode guarantees this
        entrySib->_setNextSib(firstMovingChild);
    }
    lastMovingChild->_setNextSib(exitSib);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::constructPhyloGraphAndCollapseIfNecessary(T & scaffoldNode, SupertreeContextWithSplits & sc) {
    LOG(DEBUG) << "constructPhyloGraphAndCollapseIfNecessary for " << scaffoldNode.getOttId();
    LOG(DEBUG) << "TEMP collapsing if conflict..." ;
    if (COLLAPSE_IF_CONFLICT) {
        debugNodeEmbedding(true, sc.scaffold2NodeEmbedding);
        collapseGroup(scaffoldNode, sc);
        debugNodeEmbedding(true, sc.scaffold2NodeEmbedding);
        return;
    }
    GreedyBandedForest<T, U> gpf{scaffoldNode.getOttId()};
    gpf.setPossibleMonophyletic(scaffoldNode);
    for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
        const auto laIt = loopEmbeddings.find(treeInd);
        const auto ebaIt = edgeBelowEmbeddings.find(treeInd);
        if (laIt == loopEmbeddings.end() && ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        LOG(INFO) << "      treeInd = " << treeInd;
        /* order the groupings */
        std::map<OttIdSet, PathPairing<T, U> *> mapToProvideOrder;
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
            LOG(INFO) << "        bogusGroupIndex = " << bogusGroupIndex;
            gpf.attemptToAddGrouping(d, relevantIds, static_cast<int>(treeInd), bogusGroupIndex++, &sc);
            if (!gpf.possibleMonophyleticGroupStillViable()) {
                collapseGroup(scaffoldNode, sc);
                return;
            }
            gpf.debugInvariantsCheck();
        }
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, &sc);
}

template<typename T, typename U>
OttIdSet NodeEmbedding<T, U>::getRelevantDesIds(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                std::size_t treeIndex) {
    /* find MRCA of the phylo nodes */
    auto ippV = getAllIncomingPathPairs(eForNd, treeIndex);
    OttIdSet relevantIds;
    for (auto pIt : ippV) {
        const OttIdSet otherRelevantIds = getRelevantDesIdsFromPath(*pIt);
        relevantIds.insert(otherRelevantIds.begin(), otherRelevantIds.end());
    }
    return relevantIds;
}

template<typename T, typename U>
bool NodeEmbedding<T, U>::reportIfContested(std::ostream & out,
                       const T * nd,
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
        throw OTCError("asserts are disabled, but one is not true");
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
        throw OTCError("Expecting isContested to have guarded call to reportOnConflicting");
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
