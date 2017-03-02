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
bool NodeEmbedding<T, U>::debugNodeEmbedding(const char * tag, 
                                             bool isContested,
                                             const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    for (const auto  & t2exit : edgeBelowEmbeddings) {
        const auto treeIndex = t2exit.first;
        const auto & exitPaths =  t2exit.second;
        for (const auto epref : exitPaths) {
            const auto & sNERef = sn2ne.at(epref->scaffoldDes);
            assert(contains(sNERef.edgeBelowEmbeddings, treeIndex));
            assert(contains(sNERef.edgeBelowEmbeddings.at(treeIndex), epref));
            for (auto aNode : iter_anc(*(epref->scaffoldDes))) {
                if (aNode == epref->scaffoldAnc) {
                    break;
                }
                const auto & saNERef = sn2ne.at(epref->scaffoldDes);
                assert(contains(saNERef.edgeBelowEmbeddings, treeIndex));
                assert(contains(saNERef.edgeBelowEmbeddings.at(treeIndex), epref));
            }
        }
    }
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
                    LOG(DEBUG)  << " parentless " << getDesignator(*p) << ' ' << reinterpret_cast<long>(p);
                    if (p->getParent() != nullptr) {
                        auto gp = p->getParent();
                        LOG(ERROR)  << " Failing on scaffold ott" << embeddedNode->get_ott_id() << " treeIndex " << treeIndex;
                        LOG(ERROR)  << " phylo parent node (@" << reinterpret_cast<long>(p) << ") = " << getDesignator(*p) << ' ' << reinterpret_cast<long>(p) << " is not in loop or exit edges";
                        writeNewick(std::cerr, p);
                        std::cerr << std::endl;
                        LOG(ERROR)  << " the phylo child was (@" << reinterpret_cast<long>(c2p.first) << ") = " << getDesignator(*c2p.first);
                        writeNewick(std::cerr, c2p.first);
                        LOG(ERROR)  << " phylo grandparent (@" << reinterpret_cast<long>(gp) << ") = " << getDesignator(*gp) << ' ' << reinterpret_cast<long>(gp);
                        writeNewick(std::cerr, gp);
                        std::cerr << std::endl;
                        for (auto anc : iter_anc(*embeddedNode)) {
                            auto pathPairings = sn2ne.at(anc).refersToNode(treeIndex, gp);
                            if (!pathPairings.empty()) {
                                LOG(ERROR) << " ancestor ott" << anc->get_ott_id() << " refers to the parent";
                            }
                        }
                        debugPrint(*embeddedNode, treeIndex, sn2ne);
                        LOG(ERROR) << "Failure context is " << tag;
                        assert(false);
                    }
                    assert(parentsOfExits.empty());
                }
                root = p;
            }
        }
        if (parentsOfExits.empty()) {
            // no loops, the tree is just polytomy for this subproblem
            for (auto pathPtr : childExitForThisTree) {
                const auto & rids = pathPtr->get_ott_idSet();
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
                const auto & rids = pp->get_ott_idSet();
                assert(rids.size() == 1);
                assert(rids.size() == 1);
                assert(*rids.begin() != LONG_MAX);
            }
        }
    }
    return true;
}

template<typename T, typename U>
void NodeEmbedding<T, U>::set_ott_idForExitEmbeddings(
                    T * newScaffDes,
                    long ottId,
                    std::map<const T *, NodeEmbedding<T, U> > & n2ne) {
    for (auto treeInd2eout : edgeBelowEmbeddings) {
        assert(treeInd2eout.second.size() < 2);
        for (auto eout : treeInd2eout.second) {
            LOG(DEBUG) << "for tree " << treeInd2eout.first << " set_ott_id(" << ottId<< ')';
            eout->scaffoldDes = newScaffDes;
            eout->set_ott_idSet(ottId, n2ne);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::mergeExitEmbeddingsIfMultiple() {
    std::map<std::size_t, std::set<PathPairPtr> > toCull;
    for (auto treeInd2eout : edgeBelowEmbeddings) {
        if (treeInd2eout.second.size() > 1) {
            OTC_UNREACHABLE; // now resolving earlier...
            // If an input tree has a polytomy with members of a taxon as well as its "outgroup" taxa,
            //  then the polytomy does not contest the monophyly.
            //  this function will be called after resolution of the polytomy. 
            //  TODO: We might want to add a node for the resolved node should be added to the source tree to reflect
            //      that its ambiguity has been resolved in one particular way.
            //  The crucial thing for the continuation of the export/supertree operation
            //      is the set of exit nodes from embeddedNode be reduced to 1 exit node
            U * resolvedPolyParent = nullptr;
            for (auto eout : treeInd2eout.second) {
                if (resolvedPolyParent == nullptr) {
                    resolvedPolyParent = eout->phyloParent;
                } else {
                    assert(resolvedPolyParent == eout->phyloParent);
                    toCull[treeInd2eout.first].insert(eout); //TMP just keeping the first node
                }
            }
        }
    }
    for (const auto & toCullEl : toCull) {
        for (auto doomedPath : toCullEl.second) {
            removeRefToExitPath(toCullEl.first, doomedPath);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::resolveParentInFavorOfThisNode(
                T & scaffoldNode,
                std::size_t treeIndex,
                SupertreeContextWithSplits & sc) {
    auto exitForThisTreeIt = edgeBelowEmbeddings.find(treeIndex);
    if (exitForThisTreeIt == edgeBelowEmbeddings.end()) {
        return;
    }
    PathPairSet & exitSetForThisTree{exitForThisTreeIt->second};
    const auto origNumExits = exitSetForThisTree.size();
    if (origNumExits < 2) {
        return;
    }
    const TreeMappedWithSplits * cTreePtr = sc.treesByIndex.at(treeIndex);
    // this is the only way that we modify the embedded tree. So for now, we'll
    //      tolerate the const cast, which is in very poor form...
    // @TMP @TODO
    TreeMappedWithSplits * treePtr = const_cast<TreeMappedWithSplits *>(cTreePtr);
    LOG(INFO) << "Resolving treeIndex = " << treeIndex << " name = " << treePtr->get_name() << " for OTT " << scaffoldNode.get_ott_id() << '\n';
    const std::map<U *, U *> phyloNode2PhyloPar = getExitPhyloNd2Par(treeIndex);
    // resolve the phylo tree
    U * phPar = nullptr;
    for (auto p2p : phyloNode2PhyloPar) {
        if (phPar == nullptr) {
            phPar = p2p.second;
        } else {
            assert(phPar == p2p.second); // there should only be one par.
        }
    }
    assert(phPar != nullptr);
    U * insertedNodePtr = treePtr->createNode(phPar);
    for (auto p2p : phyloNode2PhyloPar) {
        U * phChild = p2p.first;
        phChild->detachThisNode();
        insertedNodePtr->addChild(phChild);
        const OttIdSet & cd = phChild->get_data().desIds;
        insertedNodePtr->get_data().desIds.insert(cd.begin(), cd.end());
    }
    T * scaffoldAncestor = nullptr;
    for (auto epp : exitSetForThisTree) {
        if (scaffoldAncestor == nullptr) {
            scaffoldAncestor = epp->scaffoldAnc;
        } else {
            assert(scaffoldAncestor == epp->scaffoldAnc); // to be a resolution case all exits must have the same parent
        }
    }
    // create the new node pairing and path pairing objects
    sc.nodePairingsFromResolve.emplace_back(NodePairingWithSplits(&scaffoldNode, insertedNodePtr));
    NodePairingWithSplits & newNodePairing{*sc.nodePairingsFromResolve.rbegin()};
    nodeEmbeddings[treeIndex].insert(&newNodePairing);
    sc.pathPairingsFromResolve.emplace_back(PathPairingWithSplits(scaffoldAncestor, phPar, newNodePairing));
    PathPairingWithSplits & newPathPairing{*sc.pathPairingsFromResolve.rbegin()};
    // fix every exit path to treat scaffoldNode as the scaffoldAnc node.
    // If the scaffoldDes is the scaffoldNode, then this exit is becoming a loop...
    PathPairSet toMoveToLoops;
    std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold2NodeEmbedding;
    for (auto epp : exitSetForThisTree) {
        assert(epp->phyloParent == phPar);
        epp->phyloParent = insertedNodePtr;
        assert(epp->scaffoldAnc != &scaffoldNode);
        for (auto n : iter_anc(scaffoldNode)) {
            if (n == epp->scaffoldAnc) {
                break;
            }
            sn2ne.at(n).removeRefToExitPath(treeIndex, epp);
        }
        epp->scaffoldAnc = &scaffoldNode;
        if (epp->scaffoldDes == &scaffoldNode) {
            toMoveToLoops.insert(epp);
        }
    }
    // all of the former exit paths for this tree are invalid...
    exitSetForThisTree.clear();
    // add the new loops...
    if (!toMoveToLoops.empty()) {
        auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end()) {
            loopEmbeddings[treeIndex] = toMoveToLoops;
        } else {
            laIt->second.insert(begin(toMoveToLoops), end(toMoveToLoops));
        }
    }
    // insert the new (and only) exit path for this node and its ancestors...
    edgeBelowEmbeddings[treeIndex].insert(&newPathPairing);
    for (auto n : iter_anc(scaffoldNode)) {
        if (n == scaffoldAncestor) {
            break;
        }
        sn2ne.at(n).edgeBelowEmbeddings[treeIndex].insert(&newPathPairing);
    }
    assert(edgeBelowEmbeddings[treeIndex].size() == 1);
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
        LOG(DEBUG) << "  updateOttIdSetNoTraversal for " << reinterpret_cast<long>(this) << " in ";
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
        LOG(DEBUG) << "  updateOttIdSetNoTraversal for " << reinterpret_cast<long>(this);
        dbWriteOttSet("updateOttIdSetNoTraversal exit ", currChildOttIdSet);
    }
    return true;
}

template<typename T, typename U>
inline const OttIdSet & NodeEmbedding<T, U>::getRelevantDesIdsFromPath(const PathPairing<T, U> & path) {
    return path.get_ott_idSet();
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
    assert("not implemented"[0] == 'f');; // until we check for "high ranks preserve contested monophyly optimization"
}

template<typename T, typename U>
void NodeEmbedding<T, U>::collapseSourceEdgesToForceOneEntry(T & ,
                                                             PathPairSet & pps,
                                                             std::size_t treeIndex,
                                                             SupertreeContextWithSplits & sc) {
    assert("not implemented"[0] == 'f');; // until we check for "high ranks preserve contested monophyly optimization"
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
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
        const auto ebaIt = edgeBelowEmbeddings.find(treeIndex);
        if (ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        PathPairSet & pps = ebaIt->second;
        collapseSourceEdgesToForceOneEntry(scaffoldNode, pps, treeIndex, sc);
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
                                       std::size_t treeIndex,
                                       long groupIndex);
const std::string getForestDOTFilename(const std::string & prefix,
                                       const DOTFileStep step,
                                       std::size_t treeIndex,
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
        base += std::to_string(treeIndex);
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

template<typename T>
void debugPrintNd2Par(const char * p, const T & m) {
    std::cerr << "debugPrintNd2Par " << p << ": ";
    for (auto mp : m) {
        std::cerr << "  " << reinterpret_cast<long>(mp.first) << " -> " << reinterpret_cast<long>(mp.second) << '\n';
    }
}
// Only resolves cases in which there is a deeper polytomy that is compatible with 
//  the taxon - resolved to make the taxon monophyletic...
// This allows export of deeper subproblems to be performed.
// returns the name (not the full path) of the tree file written or an empty string
template<typename T, typename U>
std::string NodeEmbedding<T, U>::exportSubproblemAndResolve(
                            T & scaffoldNode,
                            const std::string & exportDir,
                            std::ostream * exportStream,
                            SupertreeContextWithSplits & sc) {
    const std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold2NodeEmbedding;
    //debugNodeEmbedding("top of export", false, sn2ne);
    //debugPrint(scaffoldNode, 215, sn2ne);
    const OttIdSet EMPTY_SET;
    const auto scaffOTTId = scaffoldNode.get_ott_id();
    std::ofstream treeFileStream;
    std::ofstream provFileStream;
    std::ostream * treeExpStream = exportStream;
    std::ostream * provExpStream = exportStream;
    std::string retStr;
    if (exportStream == nullptr) {
        std::string outFilename = exportDir;
        outFilename.append("/ott");
        outFilename += std::to_string(scaffOTTId);
        outFilename += ".tre";
        retStr.append("ott");
        retStr += std::to_string(scaffOTTId);
        retStr += ".tre";
        std::string provOutFilename = exportDir;
        provOutFilename.append("/ott");
        provOutFilename += std::to_string(scaffOTTId);
        provOutFilename += "-tree-names.txt";
        treeFileStream.open(outFilename);
        if (!treeFileStream.good()) {
            std::string m = "Could not open \"";
            m += outFilename;
            m += '\"';
            throw OTCError(m);
        }
        treeExpStream = &treeFileStream;
        provFileStream.open(provOutFilename);
        if (!provFileStream.good()) {
            std::string m = "Could not open \"";
            m += provOutFilename;
            m += '\"';
            throw OTCError(m);
        }
        provExpStream = &provFileStream;
    } else {
        *exportStream << "ott" << scaffOTTId << '\n';
    }
    
    //TMP this could be done a lot more efficiently. Writing the trees through the GBF should
    // exercise some of the code for the scaffolded supertree operation. So this should help
    //  us find bugs in that code...
    
    OttIdSet totalLeafSet;
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
        const auto * treePtr = sc.treesByIndex.at(treeIndex);
        assert(treePtr != nullptr);
        //auto prelnd2par = getLoopedPhyloNd2Par(treeIndex);
        //debugPrintNd2Par("preloops", prelnd2par);
        //auto preend2par = getExitPhyloNd2Par(treeIndex);
        //debugPrintNd2Par("preexits", preend2par);
        resolveParentInFavorOfThisNode(scaffoldNode, treeIndex, sc);
        const auto childExitForThisTree = getAllChildExitPathsForTree(scaffoldNode,
                                                                      treeIndex,
                                                                      sn2ne);
        const auto tempDBIt = edgeBelowEmbeddings.find(treeIndex);
        assert(tempDBIt == edgeBelowEmbeddings.end() || tempDBIt->second.size() < 2);
        auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end() && childExitForThisTree.empty()) {
            continue;
        }
        const auto firstLaIt = laIt;
        LOG(INFO) << "     exporting ott" << scaffOTTId << "  treeIndex = " << treeIndex << " name = " << treePtr->get_name();
        const OttIdSet relevantIds = getRelevantDesIds(sn2ne, treeIndex);
        totalLeafSet.insert(begin(relevantIds), end(relevantIds));
        auto lnd2par = getLoopedPhyloNd2Par(treeIndex);
        //debugPrintNd2Par("loops", lnd2par);
        const auto shouldHaveBeenPrunedNd2Par = getUnEmbeddedPhyloNd2Par(treeIndex);
        std::map<NodeWithSplits *, long> nd2id;
        OttIdSet shouldHaveBeenPrunedIds;
        for (auto unembeddedPair : shouldHaveBeenPrunedNd2Par) {
            auto uDes = unembeddedPair.first;
            auto uPar = unembeddedPair.second;
            assert(!contains(lnd2par, uDes));
            lnd2par[uDes] = uPar;
            auto shbpId = uDes->get_ott_id();
            nd2id[uDes] = shbpId;
            shouldHaveBeenPrunedIds.insert(shbpId);
        }
        
        auto end2par = getExitPhyloNd2Par(treeIndex);
        //debugPrintNd2Par("exits", end2par);
        U * root = nullptr;
        std::set<U *> parentsOfExits;
        for (const auto & c2p : lnd2par) {
            auto p = c2p.second;
            if (p != root && !contains(lnd2par, p)) {
                if (contains(end2par, p)) {
                    parentsOfExits.insert(end2par.at(p));
                    assert(root == nullptr || !parentsOfExits.empty());
                    assert(parentsOfExits.size() < 2);
                } else {
                    assert(p->getParent() == nullptr);
                }
                root = p;
            }
        }
        if (root == nullptr && firstLaIt == loopEmbeddings.end()) {
            // no loops, the tree is just polytomy for this subproblem
            *treeExpStream << '(';
            OttIdSet ois = shouldHaveBeenPrunedIds;
            bool first = true;
            for (auto pathPtr : childExitForThisTree) {
                if (contains(sc.prunedSubtrees[treeIndex], pathPtr->phyloChild)) {
                    continue;
                }
                auto rids = pathPtr->get_ott_idSet();
                if (rids.size() != 1) {
                    LOG(DEBUG) << "crashing while exporting OTT ID " << scaffoldNode.get_ott_id() << " for tree " << treePtr->get_name();
                    dbWriteOttSet(" rids = ", rids);
                    assert(false);
                }
                long ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                ois.insert(ottId);
            }
            totalLeafSet.insert(ois.begin(), ois.end());
            for (long ottId : ois) {
                if (!first) {
                    *treeExpStream << ',';
                }
                *treeExpStream << "ott" << ottId;
                first = false;
            }
            *treeExpStream << ");\n";
        } else {
            LOG(DEBUG) << "Before adding child exits...";
            for (auto np : lnd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << getDesignator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << getDesignator(*np.second);
            }
            LOG(DEBUG) << " root = " << reinterpret_cast<long>(root)  << "  parentsOfExits.size() = " << parentsOfExits.size();
            assert(parentsOfExits.empty() || root != nullptr);
            assert(parentsOfExits.size() < 2 );
             // if the node is uncontested, then all exit paths must have the same parent
            U * const deeperNd = (parentsOfExits.empty() ? nullptr : *begin(parentsOfExits));
            std::set<NodeWithSplits *> toDel;
            for (auto n2pEl : lnd2par) {
                if (contains(sc.prunedSubtrees[treeIndex], n2pEl.first)) {
                    toDel.insert(n2pEl.first);
                }
            }
            for (auto td : toDel) {
                lnd2par.erase(td);
            }
            for (auto pp : childExitForThisTree) {
                if (contains(sc.prunedSubtrees[treeIndex], pp->phyloChild)) {
                    continue;
                }
                assert(!contains(lnd2par, pp->phyloChild));
                if (pp->phyloParent == root) {
                    lnd2par[pp->phyloChild] = root;
                } else {
                    if (deeperNd == pp->phyloParent) {
                        lnd2par[pp->phyloChild] = root;
                    } else {
                        if (!contains(lnd2par, pp->phyloParent)) {
                            if (!contains(end2par, pp->phyloParent)) {
                                LOG(ERROR) << "crashing on treeIndex " << treeIndex;
                                LOG(ERROR) << "deeperNd is set to  " << reinterpret_cast<long>(deeperNd);
                                debugPrintPathPairing(*pp);
                                assert(false);
                            }
                            // necessary if we have node maps to scaffold node and is also
                            //  the child of a polytomy that is compatible with scaffoldnode
                            auto ep = end2par[pp->phyloParent];
                            lnd2par[pp->phyloParent] = ep;
                            while (!contains(lnd2par, ep) && contains(end2par, ep)) {
                                auto nep = end2par[ep];
                                lnd2par[ep] = nep;
                                ep = nep;
                            }
                            if (root != nullptr && !contains(lnd2par, root) && contains(end2par, root)) {
                                lnd2par[root] = end2par[root];
                            }
                        }
                        lnd2par[pp->phyloChild] = pp->phyloParent;
                    }
                }
                auto rids = pp->get_ott_idSet();
                assert(rids.size() == 1);
                long ottId = *rids.begin();
                assert(ottId != LONG_MAX);
                totalLeafSet.insert(ottId);
                nd2id[pp->phyloChild] = ottId;
            }
            RootedTreeTopologyNoData toWrite;
            LOG(DEBUG) << "After adding child exits...";
            for (auto np : lnd2par) {
                LOG(DEBUG) << "   mapping " << np.first << " " << getDesignator(*np.first);
                LOG(DEBUG) << "   to par  " << np.second << " " << getDesignator(*np.second);
            }
            try {
                copyTreeStructure(lnd2par, nd2id, toWrite);
            } catch (const OTCError & ) {
                LOG(ERROR) << "could not construct a valid tree";
                debugPrint(scaffoldNode, treeIndex, sn2ne);
                debugPrintNd2Par("lnd2par", lnd2par);
                assert(false);
            }
            sortChildOrderByLowestDesOttId(toWrite.getRoot());
            writeTreeAsNewick(*treeExpStream, toWrite);
            *treeExpStream << "\n";
        }
        *provExpStream << treePtr->get_name() << "\n";
    }
    GreedyBandedForest<T, U> gpf{scaffoldNode.get_ott_id()};
    gpf.attemptToAddGrouping(totalLeafSet, EMPTY_SET, 0, 1, &sc);
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
        for (auto snc : iter_child(scaffoldNode)) {
            assert(snc != nullptr);
        }
    }
    gpf.finishResolutionOfEmbeddedClade(scaffoldNode, this, &sc);
    if (exportStream == nullptr) {
        provFileStream.close();
        treeFileStream.close();
    }
    //debugPrint(scaffoldNode, 7, sn2ne);
    //debugNodeEmbedding("leaving exportSubproblemAndResolve", false, sn2ne);
    return retStr;
}


template<typename T, typename U>
void NodeEmbedding<T, U>::resolveGivenUncontestedMonophyly(T & scaffoldNode,
                                                           SupertreeContextWithSplits & sc) {
    const OttIdSet EMPTY_SET;
    LOG(DEBUG) << "resolveGivenUncontestedMonophyly for " << scaffoldNode.get_ott_id();
    GreedyBandedForest<T, U> gpf{scaffoldNode.get_ott_id()};
    std::set<PathPairing<T, U> *> considered;
    const auto scaffOTTId = scaffoldNode.get_ott_id();
    std::string forestDOTfile = "forestForOTT";
    forestDOTfile += std::to_string(scaffOTTId);
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
        const auto laIt = loopEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end()) {
            continue;
        }
        LOG(INFO) << "      treeIndex = " << treeIndex;
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeIndex);
        PathPairSet & pps = laIt->second;
        // leaf set of this tree for this subtree
        // for repeatability, we'll try to add groupings in reverse order of desIds sets (deeper first)
        std::map<OttIdSet, PathPairing<T, U> *> mapToProvideOrder;
        for (auto pp : pps) {
            mapToProvideOrder[pp->get_ott_idSet()] = pp;
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
                gpf.attemptToAddGrouping(d, relevantIds, static_cast<int>(treeIndex), bogusGroupIndex, &sc);
                gpf.debugInvariantsCheck();
                if (scaffOTTId == ottIDBeingDebugged) {
                    gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                    gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, INFORMATIVE_SPLIT, treeIndex, bogusGroupIndex -1));
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
            gpf.addLeaf(*inc, relevantIds, static_cast<int>(treeIndex), bogusGroupIndex++, &sc);
            if (scaffOTTId == ottIDBeingDebugged) {
                gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, TRIVIAL_SPLIT, treeIndex, bogusGroupIndex -1));
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
                appendIncludeLeafSetAsNewick("phyloStatementAttempt", pathPtr->get_ott_idSet(), pathPtr->get_ott_idSet());
            }
            gpf.attemptToAddGrouping(pathPtr->get_ott_idSet(),
                                    EMPTY_SET,
                                    static_cast<int>(bogusTreeIndex),
                                    bogusGroupIndex++,
                                    &sc);
            if (scaffOTTId == ottIDBeingDebugged) {
                gpf.dumpAcceptedPhyloStatements("acceptedPhyloStatementOut.tre");
                gpf.writeForestDOTToFN(getForestDOTFilename(forestDOTfile, NONEMBEDDED_SPLIT, bogusTreeIndex, bogusGroupIndex -1));
            }
            gpf.debugInvariantsCheck();
            considered.insert(pathPtr); // @TMP not needed
        }
    }
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
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
void debugPrintPathPairing(const PathPairing<T,U> & clp) {
    std::cerr << "  PathPairing " << reinterpret_cast<long>(&clp) << " scaffoldAnc (@" << reinterpret_cast<long>(clp.scaffoldAnc) << ") = ott" << clp.scaffoldAnc->get_ott_id() << "\n      ";
    //writeNewick(std::cerr, clp.scaffoldAnc);
    std::cerr << "            scaffoldDes (@" << reinterpret_cast<long>(clp.scaffoldDes) << ") = ott" << clp.scaffoldDes->get_ott_id() << "\n      ";
    //writeNewick(std::cerr, clp.scaffoldDes);
    std::cerr << "            phyloParent (@" << reinterpret_cast<long>(clp.phyloParent) << ") = " << getDesignator(*clp.phyloParent) << "\n      ";
    //writeNewick(std::cerr, clp.phyloParent);
    std::cerr << "            phyloChild (@" << reinterpret_cast<long>(clp.phyloChild) << ") = " << getDesignator(*clp.phyloChild) << "\n";
    //writeNewick(std::cerr, clp.phyloChild);
    
}

template<typename T, typename U>
void NodeEmbedding<T, U>::debugPrint(T & scaffoldNode,
                                     std::size_t treeIndex,
                                     const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const {
    for (auto child : iter_child(scaffoldNode)) {
        auto & cne = sn2ne.at(child);
        auto cneIt = cne.edgeBelowEmbeddings.find(treeIndex);
        if (cneIt == cne.edgeBelowEmbeddings.end()) {
            continue;
        }
        std::cerr << " Child ott" << child->get_ott_id() << "\n";
        for (const auto & cneExit : cneIt->second) {
            debugPrintPathPairing(*cneExit);
        }
    }
    auto & sne = sn2ne.at(&scaffoldNode);
    std::cerr << " Loops for parent ott" << scaffoldNode.get_ott_id() << "\n";
    auto sneIt = sne.loopEmbeddings.find(treeIndex);
    if (sneIt != sne.loopEmbeddings.end()) {
        for (auto & snExit : sneIt->second) {
            debugPrintPathPairing(*snExit);
        }
    }
    std::cerr << " Exits for parent ott" << scaffoldNode.get_ott_id() << "\n";
    sneIt = sne.edgeBelowEmbeddings.find(treeIndex);
    if (sneIt != sne.edgeBelowEmbeddings.end()) {
        for (auto & snExit : sneIt->second) {
            debugPrintPathPairing(*snExit);
        }
    }
}

template<typename T, typename U>
void NodeEmbedding<T, U>::pruneSuppressed(std::size_t treeIndex, U * phyloPar, U * phyloChild) {
    if (phyloPar->getOutDegree() == 2) {
        auto & pps = loopEmbeddings.at(treeIndex);
        PathPairSet toDel;
        for (auto pp : pps) {
            if (pp->phyloChild == phyloPar) {
                toDel.insert(pp);
            }
        }
        assert(toDel.size() < 2);
        for (auto pp : toDel) {
            pps.erase(pp);
        }
    }
    phyloPar->removeChild(phyloChild);
}

template<typename T, typename U>
void registerAsPruned(std::size_t treeIndex, T * phyloChild, U & sc) {
    auto & ps = sc.prunedSubtrees[treeIndex];
    ps.insert(phyloChild);
    auto p = phyloChild->getParent();
    for (auto c : iter_child(*p)) {
        if (!contains(ps, c)) {
            return;
        }
    }
    registerAsPruned(treeIndex, p, sc);
}
template<typename T, typename U>
void NodeEmbedding<T, U>::collapseGroup(T & scaffoldNode, SupertreeContext<T, U> & sc) {
    assert(&scaffoldNode == embeddedNode);
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
    std::map<const T *, NodeEmbedding<T, U> > & sn2ne = sc.scaffold2NodeEmbedding;
    NodeEmbedding<T, U> & parEmbedding = const_cast<NodeEmbedding<T, U> &>(sn2ne.at(p));
    LOG(DEBUG) << "TOP of collapseGroup";
    //parEmbedding.debugPrint(scaffoldNode, 7, sn2ne);
    //const auto beforePL = copyAllLoopPathPairing(p, sn2ne);
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
    std::set<std::size_t> indsOfTreesWithNewLoops;
    std::set<std::size_t> indsOfTreesMappedToInternal;
    // every exit edge for this node becomes a loop for its parent if it is not trivial
    for (auto ebai : edgeBelowEmbeddings) {
        const auto & treeIndex = ebai.first;
        std::set<PathPairPtr> pathsAblated;
        for (auto lp : ebai.second) {
            if (lp->phyloChild->isTip()
                && lp->scaffoldDes == &scaffoldNode) {
                // this only happens if a terminal was mapped to this higher level taxon
                // we don't know how to interpret this label any more, so we should drop that 
                // leaf. Unless specifically requested not to by the user.
                // The taxa will be included by other relationships (the taxonomy as
                // a last resort), so we don't need to worry about losing leaves by skipping this...
                if (sc.pruneTipsMappedToContestedTaxa) {
                    LOG(INFO) << "Ablating path leading to ott" << lp->phyloChild->get_ott_id();;
                    LOG(DEBUG) << "IGNORING scaff = " << scaffoldNode.get_ott_id() << " == phylo " << lp->phyloChild->get_ott_id();
                    assert(scaffoldNode.get_ott_id() == lp->phyloChild->get_ott_id());
                    sc.log(IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON, *lp->phyloChild);
                    OttIdSet innerOTTId;
                    innerOTTId.insert(lp->phyloChild->get_ott_id());
                    OttIdSet n = scaffoldNode.get_data().desIds; // expand the internal name to it taxonomic content
                    n.erase(lp->phyloChild->get_ott_id());
                    lp->updateDesIdsForSelfAndAnc(innerOTTId, n, sn2ne);
                    indsOfTreesMappedToInternal.insert(treeIndex);
                    pathsAblated.insert(lp);
                    registerAsPruned(treeIndex, lp->phyloChild, sc);
                } else {
                    phyloNd2ParForUnembeddedTrees[treeIndex][lp->phyloChild] = lp->phyloParent;
                }
            } else if (lp->scaffoldAnc == p) {
                //if (lp->scaffoldDes == &scaffoldNode) {
                if ((!lp->phyloChild->isTip()) && lp->scaffoldDes == &scaffoldNode) {
                    lp->scaffoldDes = p;
                    parEmbedding.loopEmbeddings[treeIndex].insert(lp);
                    indsOfTreesWithNewLoops.insert(treeIndex);
                }
                /*} else {
                    if (lp->scaffoldDes->getParent() != &scaffoldNode) {
                        LOG(ERROR) << " Anbandonding a path that ends at ott" << lp->scaffoldDes->get_ott_id();
                        LOG(ERROR) << "                    and starts at ott" << lp->scaffoldAnc->get_ott_id() << " when collapsing ott" << scaffoldNode.get_ott_id();
                        LOG(ERROR) << " lp->scaffoldDes->getParent() ott" << lp->scaffoldDes->getParent()->get_ott_id();
                        assert(false);
                    }
                }*/
            } else {
                // if the anc isn't the parent, then it must pass through scaffoldNode's par
                assert(contains(parEmbedding.edgeBelowEmbeddings[treeIndex], lp));
                if (lp->scaffoldDes == &scaffoldNode && !lp->phyloChild->isTip()) {
                    lp->scaffoldDes = p;
                }
            }
        }
        if (!pathsAblated.empty()) {
            for (auto deadPathPtr : pathsAblated) {
                debugPrintPathPairing(*deadPathPtr);
                /*T * ablatedScaffAnc = deadPathPtr->scaffoldAnc;
                U * ablatedPhyloChild = deadPathPtr->phyloChild;
                U * ablatedPhyloPar = deadPathPtr->phyloParent;
                U * curr = &scaffoldNode;
                while (curr != deadPathPtr->scaffoldAnc) {
                    sn2ne.at(curr).removeRefToExitPath(treeIndex, deadPathPtr);
                    curr = curr->getParent();
                }
                //sn2ne.at(ablatedScaffAnc).pruneSuppressed(treeIndex, ablatedPhyloPar, ablatedPhyloChild);
                */
            }
        }
    }
    for (auto child : iter_child(scaffoldNode)) {
        auto cit = sn2ne.find(child);
        if (cit == sn2ne.end()) {
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
    const auto afterPL = copyAllLoopPathPairing(p, sn2ne);

    for (auto bpl : beforePL) {
        const auto & afterVal = afterPL.at(bpl.first);
        const auto tv = parEmbedding.getAllIncomingPathPairs(sn2ne, bpl.first);
        assert(isSubset(bpl.second, afterVal));
        for (auto pp : bpl.second) {
            const PathPairing<T, U> * cpp = pp;
            LOG(DEBUG) << "Checking tree" << bpl.first << " scaff" << p->get_ott_id() << " " << cpp->scaffoldAnc->get_ott_id() << " -> " <<  cpp->scaffoldDes->get_ott_id();
            assert(vcontains(tv, cpp));
        }
    }
    */
}

template<typename T, typename U>
void NodeEmbedding<T, U>::pruneCollapsedNode(T & scaffoldNode, SupertreeContextWithSplits & sc) {
    checkAllNodePointersIter(scaffoldNode);
     LOG(DEBUG) << "collapsed paths from ott" << scaffoldNode.get_ott_id() << ", adding child to parent";
    // NOTE: it is important that we add the children of scaffoldNode the left of its location
    //  in the tree so that the postorder traversal will not iterate over them.
    assert(scaffoldNode.hasChildren());
    auto p = scaffoldNode.getParent();
    assert(p);
    if (!phyloNd2ParForUnembeddedTrees.empty()) {
        auto & scaffoldPar = sc.scaffold2NodeEmbedding.at(p);
        // We have a hacky data structure to pass along to our parent if
        //  we were asked to retain terminal nodes that are mapped to this
        //  these tips will no longer be embedded in the scaffold... Ugh.
        for (auto pmfu : phyloNd2ParForUnembeddedTrees) {
            if (pmfu.second.empty()) {
                continue;
            }
            const auto treeIndex = pmfu.first;
            auto & ppmfuForTree =  scaffoldPar.phyloNd2ParForUnembeddedTrees[treeIndex];
            for (auto dpmfuForTree : pmfu.second) {
                ppmfuForTree[dpmfuForTree.first] = dpmfuForTree.second;
            }
        }
    }
    while(scaffoldNode.hasChildren())
    {
        auto n = scaffoldNode.getFirstChild();
        n->detachThisNode();
        scaffoldNode.addSibOnLeft(n);
    }
    scaffoldNode.detachThisNode();
    sc.scaffoldTree.markAsDetached(&scaffoldNode);
    sc.detachedScaffoldNodes.insert(&scaffoldNode);
}

template<typename T, typename U>
void NodeEmbedding<T, U>::constructPhyloGraphAndCollapseIfNecessary(T & scaffoldNode, SupertreeContextWithSplits & sc) {
    assert(&scaffoldNode == embeddedNode);
    LOG(DEBUG) << "constructPhyloGraphAndCollapseIfNecessary for " << scaffoldNode.get_ott_id();
    LOG(DEBUG) << "TEMP collapsing if conflict..." ;
    if (COLLAPSE_IF_CONFLICT) {
        //debugNodeEmbedding("before collapse", true, sc.scaffold2NodeEmbedding);
        collapseGroup(scaffoldNode, sc);
        //debugNodeEmbedding("after collapse", true, sc.scaffold2NodeEmbedding);
        return;
    }
    GreedyBandedForest<T, U> gpf{scaffoldNode.get_ott_id()};
    gpf.setPossibleMonophyletic(scaffoldNode);
    for (std::size_t treeIndex = 0 ; treeIndex < sc.numTrees; ++treeIndex) {
        const auto laIt = loopEmbeddings.find(treeIndex);
        const auto ebaIt = edgeBelowEmbeddings.find(treeIndex);
        if (laIt == loopEmbeddings.end() && ebaIt == edgeBelowEmbeddings.end()) {
            continue;
        }
        LOG(INFO) << "      treeIndex = " << treeIndex;
        /* order the groupings */
        std::map<OttIdSet, PathPairing<T, U> *> mapToProvideOrder;
        for (auto pp : laIt->second) {
            mapToProvideOrder[pp->get_ott_idSet()] = pp;
        }
        for (auto pp : ebaIt->second) {
            mapToProvideOrder[pp->get_ott_idSet()] = pp;
        }
        const OttIdSet relevantIds = getRelevantDesIds(sc.scaffold2NodeEmbedding, treeIndex);
        /* try to add groups bail out when we know that the possible group is not monophyletic */
        long bogusGroupIndex = 200000; // should get this from the node!
        for (auto mpoIt : mapToProvideOrder) {
            const auto & d = mpoIt.first;
            LOG(INFO) << "        bogusGroupIndex = " << bogusGroupIndex;
            gpf.attemptToAddGrouping(d, relevantIds, static_cast<int>(treeIndex), bogusGroupIndex++, &sc);
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
        auto c = getContestingTreeIndices();
        for (auto cti : c) {
            auto ctree = treePtrByIndex.at(cti);
            const OttIdSet ls = get_ott_idSetForLeaves(*ctree);
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
    const auto scaffoldDes = set_intersection_as_set(scaffold->get_data().desIds, phyloLeafSet);
    auto epIt = begin(exitPaths);
    const PathPairing<T, U> * ep = *epIt;
    const U * phyloPar = ep->phyloParent;
    const U * deepestPhylo = nullptr;
    std::map<OttIdSet, const U *> desIdSet2NdConflicting;
    if (isProperSubset(scaffoldDes, phyloPar->get_data().desIds)) {
        deepestPhylo = phyloPar;
    } else {
        desIdSet2NdConflicting[phyloPar->get_data().desIds] = phyloPar;
        for (auto anc : iter_anc_const(*phyloPar)) {
            if (isProperSubset(scaffoldDes, anc->get_data().desIds)) {
                deepestPhylo = anc;
                break;
            }
            desIdSet2NdConflicting[anc->get_data().desIds] = anc;
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
            desIdSet2NdConflicting[anc->get_data().desIds] = anc;
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
