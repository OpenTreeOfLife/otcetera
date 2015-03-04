#include "otc/otcli.h"
#include "otc/debug.h"
using namespace otc;
constexpr bool COLLAPSE_IF_CONFLICT = true;
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &);


//currently not copying names
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &tree) {
    TreeMappedWithSplits * rawTreePtr = new TreeMappedWithSplits();
    try {
        NodeWithSplits * newRoot = rawTreePtr->createRoot();
        auto r = tree.getRoot();
        assert(r->hasOttId());
        newRoot->setOttId(r->getOttId());
        std::map<const NodeWithSplits *, NodeWithSplits *> templateToNew;
        templateToNew[r]= newRoot;
        std::map<long, NodeWithSplits *> & newMap = rawTreePtr->getData().ottIdToNode;
        rawTreePtr->getData().desIdSetsContainInternals = tree.getData().desIdSetsContainInternals;
        for (auto nd : iter_pre_const(tree)) {
            auto p = nd->getParent();
            if (p == nullptr) {
                continue;
            }
            auto t2nIt = templateToNew.find(p);
            assert(t2nIt != templateToNew.end());
            auto ntp = t2nIt->second;
            auto nn = rawTreePtr->createChild(ntp);
            assert(templateToNew.find(nd) == templateToNew.end());
            templateToNew[nd] = nn;
            if (nd->hasOttId()) {
                nn->setOttId(nd->getOttId());
                newMap[nd->getOttId()] = nn;
            } else {
                assert(false);
            }
            nn->getData().desIds = nd->getData().desIds;
        }
    } catch (...) {
        delete rawTreePtr;
        throw;
    }
    return std::unique_ptr<TreeMappedWithSplits>(rawTreePtr);
}


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

template<typename T, typename U>
class PathPairing {
    public:
    T * scaffoldDes;
    T * scaffoldAnc;
    U * phyloChild;
    U * phyloParent;
    
    PathPairing(const NodePairing<T,U> & parent, const NodePairing<T,U> & child)
        :scaffoldDes(child.scaffoldNode),
        scaffoldAnc(parent.scaffoldNode),
        phyloChild(child.phyloNode),
        phyloParent(parent.phyloNode) {
    }
    std::set<long> getPhyloChildDesID() const {
        return phyloChild->getData().desIds;
    }
};


template<typename T, typename U>
void reportOnConflicting(std::ostream & out, const std::string & prefix, const T * scaff, const std::set<PathPairing<T, U> *> & exitPaths, const std::set<long> & phyloLeafSet) {
    if (exitPaths.size() < 2) {
        assert(false);
        return;
    }
    const auto scaffDes = set_intersection_as_set(scaff->getData().desIds, phyloLeafSet);
    auto epIt = begin(exitPaths);
    const PathPairing<T, U> * ep = *epIt;
    const U * phyloPar = ep->phyloParent;
    const U * deepestPhylo = nullptr;
    std::map<std::set<long>, const U *> desIdSet2NdConflicting;
    if (isProperSubset(scaffDes, phyloPar->getData().desIds)) {
        deepestPhylo = phyloPar;
    } else {
        desIdSet2NdConflicting[phyloPar->getData().desIds] = phyloPar;
        for (auto anc : iter_anc_const(*phyloPar)) {
            if (isProperSubset(scaffDes, anc->getData().desIds)) {
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
        const std::set<long> e = set_difference_as_set(di, scaffDes);
        const std::set<long> m = set_difference_as_set(scaffDes, di);
        out << prefix;
        emitConflictDetails(out, *nd, e, m);
    }
}


template<typename T, typename U> class NodeThreading;
//template<typename T, typename U>
//using ThreadingObj = NodeThreading<T, U>;

template<typename T, typename U> class SupertreeContext;

template<typename T, typename U>
class RootedForest {
    public:
        bool empty() const {
            return roots.empty();
        }
        void addGroupToNewTree(const std::set<long> & ingroup, const std::set<long> & leafSet);
    private:
        RootedTree<T, U> nodeSrc;
        RootedTreeNode<T> * createNewRoot();
    protected:
        std::list<RootedTreeNode<T> *> roots;
        std::set<long> ottSet;
};

template<typename T, typename U>
class GreedyPhylogeneticForest: public RootedForest<RTSplits, MappedWithSplitsData> {
    public:
    void attemptToAddGrouping(PathPairing<T, U> * ppptr,
                              const std::set<long> & ingroup,
                              const std::set<long> & leafSet,
                              const SupertreeContext<T, U> &sc);
    void finalizeTree() {
        assert(false);
    }
    void setPossibleMonophyletic(U & scaffoldNode) {
        assert(false);
    }
    bool possibleMonophyleticGroupStillViable() {
        assert(false);
    }
    void resolveThreadedClade(U & scaffoldNode, NodeThreading<T, U> * );
    private:
};


template<typename T, typename U>
inline void GreedyPhylogeneticForest<T,U>::resolveThreadedClade(U & ,
                                                                NodeThreading<T, U> *) {
    assert(false);
}
enum SupertreeCtorEvent {
    COLLAPSE_TAXON,
    IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON,
};
template<typename T, typename U>
class SupertreeContext {
    public:
        using LogEvent = std::pair<SupertreeCtorEvent, std::string>;
        mutable std::list<LogEvent> events;

        const std::size_t numTrees;
        void log(SupertreeCtorEvent e, const U & node) const {
            if (e == COLLAPSE_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.getOttId())});
            } else if (e == IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.getOttId())});
            }
        }
        std::map<const NodeWithSplits *, NodeThreading<T, U> > & scaff2NodeThreading;
        SupertreeContext(std::size_t nt, std::map<const NodeWithSplits *, NodeThreading<T, U> > & taxoToAlignment)
            :numTrees(nt),
            scaff2NodeThreading(taxoToAlignment) {
        }
};

using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;
using NodePairingWithSplits = NodePairing<NodeWithSplits, NodeWithSplits>;
using PathPairingWithSplits = PathPairing<NodeWithSplits, NodeWithSplits>;
using NodeThreadingWithSplits = NodeThreading<NodeWithSplits, NodeWithSplits>;
const std::set<long> EMPTY_SET;
// does NOT add the node to roots
template<typename T, typename U>
RootedTreeNode<T> * RootedForest<T,U>::createNewRoot() {
    assert(!empty());
    auto firstRoot = *roots.begin();
    auto n = nodeSrc.createChild(firstRoot);
    // detach the node
    n->_setParent(nullptr);
    auto s = n->getPrevSib();
    if (s == nullptr) {
        firstRoot->_setLChild(nullptr);
    } else {
        s->_setNextSib(nullptr);
    }
}
template<typename T, typename U>
void RootedForest<T,U>::addGroupToNewTree(const std::set<long> & ingroup,
                                     const std::set<long> & leafSet) {
    assert(ingroup.size() < leafSet.size());
    assert(ingroup.size() > 0);
    assert(isProperSubset(ingroup, leafSet));
    assert(areDisjoint(leafSet, ottSet));
    RootedTreeNode<T> * r = nullptr;
    if (empty()) {
        r = nodeSrc.createRoot();
    } else {
        r = createNewRoot();
    }
    roots.push_back(r);
    RTreeOttIDMapping<RTSplits> & trd = nodeSrc.getData();
    auto c = nodeSrc.createChild(r);
    const auto outgroup = set_difference_as_set(leafSet, ingroup);
    assert(outgroup.size() > 0);
    for (auto i : outgroup) {
        auto n = nodeSrc.createChild(r);
        n->setOttId(i);
        trd.ottIdToNode[i] = n;
        n->getData().desIds.insert(i);
    }
    for (auto i : ingroup) {
        auto n = nodeSrc.createChild(c);
        n->setOttId(i);
        trd.ottIdToNode[i] = n;
        n->getData().desIds.insert(i);
    }
    r->getData().desIds = leafSet;
    c->getData().desIds = ingroup;
    ottSet.insert(leafSet.begin(), leafSet.end());
}

template<typename T, typename U>
inline void GreedyPhylogeneticForest<T,U>::attemptToAddGrouping(PathPairing<T, U> * ppptr,
                                                                const std::set<long> & ingroup,
                                                                const std::set<long> & leafSet,
                                                                const SupertreeContext<T,U> &sc) {
    if (this->empty()) { // first grouping, always add...
        addGroupToNewTree(ingroup, leafSet);
    }
}
template<typename T, typename U>
class NodeThreading {
    public:
    using NodePairSet = std::set<NodePairing<T, U> *>;
    using PathPairSet = std::set<PathPairing<T, U> *>;
    
    std::map<std::size_t, NodePairSet> nodeAlignments;
    std::map<std::size_t, PathPairSet> edgeBelowAlignments;
    std::map<std::size_t, PathPairSet > loopAlignments;
    std::size_t getTotalNumNodeMappings() const {
        unsigned long t = 0U;
        for (auto i : nodeAlignments) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getTotalNumLoops() const {
        unsigned long t = 0U;
        for (auto i : loopAlignments) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getTotalNumEdgeBelowTraversals() const {
        unsigned long t = 0U;
        for (auto i : edgeBelowAlignments) {
            t += i.second.size();
        }
        return t;
    }
    static bool treeContestsMonophyly(const PathPairSet & edgesBelowForTree) {
        if (edgesBelowForTree.size() > 1) {
            const T * firstSrcPar = nullptr;
            for (auto pp : edgesBelowForTree) {
                auto sp = pp->phyloParent;
                if (sp != firstSrcPar) {
                    if (firstSrcPar == nullptr) {
                        firstSrcPar = sp;
                    } else {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    bool isContested() const {
        for (auto i : edgeBelowAlignments) {
            if (treeContestsMonophyly(i.second)) {
                return true;
            }
        }
        return false;
    }
    std::list<std::size_t> getContestingTrees() const {
        std::list<std::size_t> r;
        for (auto i : edgeBelowAlignments) {
            if (treeContestsMonophyly(i.second)) {
                r.push_back(i.first);
            }
        }
        return r;
    }
    const PathPairSet & getEdgesExiting(std::size_t treeIndex) const {
        auto el = edgeBelowAlignments.find(treeIndex);
        assert(el != edgeBelowAlignments.end());
        return el->second;
    }

    // some trees contest monophyly. Return true if these trees are obviously overruled
    //   by higher ranking trees so that we can avoid the more expensive unconstrained phylo graph 
    bool highRankingTreesPreserveMonophyly(std::size_t ) {
        return false;
    }

    std::set<long> getRelevantDesIds(const PathPairSet & pps) {
        std::set<long> relevantIds;
        for (auto path : pps) {
            const auto pc = path->phyloChild;
            const auto & cdi = pc->getData().desIds;
            relevantIds.insert(begin(cdi), end(cdi));
        }
        return relevantIds;
    }

    void collapseSourceEdge(const T * phyloParent, PathPairing<T, U> *path) {
        assert(false);
    }
    // there may be 
    void collapseSourceEdgesToForceOneEntry(U & scaffoldNode, PathPairSet & pps) {
        if (pps.size() < 2) {
            return;
        }
        auto relevantIds = getRelevantDesIds(pps);
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
    void resolveGivenContestedMonophyly(U & scaffoldNode, const SupertreeContext<T, U> & sc) {
        for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
            const auto ebaIt = edgeBelowAlignments.find(treeInd);
            if (ebaIt == edgeBelowAlignments.end()) {
                continue;
            }
            PathPairSet & pps = ebaIt->second;
            collapseSourceEdgesToForceOneEntry(scaffoldNode, pps);
        }
        resolveGivenUncontestedMonophyly(scaffoldNode, sc);
    }
    std::set<PathPairing<T, U> *> getAllChildExitPaths(U & scaffoldNode, const SupertreeContext<T, U> & sc) {
        std::set<PathPairing<T, U> *> r;
        for (auto c : iter_child(scaffoldNode)) {
            const auto & thr = sc.scaff2NodeThreading.at(c);
            for (auto te : thr.edgeBelowAlignments) {
                r.insert(begin(te.second), end(te.second));
            }
        }
        return r;
    }
    void resolveGivenUncontestedMonophyly(U & scaffoldNode, const SupertreeContext<T, U> & sc) {
        LOG(DEBUG) << "resolveGivenUncontestedMonophyly for " << scaffoldNode.getOttId();
        GreedyPhylogeneticForest<T,U> gpf;
        std::set<PathPairing<T, U> *> considered;
        for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
            const auto laIt = loopAlignments.find(treeInd);
            if (laIt == loopAlignments.end()) {
                continue;
            }

            PathPairSet & pps = laIt->second;
            // leaf set of this tree for this subtree
            std::set<long> relevantIds = getRelevantDesIds(pps);
            // for repeatability, we'll try to add groupings in reverse order of desIds sets (deeper first)
            std::map<std::set<long>, PathPairing<T,U> *> mapToProvideOrder;
            for (auto pp : pps) {
                mapToProvideOrder[pp->getPhyloChildDesID()] = pp;
            }
            for (auto mpoIt : mapToProvideOrder) {
                const auto & d = mpoIt.first;
                auto ppptr = mpoIt.second;
                gpf.attemptToAddGrouping(ppptr, d, relevantIds, sc);
                considered.insert(ppptr);
            }
        }
        // we might have missed some descendants  - any child that is has
        //  "scaffoldNode" as its threaded parent, but which is not involved
        //  any loop or exiting edges.
        //  This means that we have no info on the placement of such nodes.
        //      so we'll just attach them here.
        //  First step: get the list of paths for the children.
        auto childExitPaths = getAllChildExitPaths(scaffoldNode, sc);
        for (auto pathPtr : childExitPaths) {
            if (!contains(considered, pathPtr)) {
                gpf.attemptToAddGrouping(pathPtr, pathPtr->getPhyloChildDesID(), EMPTY_SET, sc);
                considered.insert(pathPtr); // @TMP not needed
            }
        }
        for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
            for (auto sc : iter_child(scaffoldNode)) {

            }
        }
        gpf.finalizeTree();
        gpf.resolveThreadedClade(scaffoldNode, this);
    }

    void collapseGroup(U & scaffoldNode, const SupertreeContext<T,U> & sc) {
        sc.log(COLLAPSE_TAXON, scaffoldNode);
        U * p = scaffoldNode.getParent();
        assert(p != nullptr); // can't disagree with the root !
        // remap all nodes in NodePairing to parent
        for (auto nai : nodeAlignments) {
            for (auto np : nai.second) {
                np->scaffoldNode = p;
            }
        }
        NodeThreading<T, U>& parThreading = sc.scaff2NodeThreading.at(p);
        // every loop for this node becomes a loop for its parent
        for (auto lai : loopAlignments) {
            for (auto lp : lai.second) {
                if (lp->scaffoldDes) {
                    lp->scaffoldDes = p;
                }
                parThreading.loopAlignments[lai.first].insert(lp);
            }
        }
        // every exit edge for this node becomes a loop for its parent if it is not trivial
        for (auto ebai : edgeBelowAlignments) {
            for (auto lp : ebai.second) {
                if (lp->scaffoldAnc == p) {
                    if (lp->scaffoldDes == &scaffoldNode) {
                        // this only happens if a terminal was mapped to this higher level taxon
                        // we don't know how to interpret this labe any more, so we'll drop that 
                        // leaf. The taxa will be included by other relationships (the taxonomy as
                        // a last resort), so we don't need to worry about losing leaves by skipping this...
                        assert(scaffoldNode.getOttId() == lp->phyloChild->getOttId());
                        sc.log(IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON, *lp->phyloChild);
                    } else {
                        parThreading.loopAlignments[ebai.first].insert(lp);
                    }
                } else {
                    // if the anc isn't the parent, then it must pass through scaffNode's par
                    assert(contains(parThreading.edgeBelowAlignments[ebai.first], lp));
                }
            }
        }
    }
    void constructPhyloGraphAndCollapseIfNecessary(U & scaffoldNode, const SupertreeContext<T, U> & sc) {
        LOG(DEBUG) << "constructPhyloGraphAndCollapseIfNecessary for " << scaffoldNode.getOttId();
        LOG(DEBUG) << "TEMP collapsing if conflict..." ;
        if (COLLAPSE_IF_CONFLICT) {
            collapseGroup(scaffoldNode, sc);
            return;
        }
        GreedyPhylogeneticForest<T,U> gpf;
        gpf.setPossibleMonophyletic(scaffoldNode);
        for (std::size_t treeInd = 0 ; treeInd < sc.numTrees; ++treeInd) {
            const auto laIt = loopAlignments.find(treeInd);
            const auto ebaIt = edgeBelowAlignments.find(treeInd);
            if (laIt == loopAlignments.end() && ebaIt == edgeBelowAlignments.end()) {
                continue;
            }
            /* find MRCA of the phylo nodes */
            std::set<long> relevantIds;
            if (laIt != loopAlignments.end()) {
                relevantIds = getRelevantDesIds(laIt->second);
            }
            if (ebaIt != edgeBelowAlignments.end()) {
                std::set<long> otherRelevantIds = getRelevantDesIds(ebaIt->second);
                relevantIds.insert(otherRelevantIds.begin(), otherRelevantIds.end());
            }
            /* order the groupings */
            std::map<std::set<long>, PathPairing<T,U> *> mapToProvideOrder;
            for (auto pp : laIt->second) {
                mapToProvideOrder[pp->getPhyloChildDesID()] = pp;
            }
            for (auto pp : ebaIt->second) {
                mapToProvideOrder[pp->getPhyloChildDesID()] = pp;
            }
            /* try to add groups bail out when we know that the possible group is not monophyletic */
            for (auto mpoIt : mapToProvideOrder) {
                const auto & d = mpoIt.first;
                auto ppptr = mpoIt.second;
                gpf.attemptToAddGrouping(ppptr, d, relevantIds, sc);
                if (!gpf.possibleMonophyleticGroupStillViable()) {
                    collapseGroup(scaffoldNode, sc);
                    return;
                }
            }
        }
        gpf.finalizeTree();
        gpf.resolveThreadedClade(scaffoldNode, this);
    }

    bool reportIfContested(std::ostream & out,
                           const U * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const {
        if (isContested()) {
            auto c = getContestingTrees();
            for (auto cti : c) {
                auto ctree = treePtrByIndex.at(cti);
                const std::set<long> ls = getOttIdSetForLeaves(*ctree);
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
    
};




class ThreadedTree {
    protected:
    std::list<NodePairingWithSplits> nodePairings;
    std::list<PathPairingWithSplits> pathPairings;
    std::map<const NodeWithSplits *, NodeThreadingWithSplits> taxoToAlignment;
    public:

    NodePairingWithSplits * _addNodeMapping(NodeWithSplits *taxo, NodeWithSplits *nd, std::size_t treeIndex) {
        assert(taxo != nullptr);
        assert(nd != nullptr);
        nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
        auto ndPairPtr = &(*nodePairings.rbegin());
        auto & athreading = taxoToAlignment[taxo];
        athreading.nodeAlignments[treeIndex].insert(ndPairPtr);
        return ndPairPtr;
    }
    PathPairingWithSplits * _addPathMapping(NodePairingWithSplits * parentPairing,
                                            NodePairingWithSplits * childPairing,
                                            std::size_t treeIndex) {
        pathPairings.emplace_back(*parentPairing, *childPairing);
        auto pathPairPtr = &(*pathPairings.rbegin());
        // register a pointer to the path at each traversed...
        auto currTaxo = pathPairPtr->scaffoldDes;
        auto ancTaxo = pathPairPtr->scaffoldAnc;
        if (currTaxo != ancTaxo) {
            while (currTaxo != ancTaxo) {
                taxoToAlignment[currTaxo].edgeBelowAlignments[treeIndex].insert(pathPairPtr);
                currTaxo = currTaxo->getParent();
                if (currTaxo == nullptr) {
                    break;
                }
            }
        } else {
            taxoToAlignment[currTaxo].loopAlignments[treeIndex].insert(pathPairPtr);
        }
        return pathPairPtr;
    }
    void threadNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
        // do threading
        std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
        std::set<NodePairingWithSplits *> tipPairings;
        for (auto nd : iter_post(tree)) {
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            NodePairingWithSplits * ndPairPtr = nullptr;
            NodeWithSplits * taxoDes = nullptr;
            if (nd->isTip()) {
                assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
                assert(nd->hasOttId());
                auto ottId = nd->getOttId();
                taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
                assert(taxoDes != nullptr);
                ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
                for (auto former : tipPairings) {
                    if (areLinearlyRelated(taxoDes, former->scaffoldNode)) {
                        std::string m = "Repeated or nested OTT ID in tip mapping of an input tree: \"";
                        m += nd->getName();
                        m += "\" and \"";
                        m += former->phyloNode->getName();
                        m += "\" found.";
                        throw OTCError(m);
                    }
                }
                tipPairings.insert(ndPairPtr);
                currTreeNodePairings[nd] = ndPairPtr;
            } else {
                auto reuseNodePairingIt = currTreeNodePairings.find(nd);
                assert(reuseNodePairingIt != currTreeNodePairings.end());
                ndPairPtr = reuseNodePairingIt->second;
                taxoDes = ndPairPtr->scaffoldNode;
                assert(taxoDes != nullptr);
            }
            NodePairingWithSplits * parPairPtr = nullptr;
            auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
            if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
                const auto & parDesIds = par->getData().desIds;
                auto taxoAnc = searchAncForMRCAOfDesIds(taxoDes, parDesIds);
                assert(taxoAnc != nullptr);
                parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
                currTreeNodePairings[par] = parPairPtr;
            } else {
                parPairPtr = prevAddedNodePairingIt->second;
            }
            _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
        }
    }

    void threadTaxonomyClone(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
        // do threading
        std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
        for (auto nd : iter_post(tree)) {
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            NodePairingWithSplits * ndPairPtr = nullptr;
            NodeWithSplits * taxoDes = nullptr;
            if (nd->isTip()) {
                assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
                assert(nd->hasOttId());
                auto ottId = nd->getOttId();
                taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
                assert(taxoDes != nullptr);
                ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
                currTreeNodePairings[nd] = ndPairPtr;
            } else {
                auto reuseNodePairingIt = currTreeNodePairings.find(nd);
                assert(reuseNodePairingIt != currTreeNodePairings.end());
                ndPairPtr = reuseNodePairingIt->second;
                taxoDes = ndPairPtr->scaffoldNode;
                assert(taxoDes != nullptr);
            }
            NodePairingWithSplits * parPairPtr = nullptr;
            auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
            if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
                auto pottId = par->getOttId(); // since it is a taxonomy, it will have internal node labels
                auto taxoAnc = scaffoldTree.getData().getNodeForOttId(pottId);
                assert(taxoAnc != nullptr);
                parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
                currTreeNodePairings[par] = parPairPtr;
            } else {
                parPairPtr = prevAddedNodePairingIt->second;
            }
            _addPathMapping(parPairPtr, ndPairPtr, treeIndex);
        }
    }

    void writeDOTExport(std::ostream & out,
                           const NodeThreading<NodeWithSplits, NodeWithSplits> & thr,
                           const NodeWithSplits * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex) const {
        writeNewick(out, nd);
        out << "nd.ottID = " << nd->getOttId() << " --> " << (nd->getParent() ? nd->getParent()->getOttId() : 0L) << "\n";
        out << "  getTotalNumNodeMappings = " << thr.getTotalNumNodeMappings() << "\n";
        out << "  getTotalNumLoops = " << thr.getTotalNumLoops() << "\n";
        out << "  getTotalNumEdgeBelowTraversals = " << thr.getTotalNumEdgeBelowTraversals() << "\n";
        out << "  isContested = " << thr.isContested() << "\n";
    }

};


class RemapToDeepestUnlistedState
    : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits>,
    public ThreadedTree {
    public:
    int numErrors;
    std::set<long> tabooIds;
    std::map<std::unique_ptr<TreeMappedWithSplits>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedWithSplits *> treePtrByIndex;
    bool doReportAllContested;
    bool doConstructSupertree;
    std::list<long> idsListToReportOn;
    std::list<long> idListForDotExport;
    TreeMappedWithSplits * taxonomyAsSource;

    void resolveOrCollapse(NodeWithSplits * scaffNd, const SupertreeContextWithSplits & sc) {
        auto & thr = taxoToAlignment[scaffNd];
        if (thr.isContested()) {
            if (thr.highRankingTreesPreserveMonophyly(sc.numTrees)) {
                thr.resolveGivenContestedMonophyly(*scaffNd, sc);
            } else {
                thr.constructPhyloGraphAndCollapseIfNecessary(*scaffNd, sc);
            }
        } else {
            thr.resolveGivenUncontestedMonophyly(*scaffNd, sc);
        }
    }
    void constructSupertree() {
        const auto numTrees = treePtrByIndex.size();
        const SupertreeContextWithSplits sc{numTrees, taxoToAlignment};
        for (auto nd : iter_post_internal(*taxonomy)) {
            resolveOrCollapse(nd, sc);
        }
    }


    virtual ~RemapToDeepestUnlistedState(){}
    RemapToDeepestUnlistedState()
        :TaxonomyDependentTreeProcessor<TreeMappedWithSplits>(),
         numErrors(0),
         doReportAllContested(false),
         doConstructSupertree(false),
         taxonomyAsSource(nullptr) {
    }

    void reportAllConflicting(std::ostream & out, bool verbose) {
        std::map<std::size_t, unsigned long> nodeMappingDegree;
        std::map<std::size_t, unsigned long> passThroughDegree;
        std::map<std::size_t, unsigned long> loopDegree;
        unsigned long totalContested = 0;
        unsigned long redundContested = 0;
        unsigned long totalNumNodes = 0;
        for (auto nd : iter_node_internal(*taxonomy)) {
            const auto & thr = taxoToAlignment[nd];
            nodeMappingDegree[thr.getTotalNumNodeMappings()] += 1;
            passThroughDegree[thr.getTotalNumEdgeBelowTraversals()] += 1;
            loopDegree[thr.getTotalNumLoops()] += 1;
            totalNumNodes += 1;
            std::vector<NodeWithSplits *> aliasedBy = getNodesAliasedBy(nd, *taxonomy);
            if (thr.reportIfContested(out, nd, treePtrByIndex, aliasedBy, verbose)) {
                totalContested += 1;
                if (nd->getOutDegree() == 1) {
                    redundContested += 1;
                }
            }
        }
        unsigned long m = std::max(loopDegree.rbegin()->first, passThroughDegree.rbegin()->first);
        m = std::max(m, nodeMappingDegree.rbegin()->first);
        out << "Degree\tNodeMaps\tEdgeMaps\tLoops\n";
        for (unsigned long i = 0 ; i <= m; ++i) {
            out << i << '\t' << nodeMappingDegree[i]<< '\t' << passThroughDegree[i] << '\t' << loopDegree[i]<< '\n';
        }
        out << totalNumNodes << " internals\n" << totalContested << " contested\n" << (totalNumNodes - totalContested) << " uncontested\n";
        out << redundContested << " monotypic contested\n";
    }
    
    bool summarize(const OTCLI &otCLI) override {
        if (doConstructSupertree) {
            cloneTaxonomyAsASourceTree();
            constructSupertree();
        }
        std::ostream & out{otCLI.out};
        assert (taxonomy != nullptr);
        if (doReportAllContested) {
            reportAllConflicting(out, otCLI.verbose);
        } else {
            for (auto tr : idsListToReportOn) {
                auto nd = taxonomy->getData().getNodeForOttId(tr);
                if (nd == nullptr) {
                    throw OTCError(std::string("Unrecognized OTT ID in list of OTT IDs to report on: ") + std::to_string(tr));
                }
                const auto & thr = taxoToAlignment[nd];
                std::vector<NodeWithSplits *> aliasedBy = getNodesAliasedBy(nd, *taxonomy);
                thr.reportIfContested(out, nd, treePtrByIndex, aliasedBy, otCLI.verbose);
            }
        }
        for (auto tr : idListForDotExport) {
            auto nd = taxonomy->getData().getNodeForOttId(tr);
            if (nd == nullptr) {
                throw OTCError(std::string("Unrecognized OTT ID in list of OTT IDs to export to DOT: ") + std::to_string(tr));
            }
            //const auto & thr = taxoToAlignment[nd];
            //writeDOTExport(out, thr, nd, treePtrByIndex);
            for (auto n : iter_pre_n_const(nd)) {
                const auto & thr = taxoToAlignment[n];
                writeDOTExport(out, thr, n, treePtrByIndex);
            }
        }
        return true;
    }

    bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        checkTreeInvariants(*taxonomy);
        suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
        checkTreeInvariants(*taxonomy);
        for (auto nd : iter_node(*taxonomy)) {
            taxoToAlignment.emplace(nd, NodeThreadingWithSplits{});
        }
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedWithSplits * raw = treeup.get();
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        treePtrByIndex.push_back(raw);
        // Store the tree's filename
        raw->setName(otCLI.currentFilename);
        threadNewTree(*taxonomy, *raw, treeIndex);
        otCLI.out << "# pathPairings = " << pathPairings.size() << '\n';
        return true;
    }

    bool cloneTaxonomyAsASourceTree() {
        assert(taxonomy != nullptr);
        assert(taxonomyAsSource == nullptr);
        std::unique_ptr<TreeMappedWithSplits> tree = std::move(cloneTree(*taxonomy));
        taxonomyAsSource = tree.get();
        std::size_t treeIndex = inputTreesToIndex.size();
        TreeMappedWithSplits * raw = tree.get();
        inputTreesToIndex[std::move(tree)] = treeIndex;
        treePtrByIndex.push_back(taxonomyAsSource);
        // Store the tree's filename
        raw->setName("TAXONOMY");
        threadTaxonomyClone(*taxonomy, *taxonomyAsSource, treeIndex);
        return true;
    }
};

bool handleReportAllFlag(OTCLI & otCLI, const std::string &);
bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &);
bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg);
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &narg);

bool handleReportAllFlag(OTCLI & otCLI, const std::string &) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doReportAllContested = true;
    return true;
}
bool handleSuperTreeFlag(OTCLI & otCLI, const std::string &) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->doConstructSupertree = true;
    return true;
}

bool handleReportOnNodesFlag(OTCLI & otCLI, const std::string &narg) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -b argument.");
    }
    auto rs = split_string(narg, ',');
    for (auto word : rs) {
        auto ottId = ottIDFromName(word);
        if (ottId < 0) {
            throw OTCError(std::string("Expecting a list of IDs after the -b argument. Offending word: ") + word);
        }
        proc->idsListToReportOn.push_back(ottId);
    }
    return true;
}

bool handleDotNodesFlag(OTCLI & otCLI, const std::string &narg) {
    RemapToDeepestUnlistedState * proc = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -d argument.");
    }
    auto rs = split_string(narg, ',');
    for (auto word : rs) {
        auto ottId = ottIDFromName(word);
        if (ottId < 0) {
            throw OTCError(std::string("Expecting a list of IDs after the -d argument. Offending word: ") + word);
        }
        proc->idListForDotExport.push_back(ottId);
    }
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otcscaffoldedsupertree",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Crashes or emits bogus output.",
                "taxonomy.tre inp1.tre inp2.tre");
    RemapToDeepestUnlistedState proc;
    otCLI.addFlag('a',
                  "Write a report of all contested nodes",
                  handleReportAllFlag,
                  false);
    otCLI.addFlag('s',
                  "Compute a supertree",
                  handleSuperTreeFlag,
                  false);
    otCLI.addFlag('b',
                  "IDLIST should be a list of OTT IDs. A status report will be generated for those nodes",
                  handleReportOnNodesFlag,
                  true);
    otCLI.addFlag('d',
                  "IDLIST should be a list of OTT IDs. A DOT file of the nodes will be generated ",
                  handleDotNodesFlag,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
