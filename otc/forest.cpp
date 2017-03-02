#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/debug.h"

#include "otc/write_dot.h"
namespace otc {
template<typename T, typename U>
typename RootedForest<T, U>::tree_type &
RootedForest<T, U>::createNewTree() {
    std::size_t i = nextTreeId++;
    auto r = trees.emplace(std::piecewise_construct,
                           std::forward_as_tuple(i),
                           std::forward_as_tuple(nextTreeId,
                                           *this,
                                           ottIdToNodeMap));
    assert(r.second); // must be a new Tree!
    return trees.at(i);
}
   
template<typename T, typename U>
RootedForest<T, U>::RootedForest(long rootOttId)
    :nextTreeId(0U),
    ottIdToNodeMap(nodeSrc.get_data().ottIdToNode),
    rootID(rootOttId) {
}

template<typename T, typename U>
void RootedForest<T, U>::registerLeaf(long ottId) {
    if (ottId == rootID) {
        return;
    }
    auto f = ottIdToNodeMap.find(ottId);
    if (f != ottIdToNodeMap.end()) {
        return;
    }
    createLeaf(nullptr, ottId, nullptr);
}

// TMP could be faster by storing node->tree lookup
template<typename T, typename U>
bool RootedForest<T, U>::isAttached(long ottId) const {
    auto f = ottIdToNodeMap.find(ottId);
    if (f == ottIdToNodeMap.end()) {
        return false;
    }
    node_type * n = f->second;
    assert(n != nullptr);
    return (n->getParent() != nullptr);
}

template<typename T, typename U>
bool RootedForest<T, U>::nodeIsAttached(RootedTreeNode<T> & n) const {
    return (n.getParent() != nullptr);
}

template<typename T, typename U>
void RootedForest<T, U>::attachAllKnownTipsAsNewTree() {
    tree_type & t = createNewTree();
    t.root = createNode(nullptr, &t);
    for (auto & o2n : ottIdToNodeMap) {
        if (o2n.first == rootID) {
            assert(false);
            throw OTCError("root as tip");
        }
        auto nd = o2n.second;
        if (!nodeIsAttached(*nd)) {
            //t.connectedIds.insert(o2n.first);
            addAndUpdateChild(t.root, nd, t);
        }
    }
}

template<typename T, typename U>
void RootedForest<T, U>::_AndUpdateChild(RootedTreeNode<T> *p,
                                                    RootedTreeNode<T> *c,
                                                    FTree<T, U> &tree) {
    registerTreeForNode(c, &tree);
    const auto & cd = c->get_data().desIds;
    if (cd.size() == 1 && !c->hasOttId()) {
        if (*begin(cd) == c->getOttId()) {
            return;
        }
    }
    p->get_data().desIds.insert(begin(cd), end(cd));
}

template<typename T, typename U>
void RootedForest<T, U>::attachAllDetachedTips() {
    if (trees.empty()) {
        attachAllKnownTipsAsNewTree();
        return;
    }
    assert(trees.size() == 1);
    tree_type & t = trees.begin()->second;
    std::list<node_type *> excludedFromRoot;
    std::list<node_type *> attachableAtRoot;
    for (auto & o2n : ottIdToNodeMap) {
        if (o2n.first == rootID) {
            assert(false);
            throw OTCError("root as tip");
        }
        auto nd = o2n.second;
        if (!nodeIsAttached(*nd)) {
            if (t.isExcludedFromRoot(nd)) {
                excludedFromRoot.push_back(nd);
            } else {
                attachableAtRoot.push_back(nd);
            }
        }
    }
    if (!excludedFromRoot.empty()) {
        auto nr = createNode(nullptr, &t);
        addAndUpdateChild(nr, t.root, t);
        t.root = nr;
        for (auto n : excludedFromRoot) {
            addAndUpdateChild(nr, n, t);
        }
    }
    // these could be attached one node more tipward (old root), but there is no basis for that.
    for (auto n : attachableAtRoot) {
        addAndUpdateChild(t.root, n, t);
    }
    
}

template<typename T, typename U>
bool RootedForest<T, U>::addPhyloStatement(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        dbWriteOttSet(" RootedForest::addPhyloStatement\nincGroup ", ps.includeGroup);
        dbWriteOttSet(" leafSet", ps.leafSet);
    }
    ps.debugCheck();
    assert(ps.includeGroup.size() > 1);
    for (auto oid : ps.leafSet) {
        registerLeaf(oid);
    }
    if (ps.includeGroup == ps.leafSet) {
        LOG(DEBUG) << "trivial group - exiting";
        return true;
    }
    const auto incompatRedundant = checkWithPreviouslyAddedStatement(ps);
    if (incompatRedundant.first) {
        LOG(DEBUG) << "    hit incompat w/ prev added shortcircuit";
        return false;
    }
    if (incompatRedundant.second) { // we have added an identical group before
        LOG(DEBUG) << "    hit redundant w/ prev added shortcircuit";
        return true;
    }
    LOG(DEBUG) << "    checking compat w/ graph";
    if (addPhyloStatementToGraph(ps)) {
        LOG(DEBUG) << "    compat w/ graph";
        addedSplitsByLeafSet[ps.leafSet].insert(ps.includeGroup);
        return true;
    }
    LOG(DEBUG) << "    incompat w/ graph";
    return false;
}

template<typename T, typename U>
void RootedForest<T, U>::dumpAcceptedPhyloStatements(const char *fn) {
    std::ofstream phyloStatementOut;
    phyloStatementOut.open(fn);
    for (const auto & ps : novelAcceptedPSInOrder) {
        ps.writeAsNewick(phyloStatementOut);
        phyloStatementOut << '\n';
    }
    phyloStatementOut.close();
}

void appendIncludeLeafSetAsNewick(const char *fn, const OttIdSet & inc, const OttIdSet & ls) {
    PhyloStatementSource pss(-1, -1);
    const OttIdSet exc = set_difference_as_set(ls, inc);
    std::ofstream phyloStatementOut;
    PhyloStatement ps{inc, exc, ls, pss};
    phyloStatementOut.open(fn, std::ios::app);
    ps.writeAsNewick(phyloStatementOut);
    phyloStatementOut << '\n';
    phyloStatementOut.close();
}

template<typename T, typename U>
std::pair<bool, bool> RootedForest<T, U>::checkWithPreviouslyAddedStatement(const PhyloStatement &ps) const {
    if (false && debuggingOutputEnabled) {
        dbWriteOttSet(" RootedForest::conflictsWithPreviouslyAddedStatement incGroup ", ps.includeGroup);
        dbWriteOttSet(" leafSet", ps.leafSet);
    }
    for (const auto sIt : addedSplitsByLeafSet) {
        const auto & prevAddedLeafSet = sIt.first;
        const auto relLeafSet = set_intersection_as_set(prevAddedLeafSet, ps.leafSet);
        const bool exactLS = relLeafSet.size() == ps.leafSet.size();
        if (relLeafSet.size() < 3) { // no conflict is possible if the intersection is so small that no phylostatements are made
            continue;
        }
        const auto relIncGroup = set_intersection_as_set(ps.includeGroup, relLeafSet);
        const auto & setPrevInc = sIt.second;
        for (const auto & prevInc : setPrevInc) {
            if (exactLS && prevInc == ps.includeGroup) {
                return std::pair<bool, bool>(false, true);
            }
            if (culledAndCompleteIncompatWRTLeafSet(relIncGroup, prevInc, relLeafSet)) {
                return std::pair<bool, bool>(true, false);
            }
        }
    }
    return std::pair<bool, bool>(false, false);
}

template<typename T, typename U>
void consumeMapToList(std::map<T, std::list<U> > &m, std::list<U> & out) {
    for (auto & mIt : m) {
        for (auto el = begin(mIt.second) ; el != end(mIt.second); ++el) {
            out.push_back(*el);
        }
    }
}

template<typename T, typename U>
std::list<OverlapFTreePair<T, U> > RootedForest<T, U>::getSortedOverlappingTrees(const OttIdSet &inc) {
    typedef OverlapFTreePair<T, U> MyOverlapFTreePair;
    std::map<std::size_t, std::list<MyOverlapFTreePair> > byOverlapSize;
    for (auto & tpIt : trees) {
        tree_type * ftree = &(tpIt.second);
        const OttIdSet & inTree = ftree->getIncludedOttIds();
        const OttIdSet inter = set_intersection_as_set(inTree, inc);
        if (!inter.empty()) {
            const auto k = inter.size();
            auto & tsList = byOverlapSize[k];
            tsList.push_back(MyOverlapFTreePair(inter, ftree));
        }
    }
    std::list<MyOverlapFTreePair> r;
    consumeMapToList(byOverlapSize, r);
    return r;
}

template<typename T, typename U>
void RootedForest<T, U>::addIngroupDisjointPhyloStatementToGraph(const PhyloStatement &ps) {
    // this ingroup does not overlap with any ftree. find the FTree with the most overlap
    //  with the excludeGroup...
    auto byExcCardinality = getSortedOverlappingTrees(ps.excludeGroup);
    if (byExcCardinality.empty()) {
        // none of the ingroup or outgroup are attached.
        // create a new FTree...
        // this can happen if the outgroup are mentioned in exclude statements (so the 
        //  areDisjoint returns false). But sense will add all of the leaves in the 
        //  includeGroup and excludeGroup to this new tree, we don't need any new constraints
        //  so we can exit
        LOG(DEBUG) << "No exclude overlap either, using addDisjointTree";
        addDisjointTree(ps);
    } else {
        // TMP TOO GREEDY A CONNECTION - should only do this if all of the outgroup is connected...
        // we'll add the ingroup as a child of the root
        auto ftreeToAttach = byExcCardinality.begin()->second;
        LOG(DEBUG) << "Using addIncludeGroupDisjointPhyloStatement";
        ftreeToAttach->addIncludeGroupDisjointPhyloStatement(ps);
    }
    // no other trees had an includeGroup, so no need to add constraints....
}

template<typename T, typename U>
bool RootedForest<T, U>::isInABand(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (tree.isInABand(ncn)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool RootedForest<T, U>::hasNodesExcludedFromIt(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (tree.hasNodesExcludedFromIt(ncn)) {
            return true;
        }
    }
    return false;
}

// should not modify the forest if returning false
template<typename T, typename U>
bool RootedForest<T, U>::checkCanAddIngroupOverlappingPhyloStatementToGraph(
            const std::list<OverlapFTreePair<T, U> > & byIncCardinality,
            const PhyloStatement &ps,
            std::list<node_type * > & nonTrivMRCAs,
            OttIdSet & attachedElsewhere,
            std::vector<bool> & shouldResolveVec,
            std::vector<bool> & shouldCreateDeeperVec) const {
    for (const auto & incPair : byIncCardinality) {
        const auto & incGroupIntersection = incPair.first;
        attachedElsewhere.insert(incGroupIntersection.begin(), incGroupIntersection.end());
        tree_type * f = incPair.second;
        node_type * includeGroupA = nullptr;
        includeGroupA = f->getMRCA(incGroupIntersection);
        assert(includeGroupA != nullptr);
        assert(getTreeForNode(includeGroupA) == f);
        if (includeGroupA->isTip()) {
            // this can happen if the overlap is one taxon.
            includeGroupA = includeGroupA->getParent();
            assert(includeGroupA != nullptr);
            assert(getTreeForNode(includeGroupA) == f);
        }
        // If any of the ingroup are specifically excluded, then we have move deeper in the tree.
        // TMP this could be more efficient and avoid the while loop.
        while (f->anyExcludedAtNode(includeGroupA, ps.includeGroup)) {
            if (f->anyIncludedAtNode(includeGroupA, ps.excludeGroup)) {
                return false;
            }
            if (f->anyPhantomNodesAtNode(includeGroupA, ps.includeGroup)) {
                return false;
            }
            includeGroupA = includeGroupA->getParent();
            if (includeGroupA == nullptr) {
                break;
            }
            assert(getTreeForNode(includeGroupA) == f);
        }
        OttIdSet excInc;
        bool forceDeeperRoot = false;
        if (includeGroupA == nullptr) {
            includeGroupA = f->getRoot();
            forceDeeperRoot = true;
            assert(getTreeForNode(includeGroupA) == f);
        } else {
            excInc = set_intersection_as_set(includeGroupA->get_data().desIds, ps.excludeGroup);
            if (debuggingOutputEnabled) {
                LOG(DEBUG) << "     addPhyloStatementToGraph search for an ancestor of ..."; 
                dbWriteOttSet(" addPhyloStatementToGraph search for an ancestor of:  ", incGroupIntersection);
                dbWriteOttSet(" wanted to avoid =  ", ps.excludeGroup);
                dbWriteOttSet(" found a node with desIds:  ", includeGroupA->get_data().desIds);
                dbWriteOttSet(" which includes the excludegroup members:  ", excInc);
            }
            if (!canBeResolvedToDisplayIncExcGroup(includeGroupA, ps.includeGroup, excInc)) {
                return false; // the MRCA of the includeGroup had interdigitated members of the excludeGroup
            }
        }
        shouldCreateDeeperVec.push_back(forceDeeperRoot);
        shouldResolveVec.push_back(!excInc.empty());
        nonTrivMRCAs.push_back(includeGroupA);
    }
    return true;
}

 
template<typename T, typename U>
InterTreeBand<T> * RootedForest<T, U>::_createNewBand(FTree<T, U> & ,
                                                 RootedTreeNode<T> &nd,
                                                 const PhyloStatement &ps) {
    std::set<RootedTreeNode<T> *> emptySet;
    allBands.emplace_back(&nd, emptySet, ps);
    return &(*allBands.rbegin());
}

template<typename T, typename U>
bool RootedForest<T, U>::addIngroupOverlappingPhyloStatementToGraph(const std::list<OverlapFTreePair<T, U> > & byIncCardinality,
                                                                   const PhyloStatement &ps) {
    std::list<node_type * > nonTrivMRCAs;
    OttIdSet attachedElsewhere;
    std::vector<bool> shouldResolveVec;
    std::vector<bool> shouldCreateDeeperVec;
    if (!checkCanAddIngroupOverlappingPhyloStatementToGraph(byIncCardinality, ps, nonTrivMRCAs, attachedElsewhere, shouldResolveVec, shouldCreateDeeperVec)) {
        return false;
    }
    novelAcceptedPSInOrder.push_back(ps); //TMP DEBUGGING
    // all non trivial overlapping trees have approved this split...
    auto ntmIt = begin(nonTrivMRCAs);
    auto srIt = begin(shouldResolveVec);
    auto scdIt = begin(shouldCreateDeeperVec);
    unsigned i = 0;
    InterTreeBand<T> * itbp = nullptr;
    for (const auto & incPair : byIncCardinality) {
        LOG(DEBUG) << "   addIngroupOverlappingPhyloStatementToGraph mod for loop round " << ++i;
        debugInvariantsCheck();
        tree_type * f = incPair.second;
        assert(ntmIt != nonTrivMRCAs.end());
        node_type * includeGroupA = *ntmIt++;
        const bool addNode = *srIt++;
        const bool shouldCreateDeeperRoot = *scdIt;
        if (addNode) {
            includeGroupA = f->resolveToCreateCladeOfIncluded(includeGroupA, ps);
            assert(getTreeForNode(includeGroupA) == f);
            LOG(DEBUG) << "   back from resolveToCreateCladeOfIncluded for loop round " << i;
            debugInvariantsCheck();
        } else if (shouldCreateDeeperRoot) {
            f->createDeeperRoot();
            includeGroupA = f->getRoot();
            assert(getTreeForNode(includeGroupA) == f);
            LOG(DEBUG) << "   back from createDeeperRoot for loop round " << i;
            debugInvariantsCheck();
        } else {
            assert(getTreeForNode(includeGroupA) == f);
        }
        if (byIncCardinality.size() > 1 && itbp == nullptr) {
            itbp = _createNewBand(*f, *includeGroupA, ps);
        }
        auto connectedHere = f->addPhyloStatementAtNode(ps, includeGroupA, attachedElsewhere, itbp);
        if (!connectedHere.empty()) {
            attachedElsewhere.insert(begin(connectedHere), end(connectedHere));
        }
        dbWriteOttSet(" includeGroupA...desIds ", includeGroupA->get_data().desIds);
        LOG(DEBUG) << "   back from addPhyloStatementAtNode for loop round " << i;
        debugInvariantsCheck();
    }
    LOG(DEBUG) << "   addIngroupOverlappingPhyloStatementToGraph exit true " << i;
    return true;
}

template<typename T, typename U>
bool RootedForest<T, U>::addPhyloStatementToGraph(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        dbWriteOttSet(" RootedForest::addPhyloStatementToGraph incGroup ", ps.includeGroup);
        dbWriteOttSet(" leafSet", ps.leafSet);
    }
    if (ps.isTrivial()) {
        novelAcceptedPSInOrder.push_back(ps); //TMP DEBUGGING
        auto newOttIds = set_difference_as_set(ps.includeGroup, ottIdSet);
        for (auto noid : newOttIds) {
            addDetachedLeaf(noid);
        }
        debugInvariantsCheck();
        return true;
    }
    if (areDisjoint(ps.leafSet, ottIdSet)) {
        novelAcceptedPSInOrder.push_back(ps); //TMP DEBUGGING
        addDisjointTree(ps);
        debugInvariantsCheck();
        return true;
    }
    auto byIncCardinality = getSortedOverlappingTrees(ps.includeGroup);
    LOG(DEBUG) << byIncCardinality.size() << " FTree instance referred to in byIncCardinality";
    if (byIncCardinality.empty()) {
        novelAcceptedPSInOrder.push_back(ps); //TMP DEBUGGING
        for (auto o : ps.includeGroup) {
            assert(!isAttached(o));
        }
        LOG(DEBUG) << "No intersection between includeGroup of an existing FTree.";
        addIngroupDisjointPhyloStatementToGraph(ps);
        debugInvariantsCheck();
        return true;
    }
    auto & attachmentPair = *byIncCardinality.begin();
    if (attachmentPair.first.size() == 1) {
        // greedy approach is to add the rest of the ingroup as deep in this tree as possible.
        //  less greedy: make include/exclude statements at that node
        LOG(DEBUG) << "Missing opportunity to special case ingroups that only overlap with leaves of current trees.\n";
    }
    debugInvariantsCheck();
    auto rc = addIngroupOverlappingPhyloStatementToGraph(byIncCardinality, ps);
    debugInvariantsCheck();
    return rc;
}

template<typename T, typename U>
std::map<const FTree<T, U> *, const RootedTreeNode<T> *>
RootedForest<T, U>::getTreeToNodeMapForBand(const InterTreeBand<T> & itb) const {
    std::set<const node_type *> bns = itb.getBandedNodes();
    std::map<const FTree<T, U> *, const RootedTreeNode<T> *> r;
    for (auto bnp : bns) {
        auto t = getTreeForNode(bnp);
        assert(t != nullptr);
        assert(!contains(r, t));
        r[t] = bnp;
    }
    return r;
}


template<typename T, typename U>
FTree<T, U> & RootedForest<T, U>::addDisjointTree(const PhyloStatement &ps) {
    tree_type & r = createNewTree();
    r.mirrorPhyloStatement(ps);
    return r;
}

template<typename T, typename U>
void RootedForest<T, U>::writeForestDOTToFN(const std::string &fn) const {
    LOG(DEBUG) << "     creating DOT file for forest: " << fn;
    std::ofstream outf(fn);
    writeDOTForest(outf, *this);
}

#if defined(DO_DEBUG_CHECKS)

template<typename T, typename U>
void RootedForest<T, U>::debugInvariantsCheck() const {
    std::map<const node_type *, const tree_type *> root2tree;
    for (const auto & t : trees) {
        t.second.debugInvariantsCheckFT();
        auto r = t.second.getRoot();
        assert(!contains(root2tree, r));
        root2tree[r] = &(t.second);
    }
    std::map<long, const tree_type *> ottId2Tree;
    std::map<const node_type *, const tree_type *> internal2Tree;
    std::set<const node_type *> detached;
    for (auto o2n : ottIdToNodeMap) {
        auto o = o2n.first;
        assert(o != rootID);
        auto n = o2n.second;
        if (n->isTip()) {
            assert(n->getOttId() == o);
            auto d = getDeepestAnc(n);
            assert(d != nullptr);
            if (d == n) {
                assert(!contains(detached, n));
                detached.insert(n);
            } else {
                ottId2Tree[o] = root2tree.at(d);
            }
        } else {
            assert(!n->hasOttId());
            auto d = getDeepestAnc(n);
            assert(d != nullptr);
            assert(!contains(internal2Tree, n));
            if (d == n) {
                assert(contains(root2tree, n));
            }
            internal2Tree[n] = root2tree.at(d);
        }
        if (n->getParent() != nullptr) {
            assert(nullptr != getTreeForNode(n));
        }
    }
    std::set<long> connectedIdSet;
    for (const auto & t : trees) {
        for (const auto o : t.second.getConnectedOttIds()) {
            assert(!contains(connectedIdSet, o));
            //LOG(DEBUG) << "checkint ott" << o;
            assert(ottId2Tree.at(o) == &(t.second));
            connectedIdSet.insert(o);
        }
    }
    assert(connectedIdSet.size() == ottId2Tree.size()); 
}
#endif
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
