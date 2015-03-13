#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"

#include "otc/write_dot.h"
namespace otc {

bool PhyloStatement::debugCheck() const {
#ifdef DEBUGGING_PHYLO_STATEMENTS
    const OttIdSet ie = set_union_as_set(includeGroup, excludeGroup);
    if (ie != leafSet) {
        LOG(DEBUG)  << " includeGroup "; dbWriteOttSet(" ", includeGroup, " ");
        LOG(DEBUG)  << " excludeGroup "; dbWriteOttSet(" ", excludeGroup, " ");
        LOG(DEBUG)  << " leafSet "; dbWriteOttSet(" ", leafSet, " ");
        assert(false);
    }
#endif
    return true;
}

template<typename T, typename U>
void FTree<T, U>::createDeeperRoot() {
    auto nr = forest.createNode(nullptr);
    nr->addChild(root);
    root = nr;
}

template<typename T, typename U>
void FTree<T, U>::addSubtree(RootedTreeNode<T> * subtreeRoot,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedInc,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedExc,
                const std::map<node_type *, std::list<std::pair<node_type*, PhyloStatementSource> > > & otherIC,
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
        root->addChild(subtreeRoot);
        for (auto c : iter_child_const(*subtreeRoot)) {
            connectedIds.insert(c->getOttId());
        }
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
                for (auto eValNd : esn->second) {
                    includesConstraints[eValNd].push_back(GroupingConstraint{sn, bogusPSS});
                }
            }
        }
        return;
    }
    assert(false);
}

template<typename T, typename U>
typename RootedForest<T,U>::tree_type &
RootedForest<T,U>::createNewTree() {
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
RootedForest<T,U>::RootedForest(long rootOttId)
    :nextTreeId(0U),
    ottIdToNodeMap(nodeSrc.getData().ottIdToNode),
    rootID(rootOttId) {
}
template<typename T, typename U>
void RootedForest<T,U>::registerLeaf(long ottId) {
    if (ottId == rootID) {
        return;
    }
    auto f = ottIdToNodeMap.find(ottId);
    if (f != ottIdToNodeMap.end()) {
        return;
    }
    auto n = createNode(nullptr);
    n->setOttId(ottId);
    ottIdToNodeMap[ottId] = n;
}

// TMP could be faster by storing node->tree lookup
template<typename T, typename U>
bool RootedForest<T,U>::isAttached(long ottId) const {
    auto f = ottIdToNodeMap.find(ottId);
    if (f == ottIdToNodeMap.end()) {
        return false;
    }
    node_type * n = f->second;
    assert(n != nullptr);
    return (n->getParent() != nullptr);
}

template<typename T, typename U>
bool RootedForest<T,U>::nodeIsAttached(RootedTreeNode<T> & n) const {
    return (n.getParent() != nullptr);
}

template<typename T, typename U>
void RootedForest<T,U>::attachAllKnownTipsAsNewTree() {
    tree_type & t = createNewTree();
    t.root = createNode(nullptr);
    for (auto & o2n : ottIdToNodeMap) {
        if (o2n.first == rootID) {
            assert(false);
            continue;
        }
        auto nd = o2n.second;
        if (!nodeIsAttached(*nd)) {
            t.connectedIds.insert(o2n.first);
            t.root->addChild(nd);
        }
    }
}
template<typename T, typename U>
bool FTree<T,U>::isExcludedFromRoot(const node_type *n) const {
    node_type * ncn = const_cast<node_type *>(n);
    auto nit = excludesConstraints.find(ncn);
    if (nit == excludesConstraints.end()) {
        return false;
    }
    for (const auto & en : nit->second) {
        if (en.first == root) {
            return true;
        }
    }
    return false;
}
template<typename T, typename U>
void RootedForest<T,U>::attachAllDetachedTips() {
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
            continue;
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
        auto nr = createNode(nullptr);
        nr->addChild(t.root);
        t.root = nr;
        for (auto n : excludedFromRoot) {
            t.connectedIds.insert(n->getOttId());
            nr->addChild(n);
        }
    }
    for (auto n : attachableAtRoot) {
        t.connectedIds.insert(n->getOttId());
        t.root->addChild(n);
    }
    
}

template<typename T, typename U>
bool RootedForest<T,U>::addPhyloStatement(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        LOG(DEBUG) << " RootedForest::addPhyloStatement";
        LOG(DEBUG) << " incGroup "; dbWriteOttSet(" ", ps.includeGroup, " ");
        LOG(DEBUG) << " leafSet "; dbWriteOttSet(" ", ps.leafSet, " ");
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
        LOG(DEBUG) << "    hit redundandt w/ prev added shortcircuit";
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
std::pair<bool, bool> RootedForest<T,U>::checkWithPreviouslyAddedStatement(const PhyloStatement &ps) const {
    if (false && debuggingOutputEnabled) {
        LOG(DEBUG) << " RootedForest::conflictsWithPreviouslyAddedStatement";
        LOG(DEBUG) << " incGroup "; dbWriteOttSet(" ", ps.includeGroup, " ");
        LOG(DEBUG) << " leafSet "; dbWriteOttSet(" ", ps.leafSet, " ");
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
std::list<OverlapFTreePair<T, U> > RootedForest<T,U>::getSortedOverlappingTrees(const OttIdSet &inc) {
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
void RootedForest<T,U>::addIngroupDisjointPhyloStatementToGraph(const PhyloStatement &ps) {
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
        addDisjointTree(ps);
    } else {
        // TMP TOO GREEDY A CONNECTION - should only do this if all of the outgroup is connected...
        // we'll add the ingroup as a child of the root
        auto ftreeToAttach = byExcCardinality.begin()->second;
        ftreeToAttach->addIncludeGroupDisjointPhyloStatement(ps);
    }
    // no other trees had an includeGroup, so no need to add constraints....
}

template<typename T, typename U>
bool RootedForest<T,U>::isMentionedInInclude(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (contains(tree.includesConstraints, ncn)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool RootedForest<T,U>::isMentionedInExclude(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (contains(tree.excludesConstraints, ncn)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool RootedForest<T,U>::addIngroupOverlappingPhyloStatementToGraph(const std::list<OverlapFTreePair<T, U> > & byIncCardinality, const PhyloStatement &ps) {
    std::list<node_type * > nonTrivMRCAs;
    OttIdSet attachedElsewhere;
    std::vector<bool> shouldResolveVec;
    std::vector<bool> shouldCreateDeeperVec;
    for (const auto & incPair : byIncCardinality) {
        const auto & incGroupIntersection = incPair.first;
        attachedElsewhere.insert(incGroupIntersection.begin(), incGroupIntersection.end());
        tree_type * f = incPair.second;
        node_type * includeGroupA = nullptr;
        includeGroupA = f->getMRCA(incGroupIntersection);
        assert(includeGroupA != nullptr);
        if (includeGroupA->isTip()) {
            // this can happen if the overlap is one taxon.
            includeGroupA = includeGroupA->getParent();
            assert(includeGroupA != nullptr);
        }
        // If any of the ingroup are specifically excluded, then we have move deeper in the tree.
        // TMP this could be more efficient and avoid the while loop.
        while (f->anyExcludedAtNode(includeGroupA, ps.includeGroup)) {
            if (f->anyIncludedAtNode(includeGroupA, ps.excludeGroup)) {
                return false;
            }
            if (f->anyForceIncludedAtNode(includeGroupA, ps.includeGroup)) {
                return false;
            }
            includeGroupA = includeGroupA->getParent();
            if (includeGroupA == nullptr) {
                break;
            }
        }
        OttIdSet excInc;
        bool forceDeeperRoot = false;
        if (includeGroupA == nullptr) {
            includeGroupA = f->getRoot();
        } else {
            excInc = set_intersection_as_set(includeGroupA->getData().desIds, ps.excludeGroup);
            if (debuggingOutputEnabled) {
                LOG(DEBUG) << "     addPhyloStatementToGraph search for an ancestor of ..."; 
                LOG(DEBUG) << " addPhyloStatementToGraph search for an ancestor of:  "; dbWriteOttSet(" ", incGroupIntersection, " ");
                LOG(DEBUG) << "  wanted to avoid =  "; dbWriteOttSet(" ", ps.excludeGroup, " ");
                LOG(DEBUG) << "  found a node with desIds:  "; dbWriteOttSet(" ", includeGroupA->getData().desIds, " ");
                LOG(DEBUG) << "  which includes the excludegroup members:  "; dbWriteOttSet(" ", excInc, " ");
            }
            if (!canBeResolvedToDisplayExcGroup(includeGroupA, ps.includeGroup, excInc)) {
                return false; // the MRCA of the includeGroup had interdigitated members of the excludeGroup
            }
        }
        shouldCreateDeeperVec.push_back(false);
        shouldResolveVec.push_back(!excInc.empty());
        nonTrivMRCAs.push_back(includeGroupA);
    }
    // all non trivial overlapping trees have approved this split...
    auto ntmIt = begin(nonTrivMRCAs);
    auto srIt = begin(shouldResolveVec);
    auto scdIt = begin(shouldCreateDeeperVec);
    for (const auto & incPair : byIncCardinality) {
        tree_type * f = incPair.second;
        assert(ntmIt != nonTrivMRCAs.end());
        node_type * includeGroupA = *ntmIt++;
        const bool addNode = *srIt++;
        const bool shouldCreateDeepeRoot = *scdIt;
        if (addNode) {
            includeGroupA = f->resolveToCreateCladeOfIncluded(includeGroupA, ps.includeGroup);
        } else if (shouldCreateDeepeRoot) {
            f->createDeeperRoot();
            includeGroupA = f->getRoot();
        }
        auto connectedHere = f->addPhyloStatementAtNode(ps, includeGroupA, attachedElsewhere);
        if (!connectedHere.empty()) {
            attachedElsewhere.insert(begin(connectedHere), end(connectedHere));
        }
    }
    return true;
}

template<typename T, typename U>
bool RootedForest<T,U>::addPhyloStatementToGraph(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        LOG(DEBUG) << " RootedForest::addPhyloStatementToGraph";
        LOG(DEBUG) << " incGroup "; dbWriteOttSet(" ", ps.includeGroup, " ");
        LOG(DEBUG) << " leafSet "; dbWriteOttSet(" ", ps.leafSet, " ");
    }
    if (ps.isTrivial()) {
        auto newOttIds = set_difference_as_set(ps.includeGroup, ottIdSet);
        for (auto noid : newOttIds) {
            addDetachedLeaf(noid);
        }
        return true;
    }
    if (areDisjoint(ps.leafSet, ottIdSet)) {
        addDisjointTree(ps);
        return true;
    }
    auto byIncCardinality = getSortedOverlappingTrees(ps.includeGroup);
    LOG(DEBUG) << byIncCardinality.size() << " FTree instance referred to in byIncCardinality";
    if (byIncCardinality.empty()) {
        LOG(DEBUG) << "No intersection between includeGroup of an existing FTree.";
        addIngroupDisjointPhyloStatementToGraph(ps);
        return true;
    }
    auto & attachmentPair = *byIncCardinality.begin();
    if (attachmentPair.first.size() == 1) {
        // greedy approach is to add the rest of the ingroup as deep in this tree as possible.
        //  less greedy: make include/exclude statements at that node
        LOG(DEBUG) << "Missing opportunity to special case ingroups that only overlap with leaves of current trees.\n";
    }
    return addIngroupOverlappingPhyloStatementToGraph(byIncCardinality, ps); 
}

template<typename T, typename U>
bool FTree<T,U>::anyExcludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    const node_type * p = nd->getParent();
    const OttIdSet & ndi = nd->getData().desIds;
    for (auto oid : ottIdSet) {
        auto leafNd = ottIdToNodeMap.at(oid);
        auto gcIt = excludesConstraints.find(leafNd);
        if (gcIt != excludesConstraints.end()) {
            for (const auto & gc : gcIt->second) {
                node_type * en = gc.first;
                if (en == nd || en == p || isSubset(ndi, en->getData().desIds)) {
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename T, typename U>
bool FTree<T,U>::anyIncludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    if (!areDisjoint(nd->getData().desIds, ottIdSet)) {
        return true;
    }
    
    auto c = anyForceIncludedAtNode(nd, ottIdSet);
    assert(c == false);// if we are correctly updating desIds we don't need this branch.... TMP
    return c;
}

template<typename T, typename U>
bool FTree<T,U>::anyForceIncludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    for (auto oid :ottIdSet) {
        auto oidN = ottIdToNodeMap.at(oid);
        auto iclIt = includesConstraints.find(oidN);
        if (iclIt != includesConstraints.end()) {
            for (auto inNd : iclIt->second) {
                if (inNd.first == oidN) {
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId) {
    connectedIds.insert(ottId);
    return forest.createLeaf(par, ottId);
}

template<typename T, typename U>
void FTree<T,U>::addExcludeStatement(long ottId, RootedTreeNode<T> * excludedFrom, const PhyloStatementSource &pss) {
    OttIdSet x;
    x.insert(ottId);
    if (anyExcludedAtNode(excludedFrom, x)) {
        return; // already excluded from this subtree
    }
    RootedTreeNode<T> * eNode = ottIdToNodeMap.at(ottId);
    // If any of the descendants exclude this node, we can remove those exclude statements,
    //  because they'll be "dominated by this one"
    auto ecIt = excludesConstraints.find(eNode);
    if (ecIt != excludesConstraints.end()) {
        auto & listOfExc = ecIt->second;
        auto efIt = begin(listOfExc);
        for (; efIt != end(listOfExc);) {
            auto aen = efIt->first;
            if (isAncestorDesNoIter(excludedFrom, aen)) {
                efIt = listOfExc.erase(efIt);
            } else {
                ++efIt;
            }
        }
    }
    excludesConstraints[eNode].push_back(GroupingConstraint(excludedFrom, pss));
}

template<typename T, typename U>
void FTree<T,U>::addIncludeStatement(long ottId, RootedTreeNode<T> * includedIn, const PhyloStatementSource &pss) {
    assert(includedIn != nullptr);
    if (contains(includedIn->getData().desIds, ottId)) {
        return; // must be included in a des
    }
    RootedTreeNode<T> * eNode = ottIdToNodeMap.at(ottId);
    // If any of the ancestors include this node, we can remove those include statements,
    //  because they'll be "dominated by this one"
    auto icIt = includesConstraints.find(eNode);
    if (icIt != includesConstraints.end()) {
        auto & listOfInc = icIt->second;
        auto ifIt = begin(listOfInc);
        for (; ifIt != end(listOfInc);) {
            auto aen = ifIt->first;
            if (isAncestorDesNoIter(includedIn, aen)) {
                ifIt = listOfInc.erase(ifIt);
            } else {
                ++ifIt;
            }
        }
    }
    includesConstraints[eNode].push_back(GroupingConstraint(includedIn, pss));
    // Since we know that the node will be a descendant of includedIn we add its Id to desIds
    includedIn->getData().desIds.insert(ottId);
    for (auto anc : iter_anc(*includedIn)) {
        anc->getData().desIds.insert(ottId);
    }
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par, const OttIdSet & oids) {
    std::set<RootedTreeNode<T> *> cToMove;
    std::list<RootedTreeNode<T> *> orderedToMove;
    std::list<GroupingConstraint *> incToUpdate;
    for (auto oid : oids) {
        auto n = ottIdToNodeMap.at(oid);
        bool connectionFound = false;
        if (n->getParent() == par) {
            cToMove.insert(n);
            orderedToMove.push_back(n);
            connectionFound = true;
        } else {
            for (auto anc : iter_anc(*n)) {
                if (anc->getParent() == par) {
                    if (!contains(cToMove, anc)) {
                        cToMove.insert(anc);
                        orderedToMove.push_back(anc);
                        connectionFound = true;
                        break;
                    }
                }
            }
        }
        if (connectionFound) {
            continue;
        }
        auto icIt = includesConstraints.find(n);
        if (icIt != includesConstraints.end()) {
            auto & listOfConstr = icIt->second;
            auto igcIt = begin(listOfConstr);
            for (; igcIt != end(listOfConstr);) {
                auto np = igcIt->first;
                if (np == par) {
                    incToUpdate.push_back(&(*igcIt));
                    connectionFound = true;
                    break;
                }
                ++igcIt;
            }
        }
    }
    assert(cToMove.size() > 0 || incToUpdate.size() > 0);
    
    auto newNode = forest.createNode(par); // parent of includeGroup
    for (auto c : orderedToMove) {
        c->_detachThisNode();
        c->_setNextSib(nullptr);
        newNode->addChild(c);
        const auto & di = c->getData().desIds;
        newNode->getData().desIds.insert(begin(di), end(di));
    }
    for (auto gcp : incToUpdate) {
        gcp->first = newNode;
    }
    assert(!par->isOutDegreeOneNode());
    return newNode;
}
template<typename T, typename U>
OttIdSet FTree<T,U>::addPhyloStatementAtNode(const PhyloStatement & ps, 
                                             RootedTreeNode<T> * includeGroupA,
                                             const OttIdSet & attachedElsewhere) {
    OttIdSet r;
    for (auto oid : ps.includeGroup) {
        if (!ottIdIsConnected(oid)) {
            if (contains(attachedElsewhere, oid)) {
                addIncludeStatement(oid, includeGroupA, ps.provenance);
            } else {
                addLeafNoDesUpdate(includeGroupA, oid);
                r.insert(oid);
            }
        }
    }
    for (auto oid : ps.excludeGroup) {
        if (!ottIdIsConnected(oid)) {
            addExcludeStatement(oid, includeGroupA, ps.provenance);
        }
    }
    return r;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::getMRCA(const OttIdSet &ottIdSet) {
    if (ottIdSet.empty()) {
        assert(false);
        return nullptr;
    }
    const auto rel = set_intersection_as_set(ottIdSet, getConnectedOttIds());
    for (auto nextOttId : rel) {
        auto x = ottIdToNodeMap.find(nextOttId);
        assert(x != ottIdToNodeMap.end());
        node_type * aTip = x->second;
        assert(aTip != nullptr);
        if (ottIdSet.size() == 1) {
            return aTip;
        }
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    // reach here if none are connected.
    for (auto oid : ottIdSet) {
        auto x = ottIdToNodeMap.find(oid);
        assert(x != ottIdToNodeMap.end());
        node_type * aTip = x->second;
        assert(aTip != nullptr);
        if (ottIdSet.size() == 1) {
            return aTip;
        }
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    assert(false);
    return nullptr;
}

template<typename T, typename U>
FTree<T, U> & RootedForest<T,U>::addDisjointTree(const PhyloStatement &ps) {
    tree_type & r = createNewTree();
    r.mirrorPhyloStatement(ps);
    return r;
}

template<typename T, typename U>
void FTree<T, U>::mirrorPhyloStatement(const PhyloStatement &ps) {
    assert(root == nullptr);
    root = forest.createNode(nullptr);
    addPhyloStatementAsChildOfRoot(ps);
}

template<typename T, typename U>
void FTree<T, U>::addPhyloStatementAsChildOfRoot(const PhyloStatement &ps) {
    assert(root != nullptr);
    auto parOfIncGroup = forest.createNode(root); // parent of includeGroup
    supportedBy[parOfIncGroup].push_back(ps.provenance);
    assert(ps.excludeGroup.size() > 0);
    for (auto i : ps.excludeGroup) {
        addLeafNoDesUpdate(root, i);
    }
    for (auto i : ps.includeGroup) {
        addLeafNoDesUpdate(parOfIncGroup, i);
    }
    root->getData().desIds = ps.leafSet;
    parOfIncGroup->getData().desIds = ps.includeGroup;
    connectedIds.insert(begin(ps.leafSet), end(ps.leafSet));
}

template<typename T, typename U>
void RootedForest<T, U>::writeForestDOTToFN(const std::string &fn) const {
    LOG(DEBUG) << "     creating DOT file for forest: " << fn;
    std::ofstream outf(fn);
    writeDOTForest(outf, *this);
}

template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
