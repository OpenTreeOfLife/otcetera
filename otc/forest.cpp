#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
namespace otc {

bool PhyloStatement::debugCheck() const {
#ifdef DEBUGGING_PHYLO_STATEMENTS
    const OttIdSet ie = set_union_as_set(includeGroup, excludeGroup);
    if (ie != leafSet) {
        std::cerr << " includeGroup "; writeOttSet(std::cerr, " ", includeGroup, " "); std::cerr << std::endl;
        std::cerr << " excludeGroup "; writeOttSet(std::cerr, " ", excludeGroup, " "); std::cerr << std::endl;
        std::cerr << " leafSet "; writeOttSet(std::cerr, " ", leafSet, " "); std::cerr << std::endl;
        assert(false);
    }
#endif
    return true;
}

template<typename T, typename U>
typename RootedForest<T,U>::tree_type &
RootedForest<T,U>::createNewTree() {
    std::size_t i = nextTreeId++;
    auto r = trees.emplace(std::piecewise_construct,
                           std::forward_as_tuple(i),
                           std::forward_as_tuple(nextTreeId,
                                           *this,
                                           ottIdToNode));
    assert(r.second); // must be a new Tree!
    return trees.at(i);
}
   
template<typename T, typename U>
RootedForest<T,U>::RootedForest()
    :nextTreeId(0U),
    ottIdToNode(nodeSrc.getData().ottIdToNode) {
}

template<typename T, typename U>
bool RootedForest<T,U>::addPhyloStatement(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        LOG(DEBUG) << " RootedForest::addPhyloStatement";
        std::cerr << " incGroup "; writeOttSet(std::cerr, " ", ps.includeGroup, " "); std::cerr << std::endl;
        std::cerr << " leafSet "; writeOttSet(std::cerr, " ", ps.leafSet, " "); std::cerr << std::endl;
    }
    ps.debugCheck();
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
        std::cerr << " incGroup "; writeOttSet(std::cerr, " ", ps.includeGroup, " "); std::cerr << std::endl;
        std::cerr << " leafSet "; writeOttSet(std::cerr, " ", ps.leafSet, " "); std::cerr << std::endl;
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
using OverlapFTreePair = std::pair<OttIdSet, FTree<T, U> *>;

template<typename T, typename U>
void consumeMapToList(std::map<T, std::list<U> > &m, std::list<U> & out) {
    for (auto & mIt : m) {
        for (auto el = begin(mIt.second) ; el != end(mIt.second); ++el) {
            out.push_back(*el);
        }
    }
}


template<typename T, typename U>
std::list<OverlapFTreePair<T, U> > RootedForest<T,U>::getSortedOverlapping(const OttIdSet &inc) {
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
bool RootedForest<T,U>::addPhyloStatementToGraph(const PhyloStatement &ps) {
    if (debuggingOutputEnabled) {
        LOG(DEBUG) << " RootedForest::addPhyloStatementToGraph";
        std::cerr << " incGroup "; writeOttSet(std::cerr, " ", ps.includeGroup, " "); std::cerr << std::endl;
        std::cerr << " leafSet "; writeOttSet(std::cerr, " ", ps.leafSet, " "); std::cerr << std::endl;
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
    auto byIncCardinality = getSortedOverlapping(ps.includeGroup);
    LOG(DEBUG) << byIncCardinality.size() << " FTree instance referred to in byIncCardinality";
    if (byIncCardinality.empty()) {
        LOG(DEBUG) << "No intersection between includeGroup of an existing FTree.";
        // this ingroup does not overlap with any ftree. find the FTree with the most overlap
        //  with the excludeGroup...
        auto byExcCardinality = getSortedOverlapping(ps.excludeGroup);
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
        return true;
    }
    auto & attachmentPair = *byIncCardinality.begin();
    if (attachmentPair.first.size() == 1) {
        assert(false);
        // greedy approach is to add the rest of the ingroup as deep in this tree as possible.
        //  less greedy: make include/exclude statements at that node
        return true;
    } else {
        std::list<node_type * > nonTrivMRCAs;
        OttIdSet attachedElsewhere;
        std::vector<bool> shouldResolve;
        for (const auto & incPair : byIncCardinality) {
            const auto & incGroupIntersection = incPair.first;
            if (incGroupIntersection.size() == 1) {
                break; // must be compatible because we've hit the singelton nodes...
            }
            attachedElsewhere.insert(incGroupIntersection.begin(), incGroupIntersection.end());
            tree_type * f = incPair.second;
            auto includeGroupA = f->getMRCA(incGroupIntersection);
            // If any of the ingroup are specifically excluded, then we have move deeper in the tree.
            // TMP this could be more efficient and avoid the while loop.
            while (f->anyExcludedAtNode(includeGroupA, ps.includeGroup)) {
                includeGroupA = includeGroupA->getParent();
                assert(includeGroupA != nullptr);
            }
            const OttIdSet excInc = set_intersection_as_set(includeGroupA->getData().desIds, ps.excludeGroup);
            if (debuggingOutputEnabled) {
                LOG(DEBUG) << "     addPhyloStatementToGraph search for an ancestor of ..."; 
                std::cerr << " addPhyloStatementToGraph search for an ancestor of:  "; writeOttSet(std::cerr, " ", incGroupIntersection, " "); std::cerr << std::endl;
                std::cerr << "  wanted to avoid =  "; writeOttSet(std::cerr, " ", ps.excludeGroup, " "); std::cerr << std::endl;
                std::cerr << "  found a node with desIds:  "; writeOttSet(std::cerr, " ", includeGroupA->getData().desIds, " "); std::cerr << std::endl;
                std::cerr << "  which includes the excludegroup members:  "; writeOttSet(std::cerr, " ", excInc, " "); std::cerr << std::endl;
            }
            if (!canBeResolvedToDisplayExcGroup(includeGroupA, ps.includeGroup, excInc)) {
                return false; // the MRCA of the includeGroup had interdigitated members of the excludeGroup
            }
            shouldResolve.push_back(!excInc.empty());
            nonTrivMRCAs.push_back(includeGroupA);
        }
        // all non trivial overlapping trees have approved this split...
        auto ntmIt = begin(nonTrivMRCAs);
        auto srIt = begin(shouldResolve);
        for (const auto & incPair : byIncCardinality) {
            const auto & incGroupIntersection = incPair.first;
            if (incGroupIntersection.size() == 1) {
                break; // must be compatible because we've hit the singelton nodes...
            }
            tree_type * f = incPair.second;
            assert(ntmIt != nonTrivMRCAs.end());
            node_type * includeGroupA = *ntmIt++;
            const bool addNode = *srIt++;
            if (addNode) {
                includeGroupA = f->resolveToCreateCladeOfIncluded(includeGroupA, ps.includeGroup);
            }
            auto connectedHere = f->addPhyloStatementAtNode(ps, includeGroupA, attachedElsewhere);
            if (!connectedHere.empty()) {
                attachedElsewhere.insert(begin(connectedHere), end(connectedHere));
            }
        }
        return true;
    } 
}

template<typename T, typename U>
bool FTree<T,U>::anyExcludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    const node_type * p = nd->getParent();
    const OttIdSet & ndi = nd->getData().desIds;
    for (auto oid : ottIdSet) {
        auto leafNd = ottIdToNode.at(oid);
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
    RootedTreeNode<T> * eNode = ottIdToNode.at(ottId);
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
    RootedTreeNode<T> * eNode = ottIdToNode.at(ottId);
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
        auto n = ottIdToNode.at(oid);
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
        auto x = ottIdToNode.find(nextOttId);
        if (x == ottIdToNode.end()) {
            continue;
        }
        node_type * aTip = x->second;
        assert(aTip != nullptr);
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    return nullptr;
}

template<typename T, typename U>
FTree<T, U> & RootedForest<T,U>::addDisjointTree(const PhyloStatement &ps) {
    assert(areDisjoint(ps.leafSet, ottIdSet)); //
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



template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
