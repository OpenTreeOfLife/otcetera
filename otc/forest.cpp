#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
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
        FTree<T, U> * ftree = &(tpIt.second);
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
    if (byIncCardinality.empty()) {
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
    assert(false);
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
        forest.createLeaf(root, i);
    }
    for (auto i : ps.includeGroup) {
        forest.createLeaf(parOfIncGroup, i);
    }
    root->getData().desIds = ps.leafSet;
    parOfIncGroup->getData().desIds = ps.includeGroup;
    connectedIds.insert(begin(ps.leafSet), end(ps.leafSet));
}



template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
