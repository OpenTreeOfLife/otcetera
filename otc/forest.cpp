#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
namespace otc {

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
void RootedForest<T,U>::addPhyloStatement(const PhyloStatement &ps) {
    LOG(DEBUG) << " RootedForest::addPhyloStatement";
    std::cerr << " ingroup "; writeOttSet(std::cerr, " ", ps.includeGroup, " "); std::cerr << std::endl;
    std::cerr << " leafset "; writeOttSet(std::cerr, " ", ps.leafSet, " "); std::cerr << std::endl;
    ps.debugCheck();
    if (ps.isTrivial()) {
        auto newOttIds = set_difference_as_set(ps.includeGroup, ottIdSet);
        for (auto noid : newOttIds) {
            addDetachedLeaf(noid);
        }
        return;
    }
    if (areDisjoint(ps.leafSet, ottIdSet)) {
        addDisjointTree(ps);
    }
    assert(false); // not implemented
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
    auto c = forest.createNode(root); // parent of includeGroup
    supportedBy[c].push_back(ps.provenance);
    assert(ps.excludeGroup.size() > 0);
    for (auto i : ps.excludeGroup) {
        forest.createLeaf(root, i);
    }
    for (auto i : ps.includeGroup) {
        forest.createLeaf(c, i);
    }
    root->getData().desIds = ps.leafSet;
    c->getData().desIds = ps.includeGroup;
}


template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
