#ifndef OTCETERA_SUPERTREE_UTIL_H
#define OTCETERA_SUPERTREE_UTIL_H

#include <map>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include "otc/tree_data.h"
namespace otc {
std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &);

template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd,
                                const OttIdSet & oldEls,
                                const OttIdSet & newEls,
                                std::map<const T *, NodeEmbedding<T, U> > & m);


template<typename T, typename U>
void reportOnConflicting(std::ostream & out,
                         const std::string & prefix,
                         const T * scaffold,
                         const std::set<PathPairing<T, U> *> & exitPaths,
                         const OttIdSet & phyloLeafSet);

// takes 2 "includeGroups" from different PhyloStatements.
//  `culled` has been pruned down to the common leafSet, and the common leafSet is passed in as `leafSet
bool culledAndCompleteIncompatWRTLeafSet(const OttIdSet & culled, const OttIdSet & complete, const OttIdSet & leafSet);

template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly,
                                    SupertreeContextWithSplits *);
// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
bool canBeResolvedToDisplay(const T *nd, const OttIdSet & incGroup, const OttIdSet & leafSet);
template<typename T>
bool canBeResolvedToDisplayExcGroup(const T *nd, const OttIdSet & incGroup, const OttIdSet & excGroup);


// assumes that nd is the mrca of incGroup and excGroup IDs
template<typename T>
inline bool canBeResolvedToDisplay(const T *nd, const OttIdSet & incGroup, const OttIdSet & leafSet) {
    const OttIdSet excGroup = set_difference_as_set(leafSet, incGroup);
    return canBeResolvedToDisplayExcGroup(nd, incGroup, excGroup);
}

template<typename T>
void addDesIdsToNdAndAnc(T * nd, const OttIdSet & oid) {
    nd->getData().desIds.insert(begin(oid), end(oid));
    for (auto anc : iter_anc(*nd)) {
        anc->getData().desIds.insert(begin(oid), end(oid));
    }
}

template<typename T>
void removeDesIdsToNdAndAnc(T * nd, const OttIdSet & oid) {
    nd->getData().desIds = set_difference_as_set(nd->getData().desIds, oid);
    for (auto anc : iter_anc(*nd)) {
        anc->getData().desIds = set_difference_as_set(anc->getData().desIds, oid);
    }
}

// used post-suppression of monotypic taxa to create a map
//  from the alias to the original OTT ID.
template<typename T>
inline std::map<long, long> generateIdRemapping(const T & tree) {
    const auto & id2ndMap = tree.getData().ottIdToNode;
    std::map<long, long> r;
    for (const auto & idNdPair : id2ndMap) {
        const auto & inID = idNdPair.first;
        const auto outID = idNdPair.second->getOttId();
        if (outID != inID) {
            assert(!contains(r, inID));
            r[inID] = outID;
        }
    }
    return r;
}

//currently not copying names
inline std::unique_ptr<TreeMappedWithSplits> cloneTree(const TreeMappedWithSplits &tree) {
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
                throw OTCError("asserts false but not enabled");
            }
            nn->getData().desIds = nd->getData().desIds;
        }
    } catch (...) {
        delete rawTreePtr;
        throw;
    }
    return std::unique_ptr<TreeMappedWithSplits>(rawTreePtr);
}

} // namespace
#endif
