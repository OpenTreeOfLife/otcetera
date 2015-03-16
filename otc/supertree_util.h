#ifndef OTCETERA_SUPERTREE_UTIL_H
#define OTCETERA_SUPERTREE_UTIL_H

#include <map>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include "otc/tree_data.h"
namespace otc {
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
                                    SupertreeContextWithSplits &);
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

} // namespace
#endif
