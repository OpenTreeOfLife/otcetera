#ifndef OTCETERA_SUPERTREE_UTIL_H
#define OTCETERA_SUPERTREE_UTIL_H

#include <map>
#include <set>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd,
                                const OttIdSet & oldEls,
                                const OttIdSet newEls,
                                std::map<const T *, NodeEmbedding<T, U> > & m);


template<typename T, typename U>
void reportOnConflicting(std::ostream & out,
                         const std::string & prefix,
                         const T * scaffold,
                         const std::set<PathPairing<T, U> *> & exitPaths,
                         const OttIdSet & phyloLeafSet);


template<typename T, typename U>
void copyStructureToResolvePolytomy(const T * srcPoly,
                                    U & destTree,
                                    typename U::node_type * destPoly);
// assumes that nd is the mrca of ingroup and outgroup IDs
template<typename T>
bool canBeResolvedToDisplay(const T *nd, const OttIdSet & ingroup, const OttIdSet & leafSet);


} // namespace
#endif
