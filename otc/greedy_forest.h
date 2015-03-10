#ifndef OTCETERA_GREEDY_FOREST_H
#define OTCETERA_GREEDY_FOREST_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/forest.h"
#include "otc/supertree_context.h"

namespace otc {
using CouldAddResult = std::tuple<bool, NodeWithSplits *, NodeWithSplits *>;
template<typename T, typename U>
class GreedyPhylogeneticForest: public RootedForest<RTSplits, MappedWithSplitsData> {
    public:
    GreedyPhylogeneticForest() {}
    GreedyPhylogeneticForest(const GreedyPhylogeneticForest &) = delete;
    GreedyPhylogeneticForest & operator=(const GreedyPhylogeneticForest &) = delete;
    bool attemptToAddGrouping(PathPairing<T, U> * ppptr,
                              const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits &sc);
    bool addLeaf(PathPairing<T, U> * ppptr,
                              const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits &sc) {
        return attemptToAddGrouping(ppptr, incGroup, leafSet, treeIndex, groupIndex, sc);
    }
    void finalizeTree(SupertreeContextWithSplits &sc);
    void setPossibleMonophyletic(U & /*scaffoldNode*/) {
        assert(false);
    }
    bool possibleMonophyleticGroupStillViable() {
        assert(false);
    }
    void finishResolutionOfThreadedClade(U & scaffoldNode, NodeEmbedding<T, U> * , SupertreeContextWithSplits & sc);
    private:
    CouldAddResult couldAddToTree(NodeWithSplits *r, const OttIdSet & incGroup, const OttIdSet & leafSet);
    void addIngroupAtNode(NodeWithSplits *r, NodeWithSplits *ing, NodeWithSplits *outg, const OttIdSet & incGroup, const OttIdSet & leafSet);
    void graftTreesTogether(NodeWithSplits *rr,
                            NodeWithSplits *ri,
                            NodeWithSplits *delr,
                            NodeWithSplits *deli,
                            NodeWithSplits *delo,
                            const OttIdSet & incGroup,
                            const OttIdSet & leafSet);
    std::set<OttIdSet> encountered; // used so that the PhyloStatement  have refs that outlive them
    bool addGroupToNewTree(const OttIdSet & incGroup,
                           const OttIdSet & leafSet,
                           int treeIndex,
                           long groupIndex);
    std::vector<T *> getRoots();
};

using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;

} // namespace otc
#endif
