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
class MRCABandedPaths{
    public:
    NodeWithSplits * mrcaF;
    NodeWithSplits * mrcaS;
    std::set<std::pair<NodeWithSplits *, NodeWithSplits *> > bandedPairs;
};
using MergeStartInfo = std::tuple<bool, std::list<MRCABandedPaths>, std::set<NodeWithSplits *> >; // banded, unbanded

using CouldAddResult = std::tuple<bool, NodeWithSplits *, NodeWithSplits *>;
template<typename T, typename U>
class GreedyPhylogeneticForest: public RootedForest<RTSplits, MappedWithSplitsData> {
    public:
    GreedyPhylogeneticForest(long rootOttId)
      :RootedForest<RTSplits, MappedWithSplitsData>(rootOttId) {}
    GreedyPhylogeneticForest(const GreedyPhylogeneticForest &) = delete;
    GreedyPhylogeneticForest & operator=(const GreedyPhylogeneticForest &) = delete;
    bool attemptToAddGrouping(PathPairing<T, U> * ppptr,
                              const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits *sc);
    bool addLeaf(PathPairing<T, U> * ppptr,
                              const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits *sc);
    void finalizeTree(SupertreeContextWithSplits *sc);
    void setPossibleMonophyletic(U & /*scaffoldNode*/) {
        NOT_IMPLEMENTED; //
    }
    bool possibleMonophyleticGroupStillViable() {
        NOT_IMPLEMENTED; //
    }
    void finishResolutionOfEmbeddedClade(U & scaffoldNode, NodeEmbedding<T, U> * , SupertreeContextWithSplits * sc);
    private:
    void mergeBandedTrees(FTree<RTSplits, MappedWithSplitsData> & donor,
                     FTree<RTSplits, MappedWithSplitsData> & recipient,
                     const MergeStartInfo & mi,
                     SupertreeContextWithSplits *);
    void mergeForest(SupertreeContextWithSplits *);
    CouldAddResult couldAddToTree(NodeWithSplits *r, const OttIdSet & incGroup, const OttIdSet & leafSet);
    /*
    void addIngroupAtNode(NodeWithSplits *r, NodeWithSplits *ing, NodeWithSplits *outg, const OttIdSet & incGroup, const OttIdSet & leafSet);
    void graftTreesTogether(NodeWithSplits *rr,
                            NodeWithSplits *ri,
                            NodeWithSplits *delr,
                            NodeWithSplits *deli,
                            NodeWithSplits *delo,
                            const OttIdSet & incGroup,
                            const OttIdSet & leafSet);
    */
    std::set<OttIdSet> encountered; // used so that the PhyloStatement  have refs that outlive them
    bool addGroupToNewTree(const OttIdSet & incGroup,
                           const OttIdSet & leafSet,
                           int treeIndex,
                           long groupIndex);
    std::vector<T *> getRoots();
    void mergePathToNextBand(FTree<RTSplits, MappedWithSplitsData> & donor,
                             NodeWithSplits * spikeDes,
                             FTree<RTSplits, MappedWithSplitsData> & recipient, 
                             SupertreeContextWithSplits *sc);
    void transferSubtreeInForest(NodeWithSplits * des,
                             FTree<RTSplits, MappedWithSplitsData> & possDonor,
                                                    NodeWithSplits * newPar,
                                                 FTree<RTSplits, MappedWithSplitsData> & recipientTree,
                                                 FTree<RTSplits, MappedWithSplitsData> *donorTree);
    void mergeTreesToFirstPostBandHandling(SupertreeContextWithSplits *sc);
};

using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;

} // namespace otc
#endif
