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
class GreedyBandedForest: public RootedForest<RTSplits, MappedWithSplitsData> {
    public:
    GreedyBandedForest(long rootOttId)
      :RootedForest<RTSplits, MappedWithSplitsData>(rootOttId) {}
    GreedyBandedForest(const GreedyBandedForest &) = delete;
    GreedyBandedForest & operator=(const GreedyBandedForest &) = delete;
    bool attemptToAddGrouping(const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits *sc);
    bool addLeaf(const OttIdSet & incGroup, const OttIdSet & leafSet, int treeIndex, long groupIndex, SupertreeContextWithSplits *sc);
    void finalizeTree(SupertreeContextWithSplits *sc);
    void writeFirstTree(std::ostream & treeFileStream);
    void setPossibleMonophyletic(U & /*scaffoldNode*/) {
        assert("not implemented"[0] == 'f');; //
    }
    bool possibleMonophyleticGroupStillViable() {
        assert("not implemented"[0] == 'f');; //
    }
    void finishResolutionOfEmbeddedClade(U & scaffoldNode, NodeEmbedding<T, U> * , SupertreeContextWithSplits * sc);
    private:
    void mergeForest(SupertreeContextWithSplits *);
    bool mergeSingleBandedTree(
                FTree<RTSplits, MappedWithSplitsData> &donorTree,
                InterTreeBand<RTSplits> * band,
                FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                SupertreeContextWithSplits *sc);
    CouldAddResult couldAddToTree(NodeWithSplits *r, const OttIdSet & incGroup, const OttIdSet & leafSet);
    /*
    void mergeBandedTrees(FTree<RTSplits, MappedWithSplitsData> & donor,
                     FTree<RTSplits, MappedWithSplitsData> & recipient,
                     const MergeStartInfo & mi,
                     SupertreeContextWithSplits *);
    void addIngroupAtNode(NodeWithSplits *r, NodeWithSplits *ing, NodeWithSplits *outg, const OttIdSet & incGroup, const OttIdSet & leafSet);
    void graftTreesTogether(NodeWithSplits *rr,
                            NodeWithSplits *ri,
                            NodeWithSplits *delr,
                            NodeWithSplits *deli,
                            NodeWithSplits *delo,
                            const OttIdSet & incGroup,
                            const OttIdSet & leafSet);
    void mergePathToNextBand(FTree<RTSplits, MappedWithSplitsData> & donor,
                             NodeWithSplits * spikeDes,
                             FTree<RTSplits, MappedWithSplitsData> & recipient, 
                             SupertreeContextWithSplits *sc);
    */
    std::set<OttIdSet> encountered; // used so that the PhyloStatement  have refs that outlive them
    bool createAndAddPhyloStatement(const OttIdSet & incGroup,
                           const OttIdSet & leafSet,
                           int treeIndex,
                           long groupIndex);
    std::vector<T *> getRoots();
    void transferSubtreeInForest(
                      NodeWithSplits * des,
                      FTree<RTSplits, MappedWithSplitsData> & possDonor,
                      NodeWithSplits * newPar,
                      FTree<RTSplits, MappedWithSplitsData> & recipientTree,
                      FTree<RTSplits, MappedWithSplitsData> *donorTree,
                      InterTreeBand<RTSplits> * bandBeingMerged);
    void mergeTreesToFirstPostBandHandling(SupertreeContextWithSplits *sc);
    std::pair<NodeWithSplits *, NodeWithSplits*> moveAllChildren(
                      NodeWithSplits * donorParent,
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * recipientNode,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      InterTreeBand<RTSplits> * bandBeingMerged);
    bool performSingleBandMerge(
                      std::size_t treeInd,
                      InterTreeBand<RTSplits> * itb,
                      const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
                      SupertreeContextWithSplits *sc);
    bool performSetOfSingleBandMerges(
                      std::size_t treeInd,
                      std::set<InterTreeBand<typename T::data_type> *> & itbSet,
                      const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
                      SupertreeContextWithSplits *sc);
    std::pair<NodeWithSplits *, NodeWithSplits *> findGrandparentThatIsRootOrBandSharing(
                      FTree<RTSplits, MappedWithSplitsData> & donorTree,
                      NodeWithSplits * nd,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree);
    bool zipPathsFromBarrenNode(
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * donorDes,
                      NodeWithSplits * donorAnc,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      NodeWithSplits * recipientDes,
                      NodeWithSplits * recipientAnc,
                      SupertreeContextWithSplits *);
    NodeWithSplits* moveAllSibs(
                      NodeWithSplits * donorC,
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * attachPoint,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      SupertreeContextWithSplits *);
};

using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;

} // namespace otc
#endif
