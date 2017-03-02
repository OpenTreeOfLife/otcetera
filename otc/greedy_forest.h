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
    bool attempt_to_add_grouping(const OttIdSet & incGroup,
                              const OttIdSet & leafSet,
                              int treeIndex,
                              long groupIndex,
                              SupertreeContextWithSplits *sc);
    bool add_leaf(const OttIdSet & incGroup, const OttIdSet & leafSet, int treeIndex, long groupIndex, SupertreeContextWithSplits *sc);
    void finish_resolution_of_embedded_clade(U & scaffoldNode, NodeEmbedding<T, U> * , SupertreeContextWithSplits * sc);
    private:
    bool create_and_add_phylo_statement(const OttIdSet & incGroup,
                           const OttIdSet & leafSet,
                           int treeIndex,
                           long groupIndex);
    void finalize_tree(SupertreeContextWithSplits *sc);
    std::pair<NodeWithSplits *, NodeWithSplits *> find_grandparent_that_is_root_or_band_sharing(
                      FTree<RTSplits, MappedWithSplitsData> & donorTree,
                      NodeWithSplits * nd,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree);
    std::vector<T *> get_roots();
    void merge_forest(SupertreeContextWithSplits *);
    bool merge_single_banded_tree(
                FTree<RTSplits, MappedWithSplitsData> &donorTree,
                InterTreeBand<RTSplits> * band,
                FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                SupertreeContextWithSplits *sc);
    void merge_trees_to_first_post_band_handling(SupertreeContextWithSplits *sc);
    NodeWithSplits* move_all_sibs(
                      NodeWithSplits * donorC,
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * attachPoint,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      SupertreeContextWithSplits *);
    std::pair<NodeWithSplits *, NodeWithSplits*> move_all_children(
                      NodeWithSplits * donorParent,
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * recipientNode,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      InterTreeBand<RTSplits> * bandBeingMerged);
    bool perform_single_band_merge(
                      std::size_t treeInd,
                      InterTreeBand<RTSplits> * itb,
                      const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
                      SupertreeContextWithSplits *sc);
    bool perform_set_of_single_band_merges(
                      std::size_t treeInd,
                      std::set<InterTreeBand<typename T::data_type> *> & itbSet,
                      const std::vector<FTree<RTSplits, MappedWithSplitsData> *> & sortedTrees,
                      SupertreeContextWithSplits *sc);
    void transfer_subtree_in_forest(
                      NodeWithSplits * des,
                      FTree<RTSplits, MappedWithSplitsData> & possDonor,
                      NodeWithSplits * newPar,
                      FTree<RTSplits, MappedWithSplitsData> & recipientTree,
                      FTree<RTSplits, MappedWithSplitsData> *donorTree,
                      InterTreeBand<RTSplits> * bandBeingMerged);
    bool zip_paths_from_barren_node(
                      FTree<RTSplits, MappedWithSplitsData> &donorTree,
                      NodeWithSplits * donorDes,
                      NodeWithSplits * donorAnc,
                      FTree<RTSplits, MappedWithSplitsData> &recipientTree,
                      NodeWithSplits * recipientDes,
                      NodeWithSplits * recipientAnc,
                      SupertreeContextWithSplits *);
    
    std::set<OttIdSet> encountered; // used so that the PhyloStatement  have refs that outlive them

};

using SupertreeContextWithSplits = SupertreeContext<NodeWithSplits, NodeWithSplits>;

} // namespace otc
#endif
