#ifndef OTCETERA_EMBEDDING_H
#define OTCETERA_EMBEDDING_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/pairings.h"
namespace otc {
template<typename T, typename U> class SupertreeContext;

template<typename T, typename U>
inline void updateAncestralPathOttIdSet(T * nd,
                                        const OttIdSet & oldEls,
                                        const OttIdSet & newEls,
                                        std::map<const T *, NodeEmbedding<T, U> > & m) {
    auto & curr = m.at(nd);
    assert(oldEls.size() > 0);
    LOG(DEBUG) << "  " << nd->getOttId() << " calling updateAllPathsOttIdSets";
    if (!curr.updateAllPathsOttIdSets(oldEls, newEls)) {
        return;
    }
    for (auto anc : iter_anc(*nd)) {
        auto & ant = m.at(anc);
        LOG(DEBUG) << "  " << anc->getOttId() << " calling updateAllPathsOttIdSets";
        if (!ant.updateAllPathsOttIdSets(oldEls, newEls)) {
            return;
        }
    }
}

/* A NodeEmbedding object holds pointers to all of the NodePairs and PathPairs
        relevant to a node (nodeWithEmbedding) in the scaffold Tree.
    These are stored in a map in which the key is the index of the phylo tree 
        involved in the mapping.
    The relevant pairings are:
        1. all node pairings for nodeWithEmbedding,
        2. all paths that are loops for nodeWithEmbedding (nodeWithEmbedding is the scaffoldDes
            and the scaffoldAnc for the PathPair), and
        3. all paths the leave this node. ie. those that paths that:
            A. have a scaffoldDes that is nodeWithEmbedding or one of its descendants, AND
            B. have scaffoldAnc set to an ancestor of nodeWithEmbedding
    Note that if you want all of the paths that intersect with nodeWithEmbedding, you have
        to check the edgeBelowEmbeddings field of all of nodeWithEmbedding's children. This
        is because any PathPairing that has nodeWithEmbedding as the scaffoldAnc will not
        be found in the edgeBelowEmbeddings for nodeWithEmbedding. And if it is not a loop
        node for nodeWithEmbedding, the path will not be in loopEmbeddings either. These
        paths are trivial wrt resolving the tree for nodeWithEmbedding, but they do 
        contribute to the relevant leaf sets
*/
template<typename T, typename U>
class NodeEmbedding {
    public:
    using NodePairSet = std::set<NodePairing<T, U> *>;
    using PathPairSet = std::set<PathPairing<T, U> *>;
    T * nodeWithEmbedding;
    NodeEmbedding(T * scaffNode)
        :nodeWithEmbedding(scaffNode) {
    }
    std::map<std::size_t, NodePairSet> nodeEmbeddings;
    std::map<std::size_t, PathPairSet> edgeBelowEmbeddings;
    std::map<std::size_t, PathPairSet > loopEmbeddings;
    std::size_t getTotalNumNodeMappings() const {
        unsigned long t = 0U;
        for (auto i : nodeEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getNumLoopTrees() const {
        std::set<std::size_t> keys;
        for (auto i : loopEmbeddings) {
            keys.insert(i.first);
        }
        return keys.size();
    }
    std::size_t getTotalNumLoops() const {
        unsigned long t = 0U;
        for (auto i : loopEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getTotalNumEdgeBelowTraversals() const {
        unsigned long t = 0U;
        for (auto i : edgeBelowEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    static bool treeContestsMonophyly(const PathPairSet & edgesBelowForTree);
    bool isContested() const {
        for (auto i : edgeBelowEmbeddings) {
            if (treeContestsMonophyly(i.second)) {
                return true;
            }
        }
        return false;
    }
    std::list<std::size_t> getContestingTrees() const {
        std::list<std::size_t> r;
        for (auto i : edgeBelowEmbeddings) {
            if (treeContestsMonophyly(i.second)) {
                r.push_back(i.first);
            }
        }
        return r;
    }
    const PathPairSet & getEdgesExiting(std::size_t treeIndex) const {
        auto el = edgeBelowEmbeddings.find(treeIndex);
        assert(el != edgeBelowEmbeddings.end());
        return el->second;
    }

    // some trees contest monophyly. Return true if these trees are obviously overruled
    //   by higher ranking trees so that we can avoid the more expensive unconstrained phylo graph 
    bool highRankingTreesPreserveMonophyly(std::size_t ) {
        return false;
    }

    const OttIdSet & getRelevantDesIdsFromPath(const PathPairing<T, U> & pps);
    OttIdSet getRelevantDesIdsFromPathPairSet(const PathPairSet & pps);
    OttIdSet getRelevantDesIds(const std::map<const T *, NodeEmbedding<T, U> > & eForNd, std::size_t treeIndex);

    void collapseSourceEdge(const T * phyloParent, PathPairing<T, U> * path);
    void collapseSourceEdgesToForceOneEntry(U & , PathPairSet & pps, std::size_t treeIndex, SupertreeContextWithSplits &);
    void resolveGivenContestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc);
    std::set<PathPairing<T, U> *> getAllChildExitPaths(U & scaffoldNode, SupertreeContextWithSplits & sc);
    std::set<PathPairing<T, U> *> getAllChildExitPathsForTree(U & scaffoldNode, std::size_t treeIndex, SupertreeContextWithSplits & sc);
    void resolveGivenUncontestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc);
    void exportSubproblemAndFakeResolution(U & scaffoldNode, const std::string & exportDir, SupertreeContextWithSplits & sc);
    void collapseGroup(U & scaffoldNode, SupertreeContext<T,U> & sc);
    void pruneCollapsedNode(U & scaffoldNode, SupertreeContextWithSplits & sc);
    void constructPhyloGraphAndCollapseIfNecessary(U & scaffoldNode, SupertreeContextWithSplits  & sc);
    
    bool reportIfContested(std::ostream & out,
                           const U * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const;
    bool updateAllPathsOttIdSets(const OttIdSet & oldEls, const OttIdSet & newEls) {
        bool r = updateAllMappedPathsOttIdSets(loopEmbeddings, oldEls, newEls);
        return updateAllMappedPathsOttIdSets(edgeBelowEmbeddings, oldEls, newEls) || r;
    }
    std::vector<const PathPairing<T, U> *> getAllIncomingPathPairs(const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                                                                  std::size_t treeIndex) const;
};


template<typename T, typename U>
inline bool NodeEmbedding<T, U>::treeContestsMonophyly(const std::set<PathPairing<T, U> *> & edgesBelowForTree) {
    if (edgesBelowForTree.size() > 1) {
        const T * firstSrcPar = nullptr;
        for (auto pp : edgesBelowForTree) {
            auto sp = pp->phyloParent;
            if (sp != firstSrcPar) {
                if (firstSrcPar == nullptr) {
                    firstSrcPar = sp;
                } else {
                    return true;
                }
            }
        }
    }
    return false;
}


template<typename T>
inline bool updateAllMappedPathsOttIdSets(T & mPathSets, const OttIdSet & oldEls, const OttIdSet & newEls) {
    bool r = false;
    for (auto mpIt : mPathSets) {
        for (auto p : mpIt.second) {
            r = r || p->updateOttIdSetNoTraversal(oldEls, newEls);
        }
    }
    return r;
}



} // namespace
#endif
