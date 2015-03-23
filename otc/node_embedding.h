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
        relevant to a node (embeddedNode) in the scaffold Tree.
    These are stored in a map in which the key is the index of the phylo tree 
        involved in the mapping.
    The relevant pairings are:
        1. all node pairings for embeddedNode,
        2. all paths that are loops for embeddedNode (embeddedNode is the scaffoldDes
            and the scaffoldAnc for the PathPair), and
        3. all paths the leave this node. ie. those that paths that:
            A. have a scaffoldDes that is embeddedNode or one of its descendants, AND
            B. have scaffoldAnc set to an ancestor of embeddedNode
    Note that if you want all of the paths that intersect with embeddedNode, you have
        to check the edgeBelowEmbeddings field of all of embeddedNode's children. This
        is because any PathPairing that has embeddedNode as the scaffoldAnc will not
        be found in the edgeBelowEmbeddings for embeddedNode. And if it is not a loop
        node for embeddedNode, the path will not be in loopEmbeddings either. These
        paths are trivial wrt resolving the tree for embeddedNode, but they do 
        contribute to the relevant leaf sets
*/
template<typename T, typename U>
class NodeEmbedding {
    using NodePairPtr = NodePairing<T, U> *;
    using PathPairPtr = PathPairing<T, U> *;
    using NodePairSet = std::set<NodePairPtr>;
    using PathPairSet = std::set<PathPairPtr>;
    using TreeToNodePairs = std::map<std::size_t, NodePairSet>;
    using TreeToPathPairs = std::map<std::size_t, PathPairSet>;
    T * embeddedNode;
    TreeToNodePairs nodeEmbeddings;
    TreeToPathPairs edgeBelowEmbeddings;
    TreeToPathPairs loopEmbeddings;
    public:
    NodeEmbedding(T * scaffNode)
        :embeddedNode(scaffNode) {
    }
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
    OttIdSet getRelevantDesIds(const std::map<const T *,
                               NodeEmbedding<T, U> > & eForNd,
                               std::size_t treeIndex);

    void collapseSourceEdge(const T * phyloParent,
                            PathPairing<T, U> * path);
    void collapseSourceEdgesToForceOneEntry(T & ,
                                            PathPairSet & pps,
                                            std::size_t treeIndex,
                                            SupertreeContextWithSplits &);
    void resolveGivenContestedMonophyly(T & scaffoldNode,
                                        SupertreeContextWithSplits & sc);
    std::set<PathPairPtr> getAllChildExitPaths(
                            const T & scaffoldNode,
                            const std::map<const T *, NodeEmbedding<T, U> > & sc) const;
    std::set<PathPairPtr> getAllChildExitPathsForTree(
                            const T & scaffoldNode,
                            std::size_t treeIndex,
                            const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const;
    void resolveGivenUncontestedMonophyly(T & scaffoldNode,
                                          SupertreeContextWithSplits & sc);
    void exportSubproblemAndFakeResolution(T & scaffoldNode,
                                           const std::string & exportDir,
                                           std::ostream * exportStream, // nonnull to override exportdir
                                           SupertreeContextWithSplits & sc);
    void collapseGroup(T & scaffoldNode, SupertreeContext<T, U> & sc);
    void pruneCollapsedNode(T & scaffoldNode, SupertreeContextWithSplits & sc);
    void constructPhyloGraphAndCollapseIfNecessary(T & scaffoldNode, SupertreeContextWithSplits  & sc);
    
    bool reportIfContested(std::ostream & out,
                           const T * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const;
    bool updateAllPathsOttIdSets(const OttIdSet & oldEls, const OttIdSet & newEls) {
        bool r = updateAllMappedPathsOttIdSets(loopEmbeddings, oldEls, newEls);
        return updateAllMappedPathsOttIdSets(edgeBelowEmbeddings, oldEls, newEls) || r;
    }
    std::vector<const PathPairing<T, U> *> getAllIncomingPathPairs(
                        const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                        std::size_t treeIndex) const;
    bool debugNodeEmbedding(bool isUncontested,
                            const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const;
    void addNodeEmbedding(std::size_t treeIndex, NodePairPtr npp) {
        nodeEmbeddings[treeIndex].insert(npp);
    }
    void addLoopEmbedding(std::size_t treeIndex, PathPairPtr pp) {
        loopEmbeddings[treeIndex].insert(pp);
    }
    void addExitEmbedding(std::size_t treeIndex, PathPairPtr pp) {
        edgeBelowEmbeddings[treeIndex].insert(pp);
    }
    void setOttIdForExitEmbeddings(
                        T * newScaffDes,
                        long ottId,
                        std::map<const T *, NodeEmbedding<T, U> > & n2ne);
    const TreeToPathPairs & getExitEmbeddings() const {
        return edgeBelowEmbeddings;
    }
    private:
    std::map<U *, U*> getLoopedPhyloNd2Par(std::size_t treeInd) const;
    std::map<U *, U*> getExitPhyloNd2Par(std::size_t treeInd) const;
    PathPairSet refersToNode(std::size_t treeInd, U *n) const {
        PathPairSet r;
        auto pit = edgeBelowEmbeddings.find(treeInd);
        if (pit != edgeBelowEmbeddings.end()) {
            for (auto i : pit->second) {
                if (i->phyloParent == n || i->phyloChild == n) {
                    r.insert(i);
                }
            }
        }
        pit = loopEmbeddings.find(treeInd);
        if (pit != loopEmbeddings.end()) {
            for (auto i : pit->second) {
                if (i->phyloParent == n || i->phyloChild == n) {
                    r.insert(i);
                }
            }
        }
        return r;
    }

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

template<typename T, typename U, typename V>
inline std::map<U *, U*> getNd2ParForKey(const V & treeInd,
                                         const std::map<V, std::set<PathPairing<T, U> *> >& m);
template<typename T, typename U, typename V>
inline std::map<U *, U*> getNd2ParForKey(const V & treeInd,
                                         const std::map<V, std::set<PathPairing<T, U> *> >& m) {
    std::map<U *, U*> nd2par;
    if (!contains(m, treeInd)) {
        return nd2par;
    }
    for (const auto & pp : m.at(treeInd)) {
        LOG(DEBUG) << " considering the edge from the child "  << (long)pp->phyloChild << ": "; dbWriteNewick(pp->phyloChild);
        LOG(DEBUG) << "             to its parent "  << (long)pp->phyloParent << ": "; dbWriteNewick(pp->phyloParent);
        assert(!contains(nd2par, pp->phyloChild));
        nd2par[pp->phyloChild] = pp->phyloParent;
    }
    return nd2par;
}

template<typename T, typename U>
inline std::map<U *, U*> NodeEmbedding<T, U>::getLoopedPhyloNd2Par(std::size_t treeInd) const {
    return getNd2ParForKey(treeInd, loopEmbeddings);
}

template<typename T, typename U>
inline std::map<U *, U*> NodeEmbedding<T, U>::getExitPhyloNd2Par(std::size_t treeInd) const {
    return getNd2ParForKey(treeInd, edgeBelowEmbeddings);
}

} // namespace
#endif
