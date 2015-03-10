#ifndef OTCETERA_EMBEDDING_H
#define OTCETERA_EMBEDDING_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/util.h"
namespace otc {
template<typename T, typename U> class NodePairing;
template<typename T, typename U> class PathPairing;
template<typename T, typename U> class NodeEmbedding;
template<typename T, typename U> class SupertreeContext;

template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd,
                                 const OttIdSet & oldEls,
                                 const OttIdSet & newEls,
                                 std::map<const T *, NodeEmbedding<T, U> > & m);

/* a pair of aligned nodes from an embedding of a phylogeny onto a scaffold
   In NodeEmbedding objects two forms of these pairings are created:
        1. tips of the tree are mapped to the scaffoldNode assigned the same
            OTT Id,
        2. internal nodes in phylo tree are mapped to the least inclusive node
            in the scaffold tree that has all of the descendant OTT Ids from
            the phylo tree. (so the phyloNode->getData().desIds will be a subset
            of scaffoldNode->getData().desIds)
*/
template<typename T, typename U>
class NodePairing {
    public:
    T * scaffoldNode;
    U * phyloNode;
    NodePairing(T *taxo, U *phylo)
        :scaffoldNode(taxo),
        phyloNode(phylo) {
        assert(taxo != nullptr);
        assert(phylo != nullptr);
    }
};

/* Represents the mapping of an edge from phyloParent -> phyloChild onto
    a scaffold tree. The endpoints will be pairs of nodes that were aligned.
    Note that the phyloChild is a child of phyloParent, but scaffoldDes can
    be the same node as scaffoldAnc or it can be any descendant of that node.
    Thus an directed edge in the phylo tree pairs with a directed path in the
    scaffold.
*/
template<typename T, typename U>
class PathPairing {
    public:
    T * scaffoldDes;
    T * scaffoldAnc;
    U * phyloChild;
    U * phyloParent;
    OttIdSet currChildOttIdSet;
    void setOttIdSet(long oid, std::map<const U *, NodeEmbedding<T, U> > & m) {
        LOG(DEBUG) << "setOttIdSet " << oid;
        OttIdSet n;
        OttIdSet oldIds;
        std::swap(oldIds, currChildOttIdSet);
        if (contains(oldIds, oid)) {
            oldIds.erase(oid);
        } else {
            n.insert(oid);
        }
        updateAncestralPathOttIdSet(scaffoldDes, oldIds, n, m);
        currChildOttIdSet = n;
    }
    bool updateOttIdSetNoTraversal(const OttIdSet & oldEls, const OttIdSet & newEls);
    PathPairing(const NodePairing<T,U> & parent, const NodePairing<T,U> & child)
        :scaffoldDes(child.scaffoldNode),
        scaffoldAnc(parent.scaffoldNode),
        phyloChild(child.phyloNode),
        phyloParent(parent.phyloNode),
        currChildOttIdSet(child.phyloNode->getData().desIds) {
    }
    // as Paths get paired back deeper in the tree, the ID may be mapped to a higher
    // taxon. The currChildOttIdSet starts out identical to the phylogenetic node's 
    // descendant Ott Id set. But may change to reflect this remapping to the effective
    // set of IDs that include the tip.
    const OttIdSet & getOttIdSet() const {
        return currChildOttIdSet;
    }
    const OttIdSet & getPhyloChildDesID() const {
        return phyloChild->getData().desIds;
    }
};

template<typename T, typename U>
inline void updateAncestralPathOttIdSet(T * nd,
                                        const OttIdSet & oldEls,
                                        const OttIdSet & newEls,
                                        std::map<const T *, NodeEmbedding<T, U> > & m) {
    auto & curr = m.at(nd);
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
    void resolveGivenUncontestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc);

    void collapseGroup(U & scaffoldNode, SupertreeContext<T,U> & sc);
    void pruneCollapsedNode(U & scaffoldNode, SupertreeContextWithSplits & sc);
    void constructPhyloGraphAndCollapseIfNecessary(U & scaffoldNode, SupertreeContextWithSplits & sc);
    
    bool reportIfContested(std::ostream & out,
                           const U * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const;
    bool updateAllPathsOttIdSets(const OttIdSet & oldEls, const OttIdSet & newEls) {
        LOG(DEBUG) << "    updateAllPathsOttIdSets loops";
        bool r = updateAllMappedPathsOttIdSets(loopEmbeddings, oldEls, newEls);
        LOG(DEBUG) << "    updateAllPathsOttIdSets edgeBelowEmbeddings ";
        return updateAllMappedPathsOttIdSets(edgeBelowEmbeddings, oldEls, newEls) || r;
    }
    std::vector<const PathPairing<T, U> *> getAllIncomingPathPairs(const T * nd,
                                                                  const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
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
