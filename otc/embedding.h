#ifndef OTCETERA_EMBEDDING_H
#define OTCETERA_EMBEDDING_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
namespace otc {
template<typename T, typename U> class NodePairing;
template<typename T, typename U> class PathPairing;
template<typename T, typename U> class NodeThreading;
template<typename T, typename U> class SupertreeContext;

template<typename T, typename U>
void updateAncestralPathOttIdSet(T * nd,
                                 const OttIdSet & oldEls,
                                 const OttIdSet newEls,
                                 std::map<const T *, NodeThreading<T, U> > & m);

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

template<typename T, typename U>
class PathPairing {
    public:
    T * scaffoldDes;
    T * scaffoldAnc;
    U * phyloChild;
    U * phyloParent;
    OttIdSet currChildOttIdSet;
    void setOttIdSet(long oid, std::map<const U *, NodeThreading<T, U> > & m) {
        OttIdSet n;
        n.insert(oid);
        updateAncestralPathOttIdSet(scaffoldAnc, currChildOttIdSet, n, m);
        currChildOttIdSet = n;
    }
    void updateOttIdSetNoTraversal(const OttIdSet & oldEls, const OttIdSet & newEls);
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
inline void updateAncestralPathOttIdSet(T * nd, const OttIdSet & oldEls, const OttIdSet newEls, std::map<const T *, NodeThreading<T, U> > & m) {
    for (auto anc : iter_anc(*nd)) {
         auto & ant = m.at(anc);
         ant.updateAllPathsOttIdSets(oldEls, newEls);
    }
}

template<typename T, typename U>
class NodeThreading {
    public:
    using NodePairSet = std::set<NodePairing<T, U> *>;
    using PathPairSet = std::set<PathPairing<T, U> *>;
    
    std::map<std::size_t, NodePairSet> nodeAlignments;
    std::map<std::size_t, PathPairSet> edgeBelowAlignments;
    std::map<std::size_t, PathPairSet > loopAlignments;
    std::size_t getTotalNumNodeMappings() const {
        unsigned long t = 0U;
        for (auto i : nodeAlignments) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getTotalNumLoops() const {
        unsigned long t = 0U;
        for (auto i : loopAlignments) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t getTotalNumEdgeBelowTraversals() const {
        unsigned long t = 0U;
        for (auto i : edgeBelowAlignments) {
            t += i.second.size();
        }
        return t;
    }
    static bool treeContestsMonophyly(const PathPairSet & edgesBelowForTree);
    bool isContested() const {
        for (auto i : edgeBelowAlignments) {
            if (treeContestsMonophyly(i.second)) {
                return true;
            }
        }
        return false;
    }
    std::list<std::size_t> getContestingTrees() const {
        std::list<std::size_t> r;
        for (auto i : edgeBelowAlignments) {
            if (treeContestsMonophyly(i.second)) {
                r.push_back(i.first);
            }
        }
        return r;
    }
    const PathPairSet & getEdgesExiting(std::size_t treeIndex) const {
        auto el = edgeBelowAlignments.find(treeIndex);
        assert(el != edgeBelowAlignments.end());
        return el->second;
    }

    // some trees contest monophyly. Return true if these trees are obviously overruled
    //   by higher ranking trees so that we can avoid the more expensive unconstrained phylo graph 
    bool highRankingTreesPreserveMonophyly(std::size_t ) {
        return false;
    }

    OttIdSet getRelevantDesIdsFromPath(const PathPairSet & pps);

    void collapseSourceEdge(const T * phyloParent, PathPairing<T, U> * path);
    void collapseSourceEdgesToForceOneEntry(U & , PathPairSet & pps, std::size_t treeIndex);
    void resolveGivenContestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc);
    std::set<PathPairing<T, U> *> getAllChildExitPaths(U & scaffoldNode, SupertreeContextWithSplits & sc);
    void resolveGivenUncontestedMonophyly(U & scaffoldNode, SupertreeContextWithSplits & sc);

    void collapseGroup(U & scaffoldNode, SupertreeContext<T,U> & sc);
    void pruneCollapsedNode(U & scaffoldNode, SupertreeContextWithSplits & sc);
    void constructPhyloGraphAndCollapseIfNecessary(U & scaffoldNode, SupertreeContextWithSplits & sc);
    OttIdSet getRelevantDesIds(std::size_t treeIndex);

    bool reportIfContested(std::ostream & out,
                           const U * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const;
    void updateAllPathsOttIdSets(const OttIdSet & oldEls, const OttIdSet & newEls) {
        updateAllMappedPathsOttIdSets(loopAlignments, oldEls, newEls);
        updateAllMappedPathsOttIdSets(edgeBelowAlignments, oldEls, newEls);
    }
};


template<typename T, typename U>
inline bool NodeThreading<T, U>::treeContestsMonophyly(const std::set<PathPairing<T, U> *> & edgesBelowForTree) {
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
inline void updateAllMappedPathsOttIdSets(T & mPathSets, const OttIdSet & oldEls, const OttIdSet & newEls) {
    for (auto lpIt : mPathSets) {
        for (auto lp : lpIt.second) {
            lp->updateOttIdSetNoTraversal(oldEls, newEls);
        }
    }
}



} // namespace
#endif
