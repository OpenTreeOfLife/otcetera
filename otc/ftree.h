#ifndef OTCETERA_FTREE_H
#define OTCETERA_FTREE_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/util.h"
namespace otc {
template<typename T, typename U> class GreedyPhylogeneticForest;
template<typename T, typename U> class RootedForest;
template<typename T, typename U> class FTree;

enum ConstraintType {
    INCLUDES_NODE,
    EXCLUDES_NODE
};
struct PhyloStatementSource {
    PhyloStatementSource(int treeInd, long groupInd)
        :sourceTreeId(treeInd),
        cladeId(groupInd) {
    }
    const int sourceTreeId;
    const long cladeId; // ID of the node in the tree
};

struct PhyloStatement {
    /*PhyloStatement(const OttIdSet &includes, const OttIdSet & other, bool otherIsExcludes)
        :includeGroup(includes),
        excludeGroup(otherIsExcludes ? other : set_difference_as_set(other, includes)),
        leafSet(otherIsExcludes ? set_union_as_set(other, includes): other) {
    }*/
    PhyloStatement(const OttIdSet &includes,
                   const OttIdSet &excludes,
                   const OttIdSet &mentioned, 
                   PhyloStatementSource pss)
        :includeGroup(includes),
        excludeGroup(excludes),
        leafSet(mentioned),
        provenance(pss) {
        debugCheck();
    }
    void writeAsNewick(std::ofstream & out) const;
    bool debugCheck() const;
    bool isTrivial() const {
        return includeGroup.size() < 2 || includeGroup.size() == leafSet.size();
    }
    const OttIdSet & includeGroup;
    const OttIdSet & excludeGroup;
    const OttIdSet & leafSet; // just the union of the includeGroup and excludeGroup
    const PhyloStatementSource provenance;
};

// Bands connect nodes across different FTree instance when the nodes are required
//  to satisfy a grouping that has been added. The includeGroup of the ps
//  is a set of IDs that will be in the desIds of each member of the band (and
//  the ancestor nodes of the members).
//  The nodes in the ps.excludeGroup should be excluded at this node or an ancestor.
template<typename T>
class InterTreeBand {
    public:
    using node_type = RootedTreeNode<T>;
    using node_set = std::set<node_type *>;
    InterTreeBand(node_type * nd1,
                  const node_set & n1set,
                  const PhyloStatement & ps)
        :statement(ps) {
        addNode(nd1, n1set);
    }
    bool isSingleTreeBand() const {
        return (nd2phantom.size() == 1);
    }
    void insertSet(node_type * nd, const node_set & phantom) {
        for (auto p : phantom) {
            assert(p->hasOttId());
        }
        nd2phantom.at(nd).insert(begin(phantom), end(phantom));
    }
    void insert(node_type * nd, node_type * phantom) {
        assert(phantom->hasOttId());
        nd2phantom[nd].insert(phantom);
    }
    void addNode(node_type * nd, const node_set & nset) {
        for (auto p : nset) {
            assert(p->hasOttId());
        }
        assert(nd != nullptr);
        const auto s = nd2phantom.size();
        nd2phantom[nd] = nset;
        assert((s + 1)== nd2phantom.size());
    }
    bool bandPointHasAny(const node_type * nd, const node_set & ns) const {
        return !(areDisjoint(nd2phantom.at(nd), ns));
    }
    OttIdSet getPhantomIds(const node_type * nd) const {
        OttIdSet r;
        auto npIt = nd2phantom.find(nd);
        if (npIt != nd2phantom.end()) {
            for (auto np : npIt->second) {
                r.insert(np->getOttId());
            }
        }
        return r;
    }
    void reassignAttachmentNode(node_type * oldAnc, node_type *newAnc) {
        assert(!contains(nd2phantom, newAnc));
        nd2phantom[newAnc] = nd2phantom.at(oldAnc);
        nd2phantom.erase(oldAnc);
    }
    bool isTheSetOfPhantomNodes(node_type * nd, const node_set & t) const {
        return nd2phantom.at(nd) == t;
    }
    std::set<const node_type *> getBandedNodes() const {
        std::set<const node_type *> r;
        for (const auto & m : nd2phantom) {
            r.insert(m.first);
        }
        return r;
    }
    void debugInvariantsCheckITB() const;
    private:
    std::map<const node_type *, node_set> nd2phantom;
    const PhyloStatement & statement;
};

template<typename T>
class ExcludeConstraints {
    public:
    using node_type = RootedTreeNode<T>;
    using node_pair = std::pair<const node_type *, const node_type *>;
    using cnode_set = std::set<const node_type *>;
    using node2many_map = std::map<const node_type *, cnode_set >;
    bool addExcludeStatement(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    bool isExcludedFrom(const node_type * ndToCheck,
                        const node_type * potentialAttachment,
                        const std::map<long, node_type*> * ottIdToNodeMap) const;
    void debugInvariantsCheckEC() const;
    bool hasNodesExcludedFromIt(const node_type *n) const {
        return contains(byNdWithConstraints, n);
    }
    const cnode_set & getNodesExcludedFromNode(const node_type * nd) const {
        auto bc = byNdWithConstraints.find(nd);
        if (bc == byNdWithConstraints.end()) {
            return emptySet;
        }
        return bc->second;
    }
    cnode_set stealExclusions(node_type *nd) {
        auto bc = byNdWithConstraints.find(nd);
        if (bc == byNdWithConstraints.end()) {
            return cnode_set{};
        }
        cnode_set r = bc->second;
        for (auto n : r) {
            byExcludedNd[n].erase(nd);
        }
        return r;
    }
    private:
    void purgeExcludeRaw(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    void ingestExcludeRaw(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    node2many_map byExcludedNd;
    node2many_map byNdWithConstraints;
    const cnode_set emptySet;
};

template<typename T>
class InterTreeBandBookkeeping {
    public:
    using node_type = RootedTreeNode<T>;
    using band_type = InterTreeBand<T>;
    using band_set = std::set<band_type *>;
    OttIdSet getPhantomIds(const node_type * nd) const {
        OttIdSet r;
        const band_set & bs =  getBandsForNode(nd);
        for (const auto & b : bs) {
            //LOG(DEBUG) << "  reading from band " << (long) b << ' ' << std::hex << (long) b << std::dec;
            const OttIdSet bo = b->getPhantomIds(nd);
            r.insert(begin(bo), end(bo));
        }
        return r;
    }
    void reassignAttachmentNode(band_type * b, node_type * oldAnc, node_type * newAnc, const PhyloStatement & ps);
    const band_set & getBandsForNode(const node_type *n) const {
        const auto nit = node2Band.find(n);
        if (nit == node2Band.end()) {
            return emptySet;
        }
        return nit->second;
    }
    band_set stealBands(const node_type *n) {
        const auto nit = node2Band.find(n);
        if (nit == node2Band.end()) {
            return emptySet;
        }
        band_set bs = nit->second;
        node2Band.erase(n);
        for (auto bandPtr : bs) {
            band2Node.erase(bandPtr);
        }
        return bs;
    }
    bool isInABand(const node_type *n) const {
        return contains(node2Band, n);
    }
    void _addRefToBand(band_type * band, node_type * nd) {
        assert(nd != nullptr);
        assert(band != nullptr);
        band2Node[band] = nd;
        node2Band[nd].insert(band);
    }
    private:
    std::map<const band_type *, node_type *> band2Node;
    std::map<const node_type *, band_set > node2Band;
    const band_set emptySet;
};

template<typename T, typename U>
class FTree {
    public:
    using node_data_type = T;
    using node_type = RootedTreeNode<T>;
    using GroupingConstraint = std::pair<node_type*, PhyloStatementSource>;
    using NdToConstrainedAt = std::map<node_type *, std::set<node_type *> >;
    FTree(std::size_t treeID,
          RootedForest<T, U> & theForest,
          std::map<long, node_type *> & ottIdToNodeRef)
        :treeId(treeID),
         root(nullptr),
         forest(theForest),
         ottIdToNodeMap(ottIdToNodeRef) {
    }
    // const methods:
    bool isInABand(const node_type *n) const {
        return bands.isInABand(n);
    }
    bool hasNodesExcludedFromIt(const node_type *n) const {
        return exclude.hasNodesExcludedFromIt(n);
    }
    const node_type * getRoot() const {
        return root;
    }
    bool ottIdIsExcludedFromRoot(long oid) const {
        return isExcludedFromRoot(ottIdToNodeMap.at(oid));
    }
    bool isExcludedFromRoot(const node_type *n) const {
        return exclude.isExcludedFrom(n, root, &ottIdToNodeMap);
    }
    bool isExcludedFrom(const node_type * ndToCheck, const node_type * potentialAttachment) const {
        return exclude.isExcludedFrom(ndToCheck, potentialAttachment, &ottIdToNodeMap);
    }
    
    // OTT Ids of nodes on the graph only....
    const OttIdSet getConnectedOttIds() const;
    // includes OTT Ids of nodes in includesConstraints
    const OttIdSet & getIncludedOttIds() {
        return getRoot()->getData().desIds;
    }
    bool ottIdIsConnected(long ottId) const {
        return contains(getConnectedOttIds(), ottId);
    }
    const ExcludeConstraints<T> & getExclusions() const {
        return exclude;
    }
    const std::set<InterTreeBand<T> *> & getBandsForNode(const node_type *n) const {
        return bands.getBandsForNode(n);
    }
    // non-const
    node_type * getRoot() {
        return root;
    }
    void _setRoot(node_type *r) {
        root = r;
    }
    void mirrorPhyloStatement(const PhyloStatement & ps);
    void debugInvariantsCheckFT() const;
    void debugVerifyDesIdsAssumingDes(const OttIdSet &s, const RootedTreeNode<T> *nd) const;
    std::set<RootedTreeNode<T> *> ottIdSetToNodeSet(const OttIdSet &ottIdSet) const;
    bool anyExcludedAtNode(const node_type *, const OttIdSet &) const ;
    void createDeeperRoot();
    // puts a node between nd and its parent and returns the new node
    node_type * createDeeperNode(node_type *nd);
    void stealExclusionStatements(node_type * newPar,  node_type * srcNode, FTree<T, U>  & donorTree);
    OttIdSet stealInclusionStatements(node_type * newPar,  node_type * srcNode, FTree<T, U>  & donorTree);
    void registerExclusionStatementForTransferringNode(node_type * srcNode, FTree<T, U>  & donorTree);
    void registerInclusionStatementForTransferringNode(node_type * srcNode, FTree<T, U>  & donorTree);
    private:
    void addExcludeStatement(long ottId, RootedTreeNode<T> *, const PhyloStatementSource &);
    void addIncludeGroupDisjointPhyloStatement(const PhyloStatement & ps) {
        addPhyloStatementAsChildOfRoot(ps);
    }
    RootedTreeNode<T> * addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId);
    void addPhyloStatementAsChildOfRoot(const PhyloStatement &);
    // this is greedy, we should be building separate FTree instances in many cases....
    OttIdSet addPhyloStatementAtNode(const PhyloStatement & ps, 
                                     node_type * includeGroupMRCA,
                                     const OttIdSet & attachedElsewhere,
                                     InterTreeBand<T> * itbp);
    bool anyPhantomNodesAtNode(const node_type *, const OttIdSet &) const ;
    bool anyIncludedAtNode(const node_type *, const OttIdSet &) const ;
    node_type * getMRCA(const OttIdSet &id);
    bool insertIntoBandNoDesUpdate(InterTreeBand<T> * itbp, RootedTreeNode<T> * connectedNode, long phantomID);
    RootedTreeNode<T> * resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par, const PhyloStatement & ps);
    void updateToReflectResolution(node_type *oldAnc,
                                   node_type * newAnc,
                                   const std::set<node_type *> & movedTips,
                                   const PhyloStatement & ps);
    
    friend class RootedForest<T, U>;
    friend class GreedyPhylogeneticForest<T, U>;
    FTree(const FTree &) = delete;
    FTree & operator=(const FTree &) = delete;
    // data members
    const std::size_t treeId; // key for this tree in forest - used for debugging
    node_type * root;
    //OttIdSet connectedIds;
    // from excludedNode to the nodes that it is excluded from...
    ExcludeConstraints<T> exclude;
    InterTreeBandBookkeeping<T> bands;
    //std::map<node_type *, std::list<PhyloStatementSource> > supportedBy; // only for non-roots
    RootedForest<T, U> & forest;
    std::map<long, node_type *> & ottIdToNodeMap;
};

template<typename T, typename U>
inline std::set<RootedTreeNode<T> *> FTree<T,U>::ottIdSetToNodeSet(const OttIdSet &ottIdSet) const {
    std::set<RootedTreeNode<T> *> ns;
    for (auto oid :ottIdSet) {
        ns.insert(ottIdToNodeMap.at(oid));
    }
    return ns;
}

template<typename T, typename U>
inline void FTree<T,U>::addExcludeStatement(long ottId,
                                            RootedTreeNode<T> * excludedFrom,
                                            const PhyloStatementSource &) {
    auto eNode = ottIdToNodeMap.at(ottId);
    exclude.addExcludeStatement(eNode, excludedFrom);
}

} // namespace otc
#endif
