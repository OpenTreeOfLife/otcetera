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
    bool debugCheck() const;
    bool isTrivial() const {
        return includeGroup.size() < 2 || includeGroup.size() == leafSet.size();
    }
    const OttIdSet & includeGroup;
    const OttIdSet & excludeGroup;
    const OttIdSet & leafSet; // just the union of the includeGroup and excludeGroup
    const PhyloStatementSource provenance;
};

template<typename T, typename U>
class FTree {
    public:
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
    const node_type * getRoot() const {
        return root;
    }
    bool ottIdIsExcludedFromRoot(long oid) const {
        return isExcludedFromRoot(ottIdToNodeMap.at(oid));
    }
    bool isExcludedFromRoot(const node_type *) const;
    // OTT Ids of nodes on the graph only....
    const OttIdSet getConnectedOttIds() const;
    // includes OTT Ids of nodes in includesConstraints
    const OttIdSet & getIncludedOttIds() {
        return getRoot()->getData().desIds;
    }
    bool ottIdIsConnected(long ottId) const {
        return contains(getConnectedOttIds(), ottId);
    }
    const std::map<node_type *, std::list<GroupingConstraint> > & getExcluded2ConstraintMap() const {
        return excludesConstraints;
    }
    const std::map<node_type *, GroupingConstraint> & getIncluded2ConstraintMap() const {
        return includesConstraints;
    }
    // non-const
    node_type * getRoot() {
        return root;
    }
    void mirrorPhyloStatement(const PhyloStatement & ps);
    void addSubtree(RootedTreeNode<T> * subtreeRoot,
                    const NdToConstrainedAt & invIncConstr,
                    const NdToConstrainedAt & invExcConstr, 
                    const std::map<node_type *, GroupingConstraint> & incConstr,
                    const std::map<node_type *, std::list<GroupingConstraint> > & excConstr);
    void debugInvariantsCheck() const;
    void debugVerifyDesIdsAssumingDes(const OttIdSet &s, const RootedTreeNode<T> *nd) const;
    private:
    void addExcludeStatement(long ottId, RootedTreeNode<T> *, const PhyloStatementSource &);
    void addIncludeGroupDisjointPhyloStatement(const PhyloStatement & ps) {
        addPhyloStatementAsChildOfRoot(ps);
    }
    void addIncludeStatement(long ottId, RootedTreeNode<T> *, const PhyloStatementSource &);
    RootedTreeNode<T> * addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId);
    void addPhyloStatementAsChildOfRoot(const PhyloStatement &);
    // this is greedy, we should be building separate FTree instances in many cases....
    OttIdSet addPhyloStatementAtNode(const PhyloStatement & ps, 
                                     node_type * includeGroupMRCA,
                                     const OttIdSet & attachedElsewhere);
    bool anyExcludedAtNode(const node_type *, const OttIdSet &) const ;
    bool anyForceIncludedAtNode(const node_type *, const OttIdSet &) const ;
    bool anyIncludedAtNode(const node_type *, const OttIdSet &) const ;
    void createDeeperRoot();
    node_type * getMRCA(const OttIdSet &id);
    RootedTreeNode<T> * resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par, const OttIdSet & oids);
    
    friend class RootedForest<T, U>;
    friend class GreedyPhylogeneticForest<T, U>;
    FTree(const FTree &) = delete;
    FTree & operator=(const FTree &) = delete;
    // data members
    const std::size_t treeId; // key for this tree in forest - used for debugging
    node_type * root;
    //OttIdSet connectedIds;
    // from excludedNode to the nodes that it is excluded from...
    std::map<node_type *, std::list<GroupingConstraint> > excludesConstraints;
    std::map<node_type *, GroupingConstraint> includesConstraints;
    std::map<node_type *, OttIdSet> constrainedDesIds;
    std::map<node_type *, std::list<PhyloStatementSource> > supportedBy; // only for non-roots
    RootedForest<T, U> & forest;
    std::map<long, node_type *> & ottIdToNodeMap;
};

} // namespace otc
#endif
