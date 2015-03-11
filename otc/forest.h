#ifndef OTCETERA_FOREST_H
#define OTCETERA_FOREST_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/util.h"
namespace otc {
template<typename T, typename U> class RootedForest;
template<typename T, typename U> class FTree;

enum ConstraintType {
    INCLUDES_NODE,
    EXCLUDES_NODE
};



/**
NOTE: In this portion of the code, the only nodes with OTT Ids are tips, and ALL tips
    have OTT Ids. So, there is a bijection between a node and its OTT ID. This mapping might
    lead me to be a bit sloppy about writing "a set of leaves" when I mean "a set of OTT Ids"
    or vice versa.

RootedForest is a forest that was written to enable
     the GreedyAdditionOfGroupingsForest, but this base class
     provides the core methods of a forest as a graph.
Memory is managed by having one private tree that serves as a source of all nodes
     (deleting all nodes when it goes out of scope), and then having the 
     instances of the FTree class represent the connected components (the
     trees of the forest).
To enable the greedy addition of groupings, the FTree instances can store
     constraints on nodes which are not connected to them. This makes it possible
     to register the phylogenetic content of a grouping, but delay the merger of
     different trees in the forest until later (without later violating the constraints
     imposed by a grouping that we have already accepted)
An FTree is a tree in a rooted forest. The tree represented by an FTree is rooted. FTree
    instances only holds aliases to their nodes (the node life cycle is managed by the 
    RootedForest instance that owns the FTree)
Nodes referred to in the FTree can be "connected" or "nonmember".
"Connected" - When talking about an FTree, a connected node is the root of the FTree
     or a node that has a path that connects it to the root of the FTree. We will say 
     that this FTree "contains" a node that is a connected node for that FTree.
"Detached"  means not connected to any FTRee.
"Attached" - When talking about a node in the context of a forest, an "attached" node 
    means that the node is connected to 1 tree.
No node is connected with more than one FTree - if we add a grouping that would imply
     that to be the case, then the RootedForest has to engineer the merger of the 
     two FTree objects. (mathematically, trees are always connected and forests are graphs
     in when every connected component forms a tree)
"Nonmember" - When talking about an FTree, a node may not connected to this FTree, but 
    could be referred to in constraint statements of the aFTree. Such nodes are the
    nonmember nodes of the FTree. Each nonmember is either connected different FTree
    instances, or one of the "detached" nodes.
An FTree's constraint statement should not include any of the connected nodes for that
    tree. The statements can refer to "attached" nodes which are connected to a other
    FTrees, but the topological constraints on connected nodes that are relevant to an
    FTree should be encoded in the position of the node on the tree.
As with NodeWithSplits, the node's data has a desIds that lists the Ids of the nodes.
     However, to store the constraints, the FTree may also contain a constrainedDesIds
     for any of the nodes that are connected to that FTree.
Groupings make the phylogenetic statements that a group of nodes (the "includeGroup" of the
    grouping statement) share at least


We will build the forest by adding PhyloStatements
General postconditions for after we successfully add a PhyloStatement:
postcond #0: The number of attached nodes will never decrease
postcond #1: The taxon IDs will all have an Node object that represents them (registerd in
    the RootedForest's OttIdToNode map). However, note that these nodes are not necessarily
    attached.
postconditions about the includeGroup 
postcond #3: At least 1 member of the includeGroup will be attached.
postcond #4: For each FTree that contains member of the includeGroup, there will exist
    a node, the FTreeCA, that is a common ancestor of all members of the includeGroup which
    are connected to that tree. This FTreeCA node will not be the ancestor of any member of
    the excludeGroup. If any of the nodes in the includeGroup are not connected to this FTree,
    their taxon ID will be recorded in the constrainedDesIds statement of the FTreeCA node.
    They will also be mentioned in an INCLUDES_NODE constraint for the FTreeCA node, if they
    are not already mentioned in and INCLUDES_NODE constraint for one of the descendants of
    the FTreeCA node. And every member of excludeGroup will be mentioned in an EXCLUDES_NODE
    constraint for the FTreeCA node or one its ancestors.

XXX BOGUS IGNORE THIS one: postcond #X: If an FTreeCA node was new node created by the addition of the PhyloStatement,
    it will have out-degree=1. In this case, if its child is not a tip, then the edge between
    the FTreeCA and its child will be flagged as an possible-unsupported edge. This does not
    mean that the FTree removes all of the "supported by" statements for this node. Another
    PhyloStatement may have supported this grouping. However, that PhyloStatement would also
    support the branch from the FTreeCA to its parent. Later we'll post-process to remove
    unsupported nodes, and the support statements will be transferred one step closer to the root.
    These out-degree=1 
postconditions about the excludeGroup
postcondition #5: If there is no intersection between the includeGroup of a PhyloStatement and the 
    forest, then the PhyloStatement will be added as a new FTree with one internal node that corresponds
    to the includeGroupMRCA and the root that is the a parent of all member of the excludeGroup. 
    Otherwise, if a member of the excludedGroup was detached before addition of the PhyloStatement,
    then it will still be detached after the statement is added.

*/

/**
A PhlyoStatement is a rooted bipartition: a bipartition of the leafSet with one of the subsets
    designated as the "includeGroup". 
The statement claims that all members of the includeGroup share at least  one common ancestor
    (their MRCA) which is not an ancestor of any member of the excluded group (the "excludeGroup" 
    of the statement). Apart from the exclusion from that node, the statement doesn't say 
    anything about members of the excludeGroup. In particular, there may be unmentioned taxa. 
    The members of the excludeGroup are *not* necessarily more closely related to each other than
    they are to the unmentioned taxa. 
Indeed the unmentioned taxa could be placed anywhere on a tree without contradicting the PhyloStatement.
*/
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
    void mirrorPhyloStatement(const PhyloStatement & ps);
    node_type * getRoot() {
        return root;
    }
    // OTT Ids of nodes on the graph only....
    const OttIdSet & getConnectedOttIds() {
        return connectedIds;
    }
    // includes OTT Ids of nodes in includesConstraints
    const OttIdSet & getIncludedOttIds() {
        return getRoot()->getData().desIds;
    }
    FTree(std::size_t treeID,
          RootedForest<T, U> & theForest,
          std::map<long, node_type *> & ottIdToNodeRef)
        :treeId(treeID),
         root(nullptr),
         forest(theForest),
         ottIdToNode(ottIdToNodeRef) {
    }
    private:
    bool ottIdIsConnected(long ottId) const {
        return contains(connectedIds, ottId);
    }
    void addIncludeStatement(long ottId, RootedTreeNode<T> *, const PhyloStatementSource &);
    void addExcludeStatement(long ottId, RootedTreeNode<T> *, const PhyloStatementSource &);

    RootedTreeNode<T> * resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par, const OttIdSet & oids);
    RootedTreeNode<T> * addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId);
    bool anyExcludedAtNode(const node_type *, const OttIdSet &) const ;
    void addPhyloStatementAsChildOfRoot(const PhyloStatement &);
    // this is greedy, we should be building separate FTree instances in many cases....
    void addIncludeGroupDisjointPhyloStatement(const PhyloStatement & ps) {
        addPhyloStatementAsChildOfRoot(ps);
    }
    OttIdSet addPhyloStatementAtNode(const PhyloStatement & ps, 
                                     node_type * includeGroupMRCA,
                                     const OttIdSet & attachedElsewhere);
    node_type * getMRCA(const OttIdSet &id);
    
    friend class RootedForest<T, U>;
    FTree(const FTree &) = delete;
    FTree & operator=(const FTree &) = delete;
    // data members
    const std::size_t treeId; // key for this tree in forest - used for debugging
    node_type * root;
    OttIdSet connectedIds;
    // from excludedNode to the nodes that it is excluded from...
    std::map<node_type *, std::list<GroupingConstraint> > excludesConstraints;
    std::map<node_type *, std::list<GroupingConstraint> > includesConstraints;
    std::map<node_type *, OttIdSet> constrainedDesIds;
    std::map<node_type *, std::list<PhyloStatementSource> > supportedBy; // only for non-roots
    RootedForest<T, U> & forest;
    std::map<long, node_type *> & ottIdToNode;
};

// Constraint: each node w/ an OTT Id is a tip
//
template<typename T, typename U>
class RootedForest {
    public:
    using node_type = RootedTreeNode<T>;
    using tree_type = FTree<T,U>;
    RootedForest();
    RootedForest(const RootedForest &) = delete;
    RootedForest & operator=(const RootedForest &) = delete;
    bool empty() const {
        return trees.empty();
    }
    void registerLeaf(long ottId);
    bool isAttached(long ottId) const;
    bool nodeIsAttached(RootedTreeNode<T> & n) const;
    bool addPhyloStatement(const PhyloStatement &);

    // return <conflicts, redundant> pair based on cache of previously added groups.
    std::pair<bool, bool> checkWithPreviouslyAddedStatement(const PhyloStatement &ps) const;

    node_type * createNode(node_type * par);
    node_type * createLeaf(node_type * par, const OttId & oid);
    private:
    std::list<std::pair<OttIdSet, FTree<T, U> *> > getSortedOverlapping(const OttIdSet &inc);
    node_type * addDetachedLeaf(const OttId & ottId);
    tree_type & addDisjointTree(const PhyloStatement &);
    bool addPhyloStatementToGraph(const PhyloStatement &ps);
    tree_type & createNewTree();
    protected:
    void attachAllKnownTipsAsNewTree();
    RootedTree<T, U> nodeSrc; // not part of the forest, just the memory manager for the nodes
    std::map<std::size_t,  tree_type> trees;
    std::size_t nextTreeId;
    OttIdSet ottIdSet;
    std::map<OttId, node_type *> & ottIdToNode; // alias to this data field in nodeSrc for convenience
    
    // addedSplitsByLeafSet
    // TMP, store every PhlyoStatement that we have accepted. This may become too memory inefficient
    // If we can do this, we can reject a new split based on conflict with another split that was
    //  accepted (which is less cryptic than just saying that we can't add it.) Some splits can't be
    //  added because of conflict with "emergent" properties of the forest. So checking
    //  addedSplitsByLeafSet is not sufficient to know if we can keep a split
    typedef std::set<OttIdSet> SetOfOTTIdSets;
    std::map<OttIdSet, SetOfOTTIdSets> addedSplitsByLeafSet; 
};

template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T,U>::addDetachedLeaf(const OttId & ottId) {
    return createLeaf(nullptr, ottId);
}

template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T,U>::createNode(RootedTreeNode<T> * p) {
    auto r = nodeSrc.getRoot();
    if (r == nullptr) {
        auto c = nodeSrc.createRoot();
        if (p != nullptr) {
            p->addChild(c);
        }
        return c;
    }
    return nodeSrc.createNode(p);
}

// does NOT update anc desIds!
template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T,U>::createLeaf(RootedTreeNode<T> * p, const OttId & oid) {
    auto n = createNode(p);
    n->setOttId(oid);
    ottIdToNode[oid] = n;
    n->getData().desIds.insert(oid);
    ottIdSet.insert(oid);
    return n;
}

} // namespace otc
#endif
