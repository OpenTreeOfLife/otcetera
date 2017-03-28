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
#include "otc/ftree.h"
namespace otc {

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
As with NodeWithSplits, the node's data has a des_ids that lists the Ids of the nodes.
     However, to store the constraints, the FTree may also contain a constrainedDesIds
     for any of the nodes that are connected to that FTree.
Groupings make the phylogenetic statements that a group of nodes (the "include_group" of the
    grouping statement) share at least


We will build the forest by adding PhyloStatements
General postconditions for after we successfully add a PhyloStatement:
postcond #0: The number of attached nodes will never decrease
postcond #1: The taxon IDs will all have an Node object that represents them (registerd in
    the RootedForest's OttIdToNode map). However, note that these nodes are not necessarily
    attached.
postconditions about the include_group 
postcond #3: At least 1 member of the include_group will be attached.
postcond #4: For each FTree that contains member of the include_group, there will exist
    a node, the FTreeCA, that is a common ancestor of all members of the include_group which
    are connected to that tree. This FTreeCA node will not be the ancestor of any member of
    the exclude_group. If any of the nodes in the include_group are not connected to this FTree,
    their taxon ID will be recorded in the constrainedDesIds statement of the FTreeCA node.
    They will also be mentioned in an INCLUDES_NODE constraint for the FTreeCA node, if they
    are not already mentioned in and INCLUDES_NODE constraint for one of the descendants of
    the FTreeCA node. And every member of exclude_group will be mentioned in an EXCLUDES_NODE
    constraint for the FTreeCA node or one its ancestors.

XXX BOGUS IGNORE THIS one: postcond #X: If an FTreeCA node was new node created by the addition of the PhyloStatement,
    it will have out-degree=1. In this case, if its child is not a tip, then the edge between
    the FTreeCA and its child will be flagged as an possible-unsupported edge. This does not
    mean that the FTree removes all of the "supported by" statements for this node. Another
    PhyloStatement may have supported this grouping. However, that PhyloStatement would also
    support the branch from the FTreeCA to its parent. Later we'll post-process to remove
    unsupported nodes, and the support statements will be transferred one step closer to the root.
    These out-degree=1 
postconditions about the exclude_group
postcondition #5: If there is no intersection between the include_group of a PhyloStatement and the 
    forest, then the PhyloStatement will be added as a new FTree with one internal node that corresponds
    to the includeGroupMRCA and the root that is the a parent of all member of the exclude_group. 
    Otherwise, if a member of the excludedGroup was detached before addition of the PhyloStatement,
    then it will still be detached after the statement is added.

*/

/**
A PhlyoStatement is a rooted bipartition: a bipartition of the leaf_set with one of the subsets
    designated as the "include_group". 
The statement claims that all members of the include_group share at least  one common ancestor
    (their MRCA) which is not an ancestor of any member of the excluded group (the "exclude_group" 
    of the statement). Apart from the exclusion from that node, the statement doesn't say 
    anything about members of the exclude_group. In particular, there may be unmentioned taxa. 
    The members of the exclude_group are *not* necessarily more closely related to each other than
    they are to the unmentioned taxa. 
Indeed the unmentioned taxa could be placed anywhere on a tree without contradicting the PhyloStatement.
*/

template<typename T, typename U>
using OverlapFTreePair = std::pair<OttIdSet, FTree<T, U> *>;
// Constraint: each node w/ an OTT Id is a tip
//
template<typename T, typename U>
class RootedForest {
    public:
    using node_type = RootedTreeNode<T>;
    using tree_type = FTree<T, U>;
    RootedForest(OttId rootOttID);
    RootedForest(const RootedForest &) = delete;
    RootedForest & operator=(const RootedForest &) = delete;
    //accessors/queries:
    bool empty() const {
        return trees.empty();
    }
    const std::list<InterTreeBand<T> > & get_all_bands() const {
        return all_bands;
    }
    const std::map<OttId, node_type *> & get_ott_id_to_node_mapping() const {
        return ott_id_to_node_map;
    }
    const std::map<std::size_t, tree_type> & get_trees() const {
        return trees;
    }
    bool is_attached(OttId ottId) const;
    bool is_in_a_band(const node_type * ) const;
    bool has_nodes_excluded_from_it(const node_type * ) const;
    bool node_is_attached(RootedTreeNode<T> & n) const;
    std::pair<bool, bool> check_with_previously_added_statement(const PhyloStatement &ps) const;
    //modifiers
    bool add_phylo_statement(const PhyloStatement &);
    InterTreeBand<T> * _createNewBand(FTree<T, U> & ftree,
                                     RootedTreeNode<T> &nd,
                                     const PhyloStatement &ps);
    node_type * create_node(node_type * par, FTree<T, U> * ftree);
    node_type * create_leaf(node_type * par, const OttId & oid, FTree<T, U> * ftree);
    void register_tree_for_node(node_type * nd, FTree<T, U> * ftree) {
        node_to_tree[nd] = ftree;
    }
    void register_leaf(OttId ottId);
    void write_forest_dot_to_fn(const std::string &fn) const;
#if defined(DO_DEBUG_CHECKS)
    void debug_invariants_check() const;
#else
    void debug_invariants_check() const {
    }
#endif

    void dump_accepted_phylo_statements(const char *fn);
    const tree_type * get_tree_for_node(const node_type * nd) const {
        return node_to_tree.at(const_cast<node_type *>(nd));
    }
    void add_and_update_child(RootedTreeNode<T> *p, RootedTreeNode<T> *c, FTree<T, U> &tree) {
        p->add_child(c);
        _and_update_child(p, c, tree);
    }
    void _and_update_child(RootedTreeNode<T> *p, RootedTreeNode<T> *c, FTree<T, U> &tree);
    std::map<const tree_type *, const node_type *> get_tree_to_node_map_for_band(const InterTreeBand<T> &) const;
    private:
    std::list<OverlapFTreePair<T, U> > get_sorted_overlapping_trees(const OttIdSet &inc);
    node_type * add_detached_leaf(const OttId & ottId);
    tree_type & add_disjoint_tree(const PhyloStatement &);
    bool add_ingroup_overlapping_phylo_statement_to_graph(const std::list<OverlapFTreePair<T, U> > &, const PhyloStatement &);
    void add_ingroup_disjoint_phylo_statement_to_graph(const PhyloStatement &);
    bool add_phylo_statement_to_graph(const PhyloStatement &ps);
    bool check_can_add_ingroup_overlapping_phylo_statement_to_graph(
            const std::list<OverlapFTreePair<T, U> > & byIncCardinality,
            const PhyloStatement &ps,
            std::list<node_type * > & nonTrivMRCAs,
            OttIdSet & attachedElsewhere,
            std::vector<bool> & shouldResolveVec,
            std::vector<bool> & shouldCreateDeeperVec) const;
    tree_type & create_new_tree();
    protected:
    void attach_all_known_tips_as_new_tree();
    void attach_all_detached_tips();
    RootedTree<T, U> node_src; // not part of the forest, just the memory manager for the nodes
    std::map<std::size_t,  tree_type> trees;
    std::size_t next_tree_id;
    OttIdSet ott_id_set;
    std::map<OttId, node_type *> & ott_id_to_node_map; // alias to this data field in node_src for convenience
    std::map<node_type *, tree_type*> node_to_tree; 
    std::list<InterTreeBand<T> > all_bands;
    // added_splits_by_leaf_set
    // TMP, store every PhlyoStatement that we have accepted. This may become too memory inefficient
    // If we can do this, we can reject a new split based on conflict with another split that was
    //  accepted (which is less cryptic than just saying that we can't add it.) Some splits can't be
    //  added because of conflict with "emergent" properties of the forest. So checking
    //  added_splits_by_leaf_set is not sufficient to know if we can keep a split
    typedef std::set<OttIdSet> SetOfOTTIdSets;
    std::map<OttIdSet, SetOfOTTIdSets> added_splits_by_leaf_set;
    std::vector<PhyloStatement> novel_accepted_phylo_statements_in_order;// TMP debugging
    const OttId rootID;
};

template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T, U>::add_detached_leaf(const OttId & ottId) {
    return create_leaf(nullptr, ottId, nullptr);
}

template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T, U>::create_node(RootedTreeNode<T> * p, FTree<T, U> *ftree) {
    auto r = node_src.get_root();
    if (r == nullptr) {
        auto c = node_src.create_root();
        if (p != nullptr) {
            p->add_child(c);
        }
        register_tree_for_node(c, ftree);
        return c;
    }
    auto n = node_src.create_node(p);
    register_tree_for_node(n, ftree);
    return n;
}

// does NOT update anc des_ids!
template<typename T, typename U>
inline RootedTreeNode<T> * RootedForest<T, U>::create_leaf(RootedTreeNode<T> * p, const OttId & oid, FTree<T, U> *ftree) {
    assert(oid != rootID);
    auto n = create_node(p, ftree);
    assert(n->get_next_sib() == nullptr);
    n->set_ott_id(oid);
    ott_id_to_node_map[oid] = n;
    n->get_data().des_ids.insert(oid);
    ott_id_set.insert(oid);
    return n;
}

} // namespace otc
#endif
