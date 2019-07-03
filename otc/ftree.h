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
template<typename T, typename U> class GreedyBandedForest;
template<typename T, typename U> class RootedForest;
template<typename T, typename U> class FTree;

struct PhyloStatementSource {
    PhyloStatementSource(int treeInd, OttId groupInd)
        :source_tree_id(treeInd),
        clade_id(groupInd) {
    }
    const int source_tree_id;
    const OttId clade_id; // ID of the node in the tree
};

struct PhyloStatement {
    /*PhyloStatement(const OttIdSet &includes, const OttIdSet & other, bool otherIsExcludes)
        :include_group(includes),
        exclude_group(otherIsExcludes ? other : set_difference_as_set(other, includes)),
        leaf_set(otherIsExcludes ? set_union_as_set(other, includes): other) {
    }*/
    PhyloStatement(const OttIdSet &includes,
                   const OttIdSet &excludes,
                   const OttIdSet &mentioned, 
                   PhyloStatementSource pss)
        :include_group(includes),
        exclude_group(excludes),
        leaf_set(mentioned),
        provenance(pss) {
        debug_check();
    }
    void write_as_newick(std::ofstream & out) const;
    bool debug_check() const;
    bool is_trivial() const {
        return include_group.size() < 2 || include_group.size() == leaf_set.size();
    }
    const OttIdSet & include_group;
    const OttIdSet & exclude_group;
    const OttIdSet & leaf_set; // just the union of the include_group and exclude_group
    const PhyloStatementSource provenance;
};

// Bands connect nodes across different FTree instance when the nodes are required
//  to satisfy a grouping that has been added. The include_group of the ps
//  is a set of IDs that will be in the des_ids of each member of the band (and
//  the ancestor nodes of the members).
//  The nodes in the ps.exclude_group should be excluded at this node or an ancestor.
template<typename T>
class InterTreeBand {
    public:
    using node_type = RootedTreeNode<T>;
    using node_set = std::set<node_type *>;
    InterTreeBand(node_type * nd1,
                  const node_set & n1set,
                  const PhyloStatement & ps)
        :statement(ps) {
        add_node(nd1, n1set);
    }
    bool is_single_tree_band() const {
        return (node_to_phantom.size() == 1);
    }
    void insert_set(node_type * nd, const node_set & phantom) {
        for (auto p : phantom) {
            assert(p->has_ott_id());
        }
        node_to_phantom.at(nd).insert(begin(phantom), end(phantom));
    }
    void insert(node_type * nd, node_type * phantom) {
        assert(phantom->has_ott_id());
        node_to_phantom[nd].insert(phantom);
    }
    void remove_node(node_type * nd) {
        node_to_phantom.erase(nd);
    }
    void remove_from_set(node_type * nd, const node_set & nset) {
        node_set r = set_difference_as_set(node_to_phantom.at(nd), nset);
        node_to_phantom[nd] = r;
    }
    void add_node(node_type * nd, const node_set & nset) {
        for (auto p : nset) {
            assert(p->has_ott_id());
        }
        assert(nd != nullptr);
        const auto s = node_to_phantom.size();
        node_to_phantom[nd] = nset;
        assert((s + 1)== node_to_phantom.size());
    }
    bool band_point_has_any(const node_type * nd, const node_set & ns) const {
        return !(are_disjoint(node_to_phantom.at(nd), ns));
    }
    const node_set & get_phantom_nodes(const node_type * nd) const {
        auto npIt = node_to_phantom.find(nd);
        return (npIt == node_to_phantom.end() ? empty_set : npIt->second);
    }
    OttIdSet get_phantom_ids(const node_type * nd) const {
        const auto & ns = get_phantom_nodes(nd);
        OttIdSet r;
        for (auto np : ns) {
            r.insert(np->get_ott_id());
        }
        return r;
    }
    void reassign_attachment_node(node_type * oldAnc, node_type *newAnc) {
        assert(!contains(node_to_phantom, newAnc));
        node_to_phantom[newAnc] = node_to_phantom.at(oldAnc);
        node_to_phantom.erase(oldAnc);
    }
    bool is_the_set_of_phantom_nodes(node_type * nd, const node_set & t) const {
        return node_to_phantom.at(nd) == t;
    }
    bool is_a_banded_node_in_this(const node_type *q) const {
        return contains(node_to_phantom, q);
    }
    std::set<const node_type *> get_banded_nodes() const {
        std::set<const node_type *> r;
        for (const auto & m : node_to_phantom) {
            r.insert(m.first);
        }
        return r;
    }
    const std::map<const node_type *, node_set> & get_raw_map() const {
        return node_to_phantom;
    }
    private:
    std::map<const node_type *, node_set> node_to_phantom;
    const PhyloStatement & statement;
    const node_set empty_set = {};
};

template<typename T>
class ExcludeConstraints {
    public:
    using node_type = RootedTreeNode<T>;
    using node_pair = std::pair<const node_type *, const node_type *>;
    using cnode_set = std::set<const node_type *>;
    using node2many_map = std::map<const node_type *, cnode_set >;
    bool add_exclude_statement(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    bool is_excluded_from(const node_type * ndToCheck,
                        const node_type * potentialAttachment,
                        const std::map<OttId, node_type*> * ott_id_to_node_map) const;
    bool has_nodes_excluded_from_it(const node_type *n) const {
        return contains(by_node_with_constraints, n);
    }
    const cnode_set & get_nodes_excluded_from_node(const node_type * nd) const {
        auto bc = by_node_with_constraints.find(nd);
        if (bc == by_node_with_constraints.end()) {
            return empty_set;
        }
        return bc->second;
    }
    cnode_set steal_exclusions(node_type *nd) {
        auto bc = by_node_with_constraints.find(nd);
        if (bc == by_node_with_constraints.end()) {
            return cnode_set{};
        }
        cnode_set r = bc->second;
        for (auto n : r) {
            by_exclude_node[n].erase(nd);
        }
        return r;
    }
    private:
    void purge_exclude_raw(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    void ingest_exclude_raw(const node_type * nd2Exclude, const node_type * forbiddenAttach);
    node2many_map by_exclude_node;
    node2many_map by_node_with_constraints;
    const cnode_set empty_set = {};
};

template<typename T>
class InterTreeBandBookkeeping {
    public:
    using node_type = RootedTreeNode<T>;
    using band_type = InterTreeBand<T>;
    using band_set = std::set<band_type *>;
    OttIdSet get_phantom_ids(const node_type * nd) const {
        OttIdSet r;
        const band_set & bs =  get_bands_for_node(nd);
        for (const auto & b : bs) {
            const OttIdSet bo = b->get_phantom_ids(nd);
            r.insert(begin(bo), end(bo));
        }
        return r;
    }
    void reassign_attachment_node(band_type * b, node_type * oldAnc, node_type * newAnc, const PhyloStatement & ps);
    const band_set & get_bands_for_node(const node_type *n) const {
        const auto nit = node_to_band.find(n);
        if (nit == node_to_band.end()) {
            return empty_set;
        }
        return nit->second;
    }
    band_set steal_bands(const node_type *n) {
        const auto nit = node_to_band.find(n);
        if (nit == node_to_band.end()) {
            return empty_set;
        }
        band_set bs = nit->second;
        node_to_band.erase(n);
        for (auto bandPtr : bs) {
            band_to_node.erase(bandPtr);
        }
        return bs;
    }
    bool is_in_a_band(const node_type *n) const {
        return contains(node_to_band, n);
    }
    void _add_ref_to_band(band_type * band, node_type * nd) {
        assert(nd != nullptr);
        assert(band != nullptr);
        band_to_node[band] = nd;
        node_to_band[nd].insert(band);
        if (!contains(band->get_raw_map(), nd)) {
            const std::set<node_type *> emptyNodeSet;
            band->add_node(nd, emptyNodeSet);
        }
    }
    private:
    std::map<const band_type *, node_type *> band_to_node;
    std::map<const node_type *, band_set > node_to_band;
    const band_set empty_set = {};
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
          std::map<OttId, node_type *> & ottIdToNodeRef)
        :tree_id(treeID),
         root(nullptr),
         forest(theForest),
         ott_id_to_node_map(ottIdToNodeRef) {
    }
    // const methods:
    bool is_in_a_band(const node_type *n) const {
        return bands.is_in_a_band(n);
    }
    bool has_nodes_excluded_from_it(const node_type *n) const {
        return exclude.has_nodes_excluded_from_it(n);
    }
    const node_type * get_root() const {
        return root;
    }
    bool ott_id_is_excluded_from_root(OttId oid) const {
        return is_excluded_from_root(ott_id_to_node_map.at(oid));
    }
    bool is_excluded_from_root(const node_type *n) const {
        return exclude.is_excluded_from(n, root, &ott_id_to_node_map);
    }
    bool is_excluded_from(const node_type * ndToCheck, const node_type * potentialAttachment) const {
        return exclude.is_excluded_from(ndToCheck, potentialAttachment, &ott_id_to_node_map);
    }
    
    // OTT Ids of nodes on the graph only....
    const OttIdSet get_connected_ott_ids() const;
    // includes OTT Ids of nodes in includesConstraints
    const OttIdSet & get_included_ott_ids() {
        return get_root()->get_data().des_ids;
    }
    bool ott_id_is_connected(OttId ottId) const {
        return contains(get_connected_ott_ids(), ottId);
    }
    const ExcludeConstraints<T> & get_exclusions() const {
        return exclude;
    }
    const std::set<InterTreeBand<T> *> & get_bands_for_node(const node_type *n) const {
        return bands.get_bands_for_node(n);
    }
    // non-const
    node_type * get_root() {
        return root;
    }
    void _set_root(node_type *r) {
        root = r;
    }
    void mirror_phylo_statement(const PhyloStatement & ps);
#if defined(DO_DEBUG_CHECKS)
    void debug_invariants_check_ft() const;
    void debug_verify_des_ids_assuming_des(const OttIdSet &s, const RootedTreeNode<T> *nd) const;
#else
    void debug_invariants_check_ft() const{
    }
    void debug_verify_des_ids_assuming_des(const OttIdSet &, const RootedTreeNode<T> *) const {
    }
#endif
    std::set<RootedTreeNode<T> *> ott_id_set_to_node_set(const OttIdSet &ott_id_set) const;
    bool any_excluded_at_node(const node_type *, const OttIdSet &) const ;
    void create_deeper_root();
    // puts a node between nd and its parent and returns the new node
    node_type * create_deeper_node(node_type *nd);
    void steal_exclusion_statements(node_type * newPar,  node_type * srcNode, FTree<T, U>  & donorTree);
    OttIdSet steal_inclusion_statements(node_type * newPar, 
                                      node_type * srcNode,
                                      FTree<T, U>  & donorTree,
                                      InterTreeBand<T> * bandToSkip);
    void register_exclusion_statement_for_transferring_node(node_type * srcNode, FTree<T, U>  & donorTree);
    void register_inclusion_statement_for_transferring_node(node_type * srcNode, FTree<T, U>  & donorTree);
    private:
    void add_exclude_statement(OttId ottId, RootedTreeNode<T> *, const PhyloStatementSource &);
    void add_include_group_disjoint_phylo_statement(const PhyloStatement & ps) {
        add_phylo_statement_as_child_of_root(ps);
    }
    RootedTreeNode<T> * add_leaf_no_des_update(RootedTreeNode<T> * par, OttId ottId);
    void add_phylo_statement_as_child_of_root(const PhyloStatement &);
    // this is greedy, we should be building separate FTree instances in many cases....
    OttIdSet add_phylo_statement_at_node(const PhyloStatement & ps, 
                                     node_type * includeGroupMRCA,
                                     const OttIdSet & attachedElsewhere,
                                     InterTreeBand<T> * itbp);
    bool add_phantom_nodes_at_node(const node_type *, const OttIdSet &) const ;
    bool any_included_at_node(const node_type *, const OttIdSet &) const ;
    node_type * get_mrca(const OttIdSet &id);
    bool insert_into_band_no_des_update(InterTreeBand<T> * itbp, RootedTreeNode<T> * connectedNode, OttId phantomID);
    RootedTreeNode<T> * resolve_to_create_clade_of_included(RootedTreeNode<T> * par, const PhyloStatement & ps);
    void update_to_reflect_resolution(node_type *oldAnc,
                                   node_type * newAnc,
                                   const std::set<node_type *> & movedTips,
                                   const PhyloStatement & ps);
    
    friend class RootedForest<T, U>;
    friend class GreedyBandedForest<T, U>;
    FTree(const FTree &) = delete;
    FTree & operator=(const FTree &) = delete;
    // data members
    const std::size_t tree_id; // key for this tree in forest - used for debugging
    node_type * root;
    //OttIdSet connectedIds;
    // from excludedNode to the nodes that it is excluded from...
    ExcludeConstraints<T> exclude;
    InterTreeBandBookkeeping<T> bands;
    //std::map<node_type *, std::list<PhyloStatementSource> > supportedBy; // only for non-roots
    RootedForest<T, U> & forest;
    std::map<OttId, node_type *> & ott_id_to_node_map;
};

template<typename T, typename U>
inline std::set<RootedTreeNode<T> *> FTree<T, U>::ott_id_set_to_node_set(const OttIdSet &ott_id_set) const {
    std::set<RootedTreeNode<T> *> ns;
    for (auto oid :ott_id_set) {
        ns.insert(ott_id_to_node_map.at(oid));
    }
    return ns;
}

template<typename T, typename U>
inline void FTree<T, U>::add_exclude_statement(OttId ottId,
                                            RootedTreeNode<T> * excludedFrom,
                                            const PhyloStatementSource &) {
    auto eNode = ott_id_to_node_map.at(ottId);
    exclude.add_exclude_statement(eNode, excludedFrom);
}

} // namespace otc
#endif
