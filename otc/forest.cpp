#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/tree_iter.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/debug.h"

#include "otc/write_dot.h"
namespace otc {
template<typename T, typename U>
typename RootedForest<T, U>::tree_type &
RootedForest<T, U>::create_new_tree() {
    std::size_t i = next_tree_id++;
    auto r = trees.emplace(std::piecewise_construct,
                           std::forward_as_tuple(i),
                           std::forward_as_tuple(next_tree_id,
                                           *this,
                                           ott_id_to_node_map));
    assert(r.second); // must be a new Tree!
    return trees.at(i);
}
   
template<typename T, typename U>
RootedForest<T, U>::RootedForest(long rootOttId)
    :next_tree_id(0U),
    ott_id_to_node_map(node_src.get_data().ottIdToNode),
    rootID(rootOttId) {
}

template<typename T, typename U>
void RootedForest<T, U>::register_leaf(long ottId) {
    if (ottId == rootID) {
        return;
    }
    auto f = ott_id_to_node_map.find(ottId);
    if (f != ott_id_to_node_map.end()) {
        return;
    }
    create_leaf(nullptr, ottId, nullptr);
}

// TMP could be faster by storing node->tree lookup
template<typename T, typename U>
bool RootedForest<T, U>::is_attached(long ottId) const {
    auto f = ott_id_to_node_map.find(ottId);
    if (f == ott_id_to_node_map.end()) {
        return false;
    }
    node_type * n = f->second;
    assert(n != nullptr);
    return (n->getParent() != nullptr);
}

template<typename T, typename U>
bool RootedForest<T, U>::node_is_attached(RootedTreeNode<T> & n) const {
    return (n.getParent() != nullptr);
}

template<typename T, typename U>
void RootedForest<T, U>::attach_all_known_tips_as_new_tree() {
    tree_type & t = create_new_tree();
    t.root = create_node(nullptr, &t);
    for (auto & o2n : ott_id_to_node_map) {
        if (o2n.first == rootID) {
            assert(false);
            throw OTCError("root as tip");
        }
        auto nd = o2n.second;
        if (!node_is_attached(*nd)) {
            //t.connectedIds.insert(o2n.first);
            add_and_update_child(t.root, nd, t);
        }
    }
}

template<typename T, typename U>
void RootedForest<T, U>::_and_update_child(RootedTreeNode<T> *p,
                                                    RootedTreeNode<T> *c,
                                                    FTree<T, U> &tree) {
    register_tree_for_node(c, &tree);
    const auto & cd = c->get_data().desIds;
    if (cd.size() == 1 && !c->has_ott_id()) {
        if (*begin(cd) == c->get_ott_id()) {
            return;
        }
    }
    p->get_data().desIds.insert(begin(cd), end(cd));
}

template<typename T, typename U>
void RootedForest<T, U>::attach_all_detached_tips() {
    if (trees.empty()) {
        attach_all_known_tips_as_new_tree();
        return;
    }
    assert(trees.size() == 1);
    tree_type & t = trees.begin()->second;
    std::list<node_type *> excludedFromRoot;
    std::list<node_type *> attachableAtRoot;
    for (auto & o2n : ott_id_to_node_map) {
        if (o2n.first == rootID) {
            assert(false);
            throw OTCError("root as tip");
        }
        auto nd = o2n.second;
        if (!node_is_attached(*nd)) {
            if (t.is_excluded_from_root(nd)) {
                excludedFromRoot.push_back(nd);
            } else {
                attachableAtRoot.push_back(nd);
            }
        }
    }
    if (!excludedFromRoot.empty()) {
        auto nr = create_node(nullptr, &t);
        add_and_update_child(nr, t.root, t);
        t.root = nr;
        for (auto n : excludedFromRoot) {
            add_and_update_child(nr, n, t);
        }
    }
    // these could be attached one node more tipward (old root), but there is no basis for that.
    for (auto n : attachableAtRoot) {
        add_and_update_child(t.root, n, t);
    }
    
}

template<typename T, typename U>
bool RootedForest<T, U>::add_phylo_statement(const PhyloStatement &ps) {
    if (debugging_output_enabled) {
        dbWriteOttSet(" RootedForest::add_phylo_statement\nincGroup ", ps.include_group);
        dbWriteOttSet(" leaf_set", ps.leaf_set);
    }
    ps.debug_check();
    assert(ps.include_group.size() > 1);
    for (auto oid : ps.leaf_set) {
        register_leaf(oid);
    }
    if (ps.include_group == ps.leaf_set) {
        LOG(DEBUG) << "trivial group - exiting";
        return true;
    }
    const auto incompatRedundant = check_with_previously_added_statement(ps);
    if (incompatRedundant.first) {
        LOG(DEBUG) << "    hit incompat w/ prev added shortcircuit";
        return false;
    }
    if (incompatRedundant.second) { // we have added an identical group before
        LOG(DEBUG) << "    hit redundant w/ prev added shortcircuit";
        return true;
    }
    LOG(DEBUG) << "    checking compat w/ graph";
    if (add_phylo_statement_to_graph(ps)) {
        LOG(DEBUG) << "    compat w/ graph";
        added_splits_by_leaf_set[ps.leaf_set].insert(ps.include_group);
        return true;
    }
    LOG(DEBUG) << "    incompat w/ graph";
    return false;
}

template<typename T, typename U>
void RootedForest<T, U>::dump_accepted_phylo_statements(const char *fn) {
    std::ofstream phyloStatementOut;
    phyloStatementOut.open(fn);
    for (const auto & ps : novel_accepted_phylo_statements_in_order) {
        ps.write_as_newick(phyloStatementOut);
        phyloStatementOut << '\n';
    }
    phyloStatementOut.close();
}

void appendIncludeLeafSetAsNewick(const char *fn, const OttIdSet & inc, const OttIdSet & ls) {
    PhyloStatementSource pss(-1, -1);
    const OttIdSet exc = set_difference_as_set(ls, inc);
    std::ofstream phyloStatementOut;
    PhyloStatement ps{inc, exc, ls, pss};
    phyloStatementOut.open(fn, std::ios::app);
    ps.write_as_newick(phyloStatementOut);
    phyloStatementOut << '\n';
    phyloStatementOut.close();
}

template<typename T, typename U>
std::pair<bool, bool> RootedForest<T, U>::check_with_previously_added_statement(const PhyloStatement &ps) const {
    if (false && debugging_output_enabled) {
        dbWriteOttSet(" RootedForest::conflictsWithPreviouslyAddedStatement incGroup ", ps.include_group);
        dbWriteOttSet(" leaf_set", ps.leaf_set);
    }
    for (const auto sIt : added_splits_by_leaf_set) {
        const auto & prevAddedLeafSet = sIt.first;
        const auto relLeafSet = set_intersection_as_set(prevAddedLeafSet, ps.leaf_set);
        const bool exactLS = relLeafSet.size() == ps.leaf_set.size();
        if (relLeafSet.size() < 3) { // no conflict is possible if the intersection is so small that no phylostatements are made
            continue;
        }
        const auto relIncGroup = set_intersection_as_set(ps.include_group, relLeafSet);
        const auto & setPrevInc = sIt.second;
        for (const auto & prevInc : setPrevInc) {
            if (exactLS && prevInc == ps.include_group) {
                return std::pair<bool, bool>(false, true);
            }
            if (culled_and_complete_incompat_wrt_leaf_set(relIncGroup, prevInc, relLeafSet)) {
                return std::pair<bool, bool>(true, false);
            }
        }
    }
    return std::pair<bool, bool>(false, false);
}

template<typename T, typename U>
void consumeMapToList(std::map<T, std::list<U> > &m, std::list<U> & out) {
    for (auto & mIt : m) {
        for (auto el = begin(mIt.second) ; el != end(mIt.second); ++el) {
            out.push_back(*el);
        }
    }
}

template<typename T, typename U>
std::list<OverlapFTreePair<T, U> > RootedForest<T, U>::get_sorted_overlapping_trees(const OttIdSet &inc) {
    typedef OverlapFTreePair<T, U> MyOverlapFTreePair;
    std::map<std::size_t, std::list<MyOverlapFTreePair> > byOverlapSize;
    for (auto & tpIt : trees) {
        tree_type * ftree = &(tpIt.second);
        const OttIdSet & inTree = ftree->get_included_ott_ids();
        const OttIdSet inter = set_intersection_as_set(inTree, inc);
        if (!inter.empty()) {
            const auto k = inter.size();
            auto & tsList = byOverlapSize[k];
            tsList.push_back(MyOverlapFTreePair(inter, ftree));
        }
    }
    std::list<MyOverlapFTreePair> r;
    consumeMapToList(byOverlapSize, r);
    return r;
}

template<typename T, typename U>
void RootedForest<T, U>::add_ingroup_disjoint_phylo_statement_to_graph(const PhyloStatement &ps) {
    // this ingroup does not overlap with any ftree. find the FTree with the most overlap
    //  with the exclude_group...
    auto byExcCardinality = get_sorted_overlapping_trees(ps.exclude_group);
    if (byExcCardinality.empty()) {
        // none of the ingroup or outgroup are attached.
        // create a new FTree...
        // this can happen if the outgroup are mentioned in exclude statements (so the 
        //  areDisjoint returns false). But sense will add all of the leaves in the 
        //  include_group and exclude_group to this new tree, we don't need any new constraints
        //  so we can exit
        LOG(DEBUG) << "No exclude overlap either, using add_disjoint_tree";
        add_disjoint_tree(ps);
    } else {
        // TMP TOO GREEDY A CONNECTION - should only do this if all of the outgroup is connected...
        // we'll add the ingroup as a child of the root
        auto ftreeToAttach = byExcCardinality.begin()->second;
        LOG(DEBUG) << "Using add_include_group_disjoint_phylo_statement";
        ftreeToAttach->add_include_group_disjoint_phylo_statement(ps);
    }
    // no other trees had an include_group, so no need to add constraints....
}

template<typename T, typename U>
bool RootedForest<T, U>::is_in_a_band(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (tree.is_in_a_band(ncn)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool RootedForest<T, U>::has_nodes_excluded_from_it(const node_type * nd) const {
    node_type * ncn = const_cast<node_type *>(nd);
    for (const auto & tp : trees) {
        const auto & tree = tp.second;
        if (tree.has_nodes_excluded_from_it(ncn)) {
            return true;
        }
    }
    return false;
}

// should not modify the forest if returning false
template<typename T, typename U>
bool RootedForest<T, U>::check_can_add_ingroup_overlapping_phylo_statement_to_graph(
            const std::list<OverlapFTreePair<T, U> > & byIncCardinality,
            const PhyloStatement &ps,
            std::list<node_type * > & nonTrivMRCAs,
            OttIdSet & attachedElsewhere,
            std::vector<bool> & shouldResolveVec,
            std::vector<bool> & shouldCreateDeeperVec) const {
    for (const auto & incPair : byIncCardinality) {
        const auto & incGroupIntersection = incPair.first;
        attachedElsewhere.insert(incGroupIntersection.begin(), incGroupIntersection.end());
        tree_type * f = incPair.second;
        node_type * includeGroupA = nullptr;
        includeGroupA = f->get_mrca(incGroupIntersection);
        assert(includeGroupA != nullptr);
        assert(get_tree_for_node(includeGroupA) == f);
        if (includeGroupA->isTip()) {
            // this can happen if the overlap is one taxon.
            includeGroupA = includeGroupA->getParent();
            assert(includeGroupA != nullptr);
            assert(get_tree_for_node(includeGroupA) == f);
        }
        // If any of the ingroup are specifically excluded, then we have move deeper in the tree.
        // TMP this could be more efficient and avoid the while loop.
        while (f->any_excluded_at_node(includeGroupA, ps.include_group)) {
            if (f->any_included_at_node(includeGroupA, ps.exclude_group)) {
                return false;
            }
            if (f->add_phantom_nodes_at_node(includeGroupA, ps.include_group)) {
                return false;
            }
            includeGroupA = includeGroupA->getParent();
            if (includeGroupA == nullptr) {
                break;
            }
            assert(get_tree_for_node(includeGroupA) == f);
        }
        OttIdSet excInc;
        bool forceDeeperRoot = false;
        if (includeGroupA == nullptr) {
            includeGroupA = f->get_root();
            forceDeeperRoot = true;
            assert(get_tree_for_node(includeGroupA) == f);
        } else {
            excInc = set_intersection_as_set(includeGroupA->get_data().desIds, ps.exclude_group);
            if (debugging_output_enabled) {
                LOG(DEBUG) << "     add_phylo_statement_to_graph search for an ancestor of ..."; 
                dbWriteOttSet(" add_phylo_statement_to_graph search for an ancestor of:  ", incGroupIntersection);
                dbWriteOttSet(" wanted to avoid =  ", ps.exclude_group);
                dbWriteOttSet(" found a node with desIds:  ", includeGroupA->get_data().desIds);
                dbWriteOttSet(" which includes the excludegroup members:  ", excInc);
            }
            if (!can_be_resolved_to_display_inc_exc_group(includeGroupA, ps.include_group, excInc)) {
                return false; // the MRCA of the include_group had interdigitated members of the exclude_group
            }
        }
        shouldCreateDeeperVec.push_back(forceDeeperRoot);
        shouldResolveVec.push_back(!excInc.empty());
        nonTrivMRCAs.push_back(includeGroupA);
    }
    return true;
}

 
template<typename T, typename U>
InterTreeBand<T> * RootedForest<T, U>::_createNewBand(FTree<T, U> & ,
                                                 RootedTreeNode<T> &nd,
                                                 const PhyloStatement &ps) {
    std::set<RootedTreeNode<T> *> empty_set;
    all_bands.emplace_back(&nd, empty_set, ps);
    return &(*all_bands.rbegin());
}

template<typename T, typename U>
bool RootedForest<T, U>::add_ingroup_overlapping_phylo_statement_to_graph(const std::list<OverlapFTreePair<T, U> > & byIncCardinality,
                                                                   const PhyloStatement &ps) {
    std::list<node_type * > nonTrivMRCAs;
    OttIdSet attachedElsewhere;
    std::vector<bool> shouldResolveVec;
    std::vector<bool> shouldCreateDeeperVec;
    if (!check_can_add_ingroup_overlapping_phylo_statement_to_graph(byIncCardinality, ps, nonTrivMRCAs, attachedElsewhere, shouldResolveVec, shouldCreateDeeperVec)) {
        return false;
    }
    novel_accepted_phylo_statements_in_order.push_back(ps); //TMP DEBUGGING
    // all non trivial overlapping trees have approved this split...
    auto ntmIt = begin(nonTrivMRCAs);
    auto srIt = begin(shouldResolveVec);
    auto scdIt = begin(shouldCreateDeeperVec);
    unsigned i = 0;
    InterTreeBand<T> * itbp = nullptr;
    for (const auto & incPair : byIncCardinality) {
        LOG(DEBUG) << "   add_ingroup_overlapping_phylo_statement_to_graph mod for loop round " << ++i;
        debug_invariants_check();
        tree_type * f = incPair.second;
        assert(ntmIt != nonTrivMRCAs.end());
        node_type * includeGroupA = *ntmIt++;
        const bool add_node = *srIt++;
        const bool shouldCreateDeeperRoot = *scdIt;
        if (add_node) {
            includeGroupA = f->resolve_to_create_clade_of_included(includeGroupA, ps);
            assert(get_tree_for_node(includeGroupA) == f);
            LOG(DEBUG) << "   back from resolve_to_create_clade_of_included for loop round " << i;
            debug_invariants_check();
        } else if (shouldCreateDeeperRoot) {
            f->create_deeper_root();
            includeGroupA = f->get_root();
            assert(get_tree_for_node(includeGroupA) == f);
            LOG(DEBUG) << "   back from create_deeper_root for loop round " << i;
            debug_invariants_check();
        } else {
            assert(get_tree_for_node(includeGroupA) == f);
        }
        if (byIncCardinality.size() > 1 && itbp == nullptr) {
            itbp = _createNewBand(*f, *includeGroupA, ps);
        }
        auto connectedHere = f->add_phylo_statement_at_node(ps, includeGroupA, attachedElsewhere, itbp);
        if (!connectedHere.empty()) {
            attachedElsewhere.insert(begin(connectedHere), end(connectedHere));
        }
        dbWriteOttSet(" includeGroupA...desIds ", includeGroupA->get_data().desIds);
        LOG(DEBUG) << "   back from add_phylo_statement_at_node for loop round " << i;
        debug_invariants_check();
    }
    LOG(DEBUG) << "   add_ingroup_overlapping_phylo_statement_to_graph exit true " << i;
    return true;
}

template<typename T, typename U>
bool RootedForest<T, U>::add_phylo_statement_to_graph(const PhyloStatement &ps) {
    if (debugging_output_enabled) {
        dbWriteOttSet(" RootedForest::add_phylo_statement_to_graph incGroup ", ps.include_group);
        dbWriteOttSet(" leaf_set", ps.leaf_set);
    }
    if (ps.is_trivial()) {
        novel_accepted_phylo_statements_in_order.push_back(ps); //TMP DEBUGGING
        auto newOttIds = set_difference_as_set(ps.include_group, ott_id_set);
        for (auto noid : newOttIds) {
            add_detached_leaf(noid);
        }
        debug_invariants_check();
        return true;
    }
    if (areDisjoint(ps.leaf_set, ott_id_set)) {
        novel_accepted_phylo_statements_in_order.push_back(ps); //TMP DEBUGGING
        add_disjoint_tree(ps);
        debug_invariants_check();
        return true;
    }
    auto byIncCardinality = get_sorted_overlapping_trees(ps.include_group);
    LOG(DEBUG) << byIncCardinality.size() << " FTree instance referred to in byIncCardinality";
    if (byIncCardinality.empty()) {
        novel_accepted_phylo_statements_in_order.push_back(ps); //TMP DEBUGGING
        for (auto o : ps.include_group) {
            assert(!is_attached(o));
        }
        LOG(DEBUG) << "No intersection between include_group of an existing FTree.";
        add_ingroup_disjoint_phylo_statement_to_graph(ps);
        debug_invariants_check();
        return true;
    }
    auto & attachmentPair = *byIncCardinality.begin();
    if (attachmentPair.first.size() == 1) {
        // greedy approach is to add the rest of the ingroup as deep in this tree as possible.
        //  less greedy: make include/exclude statements at that node
        LOG(DEBUG) << "Missing opportunity to special case ingroups that only overlap with leaves of current trees.\n";
    }
    debug_invariants_check();
    auto rc = add_ingroup_overlapping_phylo_statement_to_graph(byIncCardinality, ps);
    debug_invariants_check();
    return rc;
}

template<typename T, typename U>
std::map<const FTree<T, U> *, const RootedTreeNode<T> *>
RootedForest<T, U>::get_tree_to_node_map_for_band(const InterTreeBand<T> & itb) const {
    std::set<const node_type *> bns = itb.get_banded_nodes();
    std::map<const FTree<T, U> *, const RootedTreeNode<T> *> r;
    for (auto bnp : bns) {
        auto t = get_tree_for_node(bnp);
        assert(t != nullptr);
        assert(!contains(r, t));
        r[t] = bnp;
    }
    return r;
}


template<typename T, typename U>
FTree<T, U> & RootedForest<T, U>::add_disjoint_tree(const PhyloStatement &ps) {
    tree_type & r = create_new_tree();
    r.mirror_phylo_statement(ps);
    return r;
}

template<typename T, typename U>
void RootedForest<T, U>::write_forest_dot_to_fn(const std::string &fn) const {
    LOG(DEBUG) << "     creating DOT file for forest: " << fn;
    std::ofstream outf(fn);
    writeDOTForest(outf, *this);
}

#if defined(DO_DEBUG_CHECKS)

template<typename T, typename U>
void RootedForest<T, U>::debug_invariants_check() const {
    std::map<const node_type *, const tree_type *> root2tree;
    for (const auto & t : trees) {
        t.second.debug_invariants_check_ft();
        auto r = t.second.get_root();
        assert(!contains(root2tree, r));
        root2tree[r] = &(t.second);
    }
    std::map<long, const tree_type *> ottId2Tree;
    std::map<const node_type *, const tree_type *> internal2Tree;
    std::set<const node_type *> detached;
    for (auto o2n : ott_id_to_node_map) {
        auto o = o2n.first;
        assert(o != rootID);
        auto n = o2n.second;
        if (n->isTip()) {
            assert(n->get_ott_id() == o);
            auto d = getDeepestAnc(n);
            assert(d != nullptr);
            if (d == n) {
                assert(!contains(detached, n));
                detached.insert(n);
            } else {
                ottId2Tree[o] = root2tree.at(d);
            }
        } else {
            assert(!n->has_ott_id());
            auto d = getDeepestAnc(n);
            assert(d != nullptr);
            assert(!contains(internal2Tree, n));
            if (d == n) {
                assert(contains(root2tree, n));
            }
            internal2Tree[n] = root2tree.at(d);
        }
        if (n->getParent() != nullptr) {
            assert(nullptr != get_tree_for_node(n));
        }
    }
    std::set<long> connectedIdSet;
    for (const auto & t : trees) {
        for (const auto o : t.second.get_connected_ott_ids()) {
            assert(!contains(connectedIdSet, o));
            //LOG(DEBUG) << "checkint ott" << o;
            assert(ottId2Tree.at(o) == &(t.second));
            connectedIdSet.insert(o);
        }
    }
    assert(connectedIdSet.size() == ottId2Tree.size()); 
}
#endif
template class RootedForest<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
