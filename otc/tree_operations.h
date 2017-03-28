#ifndef OTCETERA_TREE_OPERATIONS_H
#define OTCETERA_TREE_OPERATIONS_H
// Functions that operate on trees - may include iteration
//  over trees (contrast w/tree_util.h)
// Depends on: tree.h tree_util.h tree_iter.h 
// Depended on by: tools
#include <vector>
#include <unordered_map>
#include "otc/otc_base_includes.h"
#include "otc/tree_iter.h"
#include "otc/error.h"
#include "otc/util.h"
#include "otc/debug.h"
namespace otc {
template<typename T, typename CONTAINER>
void show_children_in_set(std::ostream & out, T nd, const CONTAINER & ancestral);
template<typename T, typename CONTAINER>
int count_children_in_set(T nd, const CONTAINER & ancestral);
template<typename T, typename CONTAINER>
bool is_monotypic_in_set(T nd, const CONTAINER & ancestral);
template<typename T>
void collapse_split_and_del_node(T* nd);
template<typename T>
std::size_t n_internal_with_ott_id(const T& tree);
template<typename T>
std::size_t n_internal(const T& tree);
template<typename T>
std::size_t n_internal_out_degree_1(const T& tree);
template<typename T>
std::size_t n_internal_out_degree_many(const T & tree);
template<typename T>
std::size_t n_nodes(const T& tree);
template<typename T>
std::string newick(const T &t);
template <typename T>
T* bisect_branch_with_new_child(T* x);


template<typename T>
unsigned int count_polytomies(const T & tree);
template<typename T>
std::size_t check_for_unknown_taxa(std::ostream & err, const T & toCheck, const T & taxonomy);
template<typename T>
typename T::node_type * find_mrca_from_id_set(T & tree, const OttIdSet & idSet, OttId trigger);

template<typename T>
void prune_and_delete(T & tree, typename T::node_type *toDel);
template<typename T, typename U>
inline void cull_refs_to_node_from_data(RootedTree<T, U> & tree, RootedTreeNode<T> *toDel);

const OttIdSet& get_des_ids(RootedTreeNode<RTSplits> & nd);
template<typename T>
std::size_t prune_tips_without_ids(T & tree);
template <typename T>
std::vector<typename T::node_type*> all_nodes(T& tree);
template <typename N>
std::vector<N *> all_children(N * node);


//// impl
template <typename N>
inline std::vector<N *> all_children(N * node) {
    std::vector<N *> r;
    for (auto c : iter_child(*node)) {
        r.push_back(c);
    }
    return r;
}

// breaks the branch from x to its parent by allocating
//  a new node (which is returned) and assigning all of
//  the children of `x` to that node.
//  The returned node will be the only child of `x` and
//  all of the original children of x will be the children
//  of the returned node.
template<typename T>
T* bisect_branch_with_new_child(T* x) {
    assert(x);
    auto xc = new T(x);
    while (x->get_first_child()) {
        auto nd = x->get_first_child();
        nd->detach_this_node();
        xc->add_child(nd);
    }
    assert(not x->has_children());
    x->add_child(xc);
    return xc;
}

// breaks the branch from x to its parent by allocating
//  a new node (which is returned) and putting it in between
//  `x` and the parent of `x`
template<typename T>
T* bisect_branch_with_new_parent(T* x) {
    assert(x);
    auto xc = new T(x);
    while (x->get_first_child()) {
        auto nd = x->get_first_child();
        nd->detach_this_node();
        xc->add_child(nd);
    }
    assert(not x->has_children());
    x->add_child(xc);
    return xc;
}
template<typename T, typename CONTAINER>
inline void show_children_in_set(std::ostream & out, T nd, const CONTAINER& ancestral) {
    out << nd->get_name() << " children: ";
    for (auto c = nd->get_first_child(); c; c = c->get_next_sib()) {
        if (ancestral.count(c)) {
            out << "'" << c->get_name() << "' ";
        }
    }
    std::cerr << std::endl;
}

template<typename T, typename CONTAINER>
inline int count_children_in_set(T nd, const CONTAINER & ancestral) {
    int count = 0;
    for(auto c = nd->get_first_child(); c; c = c->get_next_sib()) {
        if (ancestral.count(c)) {
            count++;
        }
    }
    return count;
}

template<typename T, typename CONTAINER>
inline bool is_monotypic_in_set(T nd, const CONTAINER & ancestral) {
    return (count_children_in_set(nd, ancestral) == 1);
}

template<typename T>
void collapse_split_and_del_node(T* nd) {
    while (nd->get_first_child()) {
        auto nd2 = nd->get_first_child();
        nd2->detach_this_node();
        nd->add_sib_on_left(nd2);
    }
    nd->detach_this_node();
    delete nd;
}

template<typename T>
inline std::size_t n_internal_with_ott_id(const T& tree) {
    std::size_t count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->is_tip() and nd->has_ott_id()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline std::size_t n_internal(const T& tree) {
    std::size_t count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->is_tip()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline std::size_t n_internal_out_degree_1(const T& tree) {
    std::size_t count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->is_tip() and nd->is_outdegree_one_node()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline std::size_t n_internal_out_degree_many(const T & tree) {
    std::size_t count = 0;
    for(auto nd: iter_post_const(tree)) {
        if (not nd->is_tip() and not nd->is_outdegree_one_node()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline std::size_t n_nodes(const T& tree) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    std::size_t count = 0;
    for(auto nd: iter_post_const(tree)){
        count++;
    }
    return count;
}

template <typename T>
inline std::vector<typename T::node_type*> all_nodes(T& tree) {
    std::vector<typename T::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    return nodes;
}

template <typename T>
inline std::vector<typename T::node_type*> all_internal_nodes_post(T& tree) {
    std::vector<typename T::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        if (!nd->is_tip()) {
            nodes.push_back(nd);
        }
    }
    return nodes;
}


template<typename T>
unsigned int count_polytomies(const T & tree) {
    unsigned int n = 0U;
    for (auto node : iter_post_internal_const(tree)) {
        if (node->get_out_degree() > 2) {
            n += 1;
        }
    }
    return n;
}

template<typename T>
void fill_des_ids(T & tree) {
    // assumes OttId is set for each tip
    tree.get_data().des_id_sets_contain_internals = false;
    for (auto node : iter_post(tree)) {
        auto & des_ids = node->get_data().des_ids;
        if (node->is_tip()) {
            if (node->has_ott_id()) {
                    des_ids.insert(node->get_ott_id());
            }
        } else {
            for (auto child : iter_child(*node)) {
                auto & cDesIds = child->get_data().des_ids;
                des_ids.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}
template<typename T>
void clear_and_fill_des_ids(T & tree) {
    // assumes OttId is set for each tip
    tree.get_data().des_id_sets_contain_internals = false;
    for (auto node : iter_post(tree)) {
        auto & des_ids = node->get_data().des_ids;
        des_ids.clear();
        if (node->is_tip()) {
            des_ids.insert(node->get_ott_id());
        } else {
            for (auto child : iter_child(*node)) {
                auto & cDesIds = child->get_data().des_ids;
                des_ids.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}

template<typename T>
void fill_des_ids_including_internals(T & tree) {
    // assumes OttId is set for each tip
    tree.get_data().des_id_sets_contain_internals = true;
    for (auto node : iter_post(tree)) {
        auto & des_ids = node->get_data().des_ids;
        if (node->is_tip()) {
            des_ids.insert(node->get_ott_id());
        } else {
            if (node->has_ott_id()) {
                des_ids.insert(node->get_ott_id());
            }
            for (auto child : iter_child(*node)) {
                auto & cDesIds = child->get_data().des_ids;
                des_ids.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}

// fills in the depth data member for each node.
template <typename T>
void compute_depth(T& tree) {
    tree.get_root()->get_data().depth = 1;
    for (auto nd: iter_pre(tree)) {
        if (nd->get_parent()) {
            nd->get_data().depth = nd->get_parent()->get_data().depth + 1;
        }
    }
}

// uses ottID->node mapping, but not the split sets of the nodes
template<typename T>
typename T::node_type * find_mrca_from_id_set(T & tree, const OttIdSet & idSet, OttId trigger) {
    typedef typename T::node_type NT_t;
    auto ott_id_to_node = tree.get_data().ott_id_to_node;
    std::map<NT_t *, unsigned int> n2c;
    long shortestPathLen = -1;
    NT_t * shortestPathNode = nullptr;
    for (const auto & i : idSet) {
        const auto rIt = ott_id_to_node.find(i);
        if (rIt == ott_id_to_node.end()) {
            std::string em = "tip ";
            em += std::to_string(i);
            if (trigger >= 0 && trigger != std::numeric_limits<OttId>::max()) {
                em += " a descendant of ";
                em += std::to_string(trigger); 
            }
            em += " not found.";
            throw OTCError(em);
        }
        auto nd = rIt->second;
        long currPathLen = 0;
        while (nd != nullptr) {
            n2c[nd] += 1;
            currPathLen += 1;
            nd = nd->get_parent();
        }
        if (shortestPathLen < 0 || currPathLen < shortestPathLen) {
            shortestPathLen = currPathLen;
            shortestPathNode = rIt->second;
        }
    }
    auto nTips = idSet.size();
    auto cn = shortestPathNode;
    while (true) {
        if (n2c[cn] == nTips) {
            return cn;
        }
        cn = cn->get_parent();
        assert(cn != nullptr);
    }
    assert(false);
    throw OTCError("asserts disabled, but false");
}

template<typename T>
std::size_t check_for_unknown_taxa(std::ostream & err, const T & toCheck, const T & taxonomy) {
    const auto & taxOttIds = taxonomy.get_root()->get_data().des_ids;
    const auto & toCheckOttIds = toCheck.get_root()->get_data().des_ids;
    const auto extras = set_sym_difference_as_set(toCheckOttIds, taxOttIds);
    if (!extras.empty()) {
        err << "OTT Ids found in an input tree,  but not in the taxonomy:\n";
        write_ott_id_set(err, "  ", extras, "\n");
        return extras.size();
    }
    return 0U;
}

template<typename T>
inline typename T::node_type * add_child_for_ott_id(typename T::node_type & nd, OttId ottId, T & tree) {
    auto nn = tree.create_child(&nd);
    nn->set_ott_id(ottId);
    tree.get_data()->ott_id_to_node[ottId] = nn;
    return nn;
}

template<>
inline NodeWithSplits * add_child_for_ott_id<TreeMappedWithSplits>(NodeWithSplits & nd, OttId ottId, TreeMappedWithSplits & tree) {
    auto nn = tree.create_child(&nd);
    nn->set_ott_id(ottId);
    tree.get_data().ott_id_to_node[ottId] = nn;
    return nn;
}

inline const OttIdSet & get_des_ids(RootedTreeNode<RTSplits> & nd) {
    return nd.get_data().des_ids;
}

inline void fix_des_ids(RootedTreeNode<RTSplits> & nd, const OttIdSet & ls) {
    const auto toRemove = nd.get_data().des_ids;
    assert(!toRemove.empty());
    nd.get_data().des_ids = ls;
    for (auto anc : iter_anc(nd)) {
        assert(anc != nullptr);
        assert(!anc->get_data().des_ids.empty());
        for (auto tr : toRemove) {
            assert(contains(anc->get_data().des_ids, tr));
            anc->get_data().des_ids.erase(tr);
        }
        anc->get_data().des_ids.insert(begin(ls), end(ls));
    }
}

// throws an OTCError if a tip is mapped to a non-terminal taxon
template<typename N, typename T, typename M>
void require_tips_to_be_mapped_to_terminal_taxa(const RootedTree<N,T> & toExpand, const M& ott_id_to_node) {
    for (auto nd : iter_leaf_const(toExpand)) {
        assert(nd->is_tip());
        assert(nd->has_ott_id());
        auto ottId = nd->get_ott_id();

        if (not ott_id_to_node.count(ottId))
            throw OTCError()<<"OTT Id "<<ottId<<" / Label '"<<nd->get_name()<<"' not found in taxonomy!";

        auto taxNd = ott_id_to_node.at(ottId);
        if (not taxNd->is_tip())
            throw OTCError()<<"Tips must be mapped to terminal taxa: ott"<<ottId<<" found.";
    }
}

template<typename N, typename T>
void require_tips_to_be_mapped_to_terminal_taxa(const RootedTree<N,T> & toExpand, const RootedTree<N,T>& taxonomy) {
    std::map<OttId,const RootedTreeNode<N>*> ott_id_to_node;
    for(auto nd: iter_post_const(taxonomy))
        if (nd->has_ott_id())
            ott_id_to_node[nd->get_ott_id()] = nd;
    
    require_tips_to_be_mapped_to_terminal_taxa(toExpand, ott_id_to_node);
}

template<typename N>
void require_tips_to_be_mapped_to_terminal_taxa(const RootedTree<N,RTreeOttIDMapping<N>> & toExpand, const RootedTree<N,RTreeOttIDMapping<N>> & taxonomy) {
    const auto & taxData = taxonomy.get_data();
    require_tips_to_be_mapped_to_terminal_taxa(toExpand, taxData.ott_id_to_node);
}

// throws an OTCError if node has out-degree=1
template<typename T>
void require_nonredundant_tree(const T & toExpand) {
    std::map<typename T::node_type *, OttIdSet > replaceNodes;
    for (auto nd : iter_post_internal_const(toExpand)) {
        if (nd->is_outdegree_one_node()) {
            std::string msg;
            msg = "Node with outdegree=1 found";
            if (nd->has_ott_id()) {
                msg = msg + ":  ott" + std::to_string(nd->get_ott_id()) + ".";
            }
            throw OTCError(msg);
        }
    }
}


template<typename T>
std::set<const typename T::node_type *> expand_ott_internals_which_are_leaves(T & toExpand, const T & taxonomy) {
    const auto & taxData = taxonomy.get_data();
    std::map<typename T::node_type *, OttIdSet > replaceNodes;
    for (auto nd : iter_leaf(toExpand)) {
        assert(nd->is_tip());
        assert(nd->has_ott_id());
        auto ottId = nd->get_ott_id();
        auto taxNd = taxData.get_node_by_ott_id(ottId);
        if (not taxNd)
        throw OTCError()<<"OTT Id "<<ottId<<" / Label '"<<nd->get_name()<<"' not found in taxonomy!";
        if (!taxNd->is_tip()) {
            const auto & leaf_set = get_des_ids(*taxNd);
            replaceNodes[nd] = leaf_set;
        }
    }
    std::set<const typename T::node_type *> expanded;
    for (const auto & r : replaceNodes) {
        const auto & oldNode = r.first;
        const auto & ls = r.second;
        assert(ls.size() > 0);
        expanded.insert(oldNode);
        for (auto loid : ls) {
            add_child_for_ott_id<T>(*oldNode, loid, toExpand);
        }
        fix_des_ids(*oldNode, ls);
    }
    return expanded;
}

// Starts at the node in `fullTree` identified by `ottId`
// Walks back to the root, and adds `ottId` to value of the node->OttIdSet mapping `n2m`
// Creates new entries in that table, as needed.
// Multiple calls to this by tools like check-supertree to associate every node in a 
//      full tree with an OttIdSet which contains all of the descendants for a reduced
//      leafset.
template<typename T>
void mark_path_to_root(const T & fullTree,
                    OttId ottId,
                    std::map<const typename T::node_type *, OttIdSet > &n2m){
    auto startNd = fullTree.get_data().get_node_by_ott_id(ottId);
    if (startNd == nullptr) {
        std::string m = "OTT id not found ";
        m += std::to_string(ottId);
        throw OTCError(m);
    }
    assert(startNd != nullptr);
    n2m[startNd].insert(ottId);
    for (auto nd : iter_anc(*startNd)) {
        n2m[nd].insert(ottId);
    }
}


// find most recent anc of nd with out-degree > 1
template<typename T>
inline T * find_first_forking_anc(T * nd) {
    T * anc = nd->get_parent();
    if (anc == nullptr) {
        return nullptr;
    }
    while (anc->is_outdegree_one_node()) {
        anc = anc->get_parent();
        if (anc == nullptr) {
            return nullptr;
        }
    }
    return anc;
}

// find most recent anc of nd with out-degree > 1
template<typename T>
inline T * find_first_forking_self_or_des(T * nd) {
    assert(nd);
    T * des = nd;
    while (des->is_outdegree_one_node()) {
        des = des->get_first_child();
    }
    return des;
}


template<typename T, typename U>
inline bool multiple_children_in_map(const T & nd,
                                  const std::map<const T *, U> & markedMap,
                                  const T **first) {
    assert(first);
    bool foundFirst = false;
    *first = nullptr;
    for(auto c : iter_child_const(nd)) {
        if (markedMap.count(c) > 0) {
            if (foundFirst) {
                return true;
            }
            foundFirst = true;
            *first = c;
        }
    }
    return false;
}

template<typename T>
inline const T * find_next_significant_node(const T * node, const std::map<const T *, OttIdSet > & markedMap) {
    assert(node != nullptr);
    auto currNode = node;
    for (;;) {
        const T * sc;
        if (multiple_children_in_map(*currNode, markedMap, &sc)) {
            return currNode;
        }
        if (sc == 0L) {
            const char * msg = "Failing. Node found with ottIDs marked, but no children with ottIDs marked";
            throw OTCError(msg);
        }
        currNode = sc;
    }
}

template<typename T>
inline void write_pruned_subtree_newick_for_marked_nodes(std::ostream & out,
                                            const T & srcNd,
                                            const std::map<const T *, OttIdSet > & markedMap) {
    const auto nIt = markedMap.find(&srcNd);
    assert(nIt != markedMap .end());
    const auto & ottIDSet = nIt->second;
    if (ottIDSet.size() == 1) {
        out << "ott" << *ottIDSet.begin();
    } else {
        assert(ottIDSet.size() > 1);
        auto nsn = find_next_significant_node<T>(&srcNd, markedMap);
        out << '(';
        unsigned numcwritten = 0;
        for (auto child : iter_child_const(*nsn)) {
            if (markedMap.find(child) != markedMap.end()){
                if (numcwritten > 0) {
                    out << ',';
                }
                write_pruned_subtree_newick_for_marked_nodes(out, *child, markedMap);
                ++numcwritten;
            }
        }
        assert(numcwritten > 1);
        out << ')';
    }
}



template<typename T>
inline void describe_unnamed_node(const T & nd,
                                std::ostream & out,
                                unsigned int anc,
                                bool useNdNames,
                                bool addNewLine=true) {
    if (useNdNames && !nd.get_name().empty()) {
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before \"" << nd.get_name() << "\"";
        } else {
            out << "the node \"" << nd.get_name() << "\"";
        }
    }
    else if (nd.is_tip()) {
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before the leaf \"" << nd.get_name()  << "\"";
        } else {
            out << "the leaf \"" << nd.get_name()  << "\"";
        }
    } else if (nd.is_outdegree_one_node()) {
        describe_unnamed_node(*nd.get_first_child(), out, anc + 1, useNdNames);
        return;
    } else {
        const auto & left = find_leftmost_in_subtree(&nd)->get_name();
        const auto & right = find_rightmost_in_subtree(&nd)->get_name();
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right << "\"";
        } else {
            out <<  "MRCA of \"" << left << "\" and " << "\"" << right << "\"";
        }
    }
    if (addNewLine) {
        out << '\n';
    }
}

template<typename T, typename U>
inline void cull_refs_to_node_from_data(RootedTree<T, U> & , RootedTreeNode<T> *) {
}

template<>
inline void cull_refs_to_node_from_data(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & tree, 
                                   RootedTreeNode<RTNodeNoData> *nd) {
    assert(nd != nullptr);
    if (nd->has_ott_id()) {
        tree.get_data().ott_id_to_node.erase(nd->get_ott_id());
    }
}



template<>
inline void cull_refs_to_node_from_data(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree, 
                                   RootedTreeNode<RTSplits> *nd)  {
    assert(nd != nullptr);
    if (nd->has_ott_id()) {
        tree.get_data().ott_id_to_node.erase(nd->get_ott_id());
    }
    const auto & d = nd->get_data().des_ids;
    if (d.empty()) {
        return;
    }
    for (auto a : iter_anc(*nd)) {
        auto & ad = a->get_data().des_ids;
        for (auto el : d) {
            ad.erase(el);
        }
    }
}

// does NOT correct anc des_ids
template<typename T, typename U>
inline void del_ott_id_of_internal(RootedTree<T, U> & tree, RootedTreeNode<T> *nd) {
    if (nd->has_ott_id()) {
        tree.get_data().ott_id_to_node.erase(nd->get_ott_id());
        nd->del_ott_id();
    }
}

// does NOT correct anc des_ids
template<typename T, typename U>
inline void change_ott_id_of_internal(RootedTree<T, U> & tree, RootedTreeNode<T> *nd, OttId ottId) {
    del_ott_id_of_internal(tree, nd);
    if (ottId > 0) {
        nd->set_ott_id(ottId);
        tree.get_data().ott_id_to_node[ottId] = nd;
    }
}

template<typename T, typename U>
inline void collapse_node(RootedTree<T, U> & tree, RootedTreeNode<T> *nd) {
    if (nd->has_ott_id()) {
        del_ott_id_of_internal(tree, nd);
    }
    collapse_internal_into_par(nd, tree);
}


template<typename T>
inline void prune_and_delete(T & tree, typename T::node_type *toDel) {
    cull_refs_to_node_from_data<typename T::node_data_type, typename T::data_type>(tree, toDel);
    tree.prune_and_delete(toDel);
}

template <typename T, typename U>
void insert_ancestors_to_paraphyletic_set(T * nd, U & includedNodes) {
    for (auto anc : iter_anc(*nd)) {
        if (contains(includedNodes, anc)) {
            return;
        }
        includedNodes.insert(anc);
    }
}

//@TMP recursive until we have a pre-order subtree skipping iter.
template <typename T, typename U>
void insert_descendants_of_unincluded_subtrees(T * nd, U & includedNodes) {
    for (auto c : iter_child(*nd)) {
        if (!contains(includedNodes, c)) {
            includedNodes.insert(c);
            insert_descendants_of_unincluded_subtrees(c, includedNodes);
        }
    }
}


template<typename T>
inline void write_node_as_newick_label(std::ostream & out, const T *nd) {
    if (nd->is_tip()) { //TEMP tip just like internals, but at some point we may want the label-less internals to be unlabelled, regardless of OTT ID
        const auto & n = nd->get_name();
        if (n.empty()) {
            out << "ott" << nd->get_ott_id();
        } else {
            write_escaped_for_newick(out, nd->get_name());
        }
    } else if (!nd->get_name().empty()) {
        write_escaped_for_newick(out, nd->get_name());
    } else if (nd->has_ott_id()) {
        out << "ott" << nd->get_ott_id();
    }
}

template<typename T>
inline void write_closing_newick(std::ostream & out, const T *nd, const T * r) {
    out << ')';
    auto n = nd->get_parent();
    write_node_as_newick_label(out, n);
    if (n == r) {
        return;
    }
    while (n->get_next_sib() == nullptr) {
        out << ')';
        n = n->get_parent();
        assert(n != nullptr);
        write_node_as_newick_label(out, n);
        if (n == r) {
            return;
        }
    }
    out << ',';
}

template<typename T>
inline void write_newick(std::ostream & out, const T *nd) {
    assert(nd != nullptr);
    if (nd->is_tip()) {
        write_node_as_newick_label(out, nd);
    } else {
        for (auto n : iter_pre_n_const(nd)) {
            if (n->is_tip()) {
                write_node_as_newick_label(out, n);
                if (n->get_next_sib() == nullptr) {
                    write_closing_newick<T>(out, n, nd);
                } else {
                    out << ',';
                }
            } else {
                out << '(';
            }
        }
    }
}

template<typename T, typename Y>
inline void write_newick_generic(std::ostream & out,
                               T nd,
                               Y & nodeNamer,
                               long height_limit) {
    assert(nd != nullptr);
    if (!(nd->is_tip()) && height_limit != 0) {
        out << '(';
        bool first = true;
        const long nhl = height_limit - 1;
        for (auto c : iter_child_const(*nd)) {
            if (first) {
                first = false;
            } else {
                out << ',';
            }
            write_newick_generic<T, Y>(out, c, nodeNamer, nhl);
        }
        out << ')';
    }
    write_escaped_for_newick(out, nodeNamer(nd));
}


template<typename T>
inline void db_write_newick(const T *nd) {
    if (!debugging_output_enabled) {
        return;
    }
    write_newick(std::cerr, nd);
    std::cerr << std::endl;
}


template<typename T>
inline bool is_effectively_a_tip(const T *nd, std::function<bool(const T &)> subtree_filter) {
    if (nd->is_tip()) {
        return true;
    }
    for (auto child : iter_child_const(*nd)) {
        if (subtree_filter(*child)) {
            return false;
        }
    }
    return true;
}


template<typename T>
inline bool is_effectively_last_sib(const T *nd, std::function<bool(const T &)> subtree_filter) {
    nd = nd->get_next_sib();
    while (nd != nullptr) {
        if (subtree_filter(*nd)) {
            return false;
        }
        nd = nd->get_next_sib();
    }
    return true;
}

template<typename T>
inline void write_closing_newick_filtered(std::ostream & out, const T *nd, const T * r, std::function<bool(const T &)> subtree_filter) {
    assert(nd != nullptr);
    out << ')';
    auto n = nd->get_parent();
    write_node_as_newick_label(out, n);
    if (n == r) {
        return;
    }
    while (is_effectively_last_sib(n, subtree_filter)) {
        if (n == r) {
            return;
        }
        out << ')';
        n = n->get_parent();
        assert(n != nullptr);
        write_node_as_newick_label(out, n);
    }
    out << ',';
}

template<typename T>
inline void writeNewickFiltered(std::ostream & out, const T *nd, std::function<bool(const T &)> subtree_filter) {
    assert(nd != nullptr);
    if (nd->is_tip()) {
        if (!subtree_filter(*nd)) {
            return;
        }
        write_node_as_newick_label(out, nd);
    } else {
        for (auto n : iter_pre_filter_n_const(nd, subtree_filter)) {
            if (is_effectively_a_tip(n, subtree_filter)) {
                write_node_as_newick_label(out, n);
                if (is_effectively_last_sib(n, subtree_filter)) {
                    write_closing_newick_filtered<T>(out, n, nd, subtree_filter);
                } else {
                    out << ',';
                }
            } else {
                out << '(';
            }
        }
    }
}

template<typename T>
inline void write_tree_as_newick(std::ostream & out, const T &tree) {
    write_newick<typename T::node_type>(out, tree.get_root());
    out << ';';
}

template<typename T>
inline T * search_anc_for_mrca_of_des_ids(T * nd, const OttIdSet & idSet) {
    assert(nd != nullptr);
    if (is_subset(idSet, nd->get_data().des_ids)) {
        return nd;
    }
    for (auto n : iter_anc(*nd)) {
        if (is_subset(idSet, n->get_data().des_ids)) {
            return n;
        }
    }
    return nullptr;
}

template<typename T>
inline const typename T::node_type * find_node_with_matching_des_ids(const T & tree, const OttIdSet & idSet) {
    assert(!idSet.empty());
    OttId firstId = *begin(idSet);
    auto nd = tree.get_data().ott_id_to_node.at(firstId);
    assert(nd != nullptr);
    if (nd->get_data().des_ids == idSet) {
        return nd;
    }
    for (auto n : iter_anc_const(*nd)) {
        const auto & ndi = n->get_data().des_ids;
        if (ndi == idSet) {
            return n;
        }
        if (ndi.size() > idSet.size()) {
            return nullptr;
        }
    }
    return nullptr;
}

template<typename T>
inline const typename T::node_type * find_mrca_using_des_ids(const T & tree, const OttIdSet & idSet) {
    if (idSet.empty()) {
        assert(false);
        throw OTCError("asserts disabled but false");
    }
    const OttId lowestID = *idSet.begin();
    const typename T::node_type * aTip = tree.get_data().get_node_by_ott_id(lowestID);
    if (aTip == nullptr) {
        return nullptr;
    }
    return search_anc_for_mrca_of_des_ids(aTip, idSet);
}

template<typename T>
inline void erase_mappings_to_node(typename T::node_type * nd, T & tree) {
    if (nd->has_ott_id()) {
        tree.get_data().ott_id_to_node.erase(nd->get_ott_id());
    }
    tree.mark_as_detached(nd);
}

template<typename N, typename T>
inline void replace_mappings_to_node_with_alias(typename RootedTree<N,T>::node_type * , // nd
                                           typename RootedTree<N,T>::node_type * , //alias
                                           RootedTree<N,T> & ) {
}


template<typename N>
inline void replace_mappings_to_node_with_alias(typename RootedTree<N,RTreeOttIDMapping<N>>::node_type * nd,
                                           typename RootedTree<N,RTreeOttIDMapping<N>>::node_type * alias,
                                           RootedTree<N,RTreeOttIDMapping<N>>& tree) {
    if (nd->has_ott_id()) {
        tree.get_data().ott_id_to_detached_node[nd->get_ott_id()] = nd;
        tree.get_data().ott_id_to_node[nd->get_ott_id()] = alias;
        auto & iaf = tree.get_data().is_alias_for;
        iaf[alias].insert(nd->get_ott_id());
        auto na = iaf.find(nd); // if nd was an alias for other nodes, update them too...
        if (na != iaf.end()) {
            for (auto aoi : na->second) {
                iaf[alias].insert(aoi);
                tree.get_data().ott_id_to_node[aoi] = alias;
            }
        }
        iaf.erase(nd);
    }
    tree.mark_as_detached(nd);
}

// 
// detaches child from the tree (and lets it dangle)
// assigns its children to nd (which should be the parent of child)
// note in the only usage: nd is monotypic, but child is NOT.
// this suppresses monotypy by getting rid of child in case nd has a
//  taxonomic assignment that we want to preserve in the des_ids
template<typename T>
inline void suppress_monotypy_by_stealing_grandchildren(typename T::node_type * nd,
                                                    typename T::node_type * monotypic_child,
                                                    T & tree) {
    assert(nd);
    assert(monotypic_child);
    assert(monotypic_child->get_parent() == nd);
    assert(nd->get_first_child() == monotypic_child);
    assert(monotypic_child->get_next_sib() == nullptr);
    while(monotypic_child->has_children())
    {
        auto x = monotypic_child->get_first_child();
        x->detach_this_node();
        nd->add_child(x);
    }
    monotypic_child->detach_this_node();

    replace_mappings_to_node_with_alias(monotypic_child, nd, tree);
    //    delete monotypic_child;
}

template<typename T>
inline void suppress_monotypy_by_claiming_grandparent_as_par(typename T::node_type * monotypic_nd,
                                                       typename T::node_type * child,
                                                       T & tree) {
    assert(monotypic_nd);
    assert(child);
    assert(child->get_parent() == monotypic_nd);
    assert(monotypic_nd->get_first_child() == child);
    assert(child->get_next_sib() == nullptr);
    auto gp = monotypic_nd->get_parent();
    assert(gp);
    monotypic_nd->detach_this_node();
    child->detach_this_node();
    gp->add_child(child);
    assert(not monotypic_nd->has_children());
    replace_mappings_to_node_with_alias(monotypic_nd, child, tree);
    // delete monotypic_nd
}

template <typename T>
inline void del_monotypic_node(typename T::node_type* nd, T& tree) {
    assert(nd->is_outdegree_one_node());
    auto child = nd->get_first_child();
    child->detach_this_node();

    if (nd->get_parent())
    {
        nd->add_sib_on_right(child);
        nd->detach_this_node();
    }
    else
        tree._set_root(child);

    delete nd;
 }

template <typename T>
void suppress_monotypic_fast(T& tree)
{
    std::vector<typename T::node_type*> remove;
    for(auto nd:iter_pre(tree))
        if (nd->is_outdegree_one_node())
            remove.push_back(nd);

    for(auto nd: remove)
        del_monotypic_node(nd, tree);
}
 
// nd must be unnnamed (or the aliases would not be fixed because we can't call replace_mappings_to_node_with_alias)
template<typename T>
inline void collapse_internal_into_par(typename T::node_type * nd,
                                    T & tree) {
    assert(nd);
    assert(not nd->has_ott_id());
    assert(nd->has_children());

    if (nd->get_out_degree() == 1) {
        suppress_monotypy_by_claiming_grandparent_as_par(nd, nd->get_first_child(), tree);
        return;
    }
    
    while(nd->has_children())
    {
        auto child = nd->get_first_child();
        child->detach_this_node();
        nd->add_sib_on_left(child);
    }
    nd->detach_this_node();
}

// in some context we need to preserve the deepest because the des_ids may need 
// the internal nodes IDs (even if the IDs correspond to monotypic taxa)
// returns set of nodes deleted. Their get_parent() calls will tell the caller
//  their former parent (before culling). But note that you may have to 
//  call it muliple times to find the first ancestor that has not been deleted.
template<typename T>
inline std::set<typename T::node_type *> suppress_monotypic_taxa_preserve_deepest_dangle(
                    T & tree,
                    bool reset_ott_id) {
    std::set<typename T::node_type *> monotypic;
    for (auto nd : iter_node_internal(tree)) {
        if (nd->is_outdegree_one_node()) {
            monotypic.insert(nd);
        }
    }
    std::set<typename T::node_type *> removed;
    auto toProcess = monotypic;
    while (!toProcess.empty()) {
        for (auto toPIt = toProcess.begin(); toPIt != toProcess.end();) {
            auto tpn = *toPIt;
            auto child = tpn->get_first_child();
            if (!contains(toProcess, child)) {
                suppress_monotypy_by_stealing_grandchildren<T>(tpn, child, tree);
                toProcess.erase(toPIt++);
                if (reset_ott_id && child->has_ott_id()) {
                    tpn->set_ott_id(child->get_ott_id());
                }
                removed.insert(child);
            } else {
                ++toPIt;
            }
        }
    }
    return removed;
}

// in some context we need to preserve the deepest because the des_ids may need 
// the internal nodes IDs (even if the IDs correspond to monotypic taxa)
// returns set of nodes deleted. Their get_parent() calls will tell the caller
//  their former parent (before culling). But note that you may have to 
//  call it muliple times to find the first ancestor that has not been deleted.
template<typename T>
inline std::set<typename T::node_type *> suppress_monotypic_taxa_preserve_shallow_dangle(
                    T & tree) {
    std::set<typename T::node_type *> monotypic;
    for (auto nd : iter_node_internal(tree)) {
        //assert(!nd->has_ott_id()); // not a good place for this assert... true in how we use this fund
        if (nd->is_outdegree_one_node()) {
            monotypic.insert(nd);
        }
    }
    std::set<typename T::node_type *> removed;
    auto toProcess = monotypic;
    while (!toProcess.empty()) {
        for (auto toPIt = toProcess.begin(); toPIt != toProcess.end();) {
            check_tree_invariants(tree);
            auto tpn = *toPIt;
            auto par = tpn->get_parent();
            if (!contains(toProcess, par)) {
                if (tree.is_detached(tpn)) {
                    removed.insert(tpn);
                } else if (tpn->is_tip()) {
                    tree.prune_and_delete(tpn);
                    removed.insert(tpn);
                } else {
                    auto onlyChild = tpn->get_first_child();
                    if (par) {
                        suppress_monotypy_by_claiming_grandparent_as_par<T>(tpn, onlyChild, tree);
                        removed.insert(tpn);
                    } else {
                        replace_mappings_to_node_with_alias(onlyChild, tpn, tree);
                        if (onlyChild->is_tip()) {
                            tree.prune_and_delete(onlyChild);
                        } else {
                            onlyChild->del_ott_id();
                            collapse_internal_into_par(onlyChild, tree);
                        }
                        removed.insert(onlyChild);
                    }
                }
                toProcess.erase(toPIt++);
            } else {
                ++toPIt;
            }
        }
    }
    return removed;
}


template<typename T>
inline bool is_anc_des_pair(const T * nd1, const T *nd2) {
    if (nd1 == nd2) {
        return false;
    }
    for (auto a : iter_anc_const(*nd2)) {
        if (a == nd1) {
            return true;
        }
    }
    return false;
}

// returns true if the nodes are identical or one is an ancestor of the other...
template<typename T>
inline bool are_linearly_related(const T * nd1, const T *nd2) {
    return nd1 == nd2 || is_anc_des_pair(nd1, nd2) || is_anc_des_pair(nd2, nd1);
}

template<typename T>
inline OttIdSet get_ott_id_set_for_leaves(const T &tree) {
    if (!tree.get_data().des_id_sets_contain_internals) {
        return tree.get_root()->get_data().des_ids;
    }
    OttIdSet inducingIds;
    for (auto nd : iter_leaf_const(tree)) {
        inducingIds.insert(nd->get_ott_id());
    }
    return inducingIds;
}

template<typename T, typename U>
void get_induced_informative_groupings(const T & tree1, std::set<OttIdSet > & inducedSplits, const U & tree2) {
    const auto inducingIds = get_ott_id_set_for_leaves(tree2);
    auto mrca = find_mrca_using_des_ids(tree1, inducingIds);
    std::function<bool(const typename T::node_type &)> sf = [inducingIds](const typename T::node_type &nd){
        return have_intersection(inducingIds, nd.get_data().des_ids);
    };
    for (auto n : iter_pre_filter_n_const(mrca, sf)) {
        if (n == mrca) {
            continue;
        }
        auto inducedDesIds = n->get_data().des_ids;
        const auto x = intersection_of_sets(inducingIds, inducedDesIds);
        if (x.size() > 1) {
            inducedSplits.insert(std::move(x));
        }
    }
}

template<typename T>
void get_informative_groupings(const T & tree2,
                             std::set<OttIdSet > & tree2Splits) {
    auto t2r = tree2.get_root();
    for (auto n : iter_pre_internal_const(tree2)) {
        if (n == t2r) {
            continue;
        }
        const auto & x = n->get_data().des_ids;
        if (x.size() > 2) {
            tree2Splits.insert(x);
        }
    }
}

template<typename T, typename U>
void induced_clade_sets(const T & tree1,
                      const U & tree2,
                      std::set<OttIdSet > & inducedSplits,
                      std::set<OttIdSet > & tree2Splits,
                      bool firstIsSuperset) {
    assert(firstIsSuperset);
    get_induced_informative_groupings(tree1, inducedSplits, tree2);
    if (firstIsSuperset) {
        get_informative_groupings(tree2, tree2Splits);
    } else {
        throw OTCError("Why on earth did you compile without asserts?");
    }
}

template<typename T, typename U>
unsigned long induced_rf_dist(const T & tree1, const U & tree2, bool firstIsSuperset) {
    std::set<OttIdSet > inducedSplits;
    std::set<OttIdSet > tree2Splits;
    induced_clade_sets(tree1, tree2, inducedSplits, tree2Splits, firstIsSuperset);
    return size_of_symmetric_difference(tree2Splits, inducedSplits);
}

template<typename T, typename U>
std::size_t num_induced_splits_missing_in_second(const T & tree1,
                                              const U & tree2,
                                              bool firstIsSuperset) {
    std::set<OttIdSet > inducedSplits;
    std::set<OttIdSet > tree2Splits;
    induced_clade_sets(tree1, tree2, inducedSplits, tree2Splits, firstIsSuperset);
    std::size_t nm = 0;
    for (auto ics : inducedSplits) {
        if (!contains(tree2Splits, ics)) {
            nm += 1;
        }
    }
    return nm;
}

template<typename U, typename T, typename V>
inline void copy_tree_structure(const std::map<U *, U *> & nd2par,
                       const std::map<U *, OttId> & nd2id,
                       RootedTree<T, V> & toWrite) {
    std::set<RootedTreeNode<T> *> withParents;
    std::map<U *, RootedTreeNode<T> *> other2new;
    for (auto c : nd2par) {
        auto otherChild = c.first;
        auto otherParent = c.second;
        RootedTreeNode<T> * nParent;
        auto npIt = other2new.find(otherParent);
        if (npIt == other2new.end()) {
            nParent = toWrite.create_node(nullptr);
            other2new[otherParent] = nParent;
            if (otherParent->has_ott_id())
                nParent->set_ott_id(otherParent->get_ott_id());
        } else {
            nParent = npIt->second;
        }
        RootedTreeNode<T> * nChild;
        auto ncIt = other2new.find(otherChild);
        if (ncIt == other2new.end()) {
            nChild = toWrite.create_node(nParent);
            other2new[otherChild] = nChild;
            if (otherChild->has_ott_id())
                nChild->set_ott_id(otherChild->get_ott_id());
        } else {
            nChild = ncIt->second;
            nParent->add_child(nChild);
        }
        auto idIt = nd2id.find(otherChild);
        if (idIt != nd2id.end()) {
            nChild->set_ott_id(idIt->second);
        }
        idIt = nd2id.find(otherParent);
        if (idIt != nd2id.end()) {
            nParent->set_ott_id(idIt->second);
        }
        withParents.insert(nChild);
    }
    for (auto nn : withParents) {
        if ((nn->is_tip()) && !nn->has_ott_id()) {
            for (auto o2n : other2new) {
                if (o2n.second == nn) {
                    LOG(ERROR) << "tip without label in copy_tree_structure";
                    throw OTCError("tip without label");
                }
            }
        }
    }
    assert(other2new.size() == 1 + withParents.size()); // only the root should be parentless
    for (auto o2n : other2new) {
        if (!contains(withParents, o2n.second)) {
            toWrite._set_root(o2n.second);
            return;
        }
    }
}

template<typename T>
inline std::size_t prune_tips_without_ids(T & tree) {
    std::size_t r = 0;
    std::set<typename T::node_type *> toPrune;
    std::set<typename T::node_type *> parToCheck;
    for (auto nd : iter_node(tree)) {
        if (nd->is_tip() && !nd->has_ott_id()) {
            toPrune.insert(nd);
            parToCheck.insert(nd->get_parent());
        }
    }
    while (!toPrune.empty()) {
        r += toPrune.size();
        for (auto nd : toPrune) {
            prune_and_delete(tree, nd);
        }
        toPrune.clear();
        std::set<typename T::node_type *> nextParToCheck;
        for (auto nd : parToCheck) {
            if (nd == nullptr) {
                tree._set_root(nullptr);
                continue;
            }
            if (nd->is_tip()) {
                toPrune.insert(nd);
                nextParToCheck.insert(nd->get_parent());
            }
        }
        nextParToCheck.swap(parToCheck);
    }
    if (tree.get_root() != nullptr) {
        suppress_monotypic_taxa_preserve_shallow_dangle(tree);
    }
    return r;
}


template<typename T>
inline std::string newick(const T &t) {
    std::ostringstream s;
    write_tree_as_newick(s, t);
    return s.str();
}


template <typename Tree_t>
std::unordered_map<OttId, const typename Tree_t::node_type*> get_ottid_to_const_node_map(const Tree_t& T)
{
  std::unordered_map<OttId, const typename Tree_t::node_type*> ottid;
  for(auto nd: iter_pre_const(T))
    if (nd->has_ott_id())
      ottid[nd->get_ott_id()] = nd;

  return ottid;
}

template <typename Tree_t>
std::unordered_map<OttId, typename Tree_t::node_type*> get_ottid_to_node_map(Tree_t& T)
{
    std::unordered_map<OttId, typename Tree_t::node_type*> ottid;
    for(auto nd: iter_pre(T))
        if (nd->has_ott_id())
            ottid[nd->get_ott_id()] = nd;

    return ottid;
}

template <typename N>
N* get_root(N* node) {
    while (node->get_parent()) {
        node = node->get_parent();
    }
    return node;
}

template <typename node_t>
void destroy_children(node_t* node) {
    std::vector<node_t*> nodes;
    while(auto n = node->get_first_child()) {
        n->detach_this_node();
        nodes.push_back(n);
    }
    for(std::size_t i = 0; i < nodes.size(); i++) {
        while(auto n = nodes[i]->get_first_child()) {
            n->detach_this_node();
            nodes.push_back(n);
        }
        delete nodes[i];
    }
    assert(node->is_tip());
}

// Walks the subtree rooted at `nd`. for each branch of the subtree, the rootward-most
//    OTT ID is added to `ott_id_set`
// This is useful for getting a complete list of taxa within a subtree without adding every
//    leaf label.
// If `nd` has an OTT ID, then only that ID will be added to ott_id_set.
template <typename N>
inline void accumulate_closest_ott_id_for_subtree(const N * nd, OttIdSet & ott_id_set) {
    if (nd->has_ott_id()) {
        ott_id_set.insert(nd->get_ott_id());
        return;
    }
    for (auto c : iter_child_const(*nd)) {
        accumulate_closest_ott_id_for_subtree(c, ott_id_set);
    }
}

template<typename T>
void index_nodes_by_name(T & tree) {
    auto & td = tree.get_data();
    auto & m = td.name_to_node;
    for (auto nd : iter_pre(tree)) {
        m[nd->get_name()] = nd;
    }
}

template<typename T>
void set_traversal_entry_exit(T & tree) {
    std::uint32_t ind = 0;
    for (auto nd : iter_pre(tree)) {
        nd->get_data().trav_enter = ind++;
    }
    for (auto pnd : iter_post(tree)) {
        auto fc = pnd->get_last_child();
        auto & d = pnd->get_data();
        if (fc == nullptr) {
            d.trav_exit = d.trav_enter;
        } else {
            d.trav_exit = fc->get_data().trav_exit;
        }
    }
}

template<typename T>
void set_traversal_entry_exit_and_num_tips(T & tree) {
    std::uint32_t ind = 0;
    for (auto nd : iter_pre(tree)) {
        nd->get_data().trav_enter = ind++;
    }
    for (auto pnd : iter_post(tree)) {
        auto fc = pnd->get_last_child();
        auto & d = pnd->get_data();
        if (fc == nullptr) {
            d.trav_exit = d.trav_enter;
            d.num_tips = 1;
        } else {
            d.trav_exit = fc->get_data().trav_exit;
            d.num_tips = 0;
            for (auto c : iter_child_const(*pnd)) {
                d.num_tips += c->get_data().num_tips;
            }
        }
    }
}

 
}// namespace otc
#endif
