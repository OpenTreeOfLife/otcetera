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
inline void update_ancestral_path_ott_id_set(T * nd,
                                        const OttIdSet & oldEls,
                                        const OttIdSet & newEls,
                                        std::map<const T *, NodeEmbedding<T, U> > & m) {
    auto & curr = m.at(nd);
    assert(oldEls.size() > 0);
    LOG(DEBUG) << "  " << nd->get_ott_id() << " calling update_all_paths_ott_id_sets";
    if (!curr.update_all_paths_ott_id_sets(oldEls, newEls)) {
        return;
    }
    for (auto anc : iter_anc(*nd)) {
        auto & ant = m.at(anc);
        LOG(DEBUG) << "  " << anc->get_ott_id() << " calling update_all_paths_ott_id_sets";
        if (!ant.update_all_paths_ott_id_sets(oldEls, newEls)) {
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
        2. all paths that are loops for embeddedNode (embeddedNode is the scaffold_des
            and the scaffold_anc for the PathPair), and
        3. all paths the leave this node. ie. those that paths that:
            A. have a scaffold_des that is embeddedNode or one of its descendants, AND
            B. have scaffold_anc set to an ancestor of embeddedNode
    Note that if you want all of the paths that intersect with embeddedNode, you have
        to check the edgeBelowEmbeddings field of all of embeddedNode's children. This
        is because any PathPairing that has embeddedNode as the scaffold_anc will not
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
    //
    //  If we are requested to retain nodes mapped to non-terminal nodes in the scaffold, and
    //      these scaffold nodes are contested we run into some ugliness.
    //  The code was not designed for this case. The scaffold node is pruned by
    //      of collapse_group. This means that these phyloNodes will be "unembedded".
    //  The hack for covering this case is to simply propagate pointers to such nodes (and their
    //      parent nodes) back in the scaffold every time a group is collapsed. 
    //  This allows these retained tips to be emitted as a part of the subproblem when
    //      the root of the phyloTree or an uncontested node is found.
    //  It is not clear if there are some corner cases for when this would be problematic.
    std::map<std::size_t, std::map<U *, U *> > phyloNd2ParForUnembeddedTrees;
    public:
    NodeEmbedding(T * scaffNode)
        :embeddedNode(scaffNode) {
    }
    std::size_t get_total_num_node_mappings() const {
        unsigned long t = 0U;
        for (auto i : nodeEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    std::set<U *> get_phylo_nodes(std::size_t treeIndex) const {
        auto neIt = nodeEmbeddings.find(treeIndex);
        std::set<U *> r;
        if (neIt == nodeEmbeddings.end()) {
            return r;
        }
        for (const auto & np : neIt->second) {
            r.insert(np->phylo_node);
        }
        return r;
    }
    std::size_t get_num_loop_trees() const {
        std::set<std::size_t> keys;
        for (auto i : loopEmbeddings) {
            keys.insert(i.first);
        }
        return keys.size();
    }
    std::size_t get_total_num_loops() const {
        unsigned long t = 0U;
        for (auto i : loopEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    std::size_t get_total_num_edges_below_traversals() const {
        unsigned long t = 0U;
        for (auto i : edgeBelowEmbeddings) {
            t += i.second.size();
        }
        return t;
    }
    static bool does_tree_constest_monophyly(const PathPairSet & edgesBelowForTree);
    using ContestingNodeMap = std::map<const T *, std::set<const T *> >;
    using TreeIndToContestingNodeMap = std::map<std::size_t, ContestingNodeMap >;
    static std::map<const T *, std::set<const T *> > how_tree_constests_monophyly(const PathPairSet & edgesBelowForTree);

    bool is_contested() const {
        for (auto i : edgeBelowEmbeddings) {
            if (does_tree_constest_monophyly(i.second)) {
                return true;
            }
        }
        return false;
    }
    std::list<std::size_t> get_contesting_tree_indices() const {
        std::list<std::size_t> r;
        for (auto i : edgeBelowEmbeddings) {
            if (does_tree_constest_monophyly(i.second)) {
                r.push_back(i.first);
            }
        }
        return r;
    }
    TreeIndToContestingNodeMap get_how_tree_contests_monophyly_maps() const {
        TreeIndToContestingNodeMap retMap;
        for (auto i : edgeBelowEmbeddings) {
            if (does_tree_constest_monophyly(i.second)) {
                retMap[i.first] = how_tree_constests_monophyly(i.second);
            }
        }
        return retMap;
    }
    const PathPairSet & get_edges_exiting(std::size_t treeIndex) const {
        auto el = edgeBelowEmbeddings.find(treeIndex);
        assert(el != edgeBelowEmbeddings.end());
        return el->second;
    }
    // some trees contest monophyly. Return true if these trees are obviously overruled
    //   by higher ranking trees so that we can avoid the more expensive unconstrained phylo graph 
    bool high_ranking_trees_preserve_monophyly(std::size_t ) {
        return false;
    }
    const OttIdSet & get_relevant_des_ids_from_path(const PathPairing<T, U> & pps);
    OttIdSet get_relevant_des_ids_from_path_pair_set(const PathPairSet & pps);
    OttIdSet get_relevant_des_ids(const std::map<const T *,
                               NodeEmbedding<T, U> > & eForNd,
                               std::size_t treeIndex);

    void collapse_source_edge(const T * phylo_parent,
                            PathPairing<T, U> * path);
    void collapse_source_edge_to_force_one_entry(T & ,
                                            PathPairSet & pps,
                                            std::size_t treeIndex,
                                            SupertreeContextWithSplits &);
    void resolve_given_contested_monophyly(T & scaffold_node,
                                        SupertreeContextWithSplits & sc);
    std::set<PathPairPtr> get_all_child_exit_paths(
                            const T & scaffold_node,
                            const std::map<const T *, NodeEmbedding<T, U> > & sc) const;
    std::set<PathPairPtr> get_all_child_exit_paths_for_tree(
                            const T & scaffold_node,
                            std::size_t treeIndex,
                            const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const;
    void resolve_given_uncontested_monophyly(T & scaffold_node,
                                          SupertreeContextWithSplits & sc);
    std::string export_subproblem_and_resolve(T & scaffold_node,
                                    const std::string & exportDir,
                                    std::ostream * exportStream, // nonnull to override exportdir
                                    SupertreeContextWithSplits & sc);
    void collapse_group(T & scaffold_node, SupertreeContext<T, U> & sc);
    void prune_collapsed_node(T & scaffold_node, SupertreeContextWithSplits & sc);
    
    bool report_if_contested(std::ostream & out,
                           const T * nd,
                           const std::vector<TreeMappedWithSplits *> & treePtrByIndex,
                           const std::vector<NodeWithSplits *> & aliasedBy,
                           bool verbose) const;
    bool update_all_paths_ott_id_sets(const OttIdSet & oldEls, const OttIdSet & newEls) {
        bool r = update_all_mapped_paths_ott_id_sets(loopEmbeddings, oldEls, newEls);
        return update_all_mapped_paths_ott_id_sets(edgeBelowEmbeddings, oldEls, newEls) || r;
    }
    std::vector<const PathPairing<T, U> *> get_all_incoming_path_pairs(
                        const std::map<const T *, NodeEmbedding<T, U> > & eForNd,
                        std::size_t treeIndex) const;
    bool debug_node_embeddings(const char * tag,
                            bool isUncontested,
                            const std::map<const T *, NodeEmbedding<T, U> > & sn2ne) const;
    void add_node_embeddings(std::size_t treeIndex, NodePairPtr npp) {
        nodeEmbeddings[treeIndex].insert(npp);
    }
    void add_loop_embeddings(std::size_t treeIndex, PathPairPtr pp) {
        loopEmbeddings[treeIndex].insert(pp);
    }
    void add_exit_embeddings(std::size_t treeIndex, PathPairPtr pp) {
        edgeBelowEmbeddings[treeIndex].insert(pp);
    }
    void set_ott_id_for_exit_embeddings(
                        T * newScaffDes,
                        long ottId,
                        std::map<const T *, NodeEmbedding<T, U> > & n2ne);
    void merge_exit_embeddings_if_multiple();
    void resolve_parent_in_favor_of_this_node(
                        T & scaffold_node,
                        std::size_t treeIndex,
                        SupertreeContextWithSplits & sc);
    const TreeToPathPairs & get_exit_embeddings() const {
        return edgeBelowEmbeddings;
    }
    std::map<U *, U *> get_un_embedded_phylo_node_to_par(std::size_t treeInd) const;
    void debugPrint(T & scaffold_node, std::size_t treeIndex, const std::map<const T *, NodeEmbedding<T, U> > & sc) const;

    private:
    std::map<U *, U*> get_looped_phylo_node_to_par(std::size_t treeInd) const;
    std::map<U *, U*> get_exit_phylo_node_to_par(std::size_t treeInd) const;
    PathPairSet refers_to_node(std::size_t treeInd, U *n) const {
        PathPairSet r;
        auto pit = edgeBelowEmbeddings.find(treeInd);
        if (pit != edgeBelowEmbeddings.end()) {
            for (auto i : pit->second) {
                if (i->phylo_parent == n || i->phylo_child == n) {
                    r.insert(i);
                }
            }
        }
        pit = loopEmbeddings.find(treeInd);
        if (pit != loopEmbeddings.end()) {
            for (auto i : pit->second) {
                if (i->phylo_parent == n || i->phylo_child == n) {
                    r.insert(i);
                }
            }
        }
        return r;
    }
    void remove_ref_to_exit_path(std::size_t treeIndex, PathPairPtr toDel) {
        auto & pps = edgeBelowEmbeddings.at(treeIndex);
        assert(contains(pps, toDel));
        pps.erase(toDel);
    }
    void prune_suppressed(std::size_t treeIndex, U * phyloPar, U * phylo_child);
};

template<typename T, typename U>
inline bool NodeEmbedding<T, U>::does_tree_constest_monophyly(const std::set<PathPairing<T, U> *> & edgesBelowForTree) {
    if (edgesBelowForTree.size() > 1) {
        const T * firstSrcPar = nullptr;
        for (auto pp : edgesBelowForTree) {
            auto sp = pp->phylo_parent;
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

/// Returns a mapping from phylo parent to children of that taxon that belong to this taxonomic node
// If the node is contested by this tree, there should be more than one entry in the map
template<typename T, typename U>
inline std::map<const T *, std::set<const T *> > NodeEmbedding<T, U>::how_tree_constests_monophyly(const std::set<PathPairing<T, U> *> & edgesBelowForTree) {
    std::map<const T *, std::set<const T *> > retMap;
    if (edgesBelowForTree.size() > 1) {
        for (auto pp : edgesBelowForTree) {
            auto sp = pp->phylo_parent;
            auto sc = pp->phylo_child;
            retMap[sp].insert(sc);
        }
    }
    return retMap;
}

template<typename T>
inline bool update_all_mapped_paths_ott_id_sets(T & mPathSets, const OttIdSet & oldEls, const OttIdSet & newEls) {
    bool r = false;
    for (auto mpIt : mPathSets) {
        for (auto p : mpIt.second) {
            r = r || p->update_ott_id_set_no_traversal(oldEls, newEls);
        }
    }
    return r;
}

template<typename T, typename U, typename V>
inline std::map<U *, U*> get_node_to_par_for_key(const V & treeInd,
                                         const std::map<V, std::set<PathPairing<T, U> *> >& m);
template<typename T, typename U, typename V>
inline std::map<U *, U*> get_node_to_par_for_key(const V & treeInd,
                                         const std::map<V, std::set<PathPairing<T, U> *> >& m) {
    std::map<U *, U*> nd2par;
    if (!contains(m, treeInd)) {
        return nd2par;
    }
    for (const auto & pp : m.at(treeInd)) {
        LOG(DEBUG) << " considering the edge from the child "  << reinterpret_cast<long>(pp->phylo_child) << ": "; dbWriteNewick(pp->phylo_child);
        LOG(DEBUG) << "             to its parent "  << reinterpret_cast<long>(pp->phylo_parent) << ": "; dbWriteNewick(pp->phylo_parent);
        assert(!contains(nd2par, pp->phylo_child));
        nd2par[pp->phylo_child] = pp->phylo_parent;
    }
    return nd2par;
}

template<typename T, typename U>
inline std::map<U *, U*> NodeEmbedding<T, U>::get_looped_phylo_node_to_par(std::size_t treeInd) const {
    return get_node_to_par_for_key(treeInd, loopEmbeddings);
}

template<typename T, typename U>
inline std::map<U *, U *> NodeEmbedding<T, U>::get_exit_phylo_node_to_par(std::size_t treeInd) const {
    return get_node_to_par_for_key(treeInd, edgeBelowEmbeddings);
}

template<typename T, typename U>
inline std::map<U *, U *> NodeEmbedding<T, U>::get_un_embedded_phylo_node_to_par(std::size_t treeInd) const {
    if (contains(phyloNd2ParForUnembeddedTrees, treeInd)) {
        return phyloNd2ParForUnembeddedTrees.at(treeInd);
    }
    std::map<U *, U *> r;
    return r;
}



} // namespace
#endif
