#ifndef OTCETERA_SUPERTREE_CONTEXT_H
#define OTCETERA_SUPERTREE_CONTEXT_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_data.h"

namespace otc {
template<typename T, typename U> class NodeEmbedding;
enum SupertreeCtorEvent {
    COLLAPSE_TAXON,
    IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON,
    CLADE_CREATES_TREE,
    CLADE_REJECTED,
    CLADE_ADDED_TO_TREE
};
template<typename T, typename U>
class SupertreeContext {
    public:
        SupertreeContext(const SupertreeContext &) = delete;
        SupertreeContext & operator=(const SupertreeContext &) = delete;
        using LogEvent = std::pair<SupertreeCtorEvent, std::string>;

        std::list<LogEvent> events;
        std::set<const U *> detached_scaffold_nodes;
        std::vector<const TreeMappedWithSplits *> trees_by_index;
        const std::size_t num_trees;
        std::map<const NodeWithSplits *, NodeEmbedding<T, U> > & scaffold_to_node_embedding;
        std::map<long, typename U::node_type *> & scaffold_ott_id_to_node;
        RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & scaffold_tree; // should adjust the templating to make more generic
        std::map<std::size_t, std::set<NodeWithSplits *> > pruned_subtrees; // when a tip is mapped to a non-monophyletc terminal it is pruned
        std::list<NodePairingWithSplits> node_pairings_from_resolve;
        std::list<PathPairingWithSplits> path_pairings_from_resolve;
        bool prune_tips_mapped_to_contested_taxa;
        void log(SupertreeCtorEvent e, const U & node) {
            if (e == COLLAPSE_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.get_ott_id())});
            } else if (e == IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.get_ott_id())});
            }
        }
        SupertreeContext(const std::vector<TreeMappedWithSplits *> & tv,
                         std::map<const RootedTreeNode<RTSplits> *, NodeEmbedding<T, U> > & scaffoldNdToNodeEmbedding,
                         TreeMappedWithSplits & scaffTree)
            :num_trees(tv.size()),
            scaffold_to_node_embedding(scaffoldNdToNodeEmbedding),
            scaffold_ott_id_to_node(scaffTree.get_data().ottIdToNode),
            scaffold_tree(scaffTree),
            prune_tips_mapped_to_contested_taxa(true) {
            trees_by_index.reserve(num_trees); 
            for (auto tp : tv) {
                trees_by_index.push_back(const_cast<const TreeMappedWithSplits *>(tp));
            }
        }
};

} // namespace
#endif
