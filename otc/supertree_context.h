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
        std::set<const U *> detachedScaffoldNodes;
        const std::size_t numTrees;
        std::map<const NodeWithSplits *, NodeEmbedding<T, U> > & scaffold2NodeEmbedding;
        std::map<long, typename U::node_type *> & scaffoldOttId2Node;
        RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & scaffoldTree; // should adjust the templating to make more generic

        void log(SupertreeCtorEvent e, const U & node) {
            if (e == COLLAPSE_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.getOttId())});
            } else if (e == IGNORE_TIP_MAPPED_TO_NONMONOPHYLETIC_TAXON) {
                events.emplace_back(LogEvent{e, std::string("ott") + std::to_string(node.getOttId())});
            }
        }
        SupertreeContext(std::size_t nt,
                         std::map<const RootedTreeNode<RTSplits> *, NodeEmbedding<T, U> > & taxoToEmbedding,
                         TreeMappedWithSplits & scaffTree)
            :numTrees(nt),
            scaffold2NodeEmbedding(taxoToEmbedding),
            scaffoldOttId2Node(scaffTree.getData().ottIdToNode),
            scaffoldTree(scaffTree) {
        }
};

} // namespace
#endif
