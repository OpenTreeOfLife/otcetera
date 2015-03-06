#ifndef OTCETERA_SUPERTREE_CONTEXT_H
#define OTCETERA_SUPERTREE_CONTEXT_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/tree_data.h"

namespace otc {
template<typename T, typename U> class NodeThreading;
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
        std::map<const NodeWithSplits *, NodeThreading<T, U> > & scaffold2NodeThreading;
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
                         std::map<const RootedTreeNode<RTSplits> *, NodeThreading<T, U> > & taxoToAlignment,
                         TreeMappedWithSplits & scaffTree)
            :numTrees(nt),
            scaffold2NodeThreading(taxoToAlignment),
            scaffoldOttId2Node(scaffTree.getData().ottIdToNode),
            scaffoldTree(scaffTree) {
        }
};

} // namespace
#endif
