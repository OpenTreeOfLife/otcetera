#ifndef OTCETERA_TREE_DATA_H
#define OTCETERA_TREE_DATA_H
// Classes that can serve as the template args for trees and nodes
#include <map>
#include <set>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename, typename> class RootedTree;
template<typename> class RootedTreeNode;

template<typename T>
class RTreeOttIDMapping {
	public:
		typedef RootedTreeNode<T> NodeType;
		std::map<long, NodeType *> ottIdToNode;
		// if a node is pruned, the entry in ottIdToNode will refer to an alias
		//	 if the alias is later pruned, we need to know what nodes it is aliasing
		//	 so that ottIdToNode does not point to dangling nodes.
		std::map<NodeType *, std::set<long> > isAliasFor;
		std::map<long, NodeType *> ottIdToDetachedNode;
		bool desIdSetsContainInternals;
		NodeType * getNodeForOttId(long ottId) const {
			const auto it = ottIdToNode.find(ottId);
			return (it == ottIdToNode.end() ? nullptr : it->second);
		}
};

class RTSplits {
	public:
		std::set<long> desIds;
};

typedef otc::RootedTreeNode<RTSplits> NodeWithSplits;
using NodeWithSplitsPred = std::function<bool(const NodeWithSplits &)>;
typedef otc::RTreeOttIDMapping<RTSplits> MappedWithSplitsData;
typedef otc::RTreeOttIDMapping<RTNodeNoData> MappedWithEmptyNodeData;
typedef otc::RootedTree<RTSplits, MappedWithSplitsData> TreeMappedWithSplits;
typedef otc::RootedTree<RTNodeNoData, MappedWithEmptyNodeData> TreeMappedEmptyNodes;

template<typename T>
inline void verifyOttIdMapping(const T & tree) {
	for (const auto & mp : tree.getData().ottIdToNode) {
		if (mp.first != mp.second->getOttId()) {
			assert(mp.first == mp.second->getOttId());
			throw OTCError("Cache of OTT ID->node is stale " + std::to_string(mp.first) + " != " + std::to_string(mp.second->getOttId()));
		}
	}
}

} // namespace otc
#endif
