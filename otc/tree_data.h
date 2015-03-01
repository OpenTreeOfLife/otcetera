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

} // namespace otc
#endif
