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
		NodeType * getNodeForOttId(long ottId) const {
			const auto it = ottIdToNode.find(ottId);
			return (it == ottIdToNode.end() ? nullptr : it->second);
		}
};

class RTSplits {
	public:
		std::set<long> desIds;
};

} // namespace otc
#endif
