#ifndef OTCETERA_TREE_OPERATIONS_H
#define OTCETERA_TREE_OPERATIONS_H
// Functions that operate on trees - may include iteration
//	over trees (contrast w/tree_util.h)
// Depends on: tree.h tree_util.h tree_iter.h 
// Depended on by: tools

#include "otc/otc_base_includes.h"
#include "otc/tree_iter.h"

namespace otc {

template<typename T, typename U>
unsigned int countPolytomies(const RootedTree<T, U> & tree);


template<typename T, typename U>
unsigned int countPolytomies(const RootedTree<T, U> & tree) {
	unsigned int n = 0U;
	for (auto node : ConstPostorderInternalNode<T, U>{tree}) {
		if (node->getOutDegree() > 2) {
			n += 1;
		}
	}
	return n;
}

template<typename T, typename U>
void fillDesIdSets(RootedTree<T, U> & tree) {
	// assumes OttId is set for each tip
	for (auto node : PostorderInternalNode<T, U>(tree)) {
		std::set<long> & desIds = node->getData().desIds;
		if (node->IsTip()) {
			desIds.insert(node->getOttId());
		} else {
			for (auto child : ChildIterator<T, U>(*node)) {
				std::set<long> & cDesIds = child->getData().desIds;
				desIds.insert(cDesIds.begin(), cDesIds.end());
			}
		}
	}
}
} // namespace otc
#endif

