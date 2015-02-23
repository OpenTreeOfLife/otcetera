#ifndef OTCETERA_DEBUG_H
#define OTCETERA_DEBUG_H

#include "otc/otc_base_includes.h"
#include "otc/tree_operations.h"

namespace otc {

template<typename T>
bool checkNodePointers(const T & nd) {
	bool good = true;
	auto c = nd.getFirstChild();
	while (c != nullptr) {
		if (c->getParent() != &nd) {
			good = false;
			assert(c->getParent() == &nd);
		}
		c = c->getNextSib();
	}
	return good;
}

template<typename T>
bool checkAllNodePointers(const T & tree) {
	bool good = true;
	auto ns = tree.getSetOfAllNodes();
	for (auto nd : ns) {
		if (!checkNodePointers(*nd)) {
			return false;
		}
	}
	if (tree.getRoot()->getParent() != nullptr) {
		assert(tree.getRoot()->getParent() == nullptr);
		return false;
	}
	return true;
}

template<typename T>
bool checkPreorder(const T & tree) {
	auto ns = tree.getSetOfAllNodes();
	std::set<const typename T::node_type *> visited;
	for (auto nd : ConstPreorder<T>(tree)) {
		auto p = nd->getParent();
		if (p) {
			if (!contains(visited, p)) {
				assert(contains(visited, p));
				return false;
			}
		}
		visited.insert(nd);
	}
	if (visited != ns) {
		assert(visited == ns);
		return false;
	}
	return true;
}


template<typename T>
bool checkTreeInvariants(const T & tree) {
	if (!checkAllNodePointers(tree)) {
		return false;
	}
	if (!checkPreorder(tree)) {
		return false;
	}
	return true;
}

} // namespace otc
#endif

