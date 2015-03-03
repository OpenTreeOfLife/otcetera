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
	auto ns = tree.getSetOfAllAttachedNodes();
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
	auto ns = tree.getSetOfAllAttachedNodes();
	std::set<const typename T::node_type *> visited;
	for (auto nd : iter_pre_const(tree)) {
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
bool checkPostorder(const T & tree) {
	auto ns = tree.getSetOfAllAttachedNodes();
	std::set<const typename T::node_type *> visited;
	for (auto nd : iter_post_const(tree)) {
		auto p = nd->getParent();
		if (p) {
			if (contains(visited, p)) {
				assert(!contains(visited, p));
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
bool checkChildIter(const T & tree) {
	auto ns = tree.getSetOfAllAttachedNodes();
	for (auto nd : ns) {
		auto nc = nd->getOutDegree();
		auto v = 0U;
		for (auto c : iter_child_const(*nd)) {
			assert(c->getParent() == nd);
			++v;
		}
		assert(v == nc);
	}
	return true;
}

template<typename T>
bool checkDesIds(const T & tree) {
	auto ns = tree.getSetOfAllAttachedNodes();
	for (auto nd : ns) {
		if (nd->isTip()) {
			continue;
		}
		std::set<long> d;
		auto sum = 0U;
		for (auto c : iter_child_const(*nd)) {
			const auto & cd = c->getData().desIds;
			sum += cd.size();
			assert(cd.size() > 0);
			d.insert(cd.begin(), cd.end());
		}
		assert(sum == d.size());
		if (tree.getData().desIdSetsContainInternals) {
			if (nd->hasOttId()) {
				d.insert(nd->getOttId());
			}
		}
		if (!isSubset(d, nd->getData().desIds)) {
			std::cerr << "node " << nd->getOttId() << '\n';
			writeOttSetDiff(std::cerr, " ", nd->getData().desIds, " node.desId ", d, " calc.");
			assert(isSubset(d, nd->getData().desIds));
		}
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
	if (!checkPostorder(tree)) {
		return false;
	}
	if (!checkChildIter(tree)) {
		return false;
	}
	if (!checkDesIds(tree)) {
		return false;
	}
	return true;
}

} // namespace otc
#endif

