#ifndef OTCETERA_TREE_UTIL_H
#define OTCETERA_TREE_UTIL_H
// Very simple functions that operate on trees or treenode
//	without requiring iteration (contrast w/tree_operations.h)
// Depends on: tree.h 
// Depended on by: tree_iter.h tree_operation.h 

#include "otc/otc_base_includes.h"
#include "otc/tree.h"

namespace otc {

template<typename T>
bool isInternalNode(const T & nd);
template<typename T>
bool isLeaf(const T & nd);
template<typename T>
T * findLeftmostInSubtreeM(T * nd);
template<typename T>
T * findRightmostInSubtreeM(T * nd);
template<typename T>
const T * findLeftmostInSubtree(const T * nd);
template<typename T>
const T * findRightmostInSubtree(const T * nd);

template<typename T>
inline bool isInternalNode(const T & nd) {
	return !nd.isTip();
}
template<typename T>
inline bool isLeaf(const T & nd) {
	return nd.isTip();
}

template<typename T>
inline const T * findLeftmostInSubtree(const T * nd) {
	if (nd == nullptr) {
		return nullptr;
	}
	auto next = nd->getFirstChild();
	while (next != nullptr) {
		nd = next;
		next = nd->getFirstChild();
	}
	return nd;
}
template<typename T>
inline T * findLeftmostInSubtreeM(T * nd) {
	const T * c = const_cast<const T *>(nd);
	return const_cast<T*>(findLeftmostInSubtree<T>(c));
}

template<typename T>
inline const T * findRightmostInSubtree(const T * nd) {
	if (nd == nullptr) {
		return nullptr;
	}
	auto next = nd->getLastChild();
	while (next != nullptr) {
		nd = next;
		next = nd->getLastChild();
	}
	return nd;
}

template<typename T>
inline T * findRightmostInSubtreeM(T * nd) {
	const T * c = const_cast<const T *>(nd);
	return const_cast<T*>(findRightmostInSubtree<T>(c));
}

} // namespace otc
#endif

