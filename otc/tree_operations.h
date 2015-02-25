#ifndef OTCETERA_TREE_OPERATIONS_H
#define OTCETERA_TREE_OPERATIONS_H
// Functions that operate on trees - may include iteration
//	over trees (contrast w/tree_util.h)
// Depends on: tree.h tree_util.h tree_iter.h 
// Depended on by: tools

#include "otc/otc_base_includes.h"
#include "otc/tree_iter.h"

namespace otc {

template<typename T>
unsigned int countPolytomies(const T & tree);
template<typename T>
std::size_t checkForUnknownTaxa(std::ostream & err, const T & toCheck, const T & taxonomy);
template<typename T>
typename T::node_type * findMRCAFromIDSet(T & tree, const std::set<long> & idSet, long trigger);

template<typename T>
void pruneAndDelete(T & tree, typename T::node_type *toDel);
template<typename T, typename U>
inline void cullRefsToNodeFromData(RootedTree<T, U> & tree, RootedTreeNode<T> *toDel);

std::set<long> getDesOttIds(RootedTreeNode<RTSplits> & nd);

//// impl

template<typename T>
unsigned int countPolytomies(const T & tree) {
	unsigned int n = 0U;
	for (auto node : ConstPostorderInternalIter<T>{tree}) {
		if (node->getOutDegree() > 2) {
			n += 1;
		}
	}
	return n;
}

template<typename T>
void fillDesIdSets(T & tree) {
	// assumes OttId is set for each tip
	for (auto node : PostorderIter<T>(tree)) {
		std::set<long> & desIds = node->getData().desIds;
		if (node->isTip()) {
			desIds.insert(node->getOttId());
		} else {
			for (auto child : ChildIter<typename T::node_type>(*node)) {
				std::set<long> & cDesIds = child->getData().desIds;
				desIds.insert(cDesIds.begin(), cDesIds.end());
			}
		}
	}
}

// uses ottID->node mapping, but the split sets of the nodes
template<typename T>
typename T::node_type * findMRCAFromIDSet(T & tree, const std::set<long> & idSet, long trigger) {
	typedef typename T::node_type NT_t;
	auto ottIdToNode = tree.getData().ottIdToNode;
	std::map<NT_t *, unsigned int> n2c;
	long shortestPathLen = -1;
	NT_t * shortestPathNode = nullptr;
	for (auto i : idSet) {
		auto rIt = ottIdToNode.find(i);
		if (rIt == ottIdToNode.end()) {
			std::string em = "tip ";
			em += std::to_string(i);
			if (trigger >= 0) {
				em += " a descendant of ";
				em += std::to_string(trigger); 
			}
			em += " not found.";
			throw OTCError(em);
		}
		auto nd = rIt->second;
		long currPathLen = 0;
		while (nd != nullptr) {
			n2c[nd] += 1;
			currPathLen += 1;
			nd = nd->getParent();
		}
		if (shortestPathLen < 0 || currPathLen < shortestPathLen) {
			shortestPathLen = currPathLen;
			shortestPathNode = rIt->second;
		}
	}
	auto nTips = idSet.size();
	auto cn = shortestPathNode;
	while (true) {
		if (n2c[cn] == nTips) {
			return cn;
		}
		cn = cn->getParent();
		assert(cn != nullptr);
	}
	assert(false);
	return nullptr;
}

template<typename T>
std::size_t checkForUnknownTaxa(std::ostream & err, const T & toCheck, const T & taxonomy) {
	auto taxOttIds = taxonomy.getRoot()->getData().desIds;
	auto toCheckOttIds = toCheck.getRoot()->getData().desIds;
	auto extras = set_difference_as_set(toCheckOttIds, taxOttIds);
	if (!extras.empty()) {
		err << "OTT Ids found in an input tree,  but not in the taxonomy:\n";
		writeOttSet(err, "  ", extras, "\n");
		return extras.size();
	}
	return 0U;
}

template<typename T>
inline typename T::node_type * addChildForOttId(typename T::node_type & nd, long ottId, T & tree) {
	auto nn = tree.createChild(&nd);
	nn->setOttId(ottId);
	return nn;
}

inline std::set<long> getDesOttIds(RootedTreeNode<RTSplits> & nd) {
	return nd.getData().desIds;
}

inline void fixDesIdFields(RootedTreeNode<RTSplits> & nd, const std::set<long> & ls) {
	const std::set<long> toRemove = nd.getData().desIds;
	assert(!toRemove.empty());
	nd.getData().desIds = ls;
	for (auto anc : AncIter<RootedTreeNode<RTSplits> >(&nd)) {
		assert(anc != nullptr);
		assert(!anc->getData().desIds.empty());
		for (auto tr : toRemove) {
			assert(contains(anc->getData().desIds, tr));
			anc->getData().desIds.erase(tr);
		}
		anc->getData().desIds.insert(begin(ls), end(ls));
	}
}

template<typename T>
std::vector<typename T::node_type *> expandOTTInternalsWhichAreLeaves(T & toExpand, const T & taxonomy) {
	const auto & taxData = taxonomy.getData();
	std::map<typename T::node_type *, std::set<long> > replaceNodes;
	for (auto nd : LeafIter<T>(toExpand)) {
		assert(nd->isTip());
		assert(nd->hasOttId());
		auto ottId = nd->getOttId();
		auto taxNd = taxData.getNodeForOttId(ottId);
		assert(taxNd != nullptr);
		if (!taxNd->isTip()) {
			auto leafSet = getDesOttIds(*taxNd);
			replaceNodes[nd] = leafSet;
		}
	}
	std::vector<typename T::node_type *> expanded;
	expanded.reserve(replaceNodes.size());
	for (auto r : replaceNodes) {
		auto oldNode = r.first;
		auto ls = r.second;
		assert(ls.size() > 0);
		expanded.push_back(oldNode);
		for (auto loid : ls) {
			addChildForOttId<T>(*oldNode, loid, toExpand);
		}
		fixDesIdFields(*oldNode, ls);
	}
	return expanded;
}

template<typename T>
void markPathToRoot(const T & fullTree,
					long ottId,
					std::map<const typename T::node_type *, std::set<long> > &n2m){
	auto startNd = fullTree.getData().getNodeForOttId(ottId);
	assert(startNd != nullptr);
	if (startNd == nullptr) {
		std::string m = "OTT id not found ";
		m += std::to_string(ottId);
		throw OTCError(m);
	}
	n2m[startNd].insert(ottId);
	for (auto nd : AncIter<typename T::node_type>(startNd)) {
		n2m[nd].insert(ottId);
	}
}


// find most recent anc of nd with out-degree > 1
template<typename T>
inline T * findFirstBranchingAnc(T * nd) {
	T * anc = nd->getParent();
	if (anc == nullptr) {
		return nullptr;
	}
	while (anc->isOutDegreeOneNode()) {
		anc = anc->getParent();
	}
	return anc;
}

template<typename T, typename U>
inline bool multipleChildrenInMap(const T & nd,
								  const std::map<const T *, U> & markedMap,
								  const T **first) {
	assert(first);
	bool foundFirst = false;
	*first = nullptr;
	for(auto c : ConstChildIter<T>(nd)) {
		if (markedMap.find(c) != markedMap.end()) {
			if (foundFirst) {
				return true;
			}
			foundFirst = true;
			*first = c;
		}
	}
	return false;
}

template<typename T>
inline const T * findNextSignificantNode(const T * node, const std::map<const T *, std::set<long> > & markedMap) {
	assert(node != nullptr);
	auto currNode = node;
	for (;;) {
		const T * sc;
		if (multipleChildrenInMap(*currNode, markedMap, &sc)) {
			return currNode;
		}
		if (sc == 0L) {
			const char * msg = "Failing. Node found with ottIDs marked, but no children with ottIDs marked";
			throw OTCError(msg);
		}
		currNode = sc;
	}
}

template<typename T>
inline void writePrunedSubtreeNewickForMarkedNodes(std::ostream & out,
											const T & srcNd,
											const std::map<const T *, std::set<long> > & markedMap) {
	const auto nIt = markedMap.find(&srcNd);
	assert(nIt != markedMap .end());
	const auto & ottIDSet = nIt->second;
	if (ottIDSet.size() == 1) {
		out << "ott" << *ottIDSet.begin();
	} else {
		assert(ottIDSet.size() > 1);
		auto nsn = findNextSignificantNode<T>(&srcNd, markedMap);
		out << '(';
		unsigned numcwritten = 0;
		for (auto child : ConstChildIter<T>(*nsn)) {
			if (markedMap.find(child) != markedMap.end()){
				if (numcwritten > 0) {
					out << ',';
				}
				writePrunedSubtreeNewickForMarkedNodes(out, *child, markedMap);
				++numcwritten;
			}
		}
		assert(numcwritten > 1);
		out << ')';
	}
}



template<typename T>
inline void describeUnnamedNode(const T & nd,
								std::ostream & out,
								unsigned int anc,
								bool useNdNames) {
	if (useNdNames && !nd.getName().empty()) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before \"" << nd.getName() << "\"\n";
		} else {
			out << "the node \"" << nd.getName() << "\"\n";
		}
	}
	else if (nd.isTip()) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before the leaf \"" << nd.getName()  << "\"\n";
		} else {
			out << "the leaf \"" << nd.getName()  << "\"\n";
		}
	} else if (nd.isOutDegreeOneNode()) {
		describeUnnamedNode(*nd.getFirstChild(), out, anc + 1, useNdNames);
	} else {
		auto left = findLeftmostInSubtree(&nd)->getName();
		auto right = findRightmostInSubtree(&nd)->getName();
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right << "\"\n";
		} else {
			out <<  "MRCA of \"" << left << "\" and " << "\"" << right << "\"\n";
		}
	}
}

template<typename T, typename U>
inline void cullRefsToNodeFromData(RootedTree<T, U> & , RootedTreeNode<T> *) {
	std::cout << "generic!\n";
}

template<typename T>
inline void cullRefsToNodeFromData(RootedTree<T, RTreeOttIDMapping<T> > & , RootedTreeNode<T> *) {
	std::cout << "partial!\n";
}


template<typename T>
inline void pruneAndDelete(T & tree, typename T::node_type *toDel) {
	cullRefsToNodeFromData<typename T::node_data_type, typename T::data_type>(tree, toDel);
	tree._pruneAndDelete(toDel);
}

template <typename T, typename U>
void insertAncestorsToParaphyleticSet(T * nd, U & includedNodes) {
	for (auto anc : AncIter<T>(nd)) {
		if (contains(includedNodes, anc)) {
			return;
		}
		includedNodes.insert(nd);
	}
}

//@TMP recursive until we have a pre-order subtree skipping iter.
template <typename T, typename U>
void insertDescendantsOfUnincludedSubtrees(T * nd, U & includedNodes) {
	for (auto c : ChildIter<T>(*nd)) {
		if (!contains(includedNodes, c)) {
			includedNodes.insert(c);
			insertDescendantsOfUnincludedSubtrees(c, includedNodes);
		}
	}
}


template<typename T>
inline void writeNodeAsNewickLabel(std::ostream & out, const T *nd) {
	if (nd->isTip()) {
		writeEscapedForNewick(out, nd->getName());
	} else if (!nd->getName().empty()) {
		writeEscapedForNewick(out, nd->getName());
	}
}

template<typename T>
inline void writeClosingNewick(std::ostream & out, const T *nd, const T * r) {
	out << ')';
	auto n = nd->getParent();
	writeNodeAsNewickLabel(out, n);
	if (n == r) {
		return;
	}
	while (n->getNextSib() == nullptr) {
		if (n == r) {
			return;
		}
		out << ')';
		n = n->getParent();
		assert(n != nullptr);
		writeNodeAsNewickLabel(out, n);
	}
	out << ',';
}
template<typename T>
inline void writeNewick(std::ostream & out, const T *nd) {
	assert(nd != nullptr);
	if (nd->isTip()) {
		writeNodeAsNewickLabel(out, nd);
	} else {
		for (auto n : ConstPreorderIterN<T>(nd)) {
			if (n->isTip()) {
				writeNodeAsNewickLabel(out, n);
				if (n->getNextSib() == nullptr) {
					writeClosingNewick<T>(out, n, nd);
				} else {
					out << ',';
				}
			} else {
				out << '(';
			}
		}
	}
}

template<typename T>
inline void writeTreeAsNewick(std::ostream & out, const T &tree) {
	writeNewick<typename T::node_type>(out, tree.getRoot());
	out << ';';
}


} // namespace otc
#endif

