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
std::size_t checkForUnknownTaxa(std::ostream & err, const RootedTree<T, U> & toCheck, const RootedTree<T, U> & taxonomy);
template<typename T, typename U>
RootedTreeNode<T> * findMRCAFromIDSet(RootedTree<T, U> & tree, const std::set<long> & idSet, long trigger);

std::set<long> getDesOttIds(RootedTreeNode<RTSplits> & nd);

//// impl

template<typename T, typename U>
unsigned int countPolytomies(const RootedTree<T, U> & tree) {
	unsigned int n = 0U;
	for (auto node : ConstPostorderInternalNode<RootedTree<T, U> >{tree}) {
		if (node->getOutDegree() > 2) {
			n += 1;
		}
	}
	return n;
}

template<typename T, typename U>
void fillDesIdSets(RootedTree<T, U> & tree) {
	// assumes OttId is set for each tip
	for (auto node : PostorderIter<RootedTree<T, U> >(tree)) {
		std::set<long> & desIds = node->getData().desIds;
		if (node->isTip()) {
			desIds.insert(node->getOttId());
		} else {
			for (auto child : ChildIter<RootedTreeNode<T> >(*node)) {
				std::set<long> & cDesIds = child->getData().desIds;
				desIds.insert(cDesIds.begin(), cDesIds.end());
			}
		}
	}
}

// uses ottID->node mapping, but the split sets of the nodes
template<typename T, typename U>
RootedTreeNode<T> * findMRCAFromIDSet(RootedTree<T, U> & tree, const std::set<long> & idSet, long trigger) {
	typedef RootedTreeNode<T> NT_t;
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

template<typename T, typename U>
std::size_t checkForUnknownTaxa(std::ostream & err, const RootedTree<T, U> & toCheck, const RootedTree<T, U> & taxonomy) {
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

template<typename T, typename U>
inline RootedTreeNode<T> * addChildForOttId(RootedTreeNode<T> & nd, long ottId, RootedTree<T, U> & tree) {
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
	for (auto anc : AncNodeIter<RootedTreeNode<RTSplits> >(&nd)) {
		assert(anc != nullptr);
		assert(!anc->getData().desIds.empty());
		anc->getData().desIds.erase(begin(toRemove), end(toRemove));
		anc->getData().desIds.insert(begin(ls), end(ls));
	}
}

template<typename T, typename U>
std::vector<RootedTreeNode<T> *> expandOTTInternalsWhichAreLeaves(RootedTree<T, U> & toExpand, const RootedTree<T, U> & taxonomy) {
	const U & taxData = taxonomy.getData();
	std::map<RootedTreeNode<T> *, std::set<long> > replaceNodes;
	for (auto nd : LeafIter<RootedTree<T, U> >(toExpand)) {
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
	std::vector<RootedTreeNode<T> *> expanded;
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

template<typename T, typename U>
void markPathToRoot(const RootedTree<T, U> & fullTree,
					long ottId,
					std::map<const RootedTreeNode<T> *,
					std::set<long> > &n2m){
	auto startNd = fullTree.getData().getNodeForOttId(ottId);
	assert(startNd != nullptr);
	if (startNd == nullptr) {
		std::string m = "OTT id not found ";
		m += std::to_string(ottId);
		throw OTCError(m);
	}
	n2m[startNd].insert(ottId);
	for (auto nd : AncNodeIter<RootedTreeNode<T> >(startNd)) {
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
inline bool multipleChildrenInMap(const RootedTreeNode<T> & nd,
								  std::map<const RootedTreeNode<T> *, U> markedMap,
								  const RootedTreeNode<T> **first) {
	assert(first);
	bool foundFirst = false;
	*first = nullptr;
	for(auto c : ConstChildIter<RootedTreeNode<T> >(nd)) {
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
inline const RootedTreeNode<T> * findNextSignificantNode(const RootedTreeNode<T> * node,
														 const std::map<const RootedTreeNode<T> *, std::set<long> > & markedMap) {
	assert(node != nullptr);
	auto currNode = node;
	for (;;) {
		const RootedTreeNode<T> * sc;
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
											const RootedTreeNode<T> & srcNd,
											const std::map<const RootedTreeNode<T> *, std::set<long> > & markedMap) {
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
		for (auto child : ConstChildIter<RootedTreeNode<T> >(*nsn)) {
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
inline void describeUnnamedNode(const RootedTreeNode<T> & nd,
								std::ostream & out,
								unsigned int anc,
								bool useNdNames) {
	if (useNdNames && !nd.getName().empty()) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before \"" << nd.getName() << "\"";
		} else {
			out << "the node \"" << nd.getName() << "\"";
		}
	}
	else if (nd.isTip()) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before the leaf \"" << nd.getName()  << "\"";
		} else {
			out << "the leaf \"" << nd.getName()  << "\"";
		}
	} else if (nd.isOutDegreeOneNode()) {
		describeUnnamedNode(*nd.getFirstChild(), out, anc + 1, useNdNames);
	} else {
		auto left = findLeftmostInSubtree(&nd)->getName();
		auto right = findRightmostInSubtree(&nd)->getName();
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right <<'\"';
		} else {
			out <<  "MRCA of \"" << left << "\" and " << "\"" << right <<'\"';
		}
	}
	out << std::endl;
}

} // namespace otc
#endif

