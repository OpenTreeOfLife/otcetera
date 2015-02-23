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
	for (auto node : PostorderNode<T, U>(tree)) {
		std::set<long> & desIds = node->getData().desIds;
		if (node->isTip()) {
			desIds.insert(node->getOttId());
		} else {
			for (auto child : ChildIterator<T, U>(*node)) {
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
	nd.getData().desIds = ls;
	for (auto anc : AncNodeIter<RTSplits>(&nd)) {
		anc->getData().desIds.erase(begin(toRemove), end(toRemove));
		anc->getData().desIds.insert(begin(ls), end(ls));
	}
}

template<typename T, typename U>
std::vector<RootedTreeNode<T> *> expandOTTInternalsWhichAreLeaves(RootedTree<T, U> & toExpand, const RootedTree<T, U> & taxonomy) {
	const U & taxData = taxonomy.getData();
	std::map<RootedTreeNode<T> *, std::set<long> > replaceNodes;
	for (auto nd : LeafNodeIter<T, U>(toExpand)) {
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


} // namespace otc
#endif

