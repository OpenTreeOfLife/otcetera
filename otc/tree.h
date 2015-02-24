#ifndef OTCETERA_TREE_H
#define OTCETERA_TREE_H

#include <climits>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename, typename> class RootedTree;

typedef std::string namestring_t;

template<typename T>
class RootedTreeNode {
	public:
		using node_type = RootedTreeNode<T>;
		using data_type = T;
		bool isTip() const {
			return (lChild == nullptr);
		}
		const node_type * getParent() const {
			return parent;
		}
		node_type * getParent() {
			return parent;
		}
		const node_type * getFirstChild() const {
			return lChild;
		}
		node_type * getFirstChild() {
			return lChild;
		}
		const node_type * getNextSib() const {
			return rSib;
		}
		node_type * getNextSib() {
			return rSib;
		}
		const node_type * getLastChild() const {
			if (lChild == nullptr)
				return nullptr;
			return (lChild->rSib == nullptr ? lChild : lChild->getLastSib());
		}
		node_type * getLastChild() {
			return const_cast<node_type *>(const_cast<const node_type *>(this)->getLastChild());
		}
		const node_type * getLastSib() const {
			auto currNode = this->getNextSib();
			if (currNode == nullptr) {
				return nullptr;
			}
			auto nextNd = currNode->getNextSib();
			while (nextNd != nullptr) {
				currNode = nextNd;
				nextNd = currNode->getNextSib();
			}
			return currNode;
		}
		node_type * getLastSib() {
			return const_cast<node_type *>(const_cast<const node_type *>(this)->getLastSib());
		}
		
		std::vector<const node_type *> getChildren() const {
			std::vector<const node_type *> children;
			auto * currNode = getFirstChild();
			while(currNode) {
				children.push_back(currNode);
				currNode = currNode->getNextSib();
			}
			return children;
		}
		unsigned getOutDegree() const {
			unsigned n = 0;
			auto currNode = getFirstChild();
			while(currNode) {
				n += 1;
				currNode = currNode->getNextSib();
			}
			return n;
		}
		bool hasOttId() const {
			return ottId != LONG_MAX;
		}
		long getOttId() const {
			return ottId;
		}
		void setOttId(long i) {
			ottId = i;
		}
		// non-empty only for internals that are labelled with names that are NOT taxLabels
		const namestring_t & getName() const {
			return name;
		}
		void setName(const namestring_t &n) {
			name = n;
		}
		const T & getData() const {
			return data;
		}
		T & getData() {
			return data;
		}
		RootedTreeNode<T>(RootedTreeNode<T> *par)
			:lChild(nullptr),
			rSib(nullptr),
			parent(par),
			ottId(LONG_MAX) {
		}
		void addSib(node_type *n) {
			if (rSib) {
				rSib->addSib(n);
			} else {
				rSib = n;
			}
		}
		void addChild(node_type *n) {
			if (lChild)
				lChild->addSib(n);
			else
				lChild = n;
		}
		bool removeChild(node_type *n) {
			if (n == nullptr || lChild == nullptr) {
				return false;
			}
			if (lChild == n) {
				lChild = lChild->rSib;
			} else {
				auto c = lChild;
				for (;;) {
					if (c->rSib == n) {
						c->rSib = n->rSib;
						break;
					}
					if (c->rSib == nullptr) {
						return false;
					}
				}
			}
			n->parent = nullptr;
			return true;
		}
		bool isOutDegreeOneNode() const {
			return (lChild != nullptr) && (lChild->rSib == nullptr);
		}
		bool includesOnlyOneLeaf() const {
			if (isTip()) {
				return true;
			}
			return isOutDegreeOneNode() && lChild->includesOnlyOneLeaf();
		}

	public:
		void writeAsNewick(std::ostream &out,
						   bool useLeafNames,
						   const std::map<node_type *, namestring_t> *nd2name=nullptr) const;
		void addSelfAndDesToPreorder(std::vector<const node_type *> &p) const;

		void lowLevelSetFirstChild(node_type *nd) {
			lChild = nd;
		}
		void lowLevelSetNextSib(node_type *nd) {
			rSib = nd;
		}
	private:
		node_type * lChild;
		node_type * rSib;
		node_type * parent;
		namestring_t name; // non-empty only for internals that are labelled with names that are NOT taxLabels
		long ottId; // present for every leaf. UINT_MAX for internals labeled with taxlabels
		T data;
	private:
		RootedTreeNode<T>(const RootedTreeNode<T> &) = delete;
		RootedTreeNode<T> & operator=(const RootedTreeNode<T> &) = delete;
		template<typename Y, typename Z>
		friend class RootedTree;
};

template<typename T, typename U>
class RootedTree {
	public:
		using node_data_type = T;
		using node_type = RootedTreeNode<T>;
		using data_type = U;
		
		RootedTree<T, U>()
			:root(nullptr) {
		}
		~RootedTree<T, U>() {
			clear();
		}
		std::vector<const node_type *> getPreorderTraversal() const;
		const std::vector<const node_type *> & getLeavesRef() {
			return leaves;
		}
		void writeAsNewick(std::ostream &out,
						   bool nhx,
						   bool useLeafNames,
						   const std::map<node_type *, namestring_t> *nd2name=nullptr) const {
			if (root) {
				root->writeAsNewick(out, nhx, useLeafNames, nd2name);
			}
		}
		const node_type * getRoot() const {
			return root;
		}
		node_type * getRoot() {
			return root;
		}
		node_type * createRoot() {
			if (root != nullptr) {
				clear();
			}
			this->root = this->allocNewNode(nullptr);
			return this->root;
		}
		node_type * createChild(node_type *par) {
			auto c = this->allocNewNode(par);
			par->addChild(c);
			return c;
		}
		node_type * createSib(node_type *leftSib) {
			assert(leftSib->parent != nullptr);
			auto s = this->allocNewNode(leftSib->parent);
			leftSib->addSib(s);
			return s;
		}
		U & getData() {
			return this->data;
		}
		const U & getData() const {
			return this->data;
		}
	protected:
		std::vector<node_type *> allNodes;
		std::vector<node_type *> leaves;
		node_type * root;
		U data;
	public:
		node_type * allocNewNode(node_type *p) {
			node_type * nd = new node_type(p);
			allNodes.push_back(nd);
			return nd;
		}
		void clear() {
			root = NULL;
			leaves.clear();
			for (auto nIt : allNodes) {
				delete nIt;
			}
			allNodes.clear();
		}
		std::set<const node_type *> getSetOfAllNodes() const {
			std::set<const node_type *> r;
			for (auto nd : allNodes) {
				r.insert(nd);
			}
			return r;
		}
	private:
		RootedTree<T, U>(const RootedTree<T, U> &) = delete;
		RootedTree<T, U> & operator=(const RootedTree<T, U> &) = delete;
};

struct RTNodeNoData{};
struct RTreeNoData{};

typedef RootedTree<RTNodeNoData, RTreeNoData> RootedTreeTopologyNoData;

} // namespace otc
#endif

