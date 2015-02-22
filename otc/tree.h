#ifndef OTCETERA_TREE_H
#define OTCETERA_TREE_H

#include <climits>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename, typename> class RootedTree;

typedef std::string namestring_t;

template<typename T>
class RootedTreeNode {
	public:
		const RootedTreeNode<T> * getParent() const {
			return parent;
		}
		RootedTreeNode<T> * getParent() {
			return parent;
		}
		bool IsTip() const {
			return (lChild == nullptr);
		}
		const RootedTreeNode<T> * getFirstChild() const {
			return lChild;
		}
		RootedTreeNode<T> * getFirstChild() {
			return lChild;
		}
		RootedTreeNode<T> * getNextSib() const {
			return rSib;
		}
		RootedTreeNode<T> * getLastChild() const {
			auto currNode = this->getFirstChild();
			if (currNode == nullptr)
				return nullptr;
			return currNode->getLastSib();
		}
		RootedTreeNode<T> * getLastSib() const {
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
		std::vector<const RootedTreeNode<T> *> getChildren() const {
			std::vector<const RootedTreeNode<T> *> children;
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
		void addSib(RootedTreeNode<T> *n) {
			if (rSib) {
				rSib->addSib(n);
			} else {
				rSib = n;
			}
		}
		void addChild(RootedTreeNode<T> *n) {
			if (lChild)
				lChild->addSib(n);
			else
				lChild = n;
		}
		bool removeChild(RootedTreeNode<T> *n) {
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
	public:
		void writeAsNewick(std::ostream &out,
						   bool useLeafNames,
						   const std::map<RootedTreeNode<T> *, namestring_t> *nd2name=nullptr) const;
		void addSelfAndDesToPreorder(std::vector<const RootedTreeNode<T> *> &p) const;

		void lowLevelSetFirstChild(RootedTreeNode<T> *nd) {
			lChild = nd;
		}
		void lowLevelSetNextSib(RootedTreeNode<T> *nd) {
			rSib = nd;
		}
	private:
		RootedTreeNode<T> * lChild;
		RootedTreeNode<T> * rSib;
		RootedTreeNode<T> * parent;
		namestring_t name; // non-empty only for internals that are labelled with names that are NOT taxLabels
		long ottId; // present for every leaf. UINT_MAX for internals labeled with taxlabels
		T data;
	private:
		RootedTreeNode<T>(const RootedTreeNode<T> &); //not defined.  Not copyable
		RootedTreeNode<T> & operator=(const RootedTreeNode<T> &); //not defined.  Not copyable

		template<typename Y, typename Z>
		friend class RootedTree;
};

template<typename T, typename U>
class RootedTree {
	public:
		RootedTree<T, U>()
			:root(nullptr) {
		}
		~RootedTree<T, U>() {
			clear();
		}
		std::vector<const RootedTreeNode<T> *> getPreorderTraversal() const;
		const std::vector<const RootedTreeNode<T> *> & getLeavesRef() {
			return leaves;
		}
		void writeAsNewick(std::ostream &out,
						   bool nhx,
						   bool useLeafNames,
						   const std::map<RootedTreeNode<T> *, namestring_t> *nd2name=nullptr) const {
			if (root) {
				root->writeAsNewick(out, nhx, useLeafNames, nd2name);
			}
		}
		const RootedTreeNode<T> * getRoot() const {
			return root;
		}
		RootedTreeNode<T> * getRoot() {
			return root;
		}
		RootedTreeNode<T> * createRoot() {
			if (root != nullptr) {
				clear();
			}
			this->root = this->allocNewNode(nullptr);
			return this->root;
		}
		RootedTreeNode<T> * createChild(RootedTreeNode<T> *par) {
			auto c = this->allocNewNode(par);
			par->addChild(c);
			return c;
		}
		RootedTreeNode<T> * createSib(RootedTreeNode<T> *leftSib) {
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
		std::vector<RootedTreeNode<T> *> allNodes;
		std::vector<RootedTreeNode<T> *> leaves;
		RootedTreeNode<T> * root;
		U data;
	public:
		RootedTreeNode<T> * allocNewNode(RootedTreeNode<T> *p) {
			RootedTreeNode<T> * nd = new RootedTreeNode<T>(p);
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
	private:
		RootedTree<T, U>(const RootedTree<T, U> &); //not defined.  Not copyable
		RootedTree<T, U> & operator=(const RootedTree<T, U> &); //not defined.  Not copyable
};

struct RTNodeNoData{};
struct RTreeNoData{};

typedef RootedTree<RTNodeNoData, RTreeNoData> RootedTreeTopologyNoData;

} // namespace otc
#endif

