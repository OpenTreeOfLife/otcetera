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
		const RootedTreeNode<T> * GetParent() const {
			return parent;
		}
		bool IsTip() const {
			return (lChild == nullptr);
		}
		const RootedTreeNode<T> * GetFirstChild() const {
			return lChild;
		}
		RootedTreeNode<T> * GetNextSib() const {
			return rSib;
		}
		RootedTreeNode<T> * GetLastChild() const {
			auto currNode = this->GetFirstChild();
			if (currNode == nullptr)
				return nullptr;
			return currNode->GetLastSib();
		}
		RootedTreeNode<T> * GetLastSib() const {
			auto currNode = this->GetNextSib();
			if (currNode == nullptr) {
				return nullptr;
			}
			auto nextNd = currNode->GetNextSib();
			while (nextNd != nullptr) {
				currNode = nextNd;
				nextNd = currNode->GetNextSib();
			}
			return currNode;
		}
		std::vector<const RootedTreeNode<T> *> GetChildren() const {
			std::vector<const RootedTreeNode<T> *> children;
			auto * currNode = GetFirstChild();
			while(currNode) {
				children.push_back(currNode);
				currNode = currNode->GetNextSib();
			}
			return children;
		}
		unsigned GetOutDegree() const {
			unsigned n = 0;
			auto currNode = GetFirstChild();
			while(currNode) {
				n += 1;
				currNode = currNode->GetNextSib();
			}
			return n;
		}
		long GetOTUId() const {
			return otuID;
		}
		void GetOTUId(long i) {
			otuID = i;
		}
		// non-empty only for internals that are labelled with names that are NOT taxLabels
		const namestring_t & GetName() const {
			return name;
		}
		void SetName(const namestring_t &n) {
			name = n;
		}
		const T & GetData() const {
			return data;
		}
		T & GetData() {
			return data;
		}
		RootedTreeNode<T>(RootedTreeNode<T> *par)
			:lChild(nullptr),
			rSib(nullptr),
			parent(par),
			otuID(INT_MAX) {
		}
		void AddSib(RootedTreeNode<T> *n) {
			if (rSib) {
				rSib->AddSib(n);
			} else {
				rSib = n;
			}
		}
		void AddChild(RootedTreeNode<T> *n) {
			if (lChild)
				lChild->AddSib(n);
			else
				lChild = n;
		}
		bool RemoveChild(RootedTreeNode<T> *n) {
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
		void WriteAsNewick(std::ostream &out,
						   bool useLeafNames,
						   const std::map<RootedTreeNode<T> *, namestring_t> *nd2name=nullptr) const;
		void AddSelfAndDesToPreorder(std::vector<const RootedTreeNode<T> *> &p) const;

		void LowLevelSetFirstChild(RootedTreeNode<T> *nd) {
			lChild = nd;
		}
		void LowLevelSetNextSib(RootedTreeNode<T> *nd) {
			rSib = nd;
		}
	private:
		RootedTreeNode<T> * lChild;
		RootedTreeNode<T> * rSib;
		RootedTreeNode<T> * parent;
		namestring_t name; // non-empty only for internals that are labelled with names that are NOT taxLabels
		unsigned otuID; // present for every leaf. UINT_MAX for internals labeled with taxlabels
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
			Clear();
		}
		std::vector<const RootedTreeNode<T> *> GetPreorderTraversal() const;
		const std::vector<const RootedTreeNode<T> *> & GetLeavesRef() {
			return leaves;
		}
		void WriteAsNewick(std::ostream &out,
						   bool nhx,
						   bool useLeafNames,
						   const std::map<RootedTreeNode<T> *, namestring_t> *nd2name=nullptr) const {
			if (root) {
				root->WriteAsNewick(out, nhx, useLeafNames, nd2name);
			}
		}
		const RootedTreeNode<T> * GetRoot() const {
			return root;
		}
		RootedTreeNode<T> * CreateRoot() {
			if (root != nullptr) {
				Clear();
			}
			this->root = this->AllocNewNode(nullptr);
			return this->root;
		}
		RootedTreeNode<T> * CreateChild(RootedTreeNode<T> *par) {
			auto c = this->AllocNewNode(par);
			par->AddChild(c);
			return c;
		}
		RootedTreeNode<T> * CreateSib(RootedTreeNode<T> *leftSib) {
			assert(leftSib->parent != nullptr);
			auto s = this->AllocNewNode(leftSib->parent);
			leftSib->AddSib(s);
			return s;
		}
		U & GetData() {
			return this->data;
		}
		const U & GetData() const {
			return this->data;
		}
	protected:
		std::vector<RootedTreeNode<T> *> allNodes;
		std::vector<RootedTreeNode<T> *> leaves;
		RootedTreeNode<T> * root;
		U data;
	public:
		RootedTreeNode<T> * AllocNewNode(RootedTreeNode<T> *p) {
			RootedTreeNode<T> * nd = new RootedTreeNode<T>(p);
			allNodes.push_back(nd);
			return nd;
		}
		void Clear() {
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

/** simple phylo functions */

template<typename T, typename U>
unsigned int countPolytomies(RootedTree<T, U> & tree);

template<typename T, typename U>
unsigned int countPolytomies(RootedTree<T, U> & tree) {
	return 0;
}

struct RTNodeNoData{};
struct RTreeNoData{};

} // namespace otc
#endif

