#ifndef PHARSE_TREE_H
#define PHARSE_TREE_H
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <climits>

namespace pharse {
template<typename T> class TreeNode;
template<typename T>
class TreeEdge {
	public:
		typedef class TreeNode<T> Node_t;
		const Node_t * GetParent() const {
			return parent;
		}
		const Node_t * GetChild() const {
			return child;
		}
		T & GetData() const {
			return child->GetData();
		}
	private:
		void SetParent(Node_t * p) {
			this->parent = p;
		}
		void WriteAsNewick(std::ostream &out, bool nhx) const;
		TreeEdge<T>(Node_t * par, Node_t  * des)
			:parent(par),
			child(des) {
		}
		Node_t * parent;
		Node_t * child;
		friend class TreeNode<T>;
};

template<typename T>
class TreeNode {
	public:
		typedef class TreeEdge<T> Edge_t;
		const Edge_t & GetEdgeToParent() const {
			return edgeToPar;
		}
		const TreeNode<T> * GetParent() const {
			return edgeToPar.GetParent;
		}
		bool IsTip() const {
			return (lChild == nullptr);
		}
		const TreeNode<T> * GetFirstChild() const {
			return lChild;
		}
		TreeNode<T> * GetNextSib() const {
			return rSib;
		}
		TreeNode<T> * GetLastChild() const {
			auto currNode = this->GetFirstChild();
			if (currNode == nullptr)
				return nullptr;
			return currNode->GetLastSib();
		}
		TreeNode<T> * GetLastSib() const {
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
		std::vector<const TreeNode<T> *> GetChildren() const {
			std::vector<const TreeNode<T> *> children;
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
		const std::string & GetName() const {
			return name;
		}
		void SetName(const std::string &n) {
			name = n;
		}
		T & GetData() const {
			return &data;
		}
		TreeNode<T>(TreeNode<T> *par)
			:lChild(0L),
			rSib(0L),
			edgeToPar(par, 0L),
			otuID(INT_MAX) {
			edgeToPar.child = this;
		}

		void AddSib(TreeNode<T> *n) {
			if (rSib) {
				rSib->AddSib(n);
			} else {
				rSib = n;
			}
		}
		void AddChild(TreeNode<T> *n) {
			if (lChild)
				lChild->AddSib(n);
			else
				lChild = n;
		}

		bool RemoveChild(TreeNode<T> *n) {
			if (n == 0L || lChild == 0L) {
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
					if (c->rSib == 0L) {
						return false;
					}
				}
			}
			n->edgeToPar.parent = 0L;
			return true;
		}

	public:
		void WriteAsNewick(std::ostream &out,
						   bool useLeafNames,
						   const std::map<TreeNode<T> *, std::string> *nd2name=0L) const;
		void AddSelfAndDesToPreorder(std::vector<const TreeNode<T> *> &p) const;

		void LowLevelSetFirstChild(TreeNode<T> *nd) {
			lChild = nd;
		}
		void LowLevelSetNextSib(TreeNode<T> *nd) {
			rSib = nd;
		}
	private:
		TreeNode<T> * lChild;
		TreeNode<T> * rSib;
		TreeEdge<T> edgeToPar;
		std::string name; // non-empty only for internals that are labelled with names that are NOT taxLabels
		unsigned otuID; // present for every leaf. UINT_MAX for internals labeled with taxlabels
		T data;
		//friend class NxsSimpleTree;
};

} // namespace pharse
#endif

