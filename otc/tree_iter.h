#ifndef OTCETERA_TREE_ITER_H
#define OTCETERA_TREE_ITER_H
// Iterators for traversing trees
// Depends on: tree.h tree_util.h
// Depended on by: tree_iter.h tree_operation.h 

#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/tree_util.h"

namespace otc {

template<typename T>
using NdFilterFn = std::function<bool(const RootedTreeNode<T> &)>;

/// visits ancestors before descendants
/// if filterFn is supplied, then only the nodes associated with a true
///		response will be returned (but the descendants of a "false" node
///		will still be visited.
template<typename T>
class const_preorder_iterator : std::forward_iterator_tag {
	private:
		const NdFilterFn<T> filterFn;
		const RootedTreeNode<T> * curr;
		bool movingDown;

		void _advance() {
			do {
				assert(curr != nullptr);
				if (movingDown) {
					while (movingDown) {
						if (curr->GetNextSib()) {
							curr = curr->GetNextSib();
							movingDown = false;
						} else {
							curr = curr->GetParent();
							if (curr == nullptr) {
								return;
							}
						}
					}
				} else if (curr->IsTip()) {
					if (curr->GetNextSib()) {
						curr = curr->GetNextSib();
					} else {
						movingDown = true;
						_advance();
						if (curr == nullptr) {
							return;
						}
					}
				} else {
					curr = curr->GetFirstChild();
				}
			} while ((filterFn != nullptr && !filterFn(*curr)));
		}
	public:
		const_preorder_iterator(const RootedTreeNode<T> *c)
			:filterFn{nullptr},
			curr(c),
			movingDown(false) {
		}
		const_preorder_iterator(const RootedTreeNode<T> *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(c),
			movingDown(false) {
			if (c != nullptr && filterFn && !filterFn(*curr)) {
				_advance();
			}
		}
		bool operator==(const const_preorder_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const const_preorder_iterator &other) {
			return this->curr != other.curr;
		}
		const RootedTreeNode<T> * operator*() const {
			return curr;
		}
		const_preorder_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead const_preorder_iterator");
			}
			_advance();
			return *this;
		}
};


/// descendants before ancestors, but not guaranteed to be the reverse of const_preorder_iterator
template<typename T>
class child_iterator : std::forward_iterator_tag {
	private:
		typedef RootedTreeNode<T> NodeType;
		NdFilterFn<T> filterFn;
		NodeType * curr;

		void _advance() {
			assert(curr != nullptr);
			curr = curr->GetNextSib();
		}
	public:
		child_iterator(RootedTreeNode<T> *c)
			:filterFn{nullptr},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->GetFirstChild();
			}
		}
		child_iterator(RootedTreeNode<T> *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->GetFirstChild();
			}
			if (c != nullptr && filterFn && !filterFn(*c)) {
				_advance();
			}
		}
		bool operator==(child_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(child_iterator &other) {
			return this->curr != other.curr;
		}
		RootedTreeNode<T> * operator*() const {
			return curr;
		}
		child_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead child_iterator");
			}
			_advance();
			return *this;
		}
};

/// descendants before ancestors, but not guaranteed to be the reverse of const_preorder_iterator
template<typename T>
class const_postorder_iterator : std::forward_iterator_tag {
	private:
		typedef const RootedTreeNode<T> NodeType;
		const NdFilterFn<T> filterFn;
		const NodeType * curr;
		const NodeType * lastNode;

		void _advance() {
			if (curr == lastNode) {
				curr = nullptr;
			} else {
				auto n = curr->GetNextSib();
				curr = (n == nullptr ? curr->GetParent() : findLeftmostInSubtree<NodeType>(n));
			}
		}
	public:
		const_postorder_iterator(const RootedTreeNode<T> *c)
			:filterFn{nullptr},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<NodeType>(lastNode);
			}
		}
		const_postorder_iterator(const RootedTreeNode<T> *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<NodeType>(lastNode);
			}
			if (c != nullptr && filterFn && !filterFn(*lastNode)) {
				_advance();
			}
		}
		bool operator==(const const_postorder_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const const_postorder_iterator &other) {
			return this->curr != other.curr;
		}
		const RootedTreeNode<T> * operator*() const {
			return curr;
		}
		const_postorder_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead const_postorder_iterator");
			}
			_advance();
			return *this;
		}
};

/// descendants before ancestors, but not guaranteed to be the reverse of const_preorder_iterator
template<typename T>
class postorder_iterator : std::forward_iterator_tag {
	private:
		typedef RootedTreeNode<T> NodeType;
		const NdFilterFn<T> filterFn;
		NodeType * curr;
		NodeType * lastNode;

		void _advance() {
			if (curr == lastNode) {
				curr = nullptr;
			} else {
				auto n = curr->GetNextSib();
				curr = (n == nullptr ? curr->GetParent() : findLeftmostInSubtree<NodeType>(n));
			}
		}
	public:
		postorder_iterator(RootedTreeNode<T> *c)
			:filterFn{nullptr},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<NodeType>(lastNode);
			}
		}
		postorder_iterator(RootedTreeNode<T> *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<NodeType>(lastNode);
			}
			if (c != nullptr && filterFn && !filterFn(*lastNode)) {
				_advance();
			}
		}
		bool operator==(const postorder_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const postorder_iterator &other) {
			return this->curr != other.curr;
		}
		RootedTreeNode<T> * operator*() const {
			return curr;
		}
		postorder_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead postorder_iterator");
			}
			_advance();
			return *this;
		}
};

template<typename T, typename U>
class ChildIterator {
	public:
	explicit ChildIterator(RootedTreeNode<T> & n)
		:node(n) {
	}
	child_iterator<T> begin() const {
		return std::move(child_iterator<T>{&node});
	}
	child_iterator<T> end() const {
		return child_iterator<T>{nullptr};
	}
	private:
		RootedTreeNode<T> & node;
};

template<typename T, typename U>
class ConstPreorderInternalNode {
	public:
	explicit ConstPreorderInternalNode(const RootedTree<T, U> &t)
		:tree(t){
	}
	const_preorder_iterator<T> begin() const {
		return std::move(const_preorder_iterator<T>{tree.GetRoot(), isInternalNode<T>});
	}
	const_preorder_iterator<T> end() const {
		return const_preorder_iterator<T>{nullptr};
	}
	private:
		const RootedTree<T, U> & tree;
};

template<typename T, typename U>
class ConstPostorderInternalNode {
	public:
	explicit ConstPostorderInternalNode(const RootedTree<T, U> &t)
		:tree(t){
	}
	const_postorder_iterator<T> begin() const {
		return std::move(const_postorder_iterator<T>{tree.GetRoot(), isInternalNode<T>});
	}
	const_postorder_iterator<T> end() const {
		return const_postorder_iterator<T>{nullptr};
	}
	private:
		const RootedTree<T, U> & tree;
};

template<typename T, typename U>
class PostorderInternalNode {
	public:
	explicit PostorderInternalNode(RootedTree<T, U> & t)
		:tree(t){
	}
	postorder_iterator<T> begin() const {
		return std::move(postorder_iterator<T>{tree.GetRoot(), isInternalNode<T>});
	}
	postorder_iterator<T> end() const {
		return postorder_iterator<T>{nullptr};
	}
	private:
		RootedTree<T, U> & tree;
};

} // namespace otc
#endif

