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
using NdFilterFn = std::function<bool(const T &)>;

/// visits ancestors before descendants
/// if filterFn is supplied, then only the nodes associated with a true
///		response will be returned (but the descendants of a "false" node
///		will still be visited.
template<typename T>
class const_preorder_iterator : std::forward_iterator_tag {
	private:
		const NdFilterFn<T> filterFn;
		const T * curr;
		bool movingDown;

		void _advance() {
			do {
				assert(curr != nullptr);
				if (movingDown) {
					while (movingDown) {
						if (curr->getNextSib()) {
							curr = curr->getNextSib();
							movingDown = false;
						} else {
							curr = curr->getParent();
							if (curr == nullptr) {
								return;
							}
						}
					}
				} else if (curr->isTip()) {
					if (curr->getNextSib()) {
						curr = curr->getNextSib();
					} else {
						movingDown = true;
						_advance();
						if (curr == nullptr) {
							return;
						}
					}
				} else {
					curr = curr->getFirstChild();
				}
			} while ((filterFn != nullptr && !filterFn(*curr)));
		}
	public:
		const_preorder_iterator(const T *c)
			:filterFn{nullptr},
			curr(c),
			movingDown(false) {
		}
		const_preorder_iterator(const T *c, NdFilterFn<T> f)
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
		const T * operator*() const {
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
class const_child_iterator : std::forward_iterator_tag {
	private:
		NdFilterFn<T> filterFn;
		const T * curr;

		void _advance() {
			assert(curr != nullptr);
			curr = curr->getNextSib();
		}
	public:
		const_child_iterator(const T *c)
			:filterFn{nullptr},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->getFirstChild();
			}
		}
		const_child_iterator(const T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->getFirstChild();
			}
			if (c != nullptr && filterFn && !filterFn(*c)) {
				_advance();
			}
		}
		bool operator==(const_child_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const_child_iterator &other) {
			return this->curr != other.curr;
		}
		const T * operator*() const {
			return curr;
		}
		const_child_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead const_child_iterator");
			}
			_advance();
			return *this;
		}
};
/// descendants before ancestors, but not guaranteed to be the reverse of const_preorder_iterator
template<typename T>
class child_iterator : std::forward_iterator_tag {
	private:
		NdFilterFn<T> filterFn;
		T * curr;

		void _advance() {
			assert(curr != nullptr);
			curr = curr->getNextSib();
		}
	public:
		child_iterator(T *c)
			:filterFn{nullptr},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->getFirstChild();
			}
		}
		child_iterator(T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr) {
			if (c != nullptr) {
				curr = c->getFirstChild();
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
		T * operator*() const {
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
		const NdFilterFn<T> filterFn;
		const T * curr;
		const T * lastNode;

		void _advance() {
			if (curr == lastNode) {
				curr = nullptr;
			} else {
				for (;;) {
					auto n = curr->getNextSib();
					curr = (n == nullptr ? curr->getParent() : findLeftmostInSubtree<const T>(n));
					if (filterFn == nullptr || filterFn(*curr)) {
						break;
					}
					if (curr == lastNode) {
						curr = nullptr;
						break;
					}
				}
			}
		}
	public:
		const_postorder_iterator(const T *c)
			:filterFn{nullptr},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<const T>(lastNode);
			}
		}
		const_postorder_iterator(const T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<const T>(lastNode);
			}
			if (curr != nullptr && filterFn && !filterFn(*curr)) {
				_advance();
			}
		}
		bool operator==(const const_postorder_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const const_postorder_iterator &other) {
			return this->curr != other.curr;
		}
		const T * operator*() const {
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
		const NdFilterFn<T> filterFn;
		T * curr;
		T * lastNode;

		void _advance() {
			if (curr == lastNode) {
				curr = nullptr;
			} else {
				for (;;) {
					auto n = curr->getNextSib();
					curr = (n == nullptr ? curr->getParent() : findLeftmostInSubtree<T>(n));
					if (filterFn == nullptr || filterFn(*curr)) {
						break;
					}
					if (curr == lastNode) {
						curr = nullptr;
						break;
					}
				}
			}
		}
	public:
		postorder_iterator(T *c)
			:filterFn{nullptr},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<T>(lastNode);
			}
		}
		postorder_iterator(T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtree<T>(lastNode);
			}
			if (curr != nullptr && filterFn && !filterFn(*curr)) {
				_advance();
			}
		}
		bool operator==(const postorder_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const postorder_iterator &other) {
			return this->curr != other.curr;
		}
		T * operator*() const {
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

template<typename T>
class anc_iterator : std::forward_iterator_tag {
	private:
		const NdFilterFn<T> filterFn;
		T * curr;
		
		void _advance() {
			assert(curr != nullptr);
			do {
				curr = curr->getParent();
			} while(curr != nullptr && filterFn != nullptr && !filterFn(curr));
		}
	public:
		anc_iterator(T *c)
			:filterFn{nullptr},
			curr(c) {
			if (curr != nullptr) {
				_advance();
			}
		}
		anc_iterator(T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(c) {
			if (curr != nullptr) {
				_advance();
			}
		}
		bool operator==(const anc_iterator &other) {
			return this->curr == other.curr;
		}
		bool operator!=(const anc_iterator &other) {
			return this->curr != other.curr;
		}
		T * operator*() const {
			return curr;
		}
		anc_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead anc_iterator");
			}
			_advance();
			return *this;
		}
};

template<typename T>
class ChildIter {
	public:
	explicit ChildIter(T & n)
		:node(n) {
	}
	child_iterator<T> begin() const {
		return std::move(child_iterator<T>{&node});
	}
	child_iterator<T> end() const {
		return child_iterator<T>{nullptr};
	}
	private:
		T & node;
};

template<typename T>
class ConstChildIter {
	public:
	explicit ConstChildIter(const T & n)
		:node(n) {
	}
	const_child_iterator<T> begin() const {
		return std::move(const_child_iterator<T>{&node});
	}
	const_child_iterator<T> end() const {
		return const_child_iterator<T>{nullptr};
	}
	private:
		const T & node;
};

template<typename T>
class ConstPreorderInternalNode {
	public:
	explicit ConstPreorderInternalNode(const T &t)
		:tree(t){
	}
	const_preorder_iterator<typename T::node_type> begin() const {
		return std::move(const_preorder_iterator<typename T::node_type>{tree.getRoot(), isInternalNode<typename T::node_type>});
	}
	const_preorder_iterator<typename T::node_type> end() const {
		return const_preorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		const T & tree;
};

template<typename T>
class ConstPreorder {
	public:
	explicit ConstPreorder(const T &t)
		:tree(t){
	}
	const_preorder_iterator<typename T::node_type> begin() const {
		return std::move(const_preorder_iterator<typename T::node_type>{tree.getRoot()});
	}
	const_preorder_iterator<typename T::node_type> end() const {
		return const_preorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		const T & tree;
};


template<typename T>
class ConstPostorderInternalNode {
	public:
	explicit ConstPostorderInternalNode(const T &t)
		:tree(t){
	}
	const_postorder_iterator<typename T::node_type> begin() const {
		return std::move(const_postorder_iterator<typename T::node_type>{tree.getRoot(), isInternalNode<typename T::node_type>});
	}
	const_postorder_iterator<typename T::node_type> end() const {
		return const_postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		const T & tree;
};

template<typename T>
class ConstPostorderIter {
	public:
	explicit ConstPostorderIter(const T &t)
		:tree(t){
	}
	const_postorder_iterator<typename T::node_type> begin() const {
		return std::move(const_postorder_iterator<typename T::node_type>{tree.getRoot()});
	}
	const_postorder_iterator<typename T::node_type> end() const {
		return const_postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		const T & tree;
};

template<typename T>
class PostorderInternalNode {
	public:
	explicit PostorderInternalNode(T & t)
		:tree(t){
	}
	postorder_iterator<typename T::node_type> begin() const {
		return std::move(postorder_iterator<typename T::node_type>{tree.getRoot(), isInternalNode<typename T::node_type>});
	}
	postorder_iterator<typename T::node_type> end() const {
		return postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		T & tree;
};

template<typename T>
class LeafIter {
	public:
	explicit LeafIter(T & t)
		:tree(t){
	}
	postorder_iterator<typename T::node_type> begin() const {
		return std::move(postorder_iterator<typename T::node_type>{tree.getRoot(), isLeaf<typename T::node_type>});
	}
	postorder_iterator<typename T::node_type> end() const {
		return postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		T & tree;
};

template<typename T>
class ConstLeafIter {
	public:
	explicit ConstLeafIter(const T & t)
		:tree(t){
	}
	const_postorder_iterator<typename T::node_type> begin() const {
		return std::move(const_postorder_iterator<typename T::node_type>{tree.getRoot(), isLeaf<typename T::node_type>});
	}
	const_postorder_iterator<typename T::node_type> end() const {
		return const_postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		const T & tree;
};

template<typename T>
class AncNodeIter {
	public:
	explicit AncNodeIter(T * n)
		:des(n){
	}
	anc_iterator<T> begin() const {
		return std::move(anc_iterator<T>{des});
	}
	anc_iterator<T> end() const {
		return anc_iterator<T>{nullptr};
	}
	private:
		T * des;
};

template<typename T>
class PostorderIter {
	public:
	explicit PostorderIter(T & t)
		:tree(t){
	}
	postorder_iterator<typename T::node_type> begin() const {
		return std::move(postorder_iterator<typename T::node_type>{tree.getRoot()});
	}
	postorder_iterator<typename T::node_type> end() const {
		return postorder_iterator<typename T::node_type>{nullptr};
	}
	private:
		T & tree;
};

} // namespace otc
#endif

