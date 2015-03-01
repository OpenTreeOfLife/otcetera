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
template<typename T, bool isConst>
class preorder_iterator : std::forward_iterator_tag {
	private:
		const NdFilterFn<T> filterFn;
		typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
		node_pointer curr;
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
		preorder_iterator(node_pointer c)
			:filterFn{nullptr},
			curr(c),
			movingDown(false) {
		}
		preorder_iterator(node_pointer c, NdFilterFn<T> f)
			:filterFn{f},
			curr(c),
			movingDown(false) {
			if (c != nullptr && filterFn && !filterFn(*curr)) {
				_advance();
			}
		}
		bool operator==(const preorder_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const preorder_iterator &other) const {
			return this->curr != other.curr;
		}
		node_pointer operator*() const {
			return curr;
		}
		preorder_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead preorder_iterator");
			}
			_advance();
			return *this;
		}
};


/// visits ancestors before descendants
/// if filterFn is supplied, then only the nodes associated with a true
///		response will be returned (but the descendants of a "false" node
///		will still be visited.
template<typename T>
class const_skipping_preorder_iterator : std::forward_iterator_tag {
	private:
		const NdFilterFn<T> subtreeFilterFn;
		const T * curr;
		bool movingDown;

		void _advance_down() {
			while (movingDown) {
				if (curr->getNextSib()) {
					curr = curr->getNextSib();
					if (subtreeFilterFn == nullptr || subtreeFilterFn(*curr)) {
						movingDown = false;
					}
				} else {
					curr = curr->getParent();
					if (curr == nullptr) {
						return;
					}
				}
			}
		}
		void _advance() {
			assert(curr != nullptr);
			if (movingDown) {
				_advance_down();
			} else {
				_advance_up_or_over();
			}
		}
		void _advance_up_or_over() {
			if (curr->isTip()) {
				if (curr->getNextSib() == nullptr) {
					movingDown = true;
					_advance_down();
					return;
				}
				curr = curr->getNextSib();
			} else {
				curr = curr->getFirstChild();
			}
			if (subtreeFilterFn == nullptr || subtreeFilterFn(*curr)) {
				return;
			}
			while (curr->getNextSib()) {
				curr = curr->getNextSib();
				if (subtreeFilterFn == nullptr || subtreeFilterFn(*curr)) {
					return;
				}
			}
			movingDown = true;
			_advance_down();
			return;
		}
	public:
		const_skipping_preorder_iterator(const T *c)
			:subtreeFilterFn{nullptr},
			curr(c),
			movingDown(false) {
		}
		const_skipping_preorder_iterator(const T *c, NdFilterFn<T> f)
			:subtreeFilterFn{f},
			curr(c),
			movingDown(false) {
			if (c != nullptr && subtreeFilterFn && !subtreeFilterFn(*curr)) {
				curr = nullptr;
			}
		}
		bool operator==(const const_skipping_preorder_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const const_skipping_preorder_iterator &other) const {
			return this->curr != other.curr;
		}
		const T * operator*() const {
			return curr;
		}
		const_skipping_preorder_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead const_skipping_preorder_iterator");
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
		bool operator==(const const_child_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const const_child_iterator &other) const {
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
		bool operator==(const child_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const child_iterator &other) const {
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
		bool operator==(const const_postorder_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const const_postorder_iterator &other) const {
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
					curr = (n == nullptr ? curr->getParent() : findLeftmostInSubtreeM<T>(n));
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
				curr = findLeftmostInSubtreeM<T>(lastNode);
			}
		}
		postorder_iterator(T *c, NdFilterFn<T> f)
			:filterFn{f},
			curr(nullptr),
			lastNode(c) {
			if (lastNode != nullptr) {
				curr = findLeftmostInSubtreeM<T>(lastNode);
			}
			if (curr != nullptr && filterFn && !filterFn(*curr)) {
				_advance();
			}
		}
		bool operator==(const postorder_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const postorder_iterator &other) const {
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
		T * curr;
		
		void _advance() {
			assert(curr != nullptr);
			curr = curr->getParent();
		}
	public:
		anc_iterator(T *c)
			:curr(c) {
			if (curr != nullptr) {
				_advance();
			}
		}
		bool operator==(const anc_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const anc_iterator &other) const {
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
class const_anc_iterator : std::forward_iterator_tag {
	private:
		const T * curr;
		
		void _advance() {
			assert(curr != nullptr);
			curr = curr->getParent();
		}
	public:
		const_anc_iterator(const T *c)
			:curr(c) {
			if (curr != nullptr) {
				_advance();
			}
		}
		bool operator==(const const_anc_iterator &other) const {
			return this->curr == other.curr;
		}
		bool operator!=(const const_anc_iterator &other) const {
			return this->curr != other.curr;
		}
		const T * operator*() const {
			return curr;
		}
		const_anc_iterator & operator++() {
			if (curr == nullptr) {
				throw std::out_of_range("Incremented a dead const_anc_iterator");
			}
			_advance();
			return *this;
		}
};

// children of a node
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
class ConstPreorderInternalIter {
	public:
	explicit ConstPreorderInternalIter(const T &t)
		:tree(t){
	}
	preorder_iterator<typename T::node_type, true> begin() const {
		return std::move(preorder_iterator<typename T::node_type, true>{tree.getRoot(), isInternalNode<typename T::node_type>});
	}
	preorder_iterator<typename T::node_type, true> end() const {
		return preorder_iterator<typename T::node_type, true>{nullptr};
	}
	private:
		const T & tree;
};

template<typename T>
class ConstPreorderIter {
	public:
	explicit ConstPreorderIter(const T &t)
		:tree(t){
	}
	preorder_iterator<typename T::node_type, true> begin() const {
		return std::move(preorder_iterator<typename T::node_type, true>{tree.getRoot()});
	}
	preorder_iterator<typename T::node_type, true> end() const {
		return preorder_iterator<typename T::node_type, true>{nullptr};
	}
	private:
		const T & tree;
};


// preorder constructed with a node
template<typename T>
class ConstPreorderIterN {
	public:
	explicit ConstPreorderIterN(const T * t)
		:nd(t){
	}
	preorder_iterator<T, true> begin() const {
		return std::move(preorder_iterator<T, true>{nd});
	}
	preorder_iterator<T, true> end() const {
		return preorder_iterator<T, true>{nullptr};
	}
	private:
		const T * nd;
};

// preorder constructed with a node
template<typename T>
class ConstSubtreeFilteringPreorderIterN {
	public:
	explicit ConstSubtreeFilteringPreorderIterN(const T * t, std::function<bool(const T &)> subtreeFilterFn)
		:nd(t), 
		f(subtreeFilterFn) {
	}
	const_skipping_preorder_iterator<T> begin() const {
		return std::move(const_skipping_preorder_iterator<T>{nd, f});
	}
	const_skipping_preorder_iterator<T> end() const {
		return const_skipping_preorder_iterator<T>{nullptr};
	}
	private:
		const T * nd;
		std::function<bool(const T &)> f;
};


// Only tips (aka leaves)
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
class AncIter {
	public:
	explicit AncIter(T * n)
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
class ConstAncIter {
	public:
	explicit ConstAncIter(const T * n)
		:des(n){
	}
	const_anc_iterator<T> begin() const {
		return std::move(const_anc_iterator<T>{des});
	}
	const_anc_iterator<T> end() const {
		return const_anc_iterator<T>{nullptr};
	}
	private:
		const T * des;
};

// des before anc
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

template<typename T>
class ConstPostorderInternalIter {
	public:
	explicit ConstPostorderInternalIter(const T &t)
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
class PostorderInternalIter {
	public:
	explicit PostorderInternalIter(T & t)
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



// used to express iteration over all nodes where order does not matter
template<typename T>
using NodeIter = PostorderIter<T>;
template<typename T>
using ConstNodeIter = ConstPostorderIter<T>;
template<typename T>
using InternalNodeIter = PostorderInternalIter<T>;
template<typename T>
using ConstInternalNodeIter = ConstPostorderInternalIter<T>;

} // namespace otc
#endif

