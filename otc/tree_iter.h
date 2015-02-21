#ifndef OTCETERA_TREE_ITER_H
#define OTCETERA_TREE_ITER_H
// Iterators for traversing trees
// Depends on: tree.h tree_util.h
// Depended on by: tree_iter.h tree_operation.h 

#include "otc/otc_base_includes.h"
#include "otc/tree.h"
#include "otc/tree_util.h"

namespace otc {

template<typename T, typename U>
class ConstPreorderInternalNode {
	public:
	explicit ConstPreorderInternalNode(const RootedTree<T, U> &t)
		:tree(t){
	}
	class iterator : std::forward_iterator_tag {
		private:
			std::function<bool(const RootedTreeNode<T> &)> filterFn;
			const RootedTreeNode<T> * curr;
			bool movingDown;
			iterator(const RootedTreeNode<T> *c)
				:filterFn{isInternalNode<T>},
				curr(c),
				movingDown(false) {
				if (c != nullptr) {
				}
			}

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
				} while (!filterFn(*curr));
			}
			friend class ConstPreorderInternalNode<T, U>;
		public:
			bool operator==(const iterator &other) {
				return this->curr == other.curr;
			}
			bool operator!=(const iterator &other) {
				return this->curr != other.curr;
			}
			const RootedTreeNode<T> * operator*() const {
				return curr;
			}
			iterator & operator++() {
				if (curr == nullptr) {
					throw std::out_of_range("Incremented a dead PreorderInternalNode::iterator");
				}
				_advance();
				return *this;
			}
	};
	iterator begin() const {
		return std::move(iterator{tree.GetRoot()});
	}
	iterator end() const {
		return iterator{nullptr};
	}
	private:
		const RootedTree<T, U> & tree;
};

} // namespace otc
#endif

