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
bool isAncestorDesNoIter(const T *a, const T *d) {
    auto p = d->getParent();
    while (p) {
        if (p == a) {
            return true;
        }
        p = p->getParent();
    }
    return false;
}

template<typename T>
using NdFilterFn = std::function<bool(const T &)>;

/// visits ancestors before descendants
/// if filterFn is supplied, then only the nodes associated with a true
///     response will be returned (but the descendants of a "false" node
///     will still be visited.
template<typename T, bool isConst>
class preorder_iterator : std::forward_iterator_tag {
    private:
        const NdFilterFn<T> subtreeFilterFn;
        const NdFilterFn<T> filterFn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;
        bool movingDown;
        node_pointer const exitNode;
        void _advance() {
            if (subtreeFilterFn == nullptr) {
                _unfiltered_advance();
            } else {
                _filtered_advance();
            }
            if (curr == exitNode) {
                curr = nullptr;
            }
            if (exitNode == nullptr) {
                return;
            }
#if 0
            if (curr == nullptr) {
                std::cerr << "preorder_iterator.curr = nullptr";
            } else {
                std::cerr << "preorder_iterator.curr = " << getDesignator(*curr);
                if (exitNode != nullptr) {
                    assert(isAncestorDesNoIter(exitNode, curr));
                }
            }
            if (exitNode != nullptr) {
                std::cerr << " exitNode = " << getDesignator(*exitNode) << '\n';
            } else {
                std::cerr << '\n';
            }
#endif
        }
        void _unfiltered_advance() {
            do {
                assert(curr != nullptr);
                if (movingDown) {
                    while (movingDown) {
                        if (curr->getNextSib()) {
                            curr = curr->getNextSib();
                            movingDown = false;
                        } else {
                            curr = curr->getParent();
                            if (curr == nullptr || curr == exitNode) {
                                return;
                            }
                        }
                    }
                } else if (curr->isTip()) {
                    if (curr->getNextSib()) {
                        curr = curr->getNextSib();
                    } else {
                        movingDown = true;
                        _unfiltered_advance();
                        if (curr == nullptr  || curr == exitNode) {
                            return;
                        }
                    }
                } else {
                    curr = curr->getFirstChild();
                }
            } while (filterFn != nullptr && !filterFn(*curr));
        }

        void _filtered_advance_down() {
            assert(curr != nullptr);
            while (movingDown) {
                assert(curr != nullptr);
                if (curr->getNextSib()) {
                    curr = curr->getNextSib();
                    assert(curr != nullptr);
                    if (subtreeFilterFn == nullptr || subtreeFilterFn(*curr)) {
                        movingDown = false;
                    }
                } else {
                    curr = curr->getParent();
                    if (curr == nullptr || curr == exitNode) {
                        return;
                    }
                }
            }
        }
        void _filtered_advance() {
            assert(curr != nullptr);
            do {
                if (movingDown) {
                    _filtered_advance_down();
                } else {
                    _filtered_advance_up_or_over();
                }
            } while(curr != nullptr && curr != exitNode && (filterFn != nullptr && !filterFn(*curr)));
        }
        void _filtered_advance_up_or_over() {
            if (curr->isTip()) {
                if (curr->getNextSib() == nullptr) {
                    movingDown = true;
                    _filtered_advance_down();
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
            _filtered_advance_down();
            return;
        }
    public:
        preorder_iterator(node_pointer c, node_pointer exit_node)
            :subtreeFilterFn(nullptr),
            filterFn{nullptr},
            curr(c),
            movingDown(false),
            exitNode(exit_node) {
        }
        preorder_iterator(node_pointer c, node_pointer exit_node, NdFilterFn<T> f)
            :subtreeFilterFn(nullptr),
            filterFn{f},
            curr(c),
            movingDown(false),
            exitNode(exit_node)  {
            if (c != nullptr && filterFn && !filterFn(*curr)) {
                _advance();
            }
        }
        preorder_iterator(node_pointer c, node_pointer exit_node, NdFilterFn<T> f, NdFilterFn<T> s)
            :subtreeFilterFn(s),
            filterFn{f},
            curr(c),
            movingDown(false),
            exitNode(exit_node)  {
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

/// descendants before ancestors, but not guaranteed to be the reverse of const_preorder_iterator
template<typename T, bool isConst>
class child_iterator : std::forward_iterator_tag {
    private:
        NdFilterFn<T> filterFn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;

        void _advance() {
            assert(curr != nullptr);
            curr = curr->getNextSib();
        }
    public:
        child_iterator(node_pointer c)
            :filterFn{nullptr},
            curr(nullptr) {
            if (c != nullptr) {
                curr = c->getFirstChild();
            }
        }
        child_iterator(node_pointer c, NdFilterFn<T> f)
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
        node_pointer operator*() const {
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
template<typename T, bool isConst>
class postorder_iterator : std::forward_iterator_tag {
    private:
        const NdFilterFn<T> filterFn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;
        node_pointer lastNode;

        void _advance() {
            if (curr == lastNode) {
                curr = nullptr;
            } else {
                for (;;) {
                    auto n = curr->getNextSib();
                    curr = (n == nullptr ? curr->getParent() : findLeftmostInSubtree<node_pointer>(n));
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
        postorder_iterator(node_pointer c, NdFilterFn<T> f)
            :filterFn{f},
            curr(nullptr),
            lastNode(c) {
            if (lastNode != nullptr) {
                curr = findLeftmostInSubtree<node_pointer>(lastNode);
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
        node_pointer operator*() const {
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

template<typename T, bool isConst>
class anc_iterator : std::forward_iterator_tag {
    private:
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;
        
        void _advance() {
            assert(curr != nullptr);
            curr = curr->getParent();
        }
    public:
        anc_iterator(node_pointer c)
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
        node_pointer operator*() const {
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

// children of a node
template<typename T, bool isConst>
class ChildIter {
    public:
    typedef typename std::conditional<isConst, const T &, T &>::type node_ref;
    explicit ChildIter(node_ref n)
        :node(n) {
    }
    child_iterator<T, isConst> begin() const {
        return std::move(child_iterator<T, isConst>{&node});
    }
    child_iterator<T, isConst> end() const {
        return std::move(child_iterator<T, isConst>{nullptr});
    }
    private:
        node_ref node;
};

template<typename T, bool isConst>
class PreorderIter {
    public:
    typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
    explicit PreorderIter(node_pointer t,
                          node_pointer exitPointer,
                          std::function<bool(const T &)> ndFilter,
                          std::function<bool(const T &)> subtreeFilter)
        :nd(t),
        exitNode(exitPointer),
        nodeFilter(ndFilter),
        cladeFilter(subtreeFilter) {
    }
    preorder_iterator<T, isConst> begin() const {
        return std::move(preorder_iterator<T, isConst>{nd, exitNode, nodeFilter, cladeFilter});
    }
    preorder_iterator<T, isConst> end() const {
        return std::move(preorder_iterator<T, isConst>{nullptr, nullptr});
    }
    private:
        node_pointer nd;
        node_pointer const exitNode;
        std::function<bool(const T &)> nodeFilter;
        std::function<bool(const T &)> cladeFilter;
};

template<typename T, bool isConst>
class LeafIter {
    public:
    typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
    explicit LeafIter(node_pointer t)
        :node(t){
    }
    postorder_iterator<T, isConst> begin() const {
        return std::move(postorder_iterator<T, isConst>{node, isLeaf<T>});
    }
    postorder_iterator<T, isConst> end() const {
        return std::move(postorder_iterator<T, isConst>{nullptr, nullptr});
    }
    private:
        node_pointer node;
};

template<typename T, bool isConst>
class AncIter {
    public:
    typedef typename std::conditional<isConst, const T &, T &>::type node_ref;
    explicit AncIter(node_ref n)
        :des(n){
    }
    anc_iterator<T, isConst> begin() const {
        return std::move(anc_iterator<T, isConst>{&des});
    }
    anc_iterator<T, isConst> end() const {
        return std::move(anc_iterator<T, isConst>{nullptr});
    }
    private:
        node_ref des;
};

template<typename T, bool isConst>
class PostorderIter {
    public:
    typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
    explicit PostorderIter(node_pointer n, std::function<bool(const T &)> f)
        :node(n), 
        nodeFilter(f) {
    }
    postorder_iterator<T, isConst> begin() const {
        return std::move(postorder_iterator<T, isConst>{node, nodeFilter});
    }
    postorder_iterator<T, isConst> end() const {
        return std::move(postorder_iterator<T, isConst>{nullptr, nullptr});
    }
    private:
        node_pointer node;
        std::function<bool(const T &)> nodeFilter;
};

// used to express iteration over all nodes where order does not matter
template<typename T, bool isConst>
using NodeIter = PostorderIter<T, isConst>;

// the public interface
template<typename T>
inline PostorderIter<typename T::node_type, false> iter_post(T & tree) {
    return std::move(PostorderIter<typename T::node_type, false>(tree.getRoot(), nullptr));
}

template<typename T>
inline PostorderIter<typename T::node_type, true> iter_post_const(const T & tree) {
    return std::move(PostorderIter<typename T::node_type, true>(tree.getRoot(), nullptr));
}

template<typename T>
inline PostorderIter<typename T::node_type, false> iter_post_internal(T & tree) {
    PostorderIter<typename T::node_type, false> i{tree.getRoot(), isInternalNode<typename T::node_type>};
    return std::move(i);
}

template<typename T>
inline PostorderIter<typename T::node_type, true> iter_post_internal_const(const T & tree) {
    return std::move(PostorderIter<typename T::node_type, true>(tree.getRoot(), isInternalNode<typename T::node_type>));
}

template<typename T>
inline NodeIter<typename T::node_type, false> iter_node(T & tree) {
    return std::move(NodeIter<typename T::node_type, false>(tree.getRoot(), nullptr));
}

template<typename T>
inline NodeIter<typename T::node_type, true> iter_node_const(const T & tree) {
    return std::move(NodeIter<typename T::node_type, true>(tree.getRoot(), nullptr));
}

template<typename T>
inline NodeIter<typename T::node_type, false> iter_node_internal(T & tree) {
    return std::move(NodeIter<typename T::node_type, false>(tree.getRoot(), isInternalNode<typename T::node_type>));
}

template<typename T>
inline NodeIter<typename T::node_type, true> iter_node_internal_const(const T & tree) {
    return std::move(NodeIter<typename T::node_type, true>(tree.getRoot(), isInternalNode<typename T::node_type>));
}

template<typename T>
inline PreorderIter<typename T::node_type, false> iter_pre_internal(T & tree) {
    return std::move(PreorderIter<typename T::node_type, false>(tree.getRoot(), nullptr, isInternalNode<typename T::node_type>, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, true> iter_pre_internal_const(const T & tree) {
    return std::move(PreorderIter<typename T::node_type, true>(tree.getRoot(), nullptr, isInternalNode<typename T::node_type>, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, false> iter_pre(T & tree) {
    return std::move(PreorderIter<typename T::node_type, false>(tree.getRoot(), nullptr, nullptr, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, true> iter_pre_const(const T & tree) {
    return std::move(PreorderIter<typename T::node_type, true>(tree.getRoot(), nullptr, nullptr, nullptr));
}

template<typename T>
inline PreorderIter<T, false> iter_pre_n(T * node) {
    return std::move(PreorderIter<T, false>(node, node, nullptr, nullptr));
}

template<typename T>
inline PreorderIter<T, true> iter_pre_n_const(const T * node) {
    return std::move(PreorderIter<T, true>(node, node, nullptr, nullptr));
}

template<typename T>
inline PreorderIter<T, false> iter_pre_filter_n(T * node, std::function<bool(const T &)> f) {
    return std::move(PreorderIter<T, false>(node, node, nullptr, f));
}

template<typename T>
inline PreorderIter<T, true> iter_pre_filter_n_const(const T * node, std::function<bool(const T &)> f) {
    return std::move(PreorderIter<T, true>(node, node, nullptr, f));
}

template<typename T>
inline ChildIter<T, false> iter_child(T & node) {
    return std::move(ChildIter<T, false>(node));
}

template<typename T>
inline ChildIter<T, true> iter_child_const(const T & node) {
    return std::move(ChildIter<T, true>(node));
}

template<typename T>
inline LeafIter<typename T::node_type, false> iter_leaf(T & tree) {
    return std::move(LeafIter<typename T::node_type, false>(tree.getRoot()));
}

template<typename T>
inline LeafIter<typename T::node_type, true> iter_leaf_const(const T & tree) {
    return std::move(LeafIter<typename T::node_type, true>(tree.getRoot()));
}

template<typename T>
inline AncIter<T, false> iter_anc(T & node) {
    return std::move(AncIter<T, false>(node));
}

template<typename T>
inline AncIter<T, true> iter_anc_const(const T & node) {
    return std::move(AncIter<T, true>(node));
}

} // namespace otc
#endif

