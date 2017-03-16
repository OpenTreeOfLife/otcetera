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
inline bool is_ancestor_des_no_iter(const T *a, const T *d) {
    auto p = d->get_parent();
    while (p) {
        if (p == a) {
            return true;
        }
        p = p->get_parent();
    }
    return false;
}

template<typename T>
inline T * get_deepest_anc(T *nd) {
    auto p = nd->get_parent();
    while (p != nullptr) {
        nd = p;
        p = nd->get_parent();
    }
    return nd;
}

template<typename T>
using NdFilterFn = std::function<bool(const T &)>;

/// visits ancestors before descendants
/// if filter_fn is supplied, then only the nodes associated with a true
///     response will be returned (but the descendants of a "false" node
///     will still be visited.
template<typename T, bool isConst>
class preorder_iterator : std::forward_iterator_tag {
    private:
        const NdFilterFn<T> subtree_filter_fn;
        const NdFilterFn<T> filter_fn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;
        bool moving_down;
        node_pointer const exit_node;
        void _advance() {
            if (subtree_filter_fn == nullptr) {
                _unfiltered_advance();
            } else {
                _filtered_advance();
            }
            if (curr == exit_node) {
                curr = nullptr;
            }
            if (exit_node == nullptr) {
                return;
            }
        }
        void _unfiltered_advance() {
            do {
                assert(curr != nullptr);
                if (moving_down) {
                    while (moving_down) {
                        if (curr->get_next_sib()) {
                            curr = curr->get_next_sib();
                            moving_down = false;
                        } else {
                            curr = curr->get_parent();
                            if (curr == nullptr || curr == exit_node) {
                                return;
                            }
                        }
                    }
                } else if (curr->is_tip()) {
                    if (curr->get_next_sib()) {
                        curr = curr->get_next_sib();
                    } else {
                        moving_down = true;
                        _unfiltered_advance();
                        if (curr == nullptr  || curr == exit_node) {
                            return;
                        }
                    }
                } else {
                    curr = curr->get_first_child();
                }
            } while (filter_fn != nullptr && !filter_fn(*curr));
        }

        void _filtered_advance_down() {
            assert(curr != nullptr);
            while (moving_down) {
                assert(curr != nullptr);
                if (curr->get_next_sib()) {
                    curr = curr->get_next_sib();
                    assert(curr != nullptr);
                    if (subtree_filter_fn == nullptr || subtree_filter_fn(*curr)) {
                        moving_down = false;
                    }
                } else {
                    curr = curr->get_parent();
                    if (curr == nullptr || curr == exit_node) {
                        return;
                    }
                }
            }
        }
        void _filtered_advance() {
            assert(curr != nullptr);
            do {
                if (moving_down) {
                    _filtered_advance_down();
                } else {
                    _filtered_advance_up_or_over();
                }
            } while(curr != nullptr && curr != exit_node && (filter_fn != nullptr && !filter_fn(*curr)));
        }
        void _filtered_advance_up_or_over() {
            if (curr->is_tip()) {
                if (curr->get_next_sib() == nullptr) {
                    moving_down = true;
                    _filtered_advance_down();
                    return;
                }
                curr = curr->get_next_sib();
            } else {
                curr = curr->get_first_child();
            }
            if (subtree_filter_fn == nullptr || subtree_filter_fn(*curr)) {
                return;
            }
            while (curr->get_next_sib()) {
                curr = curr->get_next_sib();
                if (subtree_filter_fn == nullptr || subtree_filter_fn(*curr)) {
                    return;
                }
            }
            moving_down = true;
            _filtered_advance_down();
            return;
        }
    public:
        preorder_iterator(node_pointer c, node_pointer exit_node)
            :subtree_filter_fn(nullptr),
            filter_fn{nullptr},
            curr(c),
            moving_down(false),
            exit_node(exit_node) {
        }
        preorder_iterator(node_pointer c, node_pointer exit_node, NdFilterFn<T> f)
            :subtree_filter_fn(nullptr),
            filter_fn{f},
            curr(c),
            moving_down(false),
            exit_node(exit_node)  {
            if (c != nullptr && filter_fn && !filter_fn(*curr)) {
                _advance();
            }
        }
        preorder_iterator(node_pointer c, node_pointer exit_node, NdFilterFn<T> f, NdFilterFn<T> s)
            :subtree_filter_fn(s),
            filter_fn{f},
            curr(c),
            moving_down(false),
            exit_node(exit_node)  {
            if (c != nullptr && filter_fn && !filter_fn(*curr)) {
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
        NdFilterFn<T> filter_fn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;

        void _advance() {
            assert(curr != nullptr);
            curr = curr->get_next_sib();
        }
    public:
        child_iterator(node_pointer c)
            :filter_fn{nullptr},
            curr(nullptr) {
            if (c != nullptr) {
                curr = c->get_first_child();
            }
        }
        child_iterator(node_pointer c, NdFilterFn<T> f)
            :filter_fn{f},
            curr(nullptr) {
            if (c != nullptr) {
                curr = c->get_first_child();
            }
            if (c != nullptr && filter_fn && !filter_fn(*c)) {
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
        const NdFilterFn<T> filter_fn;
        typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
        node_pointer curr;
        node_pointer lastNode;

        void _advance() {
            if (curr == lastNode) {
                curr = nullptr;
            } else {
                for (;;) {
                    auto n = curr->get_next_sib();
                    curr = (n == nullptr ? curr->get_parent() : find_leftmost_in_subtree<node_pointer>(n));
                    if (filter_fn == nullptr || filter_fn(*curr)) {
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
            :filter_fn{f},
            curr(nullptr),
            lastNode(c) {
            if (lastNode != nullptr) {
                curr = find_leftmost_in_subtree<node_pointer>(lastNode);
            }
            if (curr != nullptr && filter_fn && !filter_fn(*curr)) {
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
            curr = curr->get_parent();
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
                          node_pointer exit_pointer,
                          std::function<bool(const T &)> nd_filter,
                          std::function<bool(const T &)> subtree_filter)
        :nd(t),
        exit_node(exit_pointer),
        node_filter(nd_filter),
        clade_filter(subtree_filter) {
    }
    preorder_iterator<T, isConst> begin() const {
        return std::move(preorder_iterator<T, isConst>{nd, exit_node, node_filter, clade_filter});
    }
    preorder_iterator<T, isConst> end() const {
        return std::move(preorder_iterator<T, isConst>{nullptr, nullptr});
    }
    private:
        node_pointer nd;
        node_pointer const exit_node;
        std::function<bool(const T &)> node_filter;
        std::function<bool(const T &)> clade_filter;
};

template<typename T, bool isConst>
class LeafIter {
    public:
    typedef typename std::conditional<isConst, const T *, T *>::type node_pointer;
    explicit LeafIter(node_pointer t)
        :node(t){
    }
    postorder_iterator<T, isConst> begin() const {
        return std::move(postorder_iterator<T, isConst>{node, is_leaf<T>});
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
        node_filter(f) {
    }
    postorder_iterator<T, isConst> begin() const {
        return std::move(postorder_iterator<T, isConst>{node, node_filter});
    }
    postorder_iterator<T, isConst> end() const {
        return std::move(postorder_iterator<T, isConst>{nullptr, nullptr});
    }
    private:
        node_pointer node;
        std::function<bool(const T &)> node_filter;
};

// used to express iteration over all nodes where order does not matter
template<typename T, bool isConst>
using NodeIter = PostorderIter<T, isConst>;

// the public interface
template<typename T>
inline PostorderIter<typename T::node_type, false> iter_post(T & tree) {
    return std::move(PostorderIter<typename T::node_type, false>(tree.get_root(), nullptr));
}

template<typename T>
inline PostorderIter<typename T::node_type, true> iter_post_const(const T & tree) {
    return std::move(PostorderIter<typename T::node_type, true>(tree.get_root(), nullptr));
}
template<typename T>
inline PostorderIter<T, false> iter_post_n(T & node) {
    return std::move(PostorderIter<T, false>(&node, nullptr));
}

template<typename T>
inline PostorderIter<T, true> iter_post_n_const(const T & node) {
    return std::move(PostorderIter<T, true>(&node, nullptr));
}

template<typename T>
inline PostorderIter<typename T::node_type, false> iter_post_internal(T & tree) {
    PostorderIter<typename T::node_type, false> i{tree.get_root(), is_internal_node<typename T::node_type>};
    return std::move(i);
}

template<typename T>
inline PostorderIter<typename T::node_type, true> iter_post_internal_const(const T & tree) {
    return std::move(PostorderIter<typename T::node_type, true>(tree.get_root(), is_internal_node<typename T::node_type>));
}

template<typename T>
inline NodeIter<typename T::node_type, false> iter_node(T & tree) {
    return std::move(NodeIter<typename T::node_type, false>(tree.get_root(), nullptr));
}

template<typename T>
inline NodeIter<typename T::node_type, true> iter_node_const(const T & tree) {
    return std::move(NodeIter<typename T::node_type, true>(tree.get_root(), nullptr));
}

template<typename T>
inline NodeIter<typename T::node_type, false> iter_node_internal(T & tree) {
    return std::move(NodeIter<typename T::node_type, false>(tree.get_root(), is_internal_node<typename T::node_type>));
}

template<typename T>
inline NodeIter<typename T::node_type, true> iter_node_internal_const(const T & tree) {
    return std::move(NodeIter<typename T::node_type, true>(tree.get_root(), is_internal_node<typename T::node_type>));
}

template<typename T>
inline PreorderIter<typename T::node_type, false> iter_pre_internal(T & tree) {
    return std::move(PreorderIter<typename T::node_type, false>(tree.get_root(), nullptr, is_internal_node<typename T::node_type>, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, true> iter_pre_internal_const(const T & tree) {
    return std::move(PreorderIter<typename T::node_type, true>(tree.get_root(), nullptr, is_internal_node<typename T::node_type>, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, false> iter_pre(T & tree) {
    return std::move(PreorderIter<typename T::node_type, false>(tree.get_root(), nullptr, nullptr, nullptr));
}

template<typename T>
inline PreorderIter<typename T::node_type, true> iter_pre_const(const T & tree) {
    return std::move(PreorderIter<typename T::node_type, true>(tree.get_root(), nullptr, nullptr, nullptr));
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
    return std::move(LeafIter<typename T::node_type, false>(tree.get_root()));
}

template<typename T>
inline LeafIter<typename T::node_type, true> iter_leaf_const(const T & tree) {
    return std::move(LeafIter<typename T::node_type, true>(tree.get_root()));
}

template<typename T>
inline LeafIter<T, false> iter_leaf_n(T & node) {
    return std::move(LeafIter<T, false>(&node));
}

template<typename T>
inline LeafIter<T, true> iter_leaf_n_const(const T & node) {
    return std::move(LeafIter<T, true>(&node));
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

