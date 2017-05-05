#ifndef OTCETERA_DEBUG_H
#define OTCETERA_DEBUG_H

#include "otc/otc_base_includes.h"
#include "otc/tree_data.h"
#include "otc/tree_operations.h"

namespace otc {

template<typename T>
bool check_node_pointers(const T & nd) {
    bool good = true;
    auto c = nd.get_first_child();
    while (c != nullptr) {
        if (c->get_parent() != &nd) {
            good = false;
            assert(c->get_parent() == &nd);
        }
        c = c->get_next_sib();
    }
    auto ps = nd.get_prev_sib();
    assert(ps == nullptr || ps->get_next_sib() == &nd);
    return good;
}

template<typename T>
bool check_all_node_pointers(const T & tree) {
    auto ns = tree.get_all_attached_nodes();
    for (auto nd : ns) {
        if (!check_node_pointers(*nd)) {
            return false;
        }
    }
    if (tree.get_root()->get_parent() != nullptr) {
        assert(tree.get_root()->get_parent() == nullptr);
        return false;
    }
    return true;
}

template<typename T>
bool check_all_node_pointers_iter(const T & node) {
   for (auto nd : iter_pre_n_const(&node)) {
        if (!check_node_pointers(*nd)) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool check_preorder(const T & tree) {
    auto ns = tree.get_set_of_all_attached_nodes();
    std::set<const typename T::node_type *> visited;
    for (auto nd : iter_pre_const(tree)) {
        auto p = nd->get_parent();
        if (p) {
            if (!contains(visited, p)) {
                assert(contains(visited, p));
                return false;
            }
        }
        visited.insert(nd);
    }
    if (visited != ns) {
        assert(visited == ns);
        return false;
    }
    return true;
}

template<typename T>
bool check_postorder(const T & tree) {
    auto ns = tree.get_set_of_all_attached_nodes();
    std::set<const typename T::node_type *> visited;
    for (auto nd : iter_post_const(tree)) {
        auto p = nd->get_parent();
        if (p) {
            if (contains(visited, p)) {
                assert(!contains(visited, p));
                return false;
            }
        }
        visited.insert(nd);
    }
    if (visited != ns) {
        assert(visited == ns);
        return false;
    }
    return true;
}

template<typename T>
bool check_child_iter(const T & tree) {
    auto ns = tree.get_all_attached_nodes();
    for (auto nd : ns) {
        auto nc = nd->get_out_degree();
        auto v = 0U;
        for (auto c : iter_child_const(*nd)) {
            assert(c->get_parent() == nd);
            ++v;
        }
        assert(v == nc);
    }
    return true;
}

template<typename T>
bool check_des_ids(const T & tree);

template<typename T>
inline bool check_des_ids(const T & ) {
    return true;
}

template<>
inline bool check_des_ids(const TreeMappedWithSplits & tree) {
    auto ns = tree.get_set_of_all_attached_nodes();
    for (auto nd : ns) {
        if (nd->is_tip()) {
            continue;
        }
        OttIdSet d;
        auto sum = 0U;
        for (auto c : iter_child_const(*nd)) {
            const auto & cd = c->get_data().des_ids;
            sum += cd.size();
            assert(cd.size() > 0);
            d.insert(cd.begin(), cd.end());
        }
        assert(sum == d.size());
        if (tree.get_data().des_id_sets_contain_internals) {
            if (nd->has_ott_id()) {
                d.insert(nd->get_ott_id());
            }
        }
        if (!is_subset(d, nd->get_data().des_ids)) {
            std::cerr << "node " << nd->get_ott_id() << '\n';
            write_ott_id_set_diff(std::cerr, " ", nd->get_data().des_ids, " node.desId ", d, " calc.");
            assert(is_subset(d, nd->get_data().des_ids));
        }
    }
    return true;
}


template<typename T>
bool check_tree_invariants(const T & tree) {
    if (!check_all_node_pointers(tree)) {
        return false;
    }
    if (!check_preorder(tree)) {
        return false;
    }
    if (!check_postorder(tree)) {
        return false;
    }
    if (!check_child_iter(tree)) {
        return false;
    }
    if (!check_des_ids(tree)) {
        return false;
    }
    return true;
}

} // namespace otc
#endif

