#ifndef OTCETERA_TREE_H
#define OTCETERA_TREE_H

#include <climits>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "otc/otc_base_includes.h"

namespace otc {
template<typename, typename> class RootedTree;

typedef std::string namestring_t;

template<typename T>
class RootedTreeNode {
    public:
        using node_type = RootedTreeNode<T>;
        using data_type = T;

        bool is_tip() const { return (lChild == nullptr); }

        bool is_internal() const { return not is_tip(); }

        const node_type * get_parent() const { return parent; }
              node_type * get_parent()       { return parent; }
              
        const node_type * get_first_child() const { return lChild; }
              node_type * get_first_child()       { return lChild; }

        const node_type * get_last_child() const { return rChild; }
              node_type * get_last_child()       { return rChild; }

        const node_type * get_prev_sib() const { return lSib; }
              node_type * get_prev_sib()       { return lSib; }

        const node_type * get_next_sib() const { return rSib; }
              node_type * get_next_sib()       { return rSib; }

        const node_type * get_first_sib() const { assert(parent); return parent->get_first_child(); }
              node_type * get_first_sib()       { assert(parent); return parent->get_first_child(); }

        const node_type * get_last_sib() const { assert(parent); return parent->get_last_child(); }
              node_type * get_last_sib()       { assert(parent); return parent->get_last_child(); }

        unsigned get_out_degree() const {
            unsigned n = 0;
            auto currNode = get_first_child();
            while(currNode) {
                n += 1;
                currNode = currNode->get_next_sib();
            }
            return n;
        }
        bool has_ott_id() const {
            return ottId != std::numeric_limits<OttId>::max();
        }
        OttId get_ott_id() const {
            assert(has_ott_id());
            return ottId;
        }
        void set_ott_id(OttId i) {
            ottId = i;
            assert(has_ott_id());
        }
        void del_ott_id() {
            ottId = std::numeric_limits<OttId>::max();
        }
        // non-empty only for internals that are labelled with names that are NOT taxLabels
        const namestring_t & get_name() const {
            return name;
        }
        void set_name(const namestring_t &n) {
            std::string t = n;
            swap(name, t);
        }
        void set_name(namestring_t && n) {
            name = std::move(n);
        }
        const T & get_data() const {
            return data;
        }
        T & get_data() {
            return data;
        }
        RootedTreeNode<T>(RootedTreeNode<T> *par)
            :parent(par) {
        }
        void add_sib_on_left(node_type *n) {
            assert(n);

            // Connect left (from n)
            n->lSib = lSib;
            if (n->lSib)
                n->lSib->rSib = n;
            else
            {
                assert(parent->lChild == this);
                parent->lChild = n;
            }

            // Connect right (from n)
            n->rSib = this;
            lSib = n;

            // Connect up (from n)
            n->parent = parent;
        }
        void add_sib_on_right(node_type *n) {
            assert(n);

            // Connect right (from n)
            n->rSib = rSib;
            if (n->rSib) {
                n->rSib->lSib = n;
            } else {
                assert(parent->rChild == this);
                parent->rChild = n;
            }
            // Connect left (from n)
            n->lSib = this;
            rSib = n;
            // Connect up (from n)
            n->parent = parent;
        }
        void add_sib_at_end(node_type *n) {
            get_last_sib()->add_sib_on_right(n);
        }
        void add_sib(node_type *n) {
            add_sib_at_end(n);
        }
        void add_child_at_front(node_type* n) {
            if (lChild) {
                lChild->add_sib_on_left(n);
            } else {
                lChild = n;
                rChild = n;
                n->lSib = nullptr;
                n->rSib = nullptr;
                n->parent = this;
            }
        }
        void add_child(node_type *n) {
            if (rChild) {
                rChild->add_sib(n);
            } else {
                lChild = n;
                rChild = n;
                n->lSib = nullptr;
                n->rSib = nullptr;
                n->parent = this;
            }
        }
        bool has_children() const {
            assert(bool(lChild) == bool(rChild));
            return lChild;
        }
        void remove_child(node_type *n) {
            assert(n);
            assert(has_children());
            assert(n->parent == this);
            n->detach_this_node();
        }
        bool is_outdegree_one_node() const {
            return lChild and not lChild->rSib;
        }
        bool is_outdegree_two_node() const {
            return lChild and lChild->rSib and not lChild->rSib->rSib;
        }
        bool is_polytomy() const {
            return lChild and lChild->rSib and lChild->rSib->rSib;
        }
        bool includes_only_one_leaf() const {
            if (is_tip()) {
                return true;
            }
            return is_outdegree_one_node() && lChild->includes_only_one_leaf();
        }
        // takes this node out of the child array of its parent, but does not fix any other pointer.
        void detach_this_node() {
            assert(parent);
            assert(parent->has_children());
            if (lSib){
                lSib->rSib = rSib;
            } else {
                parent->lChild = rSib;
            }
            if (rSib) {
                rSib->lSib = lSib;
            } else {
                parent->rChild = lSib;
            }
            parent = nullptr;
            rSib = nullptr;
            lSib = nullptr;
        }

        // places `n` into the parent node's child array in the place
        // of `this`. None of the pointer's of `this` are modified, but
        //  `this` will no longer be findable by its parent.
        void replace_this_node(node_type *n) const {
            assert (n != nullptr);
            assert(parent);
            assert(parent->has_children());
            n->lSib = lSib;
            if (lSib){
                lSib->rSib = n;
            } else {
                parent->lChild = n;
            }
            n->rSib = rSib;
            if (rSib) {
                rSib->lSib = n;
            } else {
                parent->rChild = n;
            }
            n->parent = parent;
        }
    public:
        void write_as_newick(std::ostream &out,
                           bool useLeafNames,
                           const std::map<node_type *, namestring_t> *nd2name=nullptr) const {
            auto child = get_first_child();
            if (child) {
                out << "(";
                child->write_as_newick(out,useLeafNames,nd2name);
                child = child->get_next_sib();
                for(;child;child = child->get_next_sib()) {
                    out << ",";
                    child->write_as_newick(out,useLeafNames,nd2name);
                }
                out << ")";
            }
            out << get_name();
        }
        void add_self_and_des_to_preorder(std::vector<const node_type *> &p) const;
    private:
        node_type * lChild = nullptr;
        node_type * rChild = nullptr;
        node_type * lSib = nullptr;
        node_type * rSib = nullptr;
        node_type * parent = nullptr;
        namestring_t name;     // non-empty only for internals that are labelled with names that are NOT taxLabels
        OttId ottId = std::numeric_limits<OttId>::max(); // present for every leaf. UINT_MAX for internals labeled with taxlabels
        T data;
    private:
        RootedTreeNode<T>(const RootedTreeNode<T> &) = delete;
        RootedTreeNode<T> & operator=(const RootedTreeNode<T> &) = delete;
        template<typename Y, typename Z>
        friend class RootedTree;
};

template<typename NodeType>
inline const NodeType* find_root(const NodeType* nd) {
    while(nd->get_parent()) {
        nd = nd->get_parent();
    }
    return nd;
}
    
template<typename T, typename U>
class RootedTree {
    public:
        using node_data_type = T;
        using node_type = RootedTreeNode<T>;
        using data_type = U;
        
        RootedTree<T, U>()
            :root(nullptr) {
        }
        ~RootedTree<T, U>() {
            clear();
        }
        std::vector<const node_type *> get_preorder_traversal() const;
        void write_as_newick(std::ostream &out,
                           bool useLeafNames,
                           const std::map<node_type *, namestring_t> *nd2name=nullptr) const {
            if (root) {
                root->write_as_newick(out, useLeafNames, nd2name);
            }
        }
        const node_type * get_root() const {
            return root;
        }
        node_type * get_root() {
            return root;
        }
        void _set_root(node_type * r) {
            root = r;
        }
        node_type * create_root() {
            if (root != nullptr) {
                clear();
            }
            this->root = this->alloc_new_node(nullptr);
            return this->root;
        }
        node_type * create_child(node_type *par) {
            auto c = this->alloc_new_node(par);
            par->add_child(c);
            return c;
        }
        node_type * create_node(node_type *par) {
            auto c = this->alloc_new_node(par);
            if (par != nullptr) {
                par->add_child(c);
            }
            return c;
        }
        node_type * create_sib(node_type *leftSib) {
            assert(leftSib->parent != nullptr);
            auto s = this->alloc_new_node(leftSib->parent);
            leftSib->add_sib(s);
            return s;
        }
        U & get_data() {
            return this->data;
        }
        const U & get_data() const {
            return this->data;
        }
        void prune_and_dangle(node_type * nd) {
            assert(find_root(nd) == root);
            auto p = nd->get_parent();
            if (p == nullptr) {
                root = nullptr;
                return;
            }
            p->remove_child(nd);
        }
        void prune_and_delete(node_type * nd) {
            auto nodes = get_subtree_nodes(nd);
            prune_and_dangle(nd);
            for (auto ndi: nodes) {
                delete ndi;
            }
        }
        bool is_detached(node_type * nd) {
            return contains(detached, nd);
        }
    protected:
        node_type * root;
        U data;
        std::string name;
        std::set<node_type *> detached;
        
    public:
        void set_name(const std::string &n) {
            name.assign(n);
        }
        const std::string & get_name() const {
            return name;
        }
        node_type * alloc_new_node(node_type *p) {
            node_type * nd = new node_type(p);
            return nd;
        }
        void clear() {
            for(auto nd: get_all_attached_nodes()) {
                delete nd;
            }
            root = NULL;
        }
        std::vector<const node_type *> get_all_attached_nodes() const {
            return get_subtree_nodes(root);
        }
          
        std::vector<const node_type *> get_subtree_nodes(node_type* p) const {
            std::vector<const node_type*> nodes;
            if (p) {
                nodes.push_back(p);
            }
            for(std::size_t i = 0; i < nodes.size(); i++) {
                for(auto n = nodes[i]->lChild; n; n = n->rSib) {
                    nodes.push_back(n);
                }
            }
            return nodes;
        }
        std::set<const node_type *> get_set_of_all_attached_nodes() const {
            std::set<const node_type*> nodes;
            for(auto nd: get_all_attached_nodes())
                nodes.insert(nd);
            return nodes;
        }
        void mark_as_detached(node_type * nd) {
            detached.insert(nd);
        }
        void mark_as_attached(node_type * nd) {
          detached.erase(nd);
        }
        std::set<const node_type *> get_detached() const{
            std::set<const node_type *> r;
            for (auto d : detached) {
                r.insert(d);
            }
            return r;
        }
    private:
        RootedTree<T, U>(const RootedTree<T, U> &) = delete;
        RootedTree<T, U> & operator=(const RootedTree<T, U> &) = delete;
};

class RTNodeNoData{};
class RTreeNoData{};

template<typename Tree>
inline void add_subtree(typename Tree::node_type* par, Tree& T2) {
    auto c = T2.get_root();
    T2.prune_and_dangle(c);
    par->add_child(c);
}

template<typename Tree>
void replace_with_subtree(typename Tree::node_type* n, Tree& T2) {
    // Get the parent of the tip we are replacing
    auto p = n->get_parent();
    // Remove the data from T2 and attach it to this parent
    auto c = T2.get_root();
    T2.prune_and_dangle(c);
    p->add_child(c);
    // Remove the old child from under p
    p->remove_child(n);
    // Make the old child the root of T2, which will handle deleting it
    T2._set_root(n);
}

} // namespace otc
#endif

