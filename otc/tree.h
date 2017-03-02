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

        bool isTip() const { return (lChild == nullptr); }

        bool isInternal() const { return not isTip(); }

        const node_type * getParent() const { return parent; }
              node_type * getParent()       { return parent; }
              
        const node_type * getFirstChild() const { return lChild; }
              node_type * getFirstChild()       { return lChild; }

        const node_type * getLastChild() const { return rChild; }
              node_type * getLastChild()       { return rChild; }

        const node_type * getPrevSib() const { return lSib; }
              node_type * getPrevSib()       { return lSib; }

        const node_type * getNextSib() const { return rSib; }
              node_type * getNextSib()       { return rSib; }

        const node_type * getFirstSib() const { assert(parent); return parent->getFirstChild(); }
              node_type * getFirstSib()       { assert(parent); return parent->getFirstChild(); }

        const node_type * getLastSib() const { assert(parent); return parent->getLastChild(); }
              node_type * getLastSib()       { assert(parent); return parent->getLastChild(); }

        unsigned getOutDegree() const {
            unsigned n = 0;
            auto currNode = getFirstChild();
            while(currNode) {
                n += 1;
                currNode = currNode->getNextSib();
            }
            return n;
        }
        bool has_ott_id() const {
            return ottId != LONG_MAX;
        }
        long get_ott_id() const {
            assert(has_ott_id());
            return ottId;
        }
        void set_ott_id(long i) {
            ottId = i;
            assert(has_ott_id());
        }
        void delOttId() {
            ottId = LONG_MAX;
        }
        // non-empty only for internals that are labelled with names that are NOT taxLabels
        const namestring_t & get_name() const {
            return name;
        }
        void setName(const namestring_t &n) {
            name = n;
        }
        void setName(namestring_t && n) {
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
        void addSibOnLeft(node_type *n) {
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
        void addSibOnRight(node_type *n) {
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
        void addSibAtEnd(node_type *n) {
            getLastSib()->addSibOnRight(n);
        }
        void addSib(node_type *n) {
            addSibAtEnd(n);
        }
        void addChildAtFront(node_type* n) {
            if (lChild)
                lChild->addSibOnLeft(n);
            else
            {
                lChild = n;
                rChild = n;
                n->lSib = nullptr;
                n->rSib = nullptr;
                n->parent = this;
            }
        }
        void addChild(node_type *n) {
            if (rChild) {
                rChild->addSib(n);
            } else {
                lChild = n;
                rChild = n;
                n->lSib = nullptr;
                n->rSib = nullptr;
                n->parent = this;
            }
        }
        bool hasChildren() const {
            assert(bool(lChild) == bool(rChild));
            return lChild;
        }
        void removeChild(node_type *n) {
            assert(n);
            assert(hasChildren());
            assert(n->parent == this);
            n->detachThisNode();
        }
        bool isOutDegreeOneNode() const {
            return lChild and not lChild->rSib;
        }
        bool isOutDegreeTwoNode() const {
            return lChild and lChild->rSib and not lChild->rSib->rSib;
        }
        bool isPolytomy() const {
            return lChild and lChild->rSib and lChild->rSib->rSib;
        }
        bool includesOnlyOneLeaf() const {
            if (isTip()) {
                return true;
            }
            return isOutDegreeOneNode() && lChild->includesOnlyOneLeaf();
        }
        // takes this node out of the child array of its parent, but does not fix any other pointer.
        void detachThisNode() {
            assert(parent);
            assert(parent->hasChildren());
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
        void replaceThisNode(node_type *n) const {
            assert (n != nullptr);
            assert(parent);
            assert(parent->hasChildren());
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
            auto child = getFirstChild();
            if (child) {
                out << "(";
                child->write_as_newick(out,useLeafNames,nd2name);
                child = child->getNextSib();
                for(;child;child = child->getNextSib()) {
                    out << ",";
                    child->write_as_newick(out,useLeafNames,nd2name);
                }
                out << ")";
            }
            out << get_name();
        }
        void addSelfAndDesToPreorder(std::vector<const node_type *> &p) const;
    private:
        node_type * lChild = nullptr;
        node_type * rChild = nullptr;
        node_type * lSib = nullptr;
        node_type * rSib = nullptr;
        node_type * parent = nullptr;
        namestring_t name;     // non-empty only for internals that are labelled with names that are NOT taxLabels
        long ottId = LONG_MAX; // present for every leaf. UINT_MAX for internals labeled with taxlabels
        T data;
    private:
        RootedTreeNode<T>(const RootedTreeNode<T> &) = delete;
        RootedTreeNode<T> & operator=(const RootedTreeNode<T> &) = delete;
        template<typename Y, typename Z>
        friend class RootedTree;
};

template<typename NodeType>
const NodeType* findRoot(const NodeType* nd)
{
    while(nd->getParent())
        nd = nd->getParent();
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
        std::vector<const node_type *> getPreorderTraversal() const;
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
            this->root = this->allocNewNode(nullptr);
            return this->root;
        }
        node_type * create_child(node_type *par) {
            auto c = this->allocNewNode(par);
            par->addChild(c);
            return c;
        }
        node_type * create_node(node_type *par) {
            auto c = this->allocNewNode(par);
            if (par != nullptr) {
                par->addChild(c);
            }
            return c;
        }
        node_type * createSib(node_type *leftSib) {
            assert(leftSib->parent != nullptr);
            auto s = this->allocNewNode(leftSib->parent);
            leftSib->addSib(s);
            return s;
        }
        U & get_data() {
            return this->data;
        }
        const U & get_data() const {
            return this->data;
        }
        void pruneAndDangle(node_type * nd) {
            assert(findRoot(nd) == root);
            auto p = nd->getParent();
            if (p == nullptr) {
                root = nullptr;
                return;
            }
            p->removeChild(nd);
        }
        void pruneAndDelete(node_type * nd) {
            auto nodes = getSubtreeNodes(nd);
            pruneAndDangle(nd);
            for (auto ndi: nodes) {
                delete ndi;
            }
        }
        bool isDetached(node_type * nd) {
            return contains(detached, nd);
        }
    protected:
        node_type * root;
        U data;
        std::string name;
        std::set<node_type *> detached;
        
    public:
        void setName(const std::string &n) {
            name.assign(n);
        }
        const std::string & get_name() const {
            return name;
        }
        node_type * allocNewNode(node_type *p) {
            node_type * nd = new node_type(p);
            return nd;
        }
        void clear() {
            for(auto nd: getAllAttachedNodes()) {
                delete nd;
            }
            root = NULL;
        }
        std::vector<const node_type *> getAllAttachedNodes() const {
            return getSubtreeNodes(root);
        }
          
        std::vector<const node_type *> getSubtreeNodes(node_type* p) const {
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
        std::set<const node_type *> getSetOfAllAttachedNodes() const {
            std::set<const node_type*> nodes;
            for(auto nd: getAllAttachedNodes())
                nodes.insert(nd);
            return nodes;
        }
        void markAsDetached(node_type * nd) {
            detached.insert(nd);
        }
        void markAsAttached(node_type * nd) {
          detached.erase(nd);
        }
    private:
        RootedTree<T, U>(const RootedTree<T, U> &) = delete;
        RootedTree<T, U> & operator=(const RootedTree<T, U> &) = delete;
};

class RTNodeNoData{};
class RTreeNoData{};

template<typename Tree>
void addSubtree(typename Tree::node_type* par, Tree& T2)
{
    auto c = T2.get_root();
    T2.pruneAndDangle(c);
    par->addChild(c);
}

template<typename Tree>
void replaceWithSubtree(typename Tree::node_type* n, Tree& T2)
{
    // Get the parent of the tip we are replacing
    auto p = n->getParent();
    // Remove the data from T2 and attach it to this parent
    auto c = T2.get_root();
    T2.pruneAndDangle(c);
    p->addChild(c);
    // Remove the old child from under p
    p->removeChild(n);
    // Make the old child the root of T2, which will handle deleting it
    T2._set_root(n);
}

} // namespace otc
#endif

