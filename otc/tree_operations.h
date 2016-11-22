#ifndef OTCETERA_TREE_OPERATIONS_H
#define OTCETERA_TREE_OPERATIONS_H
// Functions that operate on trees - may include iteration
//  over trees (contrast w/tree_util.h)
// Depends on: tree.h tree_util.h tree_iter.h 
// Depended on by: tools
#include <vector>
#include <unordered_map>
#include "otc/otc_base_includes.h"
#include "otc/tree_iter.h"
#include "otc/error.h"
#include "otc/util.h"
#include "otc/debug.h"
namespace otc {
template<typename T, typename CONTAINER>
void show_children_in_set(std::ostream & out, T nd, const CONTAINER & ancestral);
template<typename T, typename CONTAINER>
int count_children_in_set(T nd, const CONTAINER & ancestral);
template<typename T, typename CONTAINER>
bool is_monotypic_in_set(T nd, const CONTAINER & ancestral);
template<typename T>
void collapse_split_and_del_node(T* nd);
template<typename T>
long n_internal_with_ott_id(const T& tree);
template<typename T>
long n_internal(const T& tree);
template<typename T>
long n_internal_out_degree_1(const T& tree);
template<typename T>
long n_internal_out_degree_many(const T & tree);
template<typename T>
long n_nodes(const T& tree);
template<typename T>
std::string newick(const T &t);
template <typename T>
T* bisect_branch_with_new_child(T* x);


template<typename T>
unsigned int countPolytomies(const T & tree);
template<typename T>
std::size_t checkForUnknownTaxa(std::ostream & err, const T & toCheck, const T & taxonomy);
template<typename T>
typename T::node_type * findMRCAFromIDSet(T & tree, const std::set<long> & idSet, long trigger);

template<typename T>
void pruneAndDelete(T & tree, typename T::node_type *toDel);
template<typename T, typename U>
inline void cullRefsToNodeFromData(RootedTree<T, U> & tree, RootedTreeNode<T> *toDel);

const std::set<long> & getDesOttIds(RootedTreeNode<RTSplits> & nd);
template<typename T>
std::size_t pruneTipsWithoutIds(T & tree);
template <typename T>
std::vector<typename T::node_type*> all_nodes(T& tree);
template <typename N>
std::vector<N *> all_children(N * node);


//// impl
template <typename N>
inline std::vector<N *> all_children(N * node) {
    std::vector<N *> r;
    for (auto c : iter_child(*node)) {
        r.push_back(c);
    }
    return r;
}

// breaks the branch from x to its parent by allocating
//  a new node (which is returned) and assigning all of
//  the children of `x` to that node.
//  The returned node will be the only child of `x` and
//  all of the original children of x will be the children
//  of the returned node.
template<typename T>
T* bisect_branch_with_new_child(T* x) {
    assert(x);
    auto xc = new T(x);
    while (x->getFirstChild()) {
        auto nd = x->getFirstChild();
        nd->detachThisNode();
        xc->addChild(nd);
    }
    assert(not x->hasChildren());
    x->addChild(xc);
    return xc;
}

// breaks the branch from x to its parent by allocating
//  a new node (which is returned) and putting it in between
//  `x` and the parent of `x`
template<typename T>
T* bisect_branch_with_new_parent(T* x) {
    assert(x);
    auto xc = new T(x);
    while (x->getFirstChild()) {
        auto nd = x->getFirstChild();
        nd->detachThisNode();
        xc->addChild(nd);
    }
    assert(not x->hasChildren());
    x->addChild(xc);
    return xc;
}
template<typename T, typename CONTAINER>
inline void show_children_in_set(std::ostream & out, T nd, const CONTAINER& ancestral) {
    out << nd->getName() << " children: ";
    for (auto c = nd->getFirstChild(); c; c = c->getNextSib()) {
        if (ancestral.count(c)) {
            out << "'" << c->getName() << "' ";
        }
    }
    std::cerr << std::endl;
}

template<typename T, typename CONTAINER>
inline int count_children_in_set(T nd, const CONTAINER & ancestral) {
    int count = 0;
    for(auto c = nd->getFirstChild(); c; c = c->getNextSib()) {
        if (ancestral.count(c)) {
            count++;
        }
    }
    return count;
}

template<typename T, typename CONTAINER>
inline bool is_monotypic_in_set(T nd, const CONTAINER & ancestral) {
    return (count_children_in_set(nd, ancestral) == 1);
}

template<typename T>
void collapse_split_and_del_node(T* nd) {
    while (nd->getFirstChild()) {
        auto nd2 = nd->getFirstChild();
        nd2->detachThisNode();
        nd->addSibOnLeft(nd2);
    }
    nd->detachThisNode();
    delete nd;
}

template<typename T>
inline long n_internal_with_ott_id(const T& tree) {
    long count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->isTip() and nd->hasOttId()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline long n_internal(const T& tree) {
    long count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->isTip()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline long n_internal_out_degree_1(const T& tree) {
    long count = 0;
    for (auto nd: iter_post_const(tree)) {
        if (not nd->isTip() and nd->isOutDegreeOneNode()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline long n_internal_out_degree_many(const T & tree) {
    long count = 0;
    for(auto nd: iter_post_const(tree)) {
        if (not nd->isTip() and not nd->isOutDegreeOneNode()) {
            count++;
        }
    }
    return count;
}

template<typename T>
inline long n_nodes(const T& tree) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    long count = 0;
    for(auto nd: iter_post_const(tree)){
        count++;
    }
    return count;
}

template <typename T>
inline std::vector<typename T::node_type*> all_nodes(T& tree) {
    std::vector<typename T::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        nodes.push_back(nd);
    }
    return nodes;
}

template <typename T>
inline std::vector<typename T::node_type*> all_internal_nodes_post(T& tree) {
    std::vector<typename T::node_type*> nodes;
    for(auto nd: iter_post(tree)) {
        if (!nd->isTip()) {
            nodes.push_back(nd);
        }
    }
    return nodes;
}


template<typename T>
unsigned int countPolytomies(const T & tree) {
    unsigned int n = 0U;
    for (auto node : iter_post_internal_const(tree)) {
        if (node->getOutDegree() > 2) {
            n += 1;
        }
    }
    return n;
}

template<typename T>
void fillDesIdSets(T & tree) {
    // assumes OttId is set for each tip
    tree.getData().desIdSetsContainInternals = false;
    for (auto node : iter_post(tree)) {
        std::set<long> & desIds = node->getData().desIds;
        if (node->isTip()) {
            if (node->hasOttId()) {
                    desIds.insert(node->getOttId());
            }
        } else {
            for (auto child : iter_child(*node)) {
                std::set<long> & cDesIds = child->getData().desIds;
                desIds.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}
template<typename T>
void clearAndfillDesIdSets(T & tree) {
    // assumes OttId is set for each tip
    tree.getData().desIdSetsContainInternals = false;
    for (auto node : iter_post(tree)) {
        std::set<long> & desIds = node->getData().desIds;
        desIds.clear();
        if (node->isTip()) {
            desIds.insert(node->getOttId());
        } else {
            for (auto child : iter_child(*node)) {
                std::set<long> & cDesIds = child->getData().desIds;
                desIds.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}

template<typename T>
void fillDesIdSetsIncludingInternals(T & tree) {
    // assumes OttId is set for each tip
    tree.getData().desIdSetsContainInternals = true;
    for (auto node : iter_post(tree)) {
        std::set<long> & desIds = node->getData().desIds;
        if (node->isTip()) {
            desIds.insert(node->getOttId());
        } else {
            if (node->hasOttId()) {
                desIds.insert(node->getOttId());
            }
            for (auto child : iter_child(*node)) {
                std::set<long> & cDesIds = child->getData().desIds;
                desIds.insert(cDesIds.begin(), cDesIds.end());
            }
        }
    }
}

// fills in the depth data member for each node.
template <typename T>
void computeDepth(T& tree) {
    tree.getRoot()->getData().depth = 1;
    for (auto nd: iter_pre(tree)) {
        if (nd->getParent()) {
            nd->getData().depth = nd->getParent()->getData().depth + 1;
        }
    }
}

// uses ottID->node mapping, but not the split sets of the nodes
template<typename T>
typename T::node_type * findMRCAFromIDSet(T & tree, const std::set<long> & idSet, long trigger) {
    typedef typename T::node_type NT_t;
    auto ottIdToNode = tree.getData().ottIdToNode;
    std::map<NT_t *, unsigned int> n2c;
    long shortestPathLen = -1;
    NT_t * shortestPathNode = nullptr;
    for (const auto & i : idSet) {
        const auto rIt = ottIdToNode.find(i);
        if (rIt == ottIdToNode.end()) {
            std::string em = "tip ";
            em += std::to_string(i);
            if (trigger >= 0) {
                em += " a descendant of ";
                em += std::to_string(trigger); 
            }
            em += " not found.";
            throw OTCError(em);
        }
        auto nd = rIt->second;
        long currPathLen = 0;
        while (nd != nullptr) {
            n2c[nd] += 1;
            currPathLen += 1;
            nd = nd->getParent();
        }
        if (shortestPathLen < 0 || currPathLen < shortestPathLen) {
            shortestPathLen = currPathLen;
            shortestPathNode = rIt->second;
        }
    }
    auto nTips = idSet.size();
    auto cn = shortestPathNode;
    while (true) {
        if (n2c[cn] == nTips) {
            return cn;
        }
        cn = cn->getParent();
        assert(cn != nullptr);
    }
    assert(false);
    throw OTCError("asserts disabled, but false");
}

template<typename T>
std::size_t checkForUnknownTaxa(std::ostream & err, const T & toCheck, const T & taxonomy) {
    const auto & taxOttIds = taxonomy.getRoot()->getData().desIds;
    const auto & toCheckOttIds = toCheck.getRoot()->getData().desIds;
    const auto extras = set_sym_difference_as_set(toCheckOttIds, taxOttIds);
    if (!extras.empty()) {
        err << "OTT Ids found in an input tree,  but not in the taxonomy:\n";
        writeOttSet(err, "  ", extras, "\n");
        return extras.size();
    }
    return 0U;
}

template<typename T>
inline typename T::node_type * addChildForOttId(typename T::node_type & nd, long ottId, T & tree) {
    auto nn = tree.createChild(&nd);
    nn->setOttId(ottId);
    tree.getData()->ottIdToNode[ottId] = nn;
    return nn;
}

template<>
inline NodeWithSplits * addChildForOttId<TreeMappedWithSplits>(NodeWithSplits & nd, long ottId, TreeMappedWithSplits & tree) {
    auto nn = tree.createChild(&nd);
    nn->setOttId(ottId);
    tree.getData().ottIdToNode[ottId] = nn;
    return nn;
}

inline const std::set<long> & getDesOttIds(RootedTreeNode<RTSplits> & nd) {
    return nd.getData().desIds;
}

inline void fixDesIdFields(RootedTreeNode<RTSplits> & nd, const std::set<long> & ls) {
    const std::set<long> toRemove = nd.getData().desIds;
    assert(!toRemove.empty());
    nd.getData().desIds = ls;
    for (auto anc : iter_anc(nd)) {
        assert(anc != nullptr);
        assert(!anc->getData().desIds.empty());
        for (auto tr : toRemove) {
            assert(contains(anc->getData().desIds, tr));
            anc->getData().desIds.erase(tr);
        }
        anc->getData().desIds.insert(begin(ls), end(ls));
    }
}

// throws an OTCError if a tip is mapped to a non-terminal taxon
template<typename N, typename T, typename M>
void requireTipsToBeMappedToTerminalTaxa(const RootedTree<N,T> & toExpand, const M& ottIdToNode) {
    for (auto nd : iter_leaf_const(toExpand)) {
        assert(nd->isTip());
        assert(nd->hasOttId());
        auto ottId = nd->getOttId();

        if (not ottIdToNode.count(ottId))
            throw OTCError()<<"OTT Id "<<ottId<<" / Label '"<<nd->getName()<<"' not found in taxonomy!";

        auto taxNd = ottIdToNode.at(ottId);
        if (not taxNd->isTip())
            throw OTCError()<<"Tips must be mapped to terminal taxa: ott"<<ottId<<" found.";
    }
}

template<typename N, typename T>
void requireTipsToBeMappedToTerminalTaxa(const RootedTree<N,T> & toExpand, const RootedTree<N,T>& taxonomy) {
    std::map<long,const RootedTreeNode<N>*> ottIdToNode;
    for(auto nd: iter_post_const(taxonomy))
        if (nd->hasOttId())
            ottIdToNode[nd->getOttId()] = nd;
    
    requireTipsToBeMappedToTerminalTaxa(toExpand, ottIdToNode);
}

template<typename N>
void requireTipsToBeMappedToTerminalTaxa(const RootedTree<N,RTreeOttIDMapping<N>> & toExpand, const RootedTree<N,RTreeOttIDMapping<N>> & taxonomy) {
    const auto & taxData = taxonomy.getData();
    requireTipsToBeMappedToTerminalTaxa(toExpand, taxData.ottIdToNode);
}

// throws an OTCError if node has out-degree=1
template<typename T>
void requireNonRedundantTree(const T & toExpand) {
    std::map<typename T::node_type *, std::set<long> > replaceNodes;
    for (auto nd : iter_post_internal_const(toExpand)) {
        if (nd->isOutDegreeOneNode()) {
            std::string msg;
            msg = "Node with outdegree=1 found";
            if (nd->hasOttId()) {
                msg = msg + ":  ott" + std::to_string(nd->getOttId()) + ".";
            }
            throw OTCError(msg);
        }
    }
}


template<typename T>
std::set<const typename T::node_type *> expandOTTInternalsWhichAreLeaves(T & toExpand, const T & taxonomy) {
    const auto & taxData = taxonomy.getData();
    std::map<typename T::node_type *, std::set<long> > replaceNodes;
    for (auto nd : iter_leaf(toExpand)) {
        assert(nd->isTip());
        assert(nd->hasOttId());
        auto ottId = nd->getOttId();
        auto taxNd = taxData.getNodeForOttId(ottId);
        if (not taxNd)
        throw OTCError()<<"OTT Id "<<ottId<<" / Label '"<<nd->getName()<<"' not found in taxonomy!";
        if (!taxNd->isTip()) {
            const auto & leafSet = getDesOttIds(*taxNd);
            replaceNodes[nd] = leafSet;
        }
    }
    std::set<const typename T::node_type *> expanded;
    for (const auto & r : replaceNodes) {
        const auto & oldNode = r.first;
        const auto & ls = r.second;
        assert(ls.size() > 0);
        expanded.insert(oldNode);
        for (auto loid : ls) {
            addChildForOttId<T>(*oldNode, loid, toExpand);
        }
        fixDesIdFields(*oldNode, ls);
    }
    return expanded;
}

// Starts at the node in `fullTree` identified by `ottId`
// Walks back to the root, and adds `ottId` to value of the node->OttIdSet mapping `n2m`
// Creates new entries in that table, as needed.
// Multiple calls to this by tools like check-supertree to associate every node in a 
//      full tree with an OttIdSet which contains all of the descendants for a reduced
//      leafset.
template<typename T>
void markPathToRoot(const T & fullTree,
                    long ottId,
                    std::map<const typename T::node_type *, OttIdSet > &n2m){
    auto startNd = fullTree.getData().getNodeForOttId(ottId);
    if (startNd == nullptr) {
        std::string m = "OTT id not found ";
        m += std::to_string(ottId);
        throw OTCError(m);
    }
    assert(startNd != nullptr);
    n2m[startNd].insert(ottId);
    for (auto nd : iter_anc(*startNd)) {
        n2m[nd].insert(ottId);
    }
}


// find most recent anc of nd with out-degree > 1
template<typename T>
inline T * findFirstForkingAnc(T * nd) {
    T * anc = nd->getParent();
    if (anc == nullptr) {
        return nullptr;
    }
    while (anc->isOutDegreeOneNode()) {
        anc = anc->getParent();
        if (anc == nullptr) {
            return nullptr;
        }
    }
    return anc;
}

// find most recent anc of nd with out-degree > 1
template<typename T>
inline T * findFirstForkingSelfOrDes(T * nd) {
    assert(nd);
    T * des = nd;
    while (des->isOutDegreeOneNode()) {
        des = des->getFirstChild();
    }
    return des;
}


template<typename T, typename U>
inline bool multipleChildrenInMap(const T & nd,
                                  const std::map<const T *, U> & markedMap,
                                  const T **first) {
    assert(first);
    bool foundFirst = false;
    *first = nullptr;
    for(auto c : iter_child_const(nd)) {
        if (markedMap.find(c) != markedMap.end()) {
            if (foundFirst) {
                return true;
            }
            foundFirst = true;
            *first = c;
        }
    }
    return false;
}

template<typename T>
inline const T * findNextSignificantNode(const T * node, const std::map<const T *, std::set<long> > & markedMap) {
    assert(node != nullptr);
    auto currNode = node;
    for (;;) {
        const T * sc;
        if (multipleChildrenInMap(*currNode, markedMap, &sc)) {
            return currNode;
        }
        if (sc == 0L) {
            const char * msg = "Failing. Node found with ottIDs marked, but no children with ottIDs marked";
            throw OTCError(msg);
        }
        currNode = sc;
    }
}

template<typename T>
inline void writePrunedSubtreeNewickForMarkedNodes(std::ostream & out,
                                            const T & srcNd,
                                            const std::map<const T *, std::set<long> > & markedMap) {
    const auto nIt = markedMap.find(&srcNd);
    assert(nIt != markedMap .end());
    const auto & ottIDSet = nIt->second;
    if (ottIDSet.size() == 1) {
        out << "ott" << *ottIDSet.begin();
    } else {
        assert(ottIDSet.size() > 1);
        auto nsn = findNextSignificantNode<T>(&srcNd, markedMap);
        out << '(';
        unsigned numcwritten = 0;
        for (auto child : iter_child_const(*nsn)) {
            if (markedMap.find(child) != markedMap.end()){
                if (numcwritten > 0) {
                    out << ',';
                }
                writePrunedSubtreeNewickForMarkedNodes(out, *child, markedMap);
                ++numcwritten;
            }
        }
        assert(numcwritten > 1);
        out << ')';
    }
}



template<typename T>
inline void describeUnnamedNode(const T & nd,
                                std::ostream & out,
                                unsigned int anc,
                                bool useNdNames,
                                bool addNewLine=true) {
    if (useNdNames && !nd.getName().empty()) {
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before \"" << nd.getName() << "\"";
        } else {
            out << "the node \"" << nd.getName() << "\"";
        }
    }
    else if (nd.isTip()) {
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before the leaf \"" << nd.getName()  << "\"";
        } else {
            out << "the leaf \"" << nd.getName()  << "\"";
        }
    } else if (nd.isOutDegreeOneNode()) {
        describeUnnamedNode(*nd.getFirstChild(), out, anc + 1, useNdNames);
        return;
    } else {
        const auto & left = findLeftmostInSubtree(&nd)->getName();
        const auto & right = findRightmostInSubtree(&nd)->getName();
        if (anc > 0) {
            out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right << "\"";
        } else {
            out <<  "MRCA of \"" << left << "\" and " << "\"" << right << "\"";
        }
    }
    if (addNewLine) {
        out << '\n';
    }
}

template<typename T, typename U>
inline void cullRefsToNodeFromData(RootedTree<T, U> & , RootedTreeNode<T> *) {
}

template<>
inline void cullRefsToNodeFromData(RootedTree<RTNodeNoData, RTreeOttIDMapping<RTNodeNoData> > & tree, 
                                   RootedTreeNode<RTNodeNoData> *nd) {
    assert(nd != nullptr);
    if (nd->hasOttId()) {
        tree.getData().ottIdToNode.erase(nd->getOttId());
    }
}



template<>
inline void cullRefsToNodeFromData(RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > & tree, 
                                   RootedTreeNode<RTSplits> *nd)  {
    assert(nd != nullptr);
    if (nd->hasOttId()) {
        tree.getData().ottIdToNode.erase(nd->getOttId());
    }
    const auto & d = nd->getData().desIds;
    if (d.empty()) {
        return;
    }
    for (auto a : iter_anc(*nd)) {
        auto & ad = a->getData().desIds;
        for (auto el : d) {
            ad.erase(el);
        }
    }
}

// does NOT correct anc desIds
template<typename T, typename U>
inline void delOttIdOfInternal(RootedTree<T, U> & tree, RootedTreeNode<T> *nd) {
    if (nd->hasOttId()) {
        tree.getData().ottIdToNode.erase(nd->getOttId());
        nd->delOttId();
    }
}

// does NOT correct anc desIds
template<typename T, typename U>
inline void changeOttIdOfInternal(RootedTree<T, U> & tree, RootedTreeNode<T> *nd, OttId ottId) {
    delOttIdOfInternal(tree, nd);
    if (ottId > 0) {
        nd->setOttId(ottId);
        tree.getData().ottIdToNode[ottId] = nd;
    }
}

template<typename T, typename U>
inline void collapseNode(RootedTree<T, U> & tree, RootedTreeNode<T> *nd) {
    if (nd->hasOttId()) {
        delOttIdOfInternal(tree, nd);
    }
    collapseInternalIntoPar(nd, tree);
}


template<typename T>
inline void pruneAndDelete(T & tree, typename T::node_type *toDel) {
    cullRefsToNodeFromData<typename T::node_data_type, typename T::data_type>(tree, toDel);
    tree.pruneAndDelete(toDel);
}

template <typename T, typename U>
void insertAncestorsToParaphyleticSet(T * nd, U & includedNodes) {
    for (auto anc : iter_anc(*nd)) {
        if (contains(includedNodes, anc)) {
            return;
        }
        includedNodes.insert(anc);
    }
}

//@TMP recursive until we have a pre-order subtree skipping iter.
template <typename T, typename U>
void insertDescendantsOfUnincludedSubtrees(T * nd, U & includedNodes) {
    for (auto c : iter_child(*nd)) {
        if (!contains(includedNodes, c)) {
            includedNodes.insert(c);
            insertDescendantsOfUnincludedSubtrees(c, includedNodes);
        }
    }
}


template<typename T>
inline void writeNodeAsNewickLabel(std::ostream & out, const T *nd) {
    if (nd->isTip()) { //TEMP tip just like internals, but at some point we may want the label-less internals to be unlabelled, regardless of OTT ID
        const auto & n = nd->getName();
        if (n.empty()) {
            out << "ott" << nd->getOttId();
        } else {
            writeEscapedForNewick(out, nd->getName());
        }
    } else if (!nd->getName().empty()) {
        writeEscapedForNewick(out, nd->getName());
    } else if (nd->hasOttId()) {
        out << "ott" << nd->getOttId();
    }
}

template<typename T>
inline void writeClosingNewick(std::ostream & out, const T *nd, const T * r) {
    out << ')';
    auto n = nd->getParent();
    writeNodeAsNewickLabel(out, n);
    if (n == r) {
        return;
    }
    while (n->getNextSib() == nullptr) {
        out << ')';
        n = n->getParent();
        assert(n != nullptr);
        writeNodeAsNewickLabel(out, n);
        if (n == r) {
            return;
        }
    }
    out << ',';
}

template<typename T>
inline void writeNewick(std::ostream & out, const T *nd) {
    assert(nd != nullptr);
    if (nd->isTip()) {
        writeNodeAsNewickLabel(out, nd);
    } else {
        for (auto n : iter_pre_n_const(nd)) {
            if (n->isTip()) {
                writeNodeAsNewickLabel(out, n);
                if (n->getNextSib() == nullptr) {
                    writeClosingNewick<T>(out, n, nd);
                } else {
                    out << ',';
                }
            } else {
                out << '(';
            }
        }
    }
}
template<typename T>
inline void dbWriteNewick(const T *nd) {
    if (!debuggingOutputEnabled) {
        return;
    }
    writeNewick(std::cerr, nd);
    std::cerr << std::endl;
}


template<typename T>
inline bool isEffectivelyATip(const T *nd, std::function<bool(const T &)> subtreeFilter) {
    if (nd->isTip()) {
        return true;
    }
    for (auto child : iter_child_const(*nd)) {
        if (subtreeFilter(*child)) {
            return false;
        }
    }
    return true;
}


template<typename T>
inline bool isEffectivelyLastSib(const T *nd, std::function<bool(const T &)> subtreeFilter) {
    nd = nd->getNextSib();
    while (nd != nullptr) {
        if (subtreeFilter(*nd)) {
            return false;
        }
        nd = nd->getNextSib();
    }
    return true;
}

template<typename T>
inline void writeClosingNewickFiltered(std::ostream & out, const T *nd, const T * r, std::function<bool(const T &)> subtreeFilter) {
    assert(nd != nullptr);
    out << ')';
    auto n = nd->getParent();
    writeNodeAsNewickLabel(out, n);
    if (n == r) {
        return;
    }
    while (isEffectivelyLastSib(n, subtreeFilter)) {
        if (n == r) {
            return;
        }
        out << ')';
        n = n->getParent();
        assert(n != nullptr);
        writeNodeAsNewickLabel(out, n);
    }
    out << ',';
}

template<typename T>
inline void writeNewickFiltered(std::ostream & out, const T *nd, std::function<bool(const T &)> subtreeFilter) {
    assert(nd != nullptr);
    if (nd->isTip()) {
        if (!subtreeFilter(*nd)) {
            return;
        }
        writeNodeAsNewickLabel(out, nd);
    } else {
        for (auto n : iter_pre_filter_n_const(nd, subtreeFilter)) {
            if (isEffectivelyATip(n, subtreeFilter)) {
                writeNodeAsNewickLabel(out, n);
                if (isEffectivelyLastSib(n, subtreeFilter)) {
                    writeClosingNewickFiltered<T>(out, n, nd, subtreeFilter);
                } else {
                    out << ',';
                }
            } else {
                out << '(';
            }
        }
    }
}

template<typename T>
inline void writeTreeAsNewick(std::ostream & out, const T &tree) {
    writeNewick<typename T::node_type>(out, tree.getRoot());
    out << ';';
}

template<typename T>
inline T * searchAncForMRCAOfDesIds(T * nd, const std::set<long> & idSet) {
    assert(nd != nullptr);
    if (isSubset(idSet, nd->getData().desIds)) {
        return nd;
    }
    for (auto n : iter_anc(*nd)) {
        if (isSubset(idSet, n->getData().desIds)) {
            return n;
        }
    }
    return nullptr;
}

template<typename T>
inline const typename T::node_type * findNodeWithMatchingDesIdSet(const T & tree, const OttIdSet & idSet) {
    assert(!idSet.empty());
    OttId firstId = *begin(idSet);
    auto nd = tree.getData().ottIdToNode.at(firstId);
    assert(nd != nullptr);
    if (nd->getData().desIds == idSet) {
        return nd;
    }
    for (auto n : iter_anc_const(*nd)) {
        const auto & ndi = n->getData().desIds;
        if (ndi == idSet) {
            return n;
        }
        if (ndi.size() > idSet.size()) {
            return nullptr;
        }
    }
    return nullptr;
}

template<typename T>
inline const typename T::node_type * findMRCAUsingDesIds(const T & tree, const std::set<long> & idSet) {
    if (idSet.empty()) {
        assert(false);
        throw OTCError("asserts disabled but false");
    }
    const long lowestID = *idSet.begin();
    const typename T::node_type * aTip = tree.getData().getNodeForOttId(lowestID);
    if (aTip == nullptr) {
        return nullptr;
    }
    return searchAncForMRCAOfDesIds(aTip, idSet);
}

template<typename T>
inline void eraseMappingsToNode(typename T::node_type * nd, T & tree) {
    if (nd->hasOttId()) {
        tree.getData().ottIdToNode.erase(nd->getOttId());
    }
    tree.markAsDetached(nd);
}

template<typename N, typename T>
inline void replaceMappingsToNodeWithAlias(typename RootedTree<N,T>::node_type * , // nd
                                           typename RootedTree<N,T>::node_type * , //alias
                                           RootedTree<N,T> & ) {
}


template<typename N>
inline void replaceMappingsToNodeWithAlias(typename RootedTree<N,RTreeOttIDMapping<N>>::node_type * nd,
                                           typename RootedTree<N,RTreeOttIDMapping<N>>::node_type * alias,
                                           RootedTree<N,RTreeOttIDMapping<N>>& tree) {
    if (nd->hasOttId()) {
        tree.getData().ottIdToDetachedNode[nd->getOttId()] = nd;
        tree.getData().ottIdToNode[nd->getOttId()] = alias;
        auto & iaf = tree.getData().isAliasFor;
        iaf[alias].insert(nd->getOttId());
        auto na = iaf.find(nd); // if nd was an alias for other nodes, update them too...
        if (na != iaf.end()) {
            for (auto aoi : na->second) {
                iaf[alias].insert(aoi);
                tree.getData().ottIdToNode[aoi] = alias;
            }
        }
        iaf.erase(nd);
    }
    tree.markAsDetached(nd);
}

// 
// detaches child from the tree (and lets it dangle)
// assigns its children to nd (which should be the parent of child)
// note in the only usage: nd is monotypic, but child is NOT.
// this suppresses monotypy by getting rid of child in case nd has a
//  taxonomic assignment that we want to preserve in the desIds
template<typename T>
inline void suppressMonotypyByStealingGrandchildren(typename T::node_type * nd,
                                                    typename T::node_type * monotypic_child,
                                                    T & tree) {
    assert(nd);
    assert(monotypic_child);
    assert(monotypic_child->getParent() == nd);
    assert(nd->getFirstChild() == monotypic_child);
    assert(monotypic_child->getNextSib() == nullptr);
    LOG(DEBUG) << "suppressing " << (nd->hasOttId() ? nd->getOttId() : long(nd))
               << " by stealing children from " << getDesignator(*monotypic_child) ;

    while(monotypic_child->hasChildren())
    {
        auto x = monotypic_child->getFirstChild();
        x->detachThisNode();
        nd->addChild(x);
    }
    monotypic_child->detachThisNode();

    replaceMappingsToNodeWithAlias(monotypic_child, nd, tree);
    //    delete monotypic_child;
}

template<typename T>
inline void suppressMonotypyByClaimingGrandparentAsPar(typename T::node_type * monotypic_nd,
                                                       typename T::node_type * child,
                                                       T & tree) {
    assert(monotypic_nd);
    assert(child);
    assert(child->getParent() == monotypic_nd);
    assert(monotypic_nd->getFirstChild() == child);
    assert(child->getNextSib() == nullptr);
    auto gp = monotypic_nd->getParent();
    assert(gp);
    LOG(DEBUG) << "suppressing " << (monotypic_nd->hasOttId() ? monotypic_nd->getOttId() : long(monotypic_nd))
                << " by claiming grandparent as parent for node " << getDesignator(*child) ;

    monotypic_nd->detachThisNode();
    child->detachThisNode();
    gp->addChild(child);
    assert(not monotypic_nd->hasChildren());
    replaceMappingsToNodeWithAlias(monotypic_nd, child, tree);
    // delete monotypic_nd
}

template <typename T>
inline void delMonotypicNode(typename T::node_type* nd, T& tree) {
    assert(nd->isOutDegreeOneNode());
    auto child = nd->getFirstChild();
    child->detachThisNode();

    if (nd->getParent())
    {
        nd->addSibOnRight(child);
        nd->detachThisNode();
    }
    else
        tree._setRoot(child);

    delete nd;
 }

template <typename T>
void suppressMonotypicFast(T& tree)
{
    std::vector<typename T::node_type*> remove;
    for(auto nd:iter_pre(tree))
        if (nd->isOutDegreeOneNode())
            remove.push_back(nd);

    for(auto nd: remove)
        delMonotypicNode(nd, tree);
}
 
// nd must be unnnamed (or the aliases would not be fixed because we can't call replaceMappingsToNodeWithAlias)
template<typename T>
inline void collapseInternalIntoPar(typename T::node_type * nd,
                                    T & tree) {
    assert(nd);
    assert(not nd->hasOttId());
    assert(nd->hasChildren());

    if (nd->getOutDegree() == 1) {
        suppressMonotypyByClaimingGrandparentAsPar(nd, nd->getFirstChild(), tree);
        return;
    }
    
    while(nd->hasChildren())
    {
        auto child = nd->getFirstChild();
        child->detachThisNode();
        nd->addSibOnLeft(child);
    }
    nd->detachThisNode();
}

// in some context we need to preserve the deepest because the desIds may need 
// the internal nodes IDs (even if the IDs correspond to monotypic taxa)
// returns set of nodes deleted. Their getParent() calls will tell the caller
//  their former parent (before culling). But note that you may have to 
//  call it muliple times to find the first ancestor that has not been deleted.
template<typename T>
inline std::set<typename T::node_type *> suppressMonotypicTaxaPreserveDeepestDangle(
                    T & tree,
                    bool resetOttId) {
    std::set<typename T::node_type *> monotypic;
    for (auto nd : iter_node_internal(tree)) {
        if (nd->isOutDegreeOneNode()) {
            monotypic.insert(nd);
        }
    }
    std::set<typename T::node_type *> removed;
    auto toProcess = monotypic;
    while (!toProcess.empty()) {
        for (auto toPIt = toProcess.begin(); toPIt != toProcess.end();) {
            auto tpn = *toPIt;
            auto child = tpn->getFirstChild();
            if (!contains(toProcess, child)) {
                suppressMonotypyByStealingGrandchildren<T>(tpn, child, tree);
                toProcess.erase(toPIt++);
                if (resetOttId && child->hasOttId()) {
                    tpn->setOttId(child->getOttId());
                }
                removed.insert(child);
            } else {
                ++toPIt;
            }
        }
    }
    return removed;
}

// in some context we need to preserve the deepest because the desIds may need 
// the internal nodes IDs (even if the IDs correspond to monotypic taxa)
// returns set of nodes deleted. Their getParent() calls will tell the caller
//  their former parent (before culling). But note that you may have to 
//  call it muliple times to find the first ancestor that has not been deleted.
template<typename T>
inline std::set<typename T::node_type *> suppressMonotypicTaxaPreserveShallowDangle(
                    T & tree) {
    std::set<typename T::node_type *> monotypic;
    for (auto nd : iter_node_internal(tree)) {
        //assert(!nd->hasOttId()); // not a good place for this assert... true in how we use this fund
        if (nd->isOutDegreeOneNode()) {
            monotypic.insert(nd);
        }
    }
    std::set<typename T::node_type *> removed;
    auto toProcess = monotypic;
    while (!toProcess.empty()) {
        for (auto toPIt = toProcess.begin(); toPIt != toProcess.end();) {
            checkTreeInvariants(tree);
            auto tpn = *toPIt;
            auto par = tpn->getParent();
            if (!contains(toProcess, par)) {
                if (tree.isDetached(tpn)) {
                    removed.insert(tpn);
                } else if (tpn->isTip()) {
                    tree.pruneAndDelete(tpn);
                    removed.insert(tpn);
                } else {
                    auto onlyChild = tpn->getFirstChild();
                    if (par) {
                        suppressMonotypyByClaimingGrandparentAsPar<T>(tpn, onlyChild, tree);
                        removed.insert(tpn);
                    } else {
                        replaceMappingsToNodeWithAlias(onlyChild, tpn, tree);
                        if (onlyChild->isTip()) {
                            tree.pruneAndDelete(onlyChild);
                        } else {
                            onlyChild->delOttId();
                            collapseInternalIntoPar(onlyChild, tree);
                        }
                        removed.insert(onlyChild);
                    }
                }
                toProcess.erase(toPIt++);
            } else {
                ++toPIt;
            }
        }
    }
    return removed;
}


template<typename T>
inline bool isAncDecPair(const T * nd1, const T *nd2) {
    if (nd1 == nd2) {
        return false;
    }
    for (auto a : iter_anc_const(*nd2)) {
        if (a == nd1) {
            return true;
        }
    }
    return false;
}

// returns true if the nodes are identical or one is an ancestor of the other...
template<typename T>
inline bool areLinearlyRelated(const T * nd1, const T *nd2) {
    return nd1 == nd2 || isAncDecPair(nd1, nd2) || isAncDecPair(nd2, nd1);
}

template<typename T>
inline std::set<long> getOttIdSetForLeaves(const T &tree) {
    if (!tree.getData().desIdSetsContainInternals) {
        return tree.getRoot()->getData().desIds;
    }
    std::set<long> inducingIds;
    for (auto nd : iter_leaf_const(tree)) {
        inducingIds.insert(nd->getOttId());
    }
    return inducingIds;
}

template<typename T, typename U>
void getInducedInformativeGroupings(const T & tree1, std::set<std::set<long> > & inducedSplits, const U & tree2) {
    const auto inducingIds = getOttIdSetForLeaves(tree2);
    auto mrca = findMRCAUsingDesIds(tree1, inducingIds);
    std::function<bool(const typename T::node_type &)> sf = [inducingIds](const typename T::node_type &nd){
        return haveIntersection(inducingIds, nd.getData().desIds);
    };
    for (auto n : iter_pre_filter_n_const(mrca, sf)) {
        if (n == mrca) {
            continue;
        }
        auto inducedDesIds = n->getData().desIds;
        const auto x = intersectionOfSets(inducingIds, inducedDesIds);
        if (x.size() > 1) {
            inducedSplits.insert(std::move(x));
        }
    }
}

template<typename T>
void getInformativeGroupings(const T & tree2,
                             std::set<std::set<long> > & tree2Splits) {
    auto t2r = tree2.getRoot();
    for (auto n : iter_pre_internal_const(tree2)) {
        if (n == t2r) {
            continue;
        }
        const std::set<long> & x = n->getData().desIds;
        if (x.size() > 2) {
            tree2Splits.insert(x);
        }
    }
}

template<typename T, typename U>
void inducedCladeSets(const T & tree1,
                      const U & tree2,
                      std::set<std::set<long> > & inducedSplits,
                      std::set<std::set<long> > & tree2Splits,
                      bool firstIsSuperset) {
    assert(firstIsSuperset);
    getInducedInformativeGroupings(tree1, inducedSplits, tree2);
    if (firstIsSuperset) {
        getInformativeGroupings(tree2, tree2Splits);
    } else {
        throw OTCError("Why on earth did you compile without asserts?");
    }
}

template<typename T, typename U>
unsigned long inducedRFDist(const T & tree1, const U & tree2, bool firstIsSuperset) {
    std::set<std::set<long> > inducedSplits;
    std::set<std::set<long> > tree2Splits;
    inducedCladeSets(tree1, tree2, inducedSplits, tree2Splits, firstIsSuperset);
    return sizeOfSymmetricDifference(tree2Splits, inducedSplits);
}

template<typename T, typename U>
unsigned long numInducedSplitsMissingInSecond(const T & tree1,
                                              const U & tree2,
                                              bool firstIsSuperset) {
    std::set<std::set<long> > inducedSplits;
    std::set<std::set<long> > tree2Splits;
    inducedCladeSets(tree1, tree2, inducedSplits, tree2Splits, firstIsSuperset);
    unsigned long nm = 0;
    for (auto ics : inducedSplits) {
        if (!contains(tree2Splits, ics)) {
            nm += 1;
        }
    }
    return nm;
}

template<typename U, typename T, typename V>
inline void copyTreeStructure(const std::map<U *, U *> & nd2par,
                       const std::map<U *, long> & nd2id,
                       RootedTree<T, V> & toWrite) {
    std::set<RootedTreeNode<T> *> withParents;
    std::map<U *, RootedTreeNode<T> *> other2new;
    for (auto c : nd2par) {
        auto otherChild = c.first;
        auto otherParent = c.second;
        RootedTreeNode<T> * nParent;
        auto npIt = other2new.find(otherParent);
        if (npIt == other2new.end()) {
            nParent = toWrite.createNode(nullptr);
            other2new[otherParent] = nParent;
            if (otherParent->hasOttId())
                nParent->setOttId(otherParent->getOttId());
        } else {
            nParent = npIt->second;
        }
        RootedTreeNode<T> * nChild;
        auto ncIt = other2new.find(otherChild);
        if (ncIt == other2new.end()) {
            nChild = toWrite.createNode(nParent);
            other2new[otherChild] = nChild;
            if (otherChild->hasOttId())
                nChild->setOttId(otherChild->getOttId());
        } else {
            nChild = ncIt->second;
            nParent->addChild(nChild);
        }
        auto idIt = nd2id.find(otherChild);
        if (idIt != nd2id.end()) {
            nChild->setOttId(idIt->second);
        }
        idIt = nd2id.find(otherParent);
        if (idIt != nd2id.end()) {
            nParent->setOttId(idIt->second);
        }
        withParents.insert(nChild);
    }
    for (auto nn : withParents) {
        if ((nn->isTip()) && !nn->hasOttId()) {
            for (auto o2n : other2new) {
                if (o2n.second == nn) {
                    LOG(ERROR) << "tip without label in copyTreeStructure";
                    throw OTCError("tip without label");
                }
            }
        }
    }
    assert(other2new.size() == 1 + withParents.size()); // only the root should be parentless
    for (auto o2n : other2new) {
        if (!contains(withParents, o2n.second)) {
            toWrite._setRoot(o2n.second);
            return;
        }
    }
}

template<typename T>
inline std::size_t pruneTipsWithoutIds(T & tree) {
    std::size_t r = 0;
    std::set<typename T::node_type *> toPrune;
    std::set<typename T::node_type *> parToCheck;
    for (auto nd : iter_node(tree)) {
        if (nd->isTip() && !nd->hasOttId()) {
            toPrune.insert(nd);
            parToCheck.insert(nd->getParent());
        }
    }
    while (!toPrune.empty()) {
        r += toPrune.size();
        for (auto nd : toPrune) {
            pruneAndDelete(tree, nd);
        }
        toPrune.clear();
        std::set<typename T::node_type *> nextParToCheck;
        for (auto nd : parToCheck) {
            if (nd == nullptr) {
                tree._setRoot(nullptr);
                continue;
            }
            if (nd->isTip()) {
                toPrune.insert(nd);
                nextParToCheck.insert(nd->getParent());
            }
        }
        nextParToCheck.swap(parToCheck);
    }
    if (tree.getRoot() != nullptr) {
        suppressMonotypicTaxaPreserveShallowDangle(tree);
    }
    return r;
}


template<typename T>
inline std::string newick(const T &t) {
    std::ostringstream s;
    writeTreeAsNewick(s, t);
    return s.str();
}


template <typename Tree_t>
std::unordered_map<long, const typename Tree_t::node_type*> get_ottid_to_const_node_map(const Tree_t& T)
{
  std::unordered_map<long, const typename Tree_t::node_type*> ottid;
  for(auto nd: iter_pre_const(T))
    if (nd->hasOttId())
      ottid[nd->getOttId()] = nd;

  return ottid;
}

template <typename Tree_t>
std::unordered_map<long, typename Tree_t::node_type*> get_ottid_to_node_map(Tree_t& T)
{
    std::unordered_map<long, typename Tree_t::node_type*> ottid;
    for(auto nd: iter_pre(T))
        if (nd->hasOttId())
            ottid[nd->getOttId()] = nd;

    return ottid;
}

template <typename N>
N* get_root(N* node)
{
    while (node->getParent()){
        node = node->getParent();
    }
    return node;
}


 
}// namespace otc

template <typename node_t>
void destroy_children(node_t* node)
{
    std::vector<node_t*> nodes;
    while(auto n = node->getFirstChild())
    {
        n->detachThisNode();
        nodes.push_back(n);
    }

    for(std::size_t i = 0; i < nodes.size(); i++) {
        while(auto n = nodes[i]->getFirstChild())
        {
            n->detachThisNode();
            nodes.push_back(n);
        }
        delete nodes[i];
    }

    assert(node->isTip());
}


#endif

