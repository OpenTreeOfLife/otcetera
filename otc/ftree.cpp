#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/debug.h"

#include "otc/write_dot.h"
namespace otc {

bool PhyloStatement::debugCheck() const {
#ifdef DEBUGGING_PHYLO_STATEMENTS
    const OttIdSet ie = set_union_as_set(includeGroup, excludeGroup);
    if (ie != leafSet) {
        dbWriteOttSet(" includeGroup ",includeGroup);
        dbWriteOttSet(" excludeGroup ", excludeGroup);
        dbWriteOttSet(" leafSet ", leafSet);
        assert(false);
    }
#endif
    return true;
}

template<typename T, typename U>
void FTree<T, U>::createDeeperRoot() {
    auto nr = forest.createNode(nullptr);
    nr->addChild(root);
    nr->getData().desIds = root->getData().desIds;
    root = nr;
}

template<typename T, typename U>
const OttIdSet FTree<T, U>::getConnectedOttIds() const {
    OttIdSet r;
    if (root == nullptr) {
        return r;
    }
    for (auto t : iter_leaf_n_const(*root)) {
        assert(isAncestorDesNoIter(root, t));
        r.insert(t->getOttId());
    }
    return r;
}


template<typename T, typename U>
void FTree<T, U>::addSubtree(RootedTreeNode<T> * subtreeRoot,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedInc,
                const std::map<node_type *, std::set<RootedTreeNode<T> *> > & otherInvertedExc,
                const std::map<node_type *, std::pair<node_type*, PhyloStatementSource> > & otherIC,
                const std::map<node_type *, std::list<std::pair<node_type*, PhyloStatementSource> > > & otherEC) {
    assert(subtreeRoot);
    auto par = subtreeRoot->getParent();
    assert(par);
    const OttIdSet & sroids = subtreeRoot->getData().desIds;
    OttIdSet overlapIds = set_intersection_as_set(root->getData().desIds, sroids);
    if (overlapIds.empty()) {
        //the new tree does not have any leaves with include statements in the current tree...
        bool needDeeperRoot = false;
        for (auto oid : sroids) {
            if (ottIdIsExcludedFromRoot(oid)) {
                needDeeperRoot = true;
                break;
            }
        }
        if (needDeeperRoot) {
            createDeeperRoot();
        }
        root->addChild(subtreeRoot);
        for (auto c : iter_leaf_n_const(*subtreeRoot)) {
            //connectedIds.insert(c->getOttId());
        }
        root->getData().desIds.insert(begin(subtreeRoot->getData().desIds), end(subtreeRoot->getData().desIds));
        const PhyloStatementSource bogusPSS{-1, -1};
        for (auto sn : iter_pre_n(subtreeRoot)) { //TMP need iter_node_n
            const auto esn = otherInvertedExc.find(sn);
            if (esn != otherInvertedExc.end()) {
                for (auto eValNd : esn->second) {
                    excludesConstraints[eValNd].push_back(GroupingConstraint{sn, bogusPSS});
                }
            }
            const auto isn = otherInvertedInc.find(sn);
            if (isn != otherInvertedInc.end()) {
                for (auto eValNd : isn->second) {
                    auto pinP = includesConstraints.find(eValNd);
                    if (pinP == includesConstraints.end() || isAncestorDesNoIter(pinP->first, eValNd)) {
                        includesConstraints.emplace(eValNd, GroupingConstraint{sn, bogusPSS});
                    }
                }
            }
        }
        return;
    }
    assert(false);
}

template<typename T, typename U>
bool FTree<T,U>::isExcludedFromRoot(const node_type *n) const {
    node_type * ncn = const_cast<node_type *>(n);
    auto nit = excludesConstraints.find(ncn);
    if (nit == excludesConstraints.end()) {
        return false;
    }
    for (const auto & en : nit->second) {
        if (en.first == root) {
            return true;
        }
    }
    return false;
}
template<typename T, typename U>
bool FTree<T,U>::anyExcludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    const node_type * p = nd->getParent();
    const OttIdSet & ndi = nd->getData().desIds;
    for (auto oid : ottIdSet) {
        auto leafNd = ottIdToNodeMap.at(oid);
        auto gcIt = excludesConstraints.find(leafNd);
        if (gcIt != excludesConstraints.end()) {
            for (const auto & gc : gcIt->second) {
                node_type * en = gc.first;
                if (en == nd || en == p || isSubset(ndi, en->getData().desIds)) {
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename T, typename U>
bool FTree<T,U>::anyIncludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    if (!areDisjoint(nd->getData().desIds, ottIdSet)) {
        return true;
    }
    
    auto c = anyForceIncludedAtNode(nd, ottIdSet);
    assert(c == false);// if we are correctly updating desIds we don't need this branch.... TMP
    return c;
}

template<typename T, typename U>
bool FTree<T,U>::anyForceIncludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    for (auto oid :ottIdSet) {
        auto oidN = ottIdToNodeMap.at(oid);
        auto iclIt = includesConstraints.find(oidN);
        if (iclIt != includesConstraints.end()) {
            if (iclIt->second.first == oidN) {
                return true;
            }
        }
    }
    return false;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId) {
    //connectedIds.insert(ottId);
    return forest.createLeaf(par, ottId);
}

template<typename T, typename U>
void FTree<T,U>::addExcludeStatement(long ottId, RootedTreeNode<T> * excludedFrom, const PhyloStatementSource &pss) {
    OttIdSet x;
    x.insert(ottId);
    if (anyExcludedAtNode(excludedFrom, x)) {
        return; // already excluded from this subtree
    }
    RootedTreeNode<T> * eNode = ottIdToNodeMap.at(ottId);
    // If any of the descendants exclude this node, we can remove those exclude statements,
    //  because they'll be "dominated by this one"
    auto ecIt = excludesConstraints.find(eNode);
    if (ecIt != excludesConstraints.end()) {
        auto & listOfExc = ecIt->second;
        auto efIt = begin(listOfExc);
        for (; efIt != end(listOfExc);) {
            auto aen = efIt->first;
            if (isAncestorDesNoIter(excludedFrom, aen)) {
                efIt = listOfExc.erase(efIt);
            } else {
                ++efIt;
            }
        }
    }
    excludesConstraints[eNode].push_back(GroupingConstraint(excludedFrom, pss));
}

template<typename T, typename U>
void FTree<T,U>::addIncludeStatement(long ottId,
                                     RootedTreeNode<T> * includedIn,
                                     const PhyloStatementSource &pss) {
    assert(includedIn != nullptr);
    if (contains(includedIn->getData().desIds, ottId)) {
        LOG(DEBUG) << "  addIncludeStatement early return";
        return; // must be included in a des
    }
    RootedTreeNode<T> * eNode = ottIdToNodeMap.at(ottId);
    // If any of the ancestors include this node, we can remove those include statements,
    //  because they'll be "dominated by this one"
    auto icIt = includesConstraints.find(eNode);
    if (icIt != includesConstraints.end()) {
        auto & gcPair = icIt->second;
        auto aen = gcPair.first;
        if (aen == includedIn || isAncestorDesNoIter(includedIn, aen)) {
            LOG(DEBUG) << "  addIncludeStatement is redundant early return";
            return; // prev constraint at least as specific
        }
        if (!isAncestorDesNoIter(aen, includedIn)) {
            LOG(ERROR) << "Cannot move an inclusion constraint to a different subtree";
            assert(false);
            throw OTCError("err in FTree");
        }
    }
    includesConstraints.emplace(eNode, GroupingConstraint(includedIn, pss));
    // Since we know that the node will be a descendant of includedIn we add its Id to desIds
    includedIn->getData().desIds.insert(ottId);
    for (auto anc : iter_anc(*includedIn)) {
        anc->getData().desIds.insert(ottId);
    }
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par, const OttIdSet & oids) {
    dbWriteOttSet("  resolveToCreateCladeOfIncluded oids = ", oids);
    dbWriteOttSet("                                 nd->getData().desIds = ", par->getData().desIds);
    std::set<RootedTreeNode<T> *> cToMove;
    std::list<RootedTreeNode<T> *> orderedToMove;
    std::list<GroupingConstraint *> incToUpdate;
    for (auto oid : oids) {
        auto n = ottIdToNodeMap.at(oid);
        bool connectionFound = false;
        if (n->getParent() == par) {
            cToMove.insert(n);
            orderedToMove.push_back(n);
            connectionFound = true;
        } else {
            for (auto anc : iter_anc(*n)) {
                if (anc->getParent() == par) {
                    if (!contains(cToMove, anc)) {
                        cToMove.insert(anc);
                        orderedToMove.push_back(anc);
                        connectionFound = true;
                        break;
                    }
                }
            }
        }
        if (connectionFound) {
            continue;
        }
        auto igcIt = includesConstraints.find(n);
        if (igcIt != includesConstraints.end()) {
            auto np = igcIt->second.first;
            if (np == par) {
                incToUpdate.push_back(&(igcIt->second));
            }
        }
    }
    assert(cToMove.size() > 0 || incToUpdate.size() > 0);
    
    auto newNode = forest.createNode(par); // parent of includeGroup
    for (auto c : orderedToMove) {
        c->_detachThisNode();
        c->_setNextSib(nullptr);
        newNode->addChild(c);
        const auto & di = c->getData().desIds;
        newNode->getData().desIds.insert(begin(di), end(di));
    }
    for (auto gcp : incToUpdate) {
        gcp->first = newNode;
    }
    assert(!par->isOutDegreeOneNode());
    debugInvariantsCheck();
    return newNode;
}
template<typename T, typename U>
OttIdSet FTree<T,U>::addPhyloStatementAtNode(const PhyloStatement & ps, 
                                             RootedTreeNode<T> * includeGroupA,
                                             const OttIdSet & attachedElsewhere) {
    dbWriteOttSet(" FTree<T,U>::addPhyloStatementAtNode inc", ps.includeGroup);
    LOG(DEBUG) << "includeGroupA = " << (long) includeGroupA;
    dbWriteOttSet("    includeGroupA->getData().desIds", includeGroupA->getData().desIds);
    OttIdSet r;
    for (auto oid : ps.includeGroup) {
        if (!ottIdIsConnected(oid)) {
            LOG(DEBUG) << " not connected " << oid;
            if (contains(attachedElsewhere, oid)) {
                LOG(DEBUG) << " attachedElsewhere " << oid;
                addIncludeStatement(oid, includeGroupA, ps.provenance);
                auto ifIt = includesConstraints.find(ottIdToNodeMap.at(oid));
                if (ifIt != includesConstraints.end()) {
                    LOG(DEBUG) << "includeGroupA = " << (long) includeGroupA << " ic = " << (long) ifIt->second.first << " isAnc" << isAncestorDesNoIter(includeGroupA, ifIt->second.first);
                }
            } else {
                LOG(DEBUG) << " adding leaf " << oid;
                addLeafNoDesUpdate(includeGroupA, oid);
                r.insert(oid);
            }
        } else {
            LOG(DEBUG) << " connected " << oid;
        }
    }
    includeGroupA->getData().desIds.insert(begin(ps.includeGroup), end(ps.includeGroup));
    dbWriteOttSet("    later includeGroupA->getData().desIds", includeGroupA->getData().desIds);
    for (auto anc : iter_anc(*includeGroupA)) {
        anc->getData().desIds.insert(begin(ps.includeGroup), end(ps.includeGroup));
    }
    for (auto oid : ps.excludeGroup) {
        if (!ottIdIsConnected(oid)) {
            addExcludeStatement(oid, includeGroupA, ps.provenance);
        }
    }
    LOG(DEBUG) << "before addPhyloStatementAtNode exit";
    debugInvariantsCheck();
    LOG(DEBUG) << "foreset level check";
    forest.debugInvariantsCheck();
    LOG(DEBUG) << "bout to addPhyloStatementAtNode exit";
    return r;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::getMRCA(const OttIdSet &ottIdSet) {
    if (ottIdSet.empty()) {
        assert(false);
        return nullptr;
    }
    checkAllNodePointersIter(*root);
    const auto con = getConnectedOttIds();
    const auto rel = set_intersection_as_set(ottIdSet, con);
    dbWriteOttSet(" getMRCA ingroup", ottIdSet);
    dbWriteOttSet(" getMRCA connected ", con);
    dbWriteOttSet(" getMRCA connected ingroup", rel);
    const auto & relCheck = root->getData().desIds;
    dbWriteOttSet(" getMRCA relCheck", relCheck);
    assert(isSubset(rel, relCheck));
    for (auto nextOttId : rel) {
        auto x = ottIdToNodeMap.find(nextOttId);
        assert(x != ottIdToNodeMap.end());
        node_type * aTip = x->second;
        if (!isAncestorDesNoIter(root, aTip)) {
            LOG(ERROR) << "aTip->getOttId() = " << aTip->getOttId();
            for (auto a : iter_anc(*aTip)) {
                LOG(ERROR) << " anc address =  " << long(a);
            }
            LOG(ERROR) << " root address = " << long(root);
            assert(false);
        }
        assert(aTip != nullptr);
        if (ottIdSet.size() == 1) {
            return aTip;
        }
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    // reach here if none are connected.
    for (auto oid : ottIdSet) {
        auto x = ottIdToNodeMap.find(oid);
        assert(x != ottIdToNodeMap.end());
        node_type * aTip = x->second;
        assert(aTip != nullptr);
        if (ottIdSet.size() == 1) {
            return aTip;
        }
        return searchAncForMRCAOfDesIds(aTip, rel);
    }
    assert(false);
    return nullptr;
}

template<typename T, typename U>
void FTree<T, U>::mirrorPhyloStatement(const PhyloStatement &ps) {
    assert(root == nullptr);
    root = forest.createNode(nullptr);
    addPhyloStatementAsChildOfRoot(ps);
}

template<typename T, typename U>
void FTree<T, U>::addPhyloStatementAsChildOfRoot(const PhyloStatement &ps) {
    dbWriteOttSet(" addPhyloStatementAsChildOfRoot", ps.includeGroup);
    assert(root != nullptr);
    if (!root->isTip()) {
        debugInvariantsCheck();
    }
    if (anyExcludedAtNode(root, ps.includeGroup)) {
        createDeeperRoot();
        assert(!root->isTip());
    }
    assert(root != nullptr);
    auto parOfIncGroup = forest.createNode(root); // parent of includeGroup
    supportedBy[parOfIncGroup].push_back(ps.provenance);
    assert(ps.excludeGroup.size() > 0);
    for (auto i : ps.excludeGroup) {
        if (!forest.isAttached(i)) { // greedy
            addLeafNoDesUpdate(root, i);
            root->getData().desIds.insert(i);
        } else {
            addExcludeStatement(i, parOfIncGroup, ps.provenance);
        }
    }
    for (auto i : ps.includeGroup) {
        assert(!forest.isAttached(i));
        addLeafNoDesUpdate(parOfIncGroup, i);
    }
    root->getData().desIds.insert(begin(ps.includeGroup), end(ps.includeGroup));
    parOfIncGroup->getData().desIds = ps.includeGroup;
    debugInvariantsCheck();
    LOG(DEBUG) << "Leaving addPhyloStatementAsChildOfRoot";
}

template<typename T>
std::set<T *> getAncSet(T *nd) {
    std::set<T *> r;
    T * p = nd->getParent();
    while (p != nullptr) {
        r.insert(p);
        p = p->getParent();
    }
    return r;
}
template<typename T, typename U>
void FTree<T, U>::debugVerifyDesIdsAssumingDes(const OttIdSet &s, const RootedTreeNode<T> *nd) const{
    OttIdSet ois;
    if (nd->isTip()) {
        ois.insert(nd->getOttId());
    } else {
        for (auto c : iter_child_const(*nd)) {
            const auto & coids = c->getData().desIds;
            ois.insert(begin(coids), end(coids));
        }
    }
    for (const auto &icnS : includesConstraints) {
        auto o = icnS.first->getOttId();
        if (!contains(ois, o)) {
            if (icnS.second.first == nd) {
                ois.insert(o);
            } else {
                assert(!isAncestorDesNoIter(nd, icnS.second.first));
            }
        }
    }
    
    if(s != ois) {
        auto tn = ottIdToNodeMap.find(612446);
        if (tn != ottIdToNodeMap.end()) {
            auto icTn = includesConstraints.find(tn->second);
            if (icTn != includesConstraints.end()) {
                LOG(DEBUG) << "nd = " << (long)nd << " tn = " << (long)tn->second << " ic = " << (long) icTn->second.first << " isAnc" << isAncestorDesNoIter(nd, icTn->second.first);
            }
        }
        dbWriteOttSet("debugVerifyDesIdsAssumingDes incoming s", s);
        dbWriteOttSet("calculated:", ois);
        assert(s == ois);
    }
}
template<typename T, typename U>
void FTree<T, U>::debugInvariantsCheck() const {
    for (auto n : iter_post_n_const(*root)) {
        OttIdSet noids;
        if (n->isTip()) {
            const auto o = n->getOttId();
            assert(ottIdToNodeMap.at(o) == n);
        } else {
            assert(!n->hasOttId());
        }
        if (n != root) {
            assert(isAncestorDesNoIter(root, n));
        }
        debugVerifyDesIdsAssumingDes(n->getData().desIds, n);
        // Make sure that our ancestors do not exclude us.
        auto exIt = excludesConstraints.find(const_cast<node_type *>(n));
        if (exIt != excludesConstraints.end()) {
            const std::set<const node_type *> ancSet = getAncSet(n);
            for (auto & ex : exIt->second) {
                const node_type * exNd = ex.first;
                assert(!contains(ancSet, exNd));
            }
        }
    }
}

template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
