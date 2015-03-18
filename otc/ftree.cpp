#include "otc/forest.h"
#include "otc/util.h"
#include "otc/tree_data.h"
#include "otc/supertree_util.h"
#include "otc/tree_operations.h"
#include "otc/debug.h"

#include "otc/write_dot.h"
namespace otc {

void PhyloStatement::writeAsNewick(std::ofstream & out) const {
    if (!excludeGroup.empty()) {
        out << '(';
    }
    out << "(";
    bool first = true;
    for (const auto & oid : includeGroup) {
        if (!first) {
            out << ',';
        }
        out << "ott" << oid;
        first = false;
    }
    out << ')';
    if (!excludeGroup.empty()) {
        for (const auto & oid : excludeGroup) {
            out << ",ott" << oid;
        }
        out << ")";
    }
    out << ';';
}
template<typename T>
bool ExcludeConstraints<T>::isExcludedFrom(const node_type * ndToCheck,
                                           const node_type * potentialAttachment,
                                           const std::map<long, node_type*> * o2n) const {
    const auto & ndi = ndToCheck->getData().desIds;
    if (ndi.size() == 1) {
        auto nit = byExcludedNd.find(ndToCheck);
        if (nit == byExcludedNd.end()) {
            return false;
        }
        return contains(nit->second, potentialAttachment);
    }
    assert(o2n != nullptr);
    for (auto oid : ndi) {
        auto tn = o2n->at(oid);
        auto nit = byExcludedNd.find(tn);
        if (nit != byExcludedNd.end() && contains(nit->second, potentialAttachment)) {
            return true;
        }
    }
    return false;
}

template<typename T>
bool ExcludeConstraints<T>::addExcludeStatement(const node_type * nd2Exclude,
                                                const node_type * forbiddenAttach) {
    if (isExcludedFrom(nd2Exclude, forbiddenAttach, nullptr)) {
        return false; // already excluded from this subtree
    }
    // If any of the descendants exclude this node, we can remove those exclude statements,
    //  because they'll be "dominated by this one"
    std::list<node_pair> toRemove;
    auto pIt = byExcludedNd.find(nd2Exclude);
    if (pIt != byExcludedNd.end()) {
        for (const auto & cp : pIt->second) {
            if (isAncestorDesNoIter(forbiddenAttach, cp)) {
                toRemove.push_back(node_pair(nd2Exclude, cp));
            }
        }
    }
    for (auto & np : toRemove) {
        purgeExcludeRaw(np.first, np.second);
    }
    ingestExcludeRaw(nd2Exclude, forbiddenAttach);
    return true;
}

template<typename T>
void ExcludeConstraints<T>::purgeExcludeRaw(const node_type * nd2Exclude,
                                            const node_type * forbiddenAttach)  {
    cnode_set & bev = byExcludedNd.at(nd2Exclude);
    assert(contains(bev, forbiddenAttach));
    bev.erase(forbiddenAttach);
    if (bev.empty()) {
        byExcludedNd.erase(nd2Exclude);
    }
    cnode_set & bcv = byNdWithConstraints.at(forbiddenAttach);
    assert(contains(bcv, nd2Exclude));
    bcv.erase(nd2Exclude);
    if (bcv.empty()) {
        byNdWithConstraints.erase(forbiddenAttach);
    }
}

template<typename T>
void ExcludeConstraints<T>::ingestExcludeRaw(const node_type * nd2Exclude,
                                             const node_type * forbiddenAttach) {
    byExcludedNd[nd2Exclude].insert(forbiddenAttach);
    byNdWithConstraints[forbiddenAttach].insert(nd2Exclude);
}

template<typename T>
void InterTreeBandBookkeeping<T>::reassignAttachmentNode(InterTreeBand<T> * b,
                                         RootedTreeNode<T> * oldAnc,
                                         RootedTreeNode<T> * newAnc,
                                         const PhyloStatement & ps) {
    assert(b != nullptr);
    assert(oldAnc != nullptr);
    assert(newAnc != nullptr);
    band2Node[b] = newAnc;
    node2Band.at(oldAnc).erase(b);
    node2Band[newAnc].insert(b);
    b->reassignAttachmentNode(oldAnc, newAnc);
}

template<typename T, typename U>
void FTree<T, U>::updateToReflectResolution(node_type *oldAnc,
                                            node_type * newAnc,
                                            const std::set<node_type *> & movedTips,
                                                            const PhyloStatement & ps) {
    auto relevantBands = bands.getBandsForNode(oldAnc);
    for (auto & b : relevantBands) {
        if (b->isTheSetOfPhantomNodes(oldAnc, movedTips)) {
            bands.reassignAttachmentNode(b, oldAnc, newAnc, ps);
        } else {
            auto nitbp = forest._createNewBand(*this, *newAnc, ps);
            assert(nitbp != nullptr);
            bands._addRefToBand(nitbp, newAnc);
            nitbp->insertSet(newAnc, movedTips);
        }
    }
}

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
    auto nr = forest.createNode(nullptr, this);
    forest.addAndUpdateChild(nr, root, *this);
    root = nr;
}

template<typename T, typename U>
const OttIdSet FTree<T, U>::getConnectedOttIds() const {
    OttIdSet r;
    // root can be a tip but not a named node, in the process of stealing
    //  children from one tree in the merging of the forests
    if (root == nullptr || (root->isTip() && !root->hasOttId())) {
        return r;
    }
    for (auto t : iter_leaf_n_const(*root)) {
        assert(isAncestorDesNoIter(root, t));
        if (t->hasOttId()) {
            r.insert(t->getOttId());
        }
    }
    return r;
}



template<typename T, typename U>
bool FTree<T,U>::anyExcludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    for (auto oid : ottIdSet) {
        if (exclude.isExcludedFrom(ottIdToNodeMap.at(oid), nd, nullptr)) {
            return true;
        }
    }
    return false;
}

template<typename T, typename U>
bool FTree<T,U>::anyIncludedAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    if (!areDisjoint(nd->getData().desIds, ottIdSet)) {
        return true;
    }
    auto c = anyPhantomNodesAtNode(nd, ottIdSet);
    assert(c == false);// if we are correctly updating desIds we don't need this branch.... TMP
    return c;
}

template<typename T, typename U>
bool FTree<T,U>::anyPhantomNodesAtNode(const node_type * nd, const OttIdSet &ottIdSet) const {
    const auto & b = bands.getBandsForNode(nd);
    if (!b.empty()) {
       const auto ns = ottIdSetToNodeSet(ottIdSet);
        for (auto ob : b) {
            if (ob->bandPointHasAny(nd, ns)) {
                return true;
            }
        }
    }
    return false;
}

template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::addLeafNoDesUpdate(RootedTreeNode<T> * par, long ottId) {
    //connectedIds.insert(ottId);
    return forest.createLeaf(par, ottId, this);
}


template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::resolveToCreateCladeOfIncluded(RootedTreeNode<T> * par,
                                                               const PhyloStatement & ps) {
    const OttIdSet & oids = ps.includeGroup;
    dbWriteOttSet("  resolveToCreateCladeOfIncluded oids = ", oids);
    dbWriteOttSet("                                 nd->getData().desIds = ", par->getData().desIds);
    std::set<RootedTreeNode<T> *> cToMove;
    std::list<RootedTreeNode<T> *> orderedToMove;
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
    }
    auto newNode = forest.createNode(par, this); // parent of includeGroup
    for (auto c : orderedToMove) {
        c->_detachThisNode();
        c->_setNextSib(nullptr);
        forest.addAndUpdateChild(newNode, c, *this);
    }
    updateToReflectResolution(par, newNode, cToMove, ps);
    assert(!par->isOutDegreeOneNode());
    forest.debugInvariantsCheck();
    return newNode;
}

template<typename T, typename U>
bool FTree<T,U>::insertIntoBandNoDesUpdate(InterTreeBand<T> * itbp,
                                    RootedTreeNode<T> * connectedNode,
                                    long phantomID) {
    assert(connectedNode != nullptr);
    assert(itbp != nullptr);
    auto noid = ottIdToNodeMap.at(phantomID);
    itbp->insert(connectedNode, noid);
    return true;
}

template<typename T>
T * rootToTipSearchByDesIds(T * nd, const OttIdSet &oids) {
    assert(nd);
    T * curr = nd;
    for (;;) {
        const auto & di = curr->getData().desIds;
        assert(isSubset(oids, di));
        if (curr->isTip()) {
            return curr;
        }
        T * nextNd = nullptr;
        for (auto c : iter_child(*curr)) {
            const auto & chdi = c->getData().desIds;
            if (!areDisjoint(oids, chdi)) {
                if (nextNd == nullptr) {
                    nextNd = &(*c);
                } else {
                    return curr;
                }
            }
        }
        if (nextNd == nullptr) {
            return curr;
        }
        curr = nextNd;
    }
}
template<typename T, typename U>
RootedTreeNode<T> * FTree<T,U>::getMRCA(const OttIdSet &ottIdSet) {
    if (ottIdSet.empty()) {
        assert(false);
        throw OTCError("empty MRCA");
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
        assert(forest.getTreeForNode(aTip) == this);
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
    const auto referredTo = set_intersection_as_set(ottIdSet, root->getData().desIds);
    assert(!referredTo.empty());
    return rootToTipSearchByDesIds(root, referredTo);
}

template<typename T, typename U>
void FTree<T, U>::mirrorPhyloStatement(const PhyloStatement &ps) {
    assert(root == nullptr);
    root = forest.createNode(nullptr, this);
    addPhyloStatementAsChildOfRoot(ps);
}

template<typename T, typename U>
void FTree<T, U>::stealExclusionStatements(node_type * newPar,
                                           node_type * srcNode,
                                           FTree<T, U>  & donorTree) {
    auto e = donorTree.exclude.stealExclusions(srcNode);
    for (auto nd : e) {
        this->exclude.addExcludeStatement(nd, newPar);
    }
}

template<typename T, typename U>
void FTree<T, U>::stealInclusionStatements(node_type * newPar,
                                           node_type * srcNode,
                                           FTree<T, U>  & donorTree) {
    auto bandSet = donorTree.bands.stealBands(srcNode);
    for (auto bandPtr : bandSet) {
        bandPtr->reassignAttachmentNode(srcNode, newPar);
        this->bands._addRefToBand(bandPtr, newPar);
    }
}

// moves the exclusion statements to a different FTree (but with the same nodes)
template<typename T, typename U>
void FTree<T, U>::registerExclusionStatementForTransferringNode(node_type * srcNode,
                                                                FTree<T, U>  & donorTree) {
    auto e = donorTree.exclude.stealExclusions(srcNode);
    for (auto nd : e) {
        this->exclude.addExcludeStatement(nd, srcNode);
    }
}
// moves the inclusion statements to a different FTree (but with the same nodes)
template<typename T, typename U>
void FTree<T, U>::registerInclusionStatementForTransferringNode(node_type * srcNode,
                                                                FTree<T, U>  & donorTree) {
    auto bandSet = donorTree.bands.stealBands(srcNode);
    for (auto bandPtr : bandSet) {
        this->bands._addRefToBand(bandPtr, srcNode);
    }
}

template<typename T, typename U>
void FTree<T, U>::addPhyloStatementAsChildOfRoot(const PhyloStatement &ps) {
    dbWriteOttSet(" addPhyloStatementAsChildOfRoot", ps.includeGroup);
    assert(root != nullptr);
    if (!root->isTip()) {
        forest.debugInvariantsCheck();
    }
    if (anyExcludedAtNode(root, ps.includeGroup)) {
        createDeeperRoot();
        assert(!root->isTip());
    }
    assert(root != nullptr);
    auto parOfIncGroup = forest.createNode(root, this); // parent of includeGroup
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
    forest.debugInvariantsCheck();
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
OttIdSet FTree<T,U>::addPhyloStatementAtNode(const PhyloStatement & ps, 
                                             RootedTreeNode<T> * includeGroupA,
                                             const OttIdSet & attachedElsewhere,
                                             InterTreeBand<T> * itbp) {
    dbWriteOttSet(" FTree<T,U>::addPhyloStatementAtNode inc", ps.includeGroup);
    LOG(DEBUG) << "includeGroupA = " << (long) includeGroupA  << " " << std::hex << (long) includeGroupA << std::dec;
    LOG(DEBUG) << "itbp = " << (long) itbp << " " << std::hex << (long) itbp << std::dec;
    LOG(DEBUG) << "FTree = " << (long) this << " " << std::hex << (long) this << std::dec;
    dbWriteOttSet("    includeGroupA->getData().desIds", includeGroupA->getData().desIds);
    assert(forest.getTreeForNode(includeGroupA) == this);
    OttIdSet r;
    if (itbp != nullptr) {
        bands._addRefToBand(itbp, includeGroupA);
    }
    for (auto oid : ps.includeGroup) {
        if (!ottIdIsConnected(oid)) {
            LOG(DEBUG) << " not connected " << oid;
            if (contains(attachedElsewhere, oid)) {
                assert(itbp != nullptr);
                LOG(DEBUG) << " attachedElsewhere " << oid;
                insertIntoBandNoDesUpdate(itbp, includeGroupA, oid);
            } else {
                LOG(DEBUG) << " adding leaf " << oid;
                addLeafNoDesUpdate(includeGroupA, oid);
                r.insert(oid);
            }
        } else {
            LOG(DEBUG) << " connected " << oid;
            assert(contains(includeGroupA->getData().desIds, oid));
        }
    }
    addDesIdsToNdAndAnc(includeGroupA, ps.includeGroup);
    dbWriteOttSet("    later includeGroupA->getData().desIds", includeGroupA->getData().desIds);
    for (auto oid : ps.excludeGroup) {
        if (!ottIdIsConnected(oid)) {
            addExcludeStatement(oid, includeGroupA, ps.provenance);
        }
    }
    LOG(DEBUG) << "before addPhyloStatementAtNode exit";
    debugInvariantsCheckFT();
    LOG(DEBUG) << "foreset level check";
    forest.debugInvariantsCheck();
    LOG(DEBUG) << "bout to addPhyloStatementAtNode exit";
    return r;
}

template<typename T, typename U>
void FTree<T, U>::debugVerifyDesIdsAssumingDes(const OttIdSet &s, const RootedTreeNode<T> *nd) const{
    OttIdSet ois;
    if (nd->isTip()) {
        if (nd->hasOttId()) {
            ois.insert(nd->getOttId());
        }
    } else {
        for (auto c : iter_child_const(*nd)) {
            const auto & coids = c->getData().desIds;
            ois.insert(begin(coids), end(coids));
        }
    }
    const auto pids = bands.getPhantomIds(nd);
    ois.insert(pids.begin(), pids.end());
    if(s != ois) {
        LOG(DEBUG) << "FTree = " << (long) this << " " << std::hex << (long) this << std::dec;
        LOG(DEBUG) << "nd = " << (long) nd << " " << std::hex << (long) nd << std::dec;
        dbWriteNewick(nd);
        dbWriteOttSet("debugVerifyDesIdsAssumingDes incoming", s);
        dbWriteOttSet("calculated:", ois);
        dbWriteOttSet("inc - calc:", set_difference_as_set(s, ois));
        dbWriteOttSet("calc - inc:", set_difference_as_set(ois, s));
        assert(s == ois);
    }
}
template<typename T, typename U>
void FTree<T, U>::debugInvariantsCheckFT() const {
    //LOG(DEBUG) << " start of debugInvariantsCheckFT for " << (long) this << " " << std::hex << (long)this << std::dec;
    //dbWriteNewick(root);
    checkAllNodePointersIter(*root);
    for (auto n : iter_post_n_const(*root)) {
        if(forest.getTreeForNode(n) != this) {
            long x = (long) forest.getTreeForNode(n);
            long p = (long) n->getParent();
            LOG(DEBUG) << getDesignator(*n) << " in the wrong tree reporting " << x << " " << std::hex << x << std::dec << " p = " << p << " " << std::hex << p << std::dec;
            if (n->getParent() != nullptr) {
                std::cerr << "the parent newick\n";
                dbWriteNewick(n->getParent());
                std::cerr << std::endl;
            }
            assert(false);
        }
        OttIdSet noids;
        if (n->isTip()) {
            if (n->hasOttId()) {
                const auto o = n->getOttId();
                assert(ottIdToNodeMap.at(o) == n);
            }
            // Make sure that our ancestors do not exclude us.
            const std::set<const node_type *> ancSet = getAncSet(n);
            for (auto a : ancSet) {
                assert(!exclude.isExcludedFrom(n, a, &ottIdToNodeMap));
            }
        } else {
            assert(!n->hasOttId());
        }
        if (n != root) {
            assert(isAncestorDesNoIter(root, n));
        }
        debugVerifyDesIdsAssumingDes(n->getData().desIds, n);
    }
}

template class ExcludeConstraints<RTSplits>; // force explicit instantiaion of this template.
template class InterTreeBand<RTSplits>; // force explicit instantiaion of this template.
template class FTree<RTSplits, MappedWithSplitsData>; // force explicit instantiaion of this template.

}// namespace
