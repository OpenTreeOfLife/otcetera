#include "otc/otcli.h"
#include "otc/supertree_util.h"
using namespace otc;

template <typename T, typename U>
U * resolveNode(T & tree, U & parent, const OttIdSet & newInc) {
    std::set<U *> childrenToMove;
    for (auto nd : iter_child(parent)) {
        if (!areDisjoint(nd->getData().desIds, newInc)) {
            childrenToMove.insert(nd);
        }
    }
    assert(!childrenToMove.empty());
    U * phPar = &parent;
   assert(phPar != nullptr);
    U * insertedNodePtr = tree.createNode(phPar);
    for (auto phChild : childrenToMove) {
        phChild->_detachThisNode();
        insertedNodePtr->addChild(phChild);
        const OttIdSet & cd = phChild->getData().desIds;
        insertedNodePtr->getData().desIds.insert(cd.begin(), cd.end());
    }
    assert(phPar->getOutDegree() > 1);
    return insertedNodePtr;
}

struct FindResolutionState;

class SupportingIDSets {
    public:
    typedef std::pair<const OttIdSet *, const OttIdSet *> incParIncPair;
    std::list<incParIncPair> supporting;
    void addSupport(const OttIdSet *i, const OttIdSet *pi) {
        supporting.push_back(std::pair<const OttIdSet *, const OttIdSet *>(i, pi));
    }
    bool attemptSplitOfSupport(OTCLI & otCLI,
                               const NodeWithSplits *nd,
                               SupportingIDSets & toAddTo,
                               FindResolutionState & frs,
                               const OttIdSet & incAdded,
                               const OttIdSet & parIncAdded);
    bool causesChildToBeUnsupported(OTCLI & otCLI, const NodeWithSplits *nd, FindResolutionState & frs);
};

struct FindResolutionState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> toCheck;
    std::size_t numIncludable;
    bool addGroups;
    int numErrors;
    std::map<const NodeWithSplits *, std::set<long> > protectedPolytomy;
    std::string protectedFilename;
    bool avoidAddingUnsupportedGroups;
    std::list<std::unique_ptr<TreeMappedWithSplits> > allInps;
    std::map<const NodeWithSplits *, SupportingIDSets> supportedNodes;
    std::list<OttIdSet> ownedIds;
    OttIdSet * aliasToNewOttIdSet() {
        ownedIds.push_back(OttIdSet{});
        return &(*(ownedIds.rbegin()));
    }
    virtual ~FindResolutionState(){}
    FindResolutionState()
        :toCheck(nullptr),
        numIncludable(0U),
        addGroups(false),
        numErrors(0),
        avoidAddingUnsupportedGroups(false) {
    }

    bool summarize(OTCLI &otCLI) override {
        if (avoidAddingUnsupportedGroups) {
            doRefinementAvoidingUnsupportedGroups(otCLI);
        }
        if (addGroups) {
            otCLI.err << numIncludable << " nodes added to the tree.\n";
            if (numIncludable == 0) {
                return false;
            }
            writeTreeAsNewick(otCLI.out, *toCheck);
            otCLI.out << '\n';
            return true;
        }
        otCLI.out << numIncludable << " clades found which could be added to the tree.\n";
        numErrors = static_cast<int>(numIncludable);
        return numErrors == 0;
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set
        clearAndfillDesIdSets(*taxonomy);
        return true;
    }
    
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (toCheck == nullptr) {
            toCheck = std::move(tree);
            if (!protectedFilename.empty()) {
                parseAndProcessMRCADesignatorsFile(protectedFilename);
            }
            return true;
        }
        expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
        clearAndfillDesIdSets(*tree);
        if (avoidAddingUnsupportedGroups) {
            allInps.emplace_back(std::move(tree));
        } else {
            attemptResolutionFromSourceTree(otCLI, *tree);
        }
        return true;
    }
    void doRefinementAvoidingUnsupportedGroups(OTCLI & otCLI) {
        for (auto & treePtr : allInps) {
            recordSupportingStatements(otCLI, *treePtr);
        }
        for (auto & treePtr : allInps) {
            attemptResolutionFromSourceTree(otCLI, *treePtr);
        }
    }
    void recordSupportingStatements(OTCLI & otCLI, TreeMappedWithSplits & tree) {
        std::map<const NodeWithSplits *, std::set<long> > restrictedDesIds;
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->getOttId();
            markPathToRoot(*toCheck, ottId, restrictedDesIds);
        }
        identifySupportedNodes(otCLI, tree, restrictedDesIds);
    }
    void identifySupportedNodes(OTCLI & otCLI,
                                const TreeMappedWithSplits & tree,
                                const std::map<const NodeWithSplits *, std::set<long> > & inducedNdToEffDesId) {
        for (auto pd : inducedNdToEffDesId) {
            checkNodeForSupport(otCLI, pd.first, pd.second, tree, inducedNdToEffDesId);
        }
    }
    void checkNodeForSupport(OTCLI & ,
                             const NodeWithSplits * nd,
                             const OttIdSet & nm,
                             const TreeMappedWithSplits & tree,
                             const std::map<const NodeWithSplits *, OttIdSet > & inducedNdToEffDesId) {
        auto par = nd->getParent();
        if (par == nullptr) {
            return;
        }
        auto firstBranchingAnc = findFirstForkingAnc<const NodeWithSplits>(nd);
        if (firstBranchingAnc == nullptr) {
            return;
        }
        auto ancIt = inducedNdToEffDesId.find(firstBranchingAnc);
        assert(ancIt != inducedNdToEffDesId.end());
        const auto & anm = ancIt->second;
        const NodeWithSplits * firstNdPtr; // just used to match call
        if (!multipleChildrenInMap(*nd, inducedNdToEffDesId, &firstNdPtr)) {
            return;
        }
        if (anm == nm) {
            return;
        }
        auto srcNode = findNodeWithMatchingDesIdSet(tree, nm);
        if (srcNode != nullptr) {
            const OttIdSet * sip = &(srcNode->getData().desIds);
            auto sna = findFirstForkingAnc<const NodeWithSplits>(srcNode);
            assert(sna != nullptr);
            const OttIdSet * spip = &(sna->getData().desIds);
            supportedNodes[nd].addSupport(sip, spip);
        }
    }

    void attemptResolutionFromSourceTree(OTCLI & otCLI, TreeMappedWithSplits & tree) {
        const OttIdSet & treeLeafSet = tree.getRoot()->getData().desIds;
        std::set<const NodeWithSplits *> nodesFromInp;
        for (auto nd : iter_pre_internal_const(tree)) {
            if (nd->getParent() != nullptr) {
                nodesFromInp.insert(nd);
            }
        }
        while (!nodesFromInp.empty()) {
            std::map<NodeWithSplits *, std::list<const NodeWithSplits *> > toCheckNodeToResolves;
            for (const auto nd : nodesFromInp) {
                const OttIdSet & incGroup = nd->getData().desIds;
                auto r = findMRCAUsingDesIds(*toCheck, incGroup);
                if (contains(protectedPolytomy, r)) {
                    continue; // skip protected
                }
                NodeWithSplits * mrca = const_cast<NodeWithSplits *>(r);
                assert(mrca);
                OttIdSet md = mrca->getData().desIds;
                OttIdSet rmd;
                for (auto i : md) {
                    if (contains(treeLeafSet, i)) {
                        rmd.insert(i);
                    }
                }
                if (incGroup == rmd) {
                    continue;
                }
                OttIdSet excGroup = treeLeafSet;
                for (auto i :incGroup) {
                    excGroup.erase(i);
                }
                if (canBeResolvedToDisplayExcGroup(mrca, incGroup, excGroup)) {
                    if (addGroups) {
                        toCheckNodeToResolves[mrca].push_back(nd);
                    } else {
                        numIncludable += 1;
                        otCLI.out << otCLI.currentFilename << " node " << getDesignator(*nd) << " could be added.\n";
                    }
                }
            }
            // rather than deal with the logic for adding to the correct node in the resolved tree, we'll just iterate
            //  over the whole loop again. NOT EFFICIENT. MTH is just being lazy...
            std::set<const NodeWithSplits *> toDealWith;
            for (auto mrcaListIncPair : toCheckNodeToResolves) {
                auto mrcaPtr = mrcaListIncPair.first;
                const NodeWithSplits * toAdd = nullptr;
                for (auto tac : mrcaListIncPair.second) {
                    if (toAdd == nullptr) {
                        toAdd = tac;
                    } else {
                        toDealWith.insert(tac);
                    }
                }
                if (toAdd != nullptr) {
                    const OttIdSet & incGroup = toAdd->getData().desIds;
                    const OttIdSet & parIncAdded = toAdd->getParent()->getData().desIds;
                    NodeWithSplits * added = resolveNode(*toCheck, *mrcaPtr, incGroup);
                    if (avoidAddingUnsupportedGroups
                        && causedUnsupported(otCLI,
                                            *toCheck,
                                            mrcaPtr,
                                            incGroup,
                                            parIncAdded,
                                            added)) {
                        collapseNode(*toCheck, added);
                    } else {
                        if (avoidAddingUnsupportedGroups) {
                            supportedNodes[added].addSupport(&incGroup, &parIncAdded);
                        }
                        numIncludable += 1;
                        if (otCLI.verbose) {
                            otCLI.err << "tree \"" << tree.getName() << "\" adding split\n  inc. ";
                            writeOttSet(otCLI.err, "  ", incGroup, " ");
                            otCLI.err << "\n  exc. ";
                            OttIdSet p = toAdd->getParent()->getData().desIds;
                            p = set_difference_as_set(p, incGroup);
                            writeOttSet(otCLI.err, "  ", p, " ");
                            otCLI.err << '\n';
                            
                        }
                    }
                }
            }
            nodesFromInp = toDealWith;
        }
    }

    bool causedUnsupported(OTCLI & otCLI,
                           TreeMappedWithSplits & ,
                           NodeWithSplits * par,
                           const OttIdSet & incAdded,
                           const OttIdSet & parIncAdded,
                           NodeWithSplits * added) {
        assert(par != nullptr);
        assert(added != nullptr);
        assert(par == added->getParent());
        if (!contains(supportedNodes, par)) {
            return false;
        }
        auto & supids = supportedNodes.at(par);
        SupportingIDSets addedSupIds;
        if (!supids.attemptSplitOfSupport(otCLI,
                                          added,
                                          addedSupIds,
                                          *this,
                                          incAdded,
                                          parIncAdded)) {
            return true;
        }
        supportedNodes[added] = addedSupIds;
        return false;
    }


    void parseAndProcessMRCADesignatorsFile(const std::string &fp) {
        assert(toCheck != nullptr);
        std::list<std::set<long> > dl = parseDesignatorsFile(fp);
        for (auto d : dl) {
            markProtectedNode(d);
        }
    }

    void markProtectedNode(const std::set<long> & designators) {
        const NodeWithSplits * mrca = findMRCAFromIDSet(*toCheck, designators, -1);
        protectedPolytomy[mrca] = designators;
    }
};

inline bool SupportingIDSets::causesChildToBeUnsupported(OTCLI & , 
                                                         const NodeWithSplits *nd,
                                                         FindResolutionState & frs) {
    std::map<const NodeWithSplits *, SupportingIDSets *> forkingDes;
    for (auto c : iter_child_const(*nd)) {
        if (c->isTip()) {
            continue;
        }
        auto fc = findFirstForkingSelfOrDes(c);
        if (fc->isTip()) {
            continue;
        }
        auto snIt = frs.supportedNodes.find(fc);
        if (snIt == frs.supportedNodes.end()) {
            if (!fc->hasOttId()) {
                return false;
            }
        } else {
            forkingDes[fc] = &(snIt->second);
        }
    }
    const OttIdSet & incAdded = nd->getData().desIds;
    std::map<const NodeWithSplits *, std::list<std::list<incParIncPair>::iterator> > toDelMap;
    for (auto fdIt : forkingDes) {
        const NodeWithSplits * fc = fdIt.first;
        SupportingIDSets & suppIds = *(fdIt.second);
        std::list<std::list<incParIncPair>::iterator> & toDel = toDelMap[fc];
        std::list<incParIncPair>::iterator xIt = begin(suppIds.supporting);
        for (; xIt != end(suppIds.supporting); ++xIt) {
            incParIncPair & ipip = *xIt;
            const auto inInter = set_intersection_as_set(*ipip.first, incAdded);
            const auto parInInter = set_intersection_as_set(*ipip.second, incAdded);
            if (inInter == parInInter) {
                toDel.push_back(xIt);
            }
        }
        if (toDel.size() == suppIds.supporting.size() && (!fc->hasOttId())) {
            return false;
        }
    }
    for (auto tdmIt : toDelMap) {
        SupportingIDSets & suppIds = *(forkingDes[tdmIt.first]);
        for (auto toDelIt : tdmIt.second) {
            suppIds.supporting.erase(toDelIt);
        }
    }
    return true;
}
inline bool SupportingIDSets::attemptSplitOfSupport(OTCLI & otCLI, 
                                                    const NodeWithSplits *nd,
                                                    SupportingIDSets & toAddTo,
                                                    FindResolutionState & frs,
                                                    const OttIdSet & incAdded,
                                                    const OttIdSet & parIncAdded) {
    //every support statement for this node, should either:
    //  1. stay in SupportingIDSets,
    //  2. move to toAddTo, or 
    //  3. be removed
    //  to toAddTo
    assert(nd != nullptr);
    auto p = nd->getParent();
    assert(p != nullptr);
    const OttIdSet & nddi = nd->getData().desIds;
    const OttIdSet & pdi = p->getData().desIds;
    std::list<std::list<incParIncPair>::iterator> toMove;
    std::list<std::list<incParIncPair>::iterator> toDel;
    std::list<incParIncPair>::iterator sIt = begin(supporting);
    for (; sIt != end(supporting); ++sIt) {
        const OttIdSet & suppInc = *(sIt->first);
        const OttIdSet & suppParInc = *(sIt->second);
        if (otCLI.verbose) {
            otCLI.err << "suppInc ";
            writeOttSet(otCLI.err, "  ", suppInc, " ");
            otCLI.err << "\nsuppParInc";
            writeOttSet(otCLI.err, "  ", suppParInc, " ");
            otCLI.err << '\n';
        }

        if (isSubset(suppInc, nddi)) {
            if (haveIntersection(suppParInc, pdi)) {
                if (otCLI.verbose) {
                    otCLI.err << "toMove\n";
                }
                toMove.push_back(sIt);
            } else {
                if (otCLI.verbose) {
                    otCLI.err << "toDel\n";
                }
                toDel.push_back(sIt);
            }
        } else {
            if (otCLI.verbose) {
                otCLI.err << "toRetain\n";
            }
        }
    }
    const auto tcs = toMove.size() + toDel.size();
    if (tcs ==  supporting.size()) {
        if (otCLI.verbose) {
            otCLI.err << "All " << tcs << " supporing statements moving or deleted \n";
        }
        return false;
    }
    if (otCLI.verbose) {
        otCLI.err << toMove.size() << " moving. " << toDel.size() << " to be del.";
        otCLI.err << supporting.size() - tcs << " to remain.";
    }
    if (causesChildToBeUnsupported(otCLI, nd, frs)) {
        return false;
    }
    for (auto i : toMove) {
        const OttIdSet & suppInc = *(i->first);
        const OttIdSet & suppParInc = *(i->second);
        OttIdSet * nsi = frs.aliasToNewOttIdSet();
        OttIdSet * npsi = frs.aliasToNewOttIdSet();
        *nsi = set_intersection_as_set(suppInc, incAdded);
        *npsi = set_intersection_as_set(suppParInc, parIncAdded);
        toAddTo.addSupport(nsi, npsi);
        supporting.erase(i);
    }

    for (auto i : toDel) {
        supporting.erase(i);
    }
    return true;
}

bool handleResolve(OTCLI & otCLI, const std::string &);
bool handleDesignator(OTCLI & otCLI, const std::string &);
bool handleAvoidUnsupported(OTCLI & otCLI, const std::string &);

bool handleAvoidUnsupported(OTCLI & otCLI, const std::string &) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->avoidAddingUnsupportedGroups = true;
    return true;
}

bool handleResolve(OTCLI & otCLI, const std::string &) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->addGroups = true;
    return true;
}

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    assert(!nextArg.empty());
    proc->protectedFilename = nextArg;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-find-resolution",
                "takes at least 3 newick file paths: a taxonomy,  a full supertree, and some number of input trees.",
                "synth.tre inp1.tre inp2.tre ...");
    FindResolutionState proc;
    otCLI.addFlag('r',
                  "Resolve the supertree rather than counting number of groups that could be added.",
                  handleResolve,
                  false);
    otCLI.addFlag('u',
                  "Do not add a resolution if it will cause an unsupported group.",
                  handleAvoidUnsupported,
                  false);
    otCLI.addFlag('p',
                  "ARG=a designators file. Each line is a list of (white-space separated) OTT ids used to designate the node that is the MRCA of them. " \
                  "If the -r option is used, these nodes will be protected from resolution.",
                  handleDesignator,
                  true);
    
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
