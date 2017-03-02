#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include <tuple>
using namespace otc;
// See http://phylo.bio.ku.edu/ot/otc-find-resolution.pdf
// which is compiled from ../doc/otc-find-resolution.tex

template <typename T, typename U>
U * resolveNode(T & tree, U & parent, const OttIdSet & newInc) {
    std::set<U *> childrenToMove;
    for (auto nd : iter_child(parent)) {
        if (!areDisjoint(nd->get_data().desIds, newInc)) {
            childrenToMove.insert(nd);
        }
    }
    assert(!childrenToMove.empty());
    U * phPar = &parent;
   assert(phPar != nullptr);
    U * insertedNodePtr = tree.create_node(phPar);
    for (auto phChild : childrenToMove) {
        phChild->detachThisNode();
        insertedNodePtr->addChild(phChild);
        const OttIdSet & cd = phChild->get_data().desIds;
        insertedNodePtr->get_data().desIds.insert(cd.begin(), cd.end());
    }
    assert(phPar->getOutDegree() > 1);
    return insertedNodePtr;
}

struct FindResolutionState;


class SupportingIDSets {
    public:
    typedef OttIdSet leafSetContainer;
    typedef std::tuple<const OttIdSet *, const leafSetContainer *, const char *> incLSTreeNameTuple;
    typedef std::list<incLSTreeNameTuple> listIncLSTreeNameTuple;
    listIncLSTreeNameTuple supporting;
    typedef std::list<listIncLSTreeNameTuple::iterator> listLIncLSTreeNameIt;
    void addSupport(const OttIdSet *i, const leafSetContainer *pi, const char * treeName) {
        supporting.push_back(incLSTreeNameTuple{i, pi, treeName});
    }
    bool attemptSplitOfSupport(OTCLI & otCLI,
                               const NodeWithSplits *nd,
                               FindResolutionState & frs,
                               const TreeMappedWithSplits & inpTree,
                               const OttIdSet & incAdded);
    bool causesChildToBeUnsupported(OTCLI & otCLI,
                                    const NodeWithSplits *added,
                                    FindResolutionState & frs);
};

struct FindResolutionState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> summaryTreeToResolve;
    std::size_t numIncludable;
    bool addGroups;
    int numErrors;
    std::map<const NodeWithSplits *, std::set<long> > protectedPolytomy;
    std::string protectedFilename;
    bool avoidAddingUnsupportedGroups;
    std::list<std::unique_ptr<TreeMappedWithSplits> > allInps;
    std::map<const NodeWithSplits *, SupportingIDSets> supportStatementsByNd; // SSS in the docs
    std::list<OttIdSet> ownedIds;
    bool treatTaxonomyAsLastTree;

    OttIdSet * aliasToNewOttIdSet() {
        ownedIds.push_back(OttIdSet{});
        return &(*(ownedIds.rbegin()));
    }

    virtual ~FindResolutionState(){}
    FindResolutionState()
        :summaryTreeToResolve(nullptr),
        numIncludable(0U),
        addGroups(false),
        numErrors(0),
        avoidAddingUnsupportedGroups(false),
        treatTaxonomyAsLastTree(false) {
        BUGEXPORTINDEX = 0;
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
            writeTreeAsNewick(otCLI.out, *summaryTreeToResolve);
            otCLI.out << '\n';
            return true;
        }
        otCLI.out << numIncludable << " clades found which could be added to the tree.\n";
        numErrors = static_cast<int>(numIncludable);
        return numErrors == 0;
    }

    virtual bool process_taxonomy_tree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::process_taxonomy_tree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set
        clearAndfillDesIdSets(*taxonomy);
        return true;
    }
    
    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (summaryTreeToResolve == nullptr) {
            requireNonRedundantTree(*tree);
            summaryTreeToResolve = std::move(tree);
            if (!protectedFilename.empty()) {
                parseAndProcessMRCADesignatorsFile(protectedFilename);
            }
            return true;
        }
        requireTipsToBeMappedToTerminalTaxa(*tree, *taxonomy);
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
        if (treatTaxonomyAsLastTree) {
            recordSupportingStatements(otCLI, *taxonomy);
        }
        for (auto & treePtr : allInps) {
            attemptResolutionFromSourceTree(otCLI, *treePtr);
        }
        if (treatTaxonomyAsLastTree) {
            attemptResolutionFromSourceTree(otCLI, *taxonomy);
        }
    }

    void recordSupportingStatements(OTCLI & otCLI, TreeMappedWithSplits & tree) {
        std::map<const NodeWithSplits *, std::set<long> > restrictedDesIds;
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->get_ott_id();
            markPathToRoot(*summaryTreeToResolve, ottId, restrictedDesIds);
        }
        identifysupportStatementsByNd(otCLI, tree, restrictedDesIds);
    }
    void identifysupportStatementsByNd(OTCLI & otCLI,
                                const TreeMappedWithSplits & tree,
                                const std::map<const NodeWithSplits *, std::set<long> > & inducedNdToEffDesId) {
        for (auto pd : inducedNdToEffDesId) {
            const auto & nm = pd.second;
            if (nm.size() > 1) {
                checkNodeForSupport(otCLI,
                                    pd.first,
                                    nm, tree,
                                    inducedNdToEffDesId);
            }
        }
    }
    // This is very similar to the IsSupporteBy check and
    //      storage in RecordSupportStatement in the docs
    // inducedNdToEffDesId is map of node in S' traversed when
    //      tracing the input tree to the intersection between
    //      the input leaf set and desIds of the node
    // `nd` is a node in the summary tree being evaluated
    // `nm` is the induced set of OttIds for this nd (the intersection
    //      between nd.desId and tree.leaf_set)
    // `tree` is the input tree
    bool checkNodeForSupport(OTCLI & ,
                             const NodeWithSplits * nd,
                             const OttIdSet & nm,
                             const TreeMappedWithSplits & tree,
                             const std::map<const NodeWithSplits *, OttIdSet > & inducedNdToEffDesId) {
        auto par = nd->getParent();
        if (par == nullptr) {
            return false;
        }
        //assert(nm.size() > 1); checked by caller

        // If no node in the input tree displays this induced desId set
        //  then the tree does not support `nd`
        auto srcNode = findNodeWithMatchingDesIdSet(tree, nm);
        if (srcNode == nullptr) {
            return false;
        }

        //
        // If only one child was traversed in creating inducedNdToEffDesId
        // then this is not a MRCA. (collapsing the edge to nd would still
        //  display srcNode via the edge to the single child).
        const NodeWithSplits * firstNdPtr; // just used to match call
        if (!multipleChildrenInMap(*nd, inducedNdToEffDesId, &firstNdPtr)) {
            return false;
        }

        // If the nearest ancestor with > 1 child has the same
        //  intersection with the leaf set of the input, then the node
        //  is not supported. (collapsing the edge to nd would still
        //  display srcNode).
        auto firstBranchingAnc = findFirstForkingAnc<const NodeWithSplits>(nd);
        assert (firstBranchingAnc == par);
        auto ancIt = inducedNdToEffDesId.find(firstBranchingAnc);
        assert(ancIt != inducedNdToEffDesId.end());
        const auto & anm = ancIt->second;
        if (anm == nm) {
            return false;
        }

        //srcNode supports nd
        const OttIdSet * sip = &(srcNode->get_data().desIds);
        auto sna = findFirstForkingAnc<const NodeWithSplits>(srcNode);
        assert(sna != nullptr);
        const auto * tls = getStableLeafSetPtr(tree);
        supportStatementsByNd[nd].addSupport(sip, tls, tree.get_name().c_str());
        return true;
    }

    std::map<const TreeMappedWithSplits *, OttIdSet> tree2LeafSet;
    const OttIdSet * getStableLeafSetPtr(const TreeMappedWithSplits & tree ) {
        const auto  tp = &tree;
        const auto tpI = tree2LeafSet.find(tp);
        if (tpI == tree2LeafSet.end()) {
            tree2LeafSet[tp] = keys(tree.get_data().ottIdToNode);
            return &(tree2LeafSet[tp]);
        }
        return &(tpI->second);
    }


    void attemptResolutionFromSourceTree(OTCLI & otCLI, TreeMappedWithSplits & tree) {
        const OttIdSet & treeLeafSet = tree.get_root()->get_data().desIds;
        std::set<const NodeWithSplits *> nodesFromInp;
        for (auto nd : iter_pre_internal_const(tree)) {
            if (nd->getParent() != nullptr) {
                nodesFromInp.insert(nd);
            }
        }
        while (!nodesFromInp.empty()) {
            std::map<NodeWithSplits *, std::list<const NodeWithSplits *> > summaryTreeToResolveNodeToResolves;
            for (const auto nd : nodesFromInp) {
                const OttIdSet & incGroup = nd->get_data().desIds;
                auto r = findMRCAUsingDesIds(*summaryTreeToResolve, incGroup);
                if (contains(protectedPolytomy, r)) {
                    continue; // skip protected
                }
                NodeWithSplits * mrca = const_cast<NodeWithSplits *>(r);
                if (!mrca->isPolytomy()) {
                    continue;
                }
                assert(mrca);
                OttIdSet md = mrca->get_data().desIds;
                const OttIdSet rmd = set_intersection_as_set(md, treeLeafSet);
                // If the polytomy does not have any members of nd's 
                //  exclude group, then nd cannot resolve it and support
                //  the new edge.
                if (incGroup == rmd) {
                    continue;
                }
                const OttIdSet excGroup = set_difference_as_set(treeLeafSet, incGroup);
                if (canBeResolvedToDisplayIncExcGroup(mrca, incGroup, excGroup)) {
                    if (addGroups) {
                        summaryTreeToResolveNodeToResolves[mrca].push_back(nd);
                    } else {
                        numIncludable += 1;
                        otCLI.out << otCLI.currentFilename << " node " << getDesignator(*nd) << " could be added.\n";
                    }
                }
            }
            // rather than deal with the logic for adding to the correct node in the resolved tree, we'll just iterate
            //  over the whole loop again. NOT EFFICIENT. MTH is just being lazy...
            std::set<const NodeWithSplits *> toDealWith;
            for (auto mrcaListIncPair : summaryTreeToResolveNodeToResolves) {
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
                    processAddableNode(otCLI, *summaryTreeToResolve, mrcaPtr, tree, toAdd);
                }
            }
            nodesFromInp = toDealWith;
        }
    }
    unsigned BUGEXPORTINDEX;
    void processAddableNode(OTCLI & otCLI,
                           TreeMappedWithSplits & toResolve,
                           NodeWithSplits * ndToResolve,
                           const TreeMappedWithSplits & inpTree,
                           const NodeWithSplits * ndToAdd) {
        const OttIdSet & incGroup = ndToAdd->get_data().desIds;
        if (otCLI.verbose) {
            otCLI.err << "processAddableNode ndToResolve = ";
            describeUnnamedNode(*ndToResolve, otCLI.err, 0, false, true);
        }
        OttId BUGGY = 913397;
        const bool DEBUGGINGTHISNODE = contains(ndToResolve->get_data().desIds, BUGGY);
        if (DEBUGGINGTHISNODE) {
            ++BUGEXPORTINDEX;
        }
        if (DEBUGGINGTHISNODE) {
            std::string fn = std::string("trace-newick-pre-") + std::to_string(BUGEXPORTINDEX) + std::string(".tre");
            std::ofstream ofn(fn.c_str());
            writeNewick(ofn, ndToResolve);
            ofn << ";\n";
        }
        NodeWithSplits * added = resolveNode(toResolve, *ndToResolve, incGroup);
        if (DEBUGGINGTHISNODE) {
            std::string fn = std::string("trace-newick-resolved-") + std::to_string(BUGEXPORTINDEX) + std::string(".tre");
            std::ofstream ofn(fn.c_str());
            writeNewick(ofn, ndToResolve);
            ofn << ";\n";
        }
        if (otCLI.verbose) {
            otCLI.err << "  from = " << inpTree.get_name() << "\n";
            otCLI.err << "                   added = ";
            describeUnnamedNode(*added, otCLI.err, 0, false, true);
        }
        if (avoidAddingUnsupportedGroups
            && causedUnsupported(otCLI,
                                toResolve,
                                ndToResolve,
                                inpTree,
                                incGroup,
                                added)) {
            if (otCLI.verbose) {
                otCLI.err << "                 Rejected.\n";
            }
            collapseNode(toResolve, added);
            if (DEBUGGINGTHISNODE) {
                std::string fn = std::string("trace-newick-collapsed-") + std::to_string(BUGEXPORTINDEX) + std::string(".tre");
                std::ofstream ofn(fn.c_str());
                writeNewick(ofn, ndToResolve);
                ofn << ";\n";
            }
        } else {
            if (otCLI.verbose) {
                otCLI.err << "                 Accepted.\n";
            }
            if (avoidAddingUnsupportedGroups) {
                supportStatementsByNd[added].addSupport(&incGroup,
                                                        getStableLeafSetPtr(inpTree),
                                                        inpTree.get_name().c_str());
            }
            numIncludable += 1;
        }
    }

    bool causedUnsupported(OTCLI & otCLI,
                           TreeMappedWithSplits & ,
                           NodeWithSplits * par,
                           const TreeMappedWithSplits & inpTree,
                           const OttIdSet & incAdded,
                           NodeWithSplits * added) {
        assert(par != nullptr);
        assert(added != nullptr);
        assert(par == added->getParent());
        if (!contains(supportStatementsByNd, par)) {
            return false;
        }
        auto & supids = supportStatementsByNd.at(par);
        return ! supids.attemptSplitOfSupport(otCLI,
                                              added,
                                              *this,
                                              inpTree,
                                              incAdded);
    }


    void parseAndProcessMRCADesignatorsFile(const std::string &fp) {
        assert(summaryTreeToResolve != nullptr);
        std::list<std::set<long> > dl = parseDesignatorsFile(fp);
        for (auto d : dl) {
            markProtectedNode(d);
        }
    }

    void markProtectedNode(const std::set<long> & designators) {
        const NodeWithSplits * mrca = findMRCAFromIDSet(*summaryTreeToResolve, designators, -1);
        protectedPolytomy[mrca] = designators;
    }
};

inline bool SupportingIDSets::causesChildToBeUnsupported(OTCLI & otCLI, 
                                                         const NodeWithSplits *added,
                                                         FindResolutionState & frs) {
    std::map<const NodeWithSplits *, SupportingIDSets *> forkingDes;
    LOG(DEBUG) << added->getOutDegree() << " children of added. ";
    for (auto c : iter_child_const(*added)) {
        if (c->isTip()) {
            continue;
        }
        auto fc = findFirstForkingSelfOrDes(c);
        if (fc->isTip()) {
            continue;
        }
        auto snIt = frs.supportStatementsByNd.find(fc);
        if (snIt == frs.supportStatementsByNd.end() && !c->has_ott_id()) {
            assert(false);
        } else {
            forkingDes[fc] = &(snIt->second);
        }
    }
    LOG(DEBUG) << forkingDes.size() << " supporting statements of children to check.";
    const OttIdSet & incAdded = added->get_data().desIds;
    std::map<const NodeWithSplits *, listLIncLSTreeNameIt > toDelMap;
    for (auto fdIt : forkingDes) {
        const NodeWithSplits * fc = fdIt.first;
        SupportingIDSets & suppIds = *(fdIt.second);
        auto & toDel = toDelMap[fc];
        auto xIt = begin(suppIds.supporting);
        for (; xIt != end(suppIds.supporting); ++xIt) {
            const incLSTreeNameTuple & ipip = *xIt;
            const auto inInter = set_intersection_as_set(*std::get<0>(ipip), incAdded);
            const auto parInInter = set_intersection_as_set(*std::get<1>(ipip), incAdded);
            if (inInter == parInInter) {
                toDel.push_back(xIt);
            } else if (otCLI.verbose) {
                otCLI.err << "child still supported by split from tree " << std::get<2>(ipip) << '\n';
            }
        }
        if ((toDel.size() == suppIds.supporting.size()) && !fc->has_ott_id()) {
            return true;
        }
    }
    for (auto tdmIt : toDelMap) {
        SupportingIDSets & suppIds = *(forkingDes[tdmIt.first]);
        for (auto toDelIt : tdmIt.second) {
            suppIds.supporting.erase(toDelIt);
        }
    }
    return false;
}
    //returns true if added should be retained
inline bool SupportingIDSets::attemptSplitOfSupport(OTCLI & otCLI, 
                                                    const NodeWithSplits *added,
                                                    FindResolutionState & frs,
                                                    const TreeMappedWithSplits & ,
                                                    const OttIdSet & ) {
    //every support statement for this node, should either:
    //  1. stay in SupportingIDSets, or
    //  2. be removed
    assert(added != nullptr);
    auto p = added->getParent(); // p is the node that is the key for `this` SupportingIDSets
    assert(p != nullptr);
    const OttIdSet & nddi = added->get_data().desIds;
    listLIncLSTreeNameIt toDel;
    auto sIt = begin(supporting);
    for (; sIt != end(supporting); ++sIt) {
        const OttIdSet & suppInc = *std::get<0>(*sIt);
        //const leafSetContainer & leafSetInp = *(sIt->second);
        //if (otCLI.verbose) {
        //    otCLI.err << "suppInc ";
        //    writeOttSet(otCLI.err, "  ", suppInc, " ");
        //    otCLI.err << "\nleafSetInp";
        //    writeOttSet(otCLI.err, "  ", leafSetInp, " ");
        //    otCLI.err << '\n';
        //}
        if (isSubset(suppInc, nddi)) {
            if (otCLI.verbose) {
                otCLI.err << "toDel\n";
            }
            toDel.push_back(sIt); // now displayed on parent(p) -> p -> nd
        } else {
            if (otCLI.verbose) {
                otCLI.err << "toRetain\n";
            }
        }
    }
     if (toDel.size() ==  supporting.size() && !p->has_ott_id()) {
        if (otCLI.verbose) {
            otCLI.err << "All " << toDel.size() << " supporting statements moving or deleted \n";
        }
        return false;
    }
    if (otCLI.verbose) {
        otCLI.err << toDel.size() << " to be deleted. ";
        otCLI.err << supporting.size() - toDel.size() << " to remain.\n";
    }
    if (causesChildToBeUnsupported(otCLI, added, frs)) {
        return false;
    }
    for (auto i : toDel) {
        supporting.erase(i);
    }
    return true;
}

bool handleResolve(OTCLI & otCLI, const std::string &);
bool handleDesignator(OTCLI & otCLI, const std::string &);
bool handleCountTaxonomy(OTCLI & otCLI, const std::string &);
bool handleAvoidUnsupported(OTCLI & otCLI, const std::string &);

bool handleAvoidUnsupported(OTCLI & otCLI, const std::string &) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->avoidAddingUnsupportedGroups = true;
    return true;
}
bool handleCountTaxonomy(OTCLI & otCLI, const std::string &) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->treatTaxonomyAsLastTree = true;
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
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input in terms of supporting groups",
                  handleCountTaxonomy,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
