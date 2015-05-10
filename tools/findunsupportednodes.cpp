#include "otc/otcli.h"
using namespace otc;

enum SupportType {
    SEEN_IN_AN_INPUT_INTERNAL = 1,
    SEEN_IN_AN_INPUT_EXPANDED = 2,
    SEEN_IN_AN_INPUT_BOTH = 3,
    NAMED_NODE = 4,
    NAMED_SEEN_IN_AN_INPUT_INTERNAL = 5,
    NAMED_SEEN_IN_AN_INPUT_EXPANDED = 6,
    NAMED_SEEN_IN_AN_INPUT_BOTH = 7,
    REDUNDANT_ND = 8,
    END_SUPPORT_TYPE_FLAG = 16
};

struct FindUnsupportedState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    class SupportSummary {
        public:
        std::map<unsigned char, std::size_t> supportCounts;
        std::size_t numUnsupportedForking; // out-degree > 1
        std::size_t numUnsupportedKnuckles; // out-degree = 1 and # leaves under = 1
        std::size_t numUnsupportedElbows; // out-degree = 1 and # leaves under > 1
        std::size_t numSupportedInternals;
        std::map<const NodeWithSplits *, const NodeWithSplits *> misnamedSupported;
        std::map<const NodeWithSplits *, const NodeWithSplits *> supportedShouldHaveName;
        std::set<const NodeWithSplits *> shouldBeUnnamed;
        bool taxoChecked;
        SupportSummary()
            :numUnsupportedForking(0U),
            numUnsupportedKnuckles(0U),
            numUnsupportedElbows(0U), 
            numSupportedInternals(0U),
            taxoChecked(false) {
            for (unsigned char i = 0U; i < END_SUPPORT_TYPE_FLAG; ++i) {
                supportCounts[i] = 0; // fill in, so we can use .at on const SupportSummary
            }
        }
    };
    
    std::unique_ptr<TreeMappedWithSplits> toCheck;
    int numErrors;
    std::map<const NodeWithSplits *, std::set<long> > aPrioriProblemNodes;
    std::map<const NodeWithSplits *, unsigned char> supportedNodes;
    std::map<const NodeWithSplits *, std::set<std::size_t> > supportedBy;
    // true is default. False means that the expansion of a tip does not create
    //  an input clade that counts as support. 
    // false is a better setting (in terms of interpreting the input trees correctly)
    // TODO: Do we want to change the default?
    bool treatExpandedTipsAsSupporting;
    bool considerNamedSupported;
    bool recordSupportingTreeIdentity;
    bool processTaxonomyAsInputTree;
    bool refreshAfterTaxonomy;
    bool inTheProcessOfAnalyzingTax;
    std::vector<std::string> supportTreeNames;
    SupportSummary taxoSummary; // only used if 
    bool fixInsteadOfReport;
        
    virtual ~FindUnsupportedState(){}
    FindUnsupportedState()
        :toCheck(nullptr),
        numErrors(0),
        treatExpandedTipsAsSupporting(true),
        considerNamedSupported(true),
        processTaxonomyAsInputTree(false),
        refreshAfterTaxonomy(false),
        inTheProcessOfAnalyzingTax(false),
        fixInsteadOfReport(false) {
    }

    // abuse of const because the fixInsteadOfReport was added late.
    mutable std::map<const NodeWithSplits * , OttId> misnameQ;
    mutable std::map<const NodeWithSplits * , OttId> toNameQ;
    mutable std::set<const NodeWithSplits * > toDelNameQ;
    mutable std::set<const NodeWithSplits * > toCollapse;
    mutable std::size_t numNamesAdded;
    mutable std::size_t numNamesChanged;
    mutable std::size_t numNamesDeleted;
    mutable std::size_t numNodesCollapsed;

    void queueFixMisnamed(const NodeWithSplits * nd, const NodeWithSplits * tax) const {
        misnameQ[nd] = tax->getOttId();
    }
    void queueFixShouldHaveBeenNamed(const NodeWithSplits * nd, const NodeWithSplits * tax) const {
        toNameQ[nd] = tax->getOttId();
    }
    void queueFixShouldBeUnnamed(const NodeWithSplits * nd) const {
        toDelNameQ.insert(nd);
    }
    void queueFixCollapse(const NodeWithSplits * nd) const {
        toCollapse.insert(nd);
    }
    // helper for performFixes
    void _changeName(const std::pair<const NodeWithSplits *, OttId> & nn, std::size_t &x) const {
            auto nd = const_cast<NodeWithSplits *>(nn.first);
            auto ottId = nn.second;
            changeOttIdOfInternal(*toCheck, nd, ottId);
            x += 1;
    }
    void performFixes() const {
        // misnamed and toName are supported by taxonomy
        for (auto nn : misnameQ) {
            _changeName(nn, numNamesChanged);
        }
        misnameQ.clear();
        for (auto nn : toNameQ) {
            _changeName(nn, numNamesAdded);
        }
        toNameQ.clear();
        // toDelName might be unsupported nodes, so we check toCollapse..
        for (auto nn : toDelNameQ) {
            auto nd = const_cast<NodeWithSplits *>(nn);
            if (!contains(toCollapse, nd)) {
                delOttIdOfInternal(*toCheck, nd);
                ++numNamesDeleted;
            }
        }
        for (auto nn : toCollapse) {
            auto nd = const_cast<NodeWithSplits *>(nn);
            collapseNode(*toCheck, nd);
        }
    }

    void extendSupportedToRedundantNodes(const TreeMappedWithSplits & tree) {
        for (auto nd : iter_post_internal_const(tree)) {
            if (nd->isOutDegreeOneNode()) {
                auto c = nd->getFirstChild();
                if (contains(supportedNodes, c)) {
                    supportedNodes[nd] = REDUNDANT_ND | supportedNodes[c];
                    if (recordSupportingTreeIdentity) {
                        const auto & cb = supportedBy[c];
                        supportedBy[nd].insert(begin(cb), end(cb));
                    }
                }
            }
        }
    }

    
    SupportSummary describeUnnamedUnsupported(std::ostream &out, const TreeMappedWithSplits & tree) const {
        SupportSummary r;
        calcSupportSummary(&out, tree, r, false);
        return r;
    }
    void calcSupportSummary(std::ostream * out,
                            const TreeMappedWithSplits & tree,
                            SupportSummary & r, 
                            bool isTaxoSummary) const {
        if (!isTaxoSummary) {
            assert(out);
        }
        r.taxoChecked = isTaxoSummary;
        auto ig = iter_pre_internal_const(tree);
        auto nIt = ig.begin();
        const auto eIt = ig.end();
        ++nIt; //skip the root
        for (; nIt != eIt; ++nIt) {
            auto nd = *nIt;
            const auto snIt = supportedNodes.find(nd);
            if (snIt != supportedNodes.end()) {
                const auto t = snIt->second;
                // not in taxoSummary, or having real (non-name-only) support
                if ((!isTaxoSummary) || (t & SEEN_IN_AN_INPUT_BOTH)) {
                    r.supportCounts[t] += 1;
                    r.numSupportedInternals += 1;
                    if (isTaxoSummary) {
                        const OttIdSet & di = nd->getData().desIds;
                        if (nd->hasOttId()) {
                            const auto ottId = nd->getOttId();
                            const auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
                            if (taxoNode->getData().desIds != di) {
                                auto matchingTaxonNode = findNodeWithMatchingDesIdSet(*taxonomy, di);
                                assert(matchingTaxonNode != nullptr);
                                r.misnamedSupported[nd] = matchingTaxonNode;
                                if (fixInsteadOfReport) {
                                    queueFixMisnamed(nd, matchingTaxonNode);
                                }
                            }
                        } else {
                            auto matchingTaxonNode = findNodeWithMatchingDesIdSet(*taxonomy, di);
                            assert(matchingTaxonNode != nullptr);
                            r.supportedShouldHaveName[nd] = matchingTaxonNode;
                            if (fixInsteadOfReport) {
                                queueFixShouldHaveBeenNamed(nd, matchingTaxonNode);
                            }
                        }
                    }
                   continue;
                }
                // if in taxo mode and no hasTaxoSupport, then fall through to unsupported
            }
            auto outDegree = nd->getOutDegree();
            if (isTaxoSummary & nd->hasOttId()) {
                r.shouldBeUnnamed.insert(nd);
                if (fixInsteadOfReport) {
                    queueFixShouldBeUnnamed(nd);
                }
            }
            if (outDegree == 1) {
                if (nd->includesOnlyOneLeaf()) {
                    r.numUnsupportedKnuckles += 1;
                } else {
                    r.numUnsupportedElbows += 1;
                }
            } else {
                if (!isTaxoSummary && !fixInsteadOfReport) {
                    if (aPrioriProblemNodes.empty()) {
                        *out << "Unsupported node ";
                    } else {
                        auto gaIt = aPrioriProblemNodes.find(nd);
                        if (gaIt == aPrioriProblemNodes.end()) {
                            *out << "Novel unsupported node ";
                        } else {
                            *out << "Confirmation of unsupported node (designators =";
                            writeOttSet(*out, "", gaIt->second, " ");
                            *out << ") ";
                        }
                    }
                    describeUnnamedNode(*nd, *out, 0, false);
                }
                r.numUnsupportedForking += 1;
                if (fixInsteadOfReport && (!isTaxoSummary)) {
                    queueFixCollapse(nd);
                }
            }
        }
        if (fixInsteadOfReport) {
            performFixes();
        }
    }

    bool summarize(OTCLI &otCLI) override {
        extendSupportedToRedundantNodes(*toCheck);
        auto & out = otCLI.out;
        const auto ss = describeUnnamedUnsupported(otCLI.out, *toCheck);
        for (auto gaIt : aPrioriProblemNodes) {
            if (supportedNodes.find(gaIt.first) != supportedNodes.end()) {
                out << "Claim of unsupported apparently refuted for designators: ";
                writeOttSet(out, "", gaIt.second, " ");
                out << ". See standard error stream for details.\n";
            }
        }
        if (processTaxonomyAsInputTree) {
            out << "Taxonomy-as-the-only-input summary:\n";
            writeSummaryStats(out, taxoSummary);
        }
        out << "Final summary:\n";
        writeSummaryStats(out, ss);
        if (numErrors < 0) {
            numErrors -= ss.numUnsupportedForking;
        } else {
            numErrors = static_cast<int>(ss.numUnsupportedForking);
        }
        return numErrors == 0;
    }

    void writeSummaryStats(std::ostream & out, const SupportSummary & ss) {
        if (ss.taxoChecked) {
            out << "Summary of checks of internal nodes names/Ids:\n";
            out << ss.misnamedSupported.size() << " named internal nodes with an incorrect name.\n";
            for (auto msp : ss.misnamedSupported) {
                out << "    Node with ID ott" << msp.first->getOttId() << " should have been named ott"<< msp.second->getOttId() << ".\n";
            }
            out << ss.supportedShouldHaveName.size() << " unamed internal nodes which should have been named.\n";
            for (auto msp : ss.supportedShouldHaveName) {
                out << "    ";
                describeUnnamedNode(*msp.first, out, 0, false, false);
                out << " should have been named ott"<< msp.second->getOttId() << ".\n";
            }
            out << ss.shouldBeUnnamed.size() << " named internal nodes which correspond to no taxon.\n";
            for (auto np : ss.shouldBeUnnamed) {
                out << "    ott" << np->getOttId() << " should not be named.\n";
            }
        }
        out << ss.numSupportedInternals << " internal nodes where flagged as being supported in some sense.\n";
        out << "    This means each of these nodes was:\n    supported by at least";
        out << " one non-taxonomic input tree";
        if (considerNamedSupported) {
            out << " or the internal node has a name";
        }
        out << ".\n";
        out << "The breakdown of the types of supported nodes:\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_INTERNAL) << " unnamed internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_EXPANDED) << " unnamed internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_BOTH) << " unnamed internal nodes supported by an input tree internal node and another tree\'s expanded tip node.\n";
        out << "    " << ss.supportCounts.at(NAMED_NODE) << " named internal nodes had no support (other than the name).\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_INTERNAL) << " named internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_EXPANDED) << " named internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_BOTH) << " named internal nodes supported by an input tree internal node and another tree\'s expanded tip node.\n";
        out << "The following counts refer to redundant internal nodes (out-degree = 1)";
        out << ". These are considered supported if their child is supported.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_INTERNAL) << " unnamed redundant internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_EXPANDED) << " unnamed redundant internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_BOTH) << " unnamed  redundant internal nodes supported by an input tree internal node and another tree\'s expanded  tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_NODE) << " redundant named internal nodes had no support (other than the name).\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_INTERNAL) << " redundant named internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_EXPANDED) << " redundant named internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_BOTH) << " redundant named internal nodes supported by an input tree internal node and another tree\'s expanded tip node.\n";
        out << ss.numUnsupportedForking << " unsupported nodes had out-degree > 1\n";
        out << ss.numUnsupportedKnuckles << " unsupported nodes had out-degree == 1 and were along a terminal path.\n";
        out << ss.numUnsupportedElbows << " unsupported nodes had out-degree == 1 and were along an internal path.\n";
        out << ss.numUnsupportedElbows + ss.numSupportedInternals  + ss.numUnsupportedKnuckles + ss.numUnsupportedForking ;
        out << " total internal nodes checked.\n";
        out << std::endl;
    }

    void parseAndProcessMRCADesignatorsFile(const std::string &fp) {
        if (toCheck == nullptr) {
            throw OTCError("Designator files (if used) must be passed in after the tree to check");
        }
        std::list<std::set<long> > dl = parseDesignatorsFile(fp);
        for (auto d : dl) {
            markSuspectNode(d);
        }
    }

    void markSuspectNode(const std::set<long> & designators) {
        const NodeWithSplits * mrca = findMRCAFromIDSet(*toCheck, designators, -1);
        aPrioriProblemNodes[mrca] = designators;
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = true;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set (rather than the full set of IDs)
        clearAndfillDesIdSets(*taxonomy);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        return true;
    }

    void reinitSupportTemps() {
        supportedNodes.clear();
        if (considerNamedSupported) {
            for (auto nd : iter_pre_internal_const(*toCheck)) {
                if (nd->hasOttId()) {
                    supportedNodes[nd] = (nd->isOutDegreeOneNode() ? (REDUNDANT_ND + NAMED_NODE) : NAMED_NODE);
                }
            }
        }
    }
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        if (toCheck == nullptr) {
            toCheck = std::move(tree);
            reinitSupportTemps();
            if (processTaxonomyAsInputTree) {
                inTheProcessOfAnalyzingTax = true;
                analyzeTreeForSupport(otCLI, *taxonomy, true);
                inTheProcessOfAnalyzingTax = false;
                extendSupportedToRedundantNodes(*toCheck);
                calcSupportSummary(nullptr, *toCheck, taxoSummary, true);
                if (refreshAfterTaxonomy) {
                   reinitSupportTemps();
                }
            }
            return true;
        }
        return analyzeTreeForSupport(otCLI, *tree, true);
    }

    bool analyzeTreeForSupport(OTCLI & otCLI, TreeMappedWithSplits & tree, bool needExpansion) {
        supportTreeNames.push_back(tree.getName());
        std::set<const NodeWithSplits *> expanded;
        if (needExpansion) {
            expanded = expandOTTInternalsWhichAreLeaves(tree, *taxonomy);
        }
        return processExpandedTree(otCLI, tree, expanded);
    }

    bool processExpandedTree(OTCLI & otCLI,
                             const TreeMappedWithSplits & tree,
                             const std::set<const NodeWithSplits *> & expandedTips) {
        assert(toCheck != nullptr);
        std::map<const NodeWithSplits *, std::set<long> > restrictedDesIds;
        if (inTheProcessOfAnalyzingTax) {
            identifySupportedNodesTaxo(tree);
        } else {
            for (auto nd : iter_leaf_const(tree)) {
                auto ottId = nd->getOttId();
                markPathToRoot(*toCheck, ottId, restrictedDesIds);
            }
            identifySupportedNodes(otCLI, tree, restrictedDesIds, expandedTips);
        }
        return true;
    }

    void identifySupportedNodes(OTCLI & otCLI,
                                const TreeMappedWithSplits & tree,
                                const std::map<const NodeWithSplits *, std::set<long> > & inducedNdToEffDesId,
                                const std::set<const NodeWithSplits *> & expandedTips) {
        for (auto pd : inducedNdToEffDesId) {
            checkNodeForSupport(otCLI, pd.first, pd.second, tree, inducedNdToEffDesId, expandedTips);
        }
    }
    void identifySupportedNodesTaxo(const TreeMappedWithSplits & tree) {
        const std::set<const NodeWithSplits *> expandedTips;
        for (auto nd : iter_pre_internal_const(*toCheck)) {
            checkNodeForSupportTaxo(nd, tree, expandedTips);
        }
    }

    void checkNodeForSupportTaxo(const NodeWithSplits *nd,
                                 const TreeMappedWithSplits & tree,
                                 const std::set<const NodeWithSplits *> & expandedTips) {
        auto par = nd->getParent();
        if (par == nullptr) {
            return;
        }
        const OttIdSet * nmp = nullptr;
        if (nd->isOutDegreeOneNode()) {
            if (!nd->hasOttId()) {
                return;
            }
            nmp = &(nd->getData().desIds); //
        } else {
            auto firstBranchingAnc = findFirstBranchingAnc<const NodeWithSplits>(nd);
            if (firstBranchingAnc == nullptr) {
                return;
            }
            nmp = &(nd->getData().desIds); //
            const auto & anm = firstBranchingAnc->getData().desIds; //
            if (anm == *nmp) {
                return;
            }
        }
        auto srcNode = findNodeWithMatchingDesIdSet(tree, *nmp);
        if (srcNode != nullptr) {
            recordInputTreeSupportForNode(nd, srcNode, tree, expandedTips);
        }
    }

    void checkNodeForSupport(OTCLI & otCLI,
                             const NodeWithSplits *nd,
                             const OttIdSet & nm,
                             const TreeMappedWithSplits & tree,
                             const std::map<const NodeWithSplits *, OttIdSet > & inducedNdToEffDesId,
                             const std::set<const NodeWithSplits *> & expandedTips) {
        auto par = nd->getParent();
        if (par == nullptr) {
            return;
        }
        auto firstBranchingAnc = findFirstBranchingAnc<const NodeWithSplits>(nd);
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
            if (aPrioriProblemNodes.find(nd) != aPrioriProblemNodes.end()) {
                auto apIt = aPrioriProblemNodes.find(nd);
                otCLI.out << "ERROR!: a priori unsupported node found. Designators were ";
                writeOttSet(otCLI.out, "", apIt->second, " ");
                otCLI.out << ". A node was found, which (when pruned to the leaf set of an input tree) contained:\n";
                writeOttSet(otCLI.out, "    ", nm, " ");
                otCLI.out << "\nThe subtree from the source was: ";
                writePrunedSubtreeNewickForMarkedNodes(otCLI.out, *srcNode, inducedNdToEffDesId);
                numErrors += 1;
            }
            recordInputTreeSupportForNode(nd, srcNode, tree, expandedTips);
        }
    }

    bool treeHasClade(const TreeMappedWithSplits & tree, const OttIdSet & oids) {
        return nullptr != findNodeWithMatchingDesIdSet(tree, oids);
    }

    void recordInputTreeSupportForNode(const NodeWithSplits * treeToCheckNode,
                                       const NodeWithSplits * srcTreeNode,
                                       const TreeMappedWithSplits & srcTree, 
                                       const std::set<const NodeWithSplits *> & expandedTips) {
        if (contains(expandedTips, srcTreeNode)) {
            supportedNodes[treeToCheckNode] |= SEEN_IN_AN_INPUT_EXPANDED;
        } else {
            supportedNodes[treeToCheckNode] |= SEEN_IN_AN_INPUT_INTERNAL;
        }
        if (recordSupportingTreeIdentity) {
            supportedBy[treeToCheckNode].insert(getInputIndex(srcTree));
        }
    }
    std::size_t getInputIndex(const TreeMappedWithSplits & ) {
        return supportTreeNames.size() - 1;
    }
};

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg);
bool handleStoreSources(OTCLI & otCLI, const std::string &nextArg);
bool handleTipsNotSupport(OTCLI & otCLI, const std::string &nextArg);
bool handleForceTaxonomy(OTCLI & otCLI, const std::string &);
bool handleForceRefreshAfterTaxonomy(OTCLI & otCLI, const std::string &);
bool handleFix(OTCLI & otCLI, const std::string &);

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    assert(!nextArg.empty());
    fusp->parseAndProcessMRCADesignatorsFile(nextArg);
    return true;
}

bool handleStoreSources(OTCLI & otCLI, const std::string &) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->recordSupportingTreeIdentity = true;
    return true;
}

bool handleTipsNotSupport(OTCLI & otCLI, const std::string &) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->treatExpandedTipsAsSupporting = false;
    return true;
}

bool handleForceTaxonomy(OTCLI & otCLI, const std::string &) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->processTaxonomyAsInputTree = true;
    return true;
}

bool handleForceRefreshAfterTaxonomy(OTCLI & otCLI, const std::string &) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->refreshAfterTaxonomy = true;
    return true;
}

bool handleFix(OTCLI & otCLI, const std::string &) {
    FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->fixInsteadOfReport = true;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-find-unsupported-nodes",
                "takes at least 2 newick file paths: a full taxonomy tree, a full supertree, and some number of input trees",
                "taxonomy.tre synth.tre inp1.tre inp2.tre");
    FindUnsupportedState proc;
    otCLI.addFlag('m',
                  "ARG=a designators file. Each line is a list of (white-space separated) OTT ids used to designate the node that is the MRCA of them.",
                  handleDesignator,
                  true);
    otCLI.addFlag('t',
                  "tips mapped to non-terminal taxa should NOT be counted as support",
                  handleTipsNotSupport,
                  false);
    otCLI.addFlag('x',
                  "Automatically treat the taxonomy as an input",
                  handleForceTaxonomy,
                  false);
    otCLI.addFlag('r',
                  "Refresh support stats after analyzing taxonomy so that final summary is only based on phylo inputs",
                  handleForceRefreshAfterTaxonomy,
                  false);
    otCLI.addFlag('c',
                  "Fix the problems (clean the tree) rather than reporting on them.",
                  handleFix,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
