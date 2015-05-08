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
    std::vector<std::string> supportTreeNames;

    virtual ~FindUnsupportedState(){}
    FindUnsupportedState()
        :toCheck(nullptr),
        numErrors(0),
        treatExpandedTipsAsSupporting(true),
        considerNamedSupported(true) {
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

    class SupportSummary {
        public:
        std::map<unsigned char, std::size_t> supportCounts;
        std::size_t numUnsupportedForking; // out-degree > 1
        std::size_t numUnsupportedKnuckles; // out-degree = 1 and # leaves under = 1
        std::size_t numUnsupportedElbows; // out-degree = 1 and # leaves under > 1
        std::size_t numSupportedInternals;
        SupportSummary()
            :numUnsupportedForking(0U),
            numUnsupportedKnuckles(0U),
            numUnsupportedElbows(0U), 
            numSupportedInternals(0U) {
            for (unsigned char i = 0U; i < END_SUPPORT_TYPE_FLAG; ++i) {
                supportCounts[i] = 0; // fill in, so we can use .at on const SupportSummary
            }
        }
    };

    SupportSummary describeUnnamedUnsupported(std::ostream &out, const TreeMappedWithSplits & tree) const {
        SupportSummary r;
        auto ig = iter_pre_internal_const(tree);
        auto nIt = ig.begin();
        const auto eIt = ig.end();
        ++nIt; //skip the root
        for (; nIt != eIt; ++nIt) {
            auto nd = *nIt;
            const auto snIt = supportedNodes.find(nd);
            if (snIt != supportedNodes.end()) {
                r.supportCounts[snIt->second] += 1;
                r.numSupportedInternals += 1;
                continue;
            }
            auto outDegree = nd->getOutDegree();
            if (outDegree == 1) {
                if (nd->includesOnlyOneLeaf()) {
                    r.numUnsupportedKnuckles += 1;
                } else {
                    r.numUnsupportedElbows += 1;
                }
            } else {
                if (aPrioriProblemNodes.empty()) {
                    out << "Unsupported node ";
                } else {
                    auto gaIt = aPrioriProblemNodes.find(nd);
                    if (gaIt == aPrioriProblemNodes.end()) {
                        out << "Novel unsupported node ";
                    } else {
                        out << "Confirmation of unsupported node (designators =";
                        writeOttSet(out, "", gaIt->second, " ");
                        out << ") ";
                    }
                }
                describeUnnamedNode(*nd, out, 0, false);
                r.numUnsupportedForking += 1;
            }
        }
        return r;
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
        out << "Final summary:\n";
        out << ss.numSupportedInternals << " internal nodes where flagged as being supported in some sense.\n";
        out << "    This means each of these nodes was:\n    supported by at least";
        out << " one non-taxonomic input tree ";
        if (considerNamedSupported) {
            out << " or the internal node has a name";
        }
        out << ".\n";
        out << "The breakdown of the types of supported nodes:\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_INTERNAL) << " unnamed internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_EXPANDED) << " unnamed internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(SEEN_IN_AN_INPUT_BOTH) << " unnamed  internal nodes supported by an input tree internal node and another tree\'s expanded  tip node.\n";
        out << "    " << ss.supportCounts.at(NAMED_NODE) << " named internal nodes had no support (other than the name).\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_INTERNAL) << " named internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_EXPANDED) << " named internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(NAMED_SEEN_IN_AN_INPUT_BOTH) << " named internal nodes supported by an input tree internal node and another tree\'s expanded  tip node.\n";
        out << "The following counts refer to redundant internal nodes (out-degree = 1)";
        out << ". These are considered supported if their child is supported.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_INTERNAL) << " unnamed redundant internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_EXPANDED) << " unnamed redundant  internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + SEEN_IN_AN_INPUT_BOTH) << " unnamed  redundant  internal nodes supported by an input tree internal node and another tree\'s expanded  tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_NODE) << "  redundant named internal nodes had no support (other than the name).\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_INTERNAL) << "  redundant named internal nodes supported by an input tree internal node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_EXPANDED) << "  redundant named internal nodes supported by an input tree expanded tip node.\n";
        out << "    " << ss.supportCounts.at(REDUNDANT_ND + NAMED_SEEN_IN_AN_INPUT_BOTH) << "  redundant named internal nodes supported by an input tree internal node and another tree\'s expanded  tip node.\n";
        out << ss.numUnsupportedForking << " unsupported nodes had out-degree > 1\n";
        out << ss.numUnsupportedKnuckles << " unsupported nodes had out-degree == 1 and were along a terminal path.\n";
        out << ss.numUnsupportedElbows << " unsupported nodes had out-degree == 1 and were along an internal path.\n";
        out << ss.numUnsupportedElbows + ss.numSupportedInternals  + ss.numUnsupportedKnuckles + ss.numUnsupportedForking ;
        out << " total internal nodes checked.\n";
        out << std::endl;
        if (numErrors < 0) {
            numErrors -= ss.numUnsupportedForking;
        } else {
            numErrors = static_cast<int>(ss.numUnsupportedForking);
        }
        return numErrors == 0;
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
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        if (toCheck == nullptr) {
            toCheck = std::move(tree);
            if (considerNamedSupported) {
                for (auto nd : iter_pre_internal_const(*toCheck)) {
                    if (nd->hasOttId()) {
                        supportedNodes[nd] = NAMED_NODE;
                    }
                }
            }
            return true;
        }
        supportTreeNames.push_back(tree->getName());
        const auto expanded = expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
        auto r = processExpandedTree(otCLI, *tree, expanded);
        return r;
    }

    bool processExpandedTree(OTCLI & otCLI,
                             const TreeMappedWithSplits & tree,
                             const std::set<const NodeWithSplits *> & expandedTips) {
        assert(toCheck != nullptr);
        std::map<const NodeWithSplits *, std::set<long> > restrictedDesIds;
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->getOttId();
            markPathToRoot(*toCheck, ottId, restrictedDesIds);
        }
        std::map<std::set<long>, const NodeWithSplits *> sourceClades;
        for (auto nd : iter_post_internal_const(tree)) {
            if (nd->getParent() != nullptr && !nd->isTip()) {
                if (treatExpandedTipsAsSupporting || contains(expandedTips, nd)) {
                    sourceClades[nd->getData().desIds] = nd;
                }
            }
        }
        identifySupportedNodes(otCLI, tree, restrictedDesIds, sourceClades, expandedTips);
        return true;
    }

    void identifySupportedNodes(OTCLI & otCLI,
                                const TreeMappedWithSplits & tree,
                                const std::map<const NodeWithSplits *, std::set<long> > & inducedNdToEffDesId,
                                const std::map<std::set<long>, const NodeWithSplits *> & sourceClades,
                                const std::set<const NodeWithSplits *> & expandedTips) {
        for (auto pd : inducedNdToEffDesId) {
            auto nd = pd.first;
            auto par = nd->getParent();
            if (par == nullptr) {
                continue;
            }
            auto firstBranchingAnc = findFirstBranchingAnc<const NodeWithSplits>(nd);
            if (firstBranchingAnc == nullptr) {
                //otCLI.out << "  firstBranchingAnc null\n";
                continue;
            }
            auto nm = pd.second;
            auto ancIt = inducedNdToEffDesId.find(firstBranchingAnc);
            assert(ancIt != inducedNdToEffDesId.end());
            auto anm = ancIt->second;
            const NodeWithSplits * firstNdPtr; // just used to match call
            if (!multipleChildrenInMap(*nd, inducedNdToEffDesId, &firstNdPtr)) {
                continue;
            }
            if (anm == nm) {
                continue;
            }
            auto scIt = sourceClades.find(nm);
            if (scIt != sourceClades.end()) {
                if (aPrioriProblemNodes.find(nd) != aPrioriProblemNodes.end()) {
                    auto apIt = aPrioriProblemNodes.find(nd);
                    otCLI.out << "ERROR!: a priori unsupported node found. Designators were ";
                    writeOttSet(otCLI.out, "", apIt->second, " ");
                    otCLI.out << ". A node was found, which (when pruned to the leaf set of an input tree) contained:\n";
                    writeOttSet(otCLI.out, "    ", nm, " ");
                    otCLI.out << "\nThe subtree from the source was: ";
                    auto srcNd = scIt->second;
                    writePrunedSubtreeNewickForMarkedNodes(otCLI.out, *srcNd, inducedNdToEffDesId);
                    numErrors += 1;
                }
                recordInputTreeSupportForNode(nd, scIt->second, tree, expandedTips);
            }
        }
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
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
