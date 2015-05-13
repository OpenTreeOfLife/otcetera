#include "otc/otcli.h"
using namespace otc;

struct PruneTaxonomyState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    bool reportStats;
    int numErrors;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    std::set<const RootedTreeNodeNoData *> directlyIncludedNodes;
    virtual ~PruneTaxonomyState(){}

    PruneTaxonomyState()
        :reportStats(false),
        numErrors(0) {
    }

    bool summarize(OTCLI &otCLI) override {
        if (reportStats) {
            std::size_t numNonTerminals = 0;
            for (auto tn : directlyIncludedNodes) {
                if (!tn->isTip()) {
                    numNonTerminals++;
                }
            }
            otCLI.out << numNonTerminals << " non-terminal taxa in OTT that are mapped by at least 1 input.\n";
            otCLI.out << (directlyIncludedNodes.size() - numNonTerminals) << " terminal taxa in OTT that are mapped by at least 1 input\n";
            otCLI.out << directlyIncludedNodes.size() << " total taxa in OTT that are mapped by at least 1 input.\n";
            return true;
        }
        assert(taxonomy != nullptr && !includedNodes.empty());
        std::set<RootedTreeNodeNoData *> toPrune;
        std::size_t numLeavesPruned = 0;
        std::size_t numInternalsPruned = 0;
        for (auto nd : iter_node(*taxonomy)) {
            const RootedTreeNodeNoData *  c = const_cast<const RootedTreeNodeNoData *>(nd);
            if (!contains(includedNodes, c)) {
                if (contains(includedNodes, c->getParent())) {
                    toPrune.insert(nd);
                }
                if (c->isTip()) {
                    numLeavesPruned += 1;
                } else {
                    numInternalsPruned += 1;
                }
            }
        }
        for (auto nd : toPrune) {
            pruneAndDelete(*taxonomy, nd);
        }
        writeTreeAsNewick(otCLI.out, *taxonomy);
        otCLI.out << '\n';
        otCLI.err << numLeavesPruned << " terminal taxa pruned\n";
        otCLI.err << numInternalsPruned << " non-terminal taxa pruned\n";
        return true;
    }

    bool processSourceTree(OTCLI & , const std::unique_ptr<TreeMappedEmptyNodes> treePtr) override {
        assert(taxonomy != nullptr);
        std::map<const RootedTreeNodeNoData *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(*treePtr)) {
            auto ottId = nd->getOttId();
            auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
            assert(taxoNode != nullptr);
            if (reportStats) {
                directlyIncludedNodes.insert(taxoNode);
            } else {
                if (!contains(includedNodes, taxoNode)) {
                    includedNodes.insert(taxoNode);
                    insertAncestorsToParaphyleticSet(taxoNode, includedNodes);
                }
                insertDescendantsOfUnincludedSubtrees(taxoNode, includedNodes);
            }
        }
        return true;
    }
};

bool handleReport(OTCLI & otCLI, const std::string &);

bool handleReport(OTCLI & otCLI, const std::string &) {
    PruneTaxonomyState * proc = static_cast<PruneTaxonomyState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->reportStats  = true;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-prune-taxonomy",
                "takes at least 2 newick file paths: a full taxonomy tree some number of input trees. Prune subtrees from the taxonomy if they are not represented in the inputs",
                "taxonomy.tre inp1.tre inp2.tre");
    PruneTaxonomyState proc;
    otCLI.addFlag('r',
                  "Just report stats on how many tips are included in the inputs.",
                  handleReport,
                  false);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, false);
}
