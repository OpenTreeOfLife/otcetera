#include "otc/otcli.h"
using namespace otc;

struct PruneSynthToSubproblem : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> synthTree;
    std::set<const NodeWithSplits *> includedNodes;
    OttIdSet subproblemTipIds;
    int numErrors;
    virtual ~PruneSynthToSubproblem(){}
    PruneSynthToSubproblem()
        :synthTree(nullptr),
         numErrors(0) {
    }

    bool summarize(OTCLI & otCLI) override {
        assert(synthTree != nullptr && !includedNodes.empty());
        std::set<NodeWithSplits *> toPrune;
        for (auto nd : iter_node(*synthTree)) {
            const NodeWithSplits *  c = const_cast<const NodeWithSplits *>(nd);
            if ((!contains(includedNodes, c)) && contains(includedNodes, c->getParent())) {
                toPrune.insert(nd);
            }
        }
        for (auto nd : toPrune) {
            pruneAndDelete(*synthTree, nd);
        }
        writeTreeAsNewick(otCLI.out, *synthTree);
        otCLI.out << '\n';

        return numErrors == 0;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (synthTree == nullptr) {
            synthTree = std::move(tree);
            otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
            return true;
        }
        return processSubproblemTree(otCLI, *tree);
    }

    bool processSubproblemTree(OTCLI & otCLI, const TreeMappedWithSplits & tree) {
        std::cout << "processSubproblemTree from " << otCLI.currentFilename << "\n";
        for (const auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->getOttId();
            if (nd->hasOttId()) {
                subproblemTipIds.insert(ottId);
            }
            auto synthNode = synthTree->getData().getNodeForOttId(ottId);
            if (synthNode != nullptr) {
                if (!contains(includedNodes, synthNode)) {
                    includedNodes.insert(synthNode);
                    insertAncestorsToParaphyleticSet(synthNode, includedNodes);
                }
                insertDescendantsOfUnincludedSubtrees(synthNode, includedNodes);
            } else {
                auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
                assert(taxoNode != nullptr);
                assert(!taxoNode->isTip());
                otCLI.err << "Warning ott" << ottId << " was is an internal node that was a tip in the subproblem, but is not found in the tree being pruned.\n";
            }
        }
        return true;
    }

};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-synth-to-subproblem",
                "takes at 3 newick file paths: a full taxonomy tree, a full supertree, and a subproblem tree file. Prints a version of the synth tree restricted",
                "taxonomy.tre synth.tre step_7_scratch/export-sub-temp/ott712383.tre");
    PruneSynthToSubproblem proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
