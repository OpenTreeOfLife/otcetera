#include "otc/otcli.h"
using namespace otc;

struct PruneSynthToSubproblem : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> synthTree;
    std::set<const NodeWithSplits *> includedNodes;
    OttIdSet subproblemTipIds;
    int numErrors;
    bool pruneInpTreesNotSynth;
    std::string outDir;
    virtual ~PruneSynthToSubproblem(){}
    PruneSynthToSubproblem()
        :synthTree(nullptr),
         numErrors(0),
         pruneInpTreesNotSynth(false) {
    }

    bool summarize(OTCLI & otCLI) override {
        if (pruneInpTreesNotSynth) {
            return true;
        }
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

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (synthTree == nullptr) {
            synthTree = std::move(tree);
            otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
            return true;
        }
        return processSubproblemTree(otCLI, *tree);
    }

    bool processSubproblemTree(OTCLI & otCLI, const TreeMappedWithSplits & tree) {
        if (pruneInpTreesNotSynth) {
            std::string path = outDir + std::string("/") + otCLI.currentFilename;
            std::ofstream outstream(path.c_str());
            writeTreeAsNewick(outstream, tree);
            return true;
        }
        for (const auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->get_ott_id();
            if (nd->has_ott_id()) {
                subproblemTipIds.insert(ottId);
            }
            auto synthNode = synthTree->get_data().getNodeForOttId(ottId);
            if (synthNode != nullptr) {
                if (!contains(includedNodes, synthNode)) {
                    includedNodes.insert(synthNode);
                    insertAncestorsToParaphyleticSet(synthNode, includedNodes);
                }
            } else {
                auto taxoNode = taxonomy->get_data().getNodeForOttId(ottId);
                assert(taxoNode != nullptr);
                assert(!taxoNode->isTip());
                otCLI.err << "Warning ott" << ottId << " was is an internal node that was a tip in the subproblem, but is not found in the tree being pruned.\n";
            }
        }
        return true;
    }

};

bool handleReverseToDir(OTCLI & otCLI, const std::string &);

bool handleReverseToDir(OTCLI & otCLI, const std::string &nextArg) {
    PruneSynthToSubproblem * proc = static_cast<PruneSynthToSubproblem *>(otCLI.blob);
    assert(proc != nullptr);
    assert(!nextArg.empty());
    proc->outDir = nextArg;
    proc->pruneInpTreesNotSynth = true;
    otCLI.getParsingRules().pruneUnrecognizedInputTips = true;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-synth-to-subproblem",
                "takes at 3 newick file paths: a full taxonomy tree, a full supertree, and a subproblem tree file. Prints a version of the synth tree restricted",
                "taxonomy.tre synth.tre step_7_scratch/export-sub-temp/ott712383.tre");
    PruneSynthToSubproblem proc;
    otCLI.addFlag('r',
                  "ARG=output directory. Reverse the operation: prune the inputs instead of the synth tree. Only works if the synth tree is passed in as the first AND second argument.",
                  handleReverseToDir,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
