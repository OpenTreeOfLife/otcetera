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
            if ((!contains(includedNodes, c)) && contains(includedNodes, c->get_parent())) {
                toPrune.insert(nd);
            }
        }
        for (auto nd : toPrune) {
            prune_and_delete(*synthTree, nd);
        }
        write_tree_as_newick(otCLI.out, *synthTree);
        otCLI.out << '\n';

        return numErrors == 0;
    }

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (synthTree == nullptr) {
            synthTree = std::move(tree);
            otCLI.get_parsing_rules().include_internal_nodes_in_des_id_sets = false;
            return true;
        }
        return processSubproblemTree(otCLI, *tree);
    }

    bool processSubproblemTree(OTCLI & otCLI, const TreeMappedWithSplits & tree) {
        if (pruneInpTreesNotSynth) {
            std::string path = outDir + std::string("/") + otCLI.currentFilename;
            std::ofstream outstream(path.c_str());
            write_tree_as_newick(outstream, tree);
            return true;
        }
        for (const auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->get_ott_id();
            if (nd->has_ott_id()) {
                subproblemTipIds.insert(ottId);
            }
            auto synthNode = synthTree->get_data().get_node_by_ott_id(ottId);
            if (synthNode != nullptr) {
                if (!contains(includedNodes, synthNode)) {
                    includedNodes.insert(synthNode);
                    insert_ancestors_to_paraphyletic_set(synthNode, includedNodes);
                }
            } else {
                auto taxoNode = taxonomy->get_data().get_node_by_ott_id(ottId);
                assert(taxoNode != nullptr);
                assert(!taxoNode->is_tip());
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
    otCLI.get_parsing_rules().prune_unrecognized_input_tips = true;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-synth-to-subproblem",
                "takes at 3 newick file paths: a full taxonomy tree, a full supertree, and a subproblem tree file. Prints a version of the synth tree restricted",
                "taxonomy.tre synth.tre step_7_scratch/export-sub-temp/ott712383.tre");
    PruneSynthToSubproblem proc;
    otCLI.add_flag('r',
                  "ARG=output directory. Reverse the operation: prune the inputs instead of the synth tree. Only works if the synth tree is passed in as the first AND second argument.",
                  handleReverseToDir,
                  true);
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 3, true);
}
