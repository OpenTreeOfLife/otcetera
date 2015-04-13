#include "otc/otcli.h"
using namespace otc;

struct NonTerminalsToExemplarsState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    int numErrors;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    std::map<std::unique_ptr<TreeMappedEmptyNodes>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedEmptyNodes *> treePtrByIndex;
    std::string exportDir;
    virtual ~NonTerminalsToExemplarsState(){}
    NonTerminalsToExemplarsState()
        :numErrors(0) {
    }

    bool summarize(OTCLI &otCLI) override {
        assert(taxonomy != nullptr && !includedNodes.empty());
        std::set<RootedTreeNodeNoData *> toPrune;
        for (auto nd : iter_node(*taxonomy)) {
            const RootedTreeNodeNoData *  c = const_cast<const RootedTreeNodeNoData *>(nd);
            if ((!contains(includedNodes, c)) && contains(includedNodes, c->getParent())) {
                toPrune.insert(nd);
            }
        }
        for (auto nd : toPrune) {
            pruneAndDelete(*taxonomy, nd);
        }
        writeTreeAsNewick(otCLI.out, *taxonomy);
        otCLI.out << '\n';
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedEmptyNodes> treeup) override {
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedEmptyNodes * raw = treeup.get();
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        treePtrByIndex.push_back(raw);
        // Store the tree's filename
        raw->setName(otCLI.currentFilename);
        std::map<const RootedTreeNodeNoData *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(*raw)) {
            auto ottId = nd->getOttId();
            auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
            assert(taxoNode != nullptr);
            if (!contains(includedNodes, taxoNode)) {
                includedNodes.insert(taxoNode);
                insertAncestorsToParaphyleticSet(taxoNode, includedNodes);
            }
            if (!taxoNode->isTip()) {
                assert(false); // expecting this to trip...
            }
        }
        return true;
    }
};

bool handleExportModified(OTCLI & otCLI, const std::string &narg);

bool handleExportModified(OTCLI & otCLI, const std::string &narg) {
    NonTerminalsToExemplarsState * proc = static_cast<NonTerminalsToExemplarsState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -b argument.");
    }
    proc->exportDir = narg;
    return true;
}


int main(int argc, char *argv[]) {
    const char * helpMsg = "takes an -e flag specifying an export diretory and at least 2 newick file paths: " \
        "a full taxonomy tree some number of input trees. Any tip in non-taxonomic input that is mapped to " \
        "non-terminal taoxn will be remapped such that the parent of the non-terminal tip will hold all of " \
        "the expanded exemplars. The exemplars will be the union of tips that (a) occur below this non-terminal " \
        "taxon in the taxonomy and (b) occur, or are used as an exemplar, in another input tree. The modified " \
        "version of each input will be written in the export directory. Trees with no non-terminal tips should " \
        "be unaltered. The taxonomy written out will be the taxonomy restricted to the set of leaves that are " \
        "leaves of the exported trees.";
    OTCLI otCLI("otc-nonterminals-to-exemplars",
                helpMsg,
                "taxonomy.tre inp1.tre inp2.tre");
    NonTerminalsToExemplarsState proc;
    otCLI.addFlag('e',
                  "ARG should be the name of a directory. A .tre file will be written to this directory for each input tree",
                  handleExportModified,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, false);
}
