#include "otc/otcli.h"
using namespace otc;

struct PruneTaxonomyState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    int numErrors;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    virtual ~PruneTaxonomyState(){}

    PruneTaxonomyState()
        :numErrors(0) {
    }

    bool summarize(const OTCLI &otCLI) override {
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

    bool processSourceTree(OTCLI & , const std::unique_ptr<TreeMappedEmptyNodes> treePtr) override {
        assert(taxonomy != nullptr);
        std::map<const RootedTreeNodeNoData *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(*treePtr)) {
            auto ottId = nd->getOttId();
            auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
            assert(taxoNode != nullptr);
            if (!contains(includedNodes, taxoNode)) {
                includedNodes.insert(taxoNode);
                insertAncestorsToParaphyleticSet(taxoNode, includedNodes);
            }
            insertDescendantsOfUnincludedSubtrees(taxoNode, includedNodes);
        }
        return true;
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcprunetaxonomy",
                "takes at least 2 newick file paths: a full taxonomy tree some number of input trees. Prune subtrees from the taxonomy if they are not represented in the inputs.",
                "taxonomy.tre inp1.tre inp2.tre");
    PruneTaxonomyState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, false);
}
