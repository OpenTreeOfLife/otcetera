#include "otc/otcli.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;

struct RFState : public TaxonomyDependentTreeProcessor<Tree_t> {
    virtual ~RFState(){}
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        otCLI.out << numInducedSplitsMissingInSecond(*taxonomy, *tree, true) << '\n';
        return true;
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcmissingsplits",
                "takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the number of splits in the induced full tree that are missing from the each other tree.",
                "taxonomy.tre inp1.tre inp2.tre");
    RFState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
