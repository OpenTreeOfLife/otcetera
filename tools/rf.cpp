#include "otc/otcli.h"
using namespace otc;
struct RFState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    virtual ~RFState(){}
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        otCLI.out << inducedRFDist(*taxonomy, *tree, true) << '\n';
        return true;
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcinducedsubtree",
                "takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes)",
                "taxonomy.tre inp1.tre inp2.tre");
    RFState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
