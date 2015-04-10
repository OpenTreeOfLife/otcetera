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
    OTCLI otCLI("otc-rf",
                "takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the Robinson-Foulds symmetric difference between the induced tree from the full tree to each input tree (one RF distance per line).",
                "synth.tre inp1.tre inp2.tre");
    RFState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
