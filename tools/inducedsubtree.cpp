#include "otc/otcli.h"
using namespace otc;

struct InducedSubtreeState
  : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::set<long> inducingIds;
    virtual ~InducedSubtreeState(){}

    // write the induced tree to the output stream
    virtual bool summarize(OTCLI &otCLI) override {
        // find the LIA of the induced taxa
        auto mrca = findMRCAUsingDesIds(*taxonomy, inducingIds);
        // do a filtered pre-order traversal
        NodeWithSplitsPred sf = [this](const NodeWithSplits &nd){
            return haveIntersection(this->inducingIds, nd.getData().desIds);
        };
        writeNewickFiltered(otCLI.out, mrca, sf);
        otCLI.out << ";\n";
        return true;
    }

    // accumulate the set of leaves to include in inducingIds
    virtual bool processSourceTree(OTCLI &,
                                   std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        auto ls = getOttIdSetForLeaves(*tree);
        inducingIds.insert(ls.begin(), ls.end());
        return true;
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcinducedsubtree",
                "takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes)",
                "taxonomy.tre inp1.tre inp2.tre");
    InducedSubtreeState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}


