#include "otc/otcli.h"
#include "otc/util.h"
using namespace otc;
typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

struct UnionOfIdsState {
    std::set<long> idsEncountered;
    int numErrors;
    UnionOfIdsState()
        :numErrors(0) {
    }
    void summarize(OTCLI &otCLI) {
        for (const auto & oid : idsEncountered) {
            otCLI.out << oid << '\n';
        }
    }
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
    UnionOfIdsState * ctsp = static_cast<UnionOfIdsState *>(otCLI.blob);
    assert(ctsp != nullptr);
    std::set<long> & ier = ctsp->idsEncountered;
    for (const auto nd : iter_node_const(*tree)) {
        if (nd->hasOttId()) {
            ier.insert(nd->getOttId());
        }
    }
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-union-of-ids",
                "takes a series of newick file paths and writes the union of the OTT Ids found in all trees",
                "some.tre");
    UnionOfIdsState cts;
    otCLI.blob = static_cast<void *>(&cts);
    auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 1);
    if (rc == 0) {
        cts.summarize(otCLI);
        return cts.numErrors;
    }
    return rc;
}
