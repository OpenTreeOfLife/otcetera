#include "otc/otcli.h"
#include "otc/debug.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef RTreeOttIDMapping<RTSplits> RootedTreeForNodeType;
typedef otc::RootedTree<RTSplits, RootedTreeForNodeType> Tree_t;

bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
bool processNextTree(OTCLI & , std::unique_ptr<Tree_t> tree) {
    check_tree_invariants(*tree);
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-assert-invariants",
                "takes a tree file, parses it and runs a series of checks of invariants of the tree operations",
                "test.tre");
    return tree_processing_main<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 1);
}
