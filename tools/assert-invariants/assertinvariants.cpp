#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
#include "otc/debug.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef RTreeOttIDMapping<RTSplits> RootedTreeForNodeType;
typedef otc::RootedTree<RTSplits, RootedTreeForNodeType> Tree_t;

bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
extern const char * badNTreesMessage;
const char * badNTreesMessage = "Expecting a full tree, a full taxonomy, and some number of input trees";

inline bool processNextTree(OTCLI & , std::unique_ptr<Tree_t> tree) {
	checkTreeInvariants(*tree);
	return true;
}


int main(int argc, char *argv[]) {
	OTCLI otCLI("otcassertinvariants",
				"takes a tree file, parses it and runs a series of checks of invariants of the tree operations.",
				"test.tre");
	auto rc = treeProcessingMain<Tree_t>(otCLI,
																 argc,
																 argv,
																 processNextTree,
																 nullptr);
	return rc;
}
