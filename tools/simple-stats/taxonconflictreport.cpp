#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

struct RFState {
	std::unique_ptr<Tree_t> fulltree;
	std::set<long> ottIds;
	RFState()
		:fulltree(nullptr) {
	}
	bool processTaxonomyTree(OTCLI & otCLI) {
		ottIds = fulltree->getRoot()->getData().desIds;
		otCLI.getParsingRules().ottIdValidator = &ottIds;
		otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
		return true;
	}

	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
		assert(tree != nullptr);
		assert(fulltree != nullptr);
		reportOnInducedConflicts(otCLI.out, *fulltree, *tree, true);
		return true;
	}
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	RFState * ctsp = static_cast<RFState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->fulltree == nullptr) {
		ctsp->fulltree = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	return ctsp->processSourceTree(otCLI, std::move(tree));
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcmissingsplits",
				"takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the number of splits in the induced full tree that are missing from the each other tree.",
				"taxonomy.tre inp1.tre inp2.tre");
	RFState fus;
	otCLI.blob = static_cast<void *>(&fus);
	otCLI.getParsingRules().includeInternalNodesInDesIdSets = true;
	auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
	return rc;
}
