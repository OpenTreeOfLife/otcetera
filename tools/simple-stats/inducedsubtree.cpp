#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

struct InducedSubtreeState {
	std::unique_ptr<Tree_t> fulltree;
	std::set<long> ottIds;
	std::set<long> inducingLabels;
	
	InducedSubtreeState()
		:fulltree(nullptr) {
	}

	void summarize(const OTCLI &otCLI) {
		auto mrca = findMRCAUsingDesIds(*fulltree, inducingLabels);
		std::function<bool(const Node_t &)> sf = [this](const Node_t &nd){
			return haveIntersection(this->inducingLabels, nd.getData().desIds);
		};
		writeNewickFiltered(otCLI.out, mrca, sf);
		otCLI.out << ";\n";
	}

	bool processTaxonomyTree(OTCLI & otCLI) {
		ottIds = fulltree->getRoot()->getData().desIds;
		otCLI.getParsingRules().ottIdValidator = &ottIds;
		otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
		return true;
	}

	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> treeup) {
		assert(treeup != nullptr);
		assert(fulltree != nullptr);
		Tree_t & tree{*treeup};
		for (auto nd : ConstLeafIter<Tree_t>(tree)) {
			inducingLabels.insert(nd->getOttId());
		}
		return true;
	}
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	InducedSubtreeState * ctsp = static_cast<InducedSubtreeState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->fulltree == nullptr) {
		ctsp->fulltree = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	return ctsp->processSourceTree(otCLI, std::move(tree));
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcinducedsubtree",
				"takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes)",
				"taxonomy.tre inp1.tre inp2.tre");
	InducedSubtreeState fus;
	otCLI.blob = static_cast<void *>(&fus);
	otCLI.getParsingRules().includeInternalNodesInDesIdSets = true;
	auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
	if (rc == 0) {
		fus.summarize(otCLI);
		return 0;
	}
	return rc;
}
