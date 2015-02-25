#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTNodeNoData> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
struct PruneTaxonomyState {
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;
	std::set<long> ottIds;
	std::set<const Node_t *> includedNodes;

	PruneTaxonomyState()
		:taxonomy(nullptr),
		 numErrors(0) {
		}

	void summarize(const OTCLI &otCLI) {
		assert(taxonomy != nullptr && !includedNodes.empty());
		std::set<Node_t *> toPrune;
		for (auto nd : NodeIter<Tree_t>(*taxonomy)) {
			const Node_t *  c = const_cast<const Node_t *>(nd);
			if ((!contains(includedNodes, c)) && contains(includedNodes, c->getParent())) {
				toPrune.insert(nd);
			}
		}
		for (auto nd : toPrune) {
			pruneAndDelete(*taxonomy, nd);
		}
		writeTreeAsNewick(otCLI.out, *taxonomy);
		otCLI.out << '\n';
	}

	bool processTaxonomyTree(OTCLI & otCLI) {
		ottIds = keys(taxonomy->getData().ottIdToNode);
		otCLI.getParsingRules().ottIdValidator = &ottIds;
		return true;
	}
	bool processSourceTree(OTCLI & , const Tree_t & tree) {
		assert(taxonomy != nullptr);
		std::map<const Node_t *, std::set<long> > prunedDesId;
		for (auto nd : ConstLeafIter<Tree_t>(tree)) {
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

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	PruneTaxonomyState * ctsp = static_cast<PruneTaxonomyState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->taxonomy == nullptr) {
		ctsp->taxonomy = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	return ctsp->processSourceTree(otCLI, *tree);
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcprunetaxonomy",
				"takes at least 2 newick file paths: a full taxonomy tree some number of input trees. Prune subtrees from the taxonomy if they are not represented in the inputs.",
				"taxonomy.tre inp1.tre inp2.tre");
	PruneTaxonomyState fus;
	otCLI.blob = static_cast<void *>(&fus);
	auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
	if (rc == 0) {
		fus.summarize(otCLI);
		return fus.numErrors;
	}
	return rc;
}
