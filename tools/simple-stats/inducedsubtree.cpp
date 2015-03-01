#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;

struct InducedSubtreeState : public TaxonomyDependentTreeProcessor<Tree_t> {
	std::set<long> inducingLabels;
	virtual ~InducedSubtreeState(){}

	virtual bool summarize(const OTCLI &otCLI) override {
		auto mrca = findMRCAUsingDesIds(*taxonomy, inducingLabels);
		std::function<bool(const Node_t &)> sf = [this](const Node_t &nd){
			return haveIntersection(this->inducingLabels, nd.getData().desIds);
		};
		writeNewickFiltered(otCLI.out, mrca, sf);
		otCLI.out << ";\n";
		return true;
	}

	virtual bool processSourceTree(OTCLI &,
								   std::unique_ptr<Tree_t> tree) override {
		assert(tree != nullptr);
		assert(taxonomy != nullptr);
		auto ls = getOttIdSetForLeaves(*tree);
		inducingLabels.insert(ls.begin(), ls.end());
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


