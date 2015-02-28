#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

template<typename T>
class TaxonomyDependentTreeProcessor {
	public:
		using node_type = typename T::node_type;
		using node_data_type = typename T::node_type::data_type;
		using tree_data_type = typename T::data_type;
		using tree_type = T;

		std::unique_ptr<Tree_t> taxonomy;
		std::set<long> ottIds;

		virtual bool processTaxonomyTree(OTCLI & otCLI) {
			ottIds = taxonomy->getRoot()->getData().desIds;
			otCLI.getParsingRules().ottIdValidator = &ottIds;
			otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
			return true;
		}
		virtual bool processSourceTree(OTCLI & , std::unique_ptr<Tree_t> tree) {
			assert(tree != nullptr);
			assert(taxonomy != nullptr);
			return true;
		}
		virtual bool summarize(const OTCLI &) {
			return true;
		}
		TaxonomyDependentTreeProcessor()
			:taxonomy(nullptr) {
		}
		virtual ~TaxonomyDependentTreeProcessor(){}
	
};

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

template<typename T>
inline bool taxDependentProcessNextTree(OTCLI & otCLI, std::unique_ptr<T> tree) {
	TaxonomyDependentTreeProcessor<T> * tdtp = static_cast<TaxonomyDependentTreeProcessor<T> *>(otCLI.blob);
	assert(tdtp != nullptr);
	assert(tree != nullptr);
	if (tdtp->taxonomy == nullptr) {
		tdtp->taxonomy = std::move(tree);
		return tdtp->processTaxonomyTree(otCLI);
	}
	return tdtp->processSourceTree(otCLI, std::move(tree));
}

template<typename T>
int taxDependentTreeProcessingMain(OTCLI & otCLI,
								   int argc,
								   char *argv[],
								   TaxonomyDependentTreeProcessor<T> & proc,
								   unsigned int numTrees,
								   bool includeInternalNodesInDesIdSets) {
	assert(otCLI.blob == nullptr);
	otCLI.blob = static_cast<void *>(&proc);
	otCLI.getParsingRules().includeInternalNodesInDesIdSets = includeInternalNodesInDesIdSets;
	std::function<bool (OTCLI &, std::unique_ptr<T>)> pcb = taxDependentProcessNextTree<T>;
	auto rc = treeProcessingMain<T>(otCLI, argc, argv, pcb, nullptr, numTrees);
	if (rc == 0) {
		return (proc.summarize(otCLI) ? 0 : 1);
	}
	return rc;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcinducedsubtree",
				"takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes)",
				"taxonomy.tre inp1.tre inp2.tre");
	InducedSubtreeState proc;
	return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}


