#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

template<typename T>
bool desIdSetsConflict(const T & ns, const T &scs) {
	if (ns.size() < 2 || scs.size() < 2) {
		return false;
	}
	T sinter;
	std::set_intersection(begin(ns), end(ns), begin(scs), end(scs), std::inserter(sinter, begin(sinter)));
	const auto soi = sinter.size();
	return !(soi == 0 || soi == ns.size() || soi == scs.size());
}

struct DetectContestedState {
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;
	std::set<long> ottIds;
	std::set<const Node_t *> contestedNodes;

	DetectContestedState()
		:taxonomy(nullptr),
		 numErrors(0) {
		}

	void summarize(const OTCLI &otCLI) {
		assert (taxonomy != nullptr);
		for (auto nd : contestedNodes) {
			otCLI.out << nd->getName() << '\n';
		}
	}

	bool processTaxonomyTree(OTCLI & otCLI) {
		ottIds = keys(taxonomy->getData().ottIdToNode);
		otCLI.getParsingRules().ottIdValidator = &ottIds;
		return true;
	}

	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
		assert(taxonomy != nullptr);
		assert(tree != nullptr);
		expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
		return processExpandedTree(otCLI, *tree);
	}

	bool processExpandedTree(OTCLI & otCLI, Tree_t & tree) {
		std::map<const Node_t *, std::set<long> > prunedDesId;
		for (auto nd : ConstLeafIter<Tree_t>(tree)) {
			auto ottId = nd->getOttId();
			markPathToRoot(*taxonomy, ottId, prunedDesId);
		}
		std::set<std::set<long> > sourceClades;
		for (auto nd : PostorderInternalIter<Tree_t>(tree)) {
			if (nd->getParent() != nullptr && !nd->isTip()) {
				sourceClades.insert(std::move(nd->getData().desIds));
			}
		}
		auto numLeaves = tree.getRoot()->getData().desIds.size();
		recordContested(prunedDesId, sourceClades, contestedNodes, numLeaves);
		return true;
	}

	void recordContested(const std::map<const Node_t *, std::set<long> > & prunedDesId,
						 const std::set<std::set<long> > & sourceClades,
						 std::set<const Node_t *> & contestedSet,
						 std::size_t numLeaves) {
		for (auto pd : prunedDesId) {
			auto nd = pd.first;
			if (contains(contestedSet, nd)) {
				continue;
			}
			const auto & ns = nd->getData().desIds;
			const auto nss = ns.size();
			if (nss == 1 || nss == numLeaves) {
				continue;
			}
			for (auto sc : sourceClades) {
				if (desIdSetsConflict(ns, sc)) {
					contestedSet.insert(nd);
				}
			}
		}
	}
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	DetectContestedState * ctsp = static_cast<DetectContestedState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->taxonomy == nullptr) {
		ctsp->taxonomy = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	return ctsp->processSourceTree(otCLI, std::move(tree));
}


int main(int argc, char *argv[]) {
	OTCLI otCLI("otcdetectcontested",
				"takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Writes the OTT IDs of clades in the taxonomy whose monophyly is questioned by at least one input",
				"taxonomy.tre inp1.tre inp2.tre");
	DetectContestedState fus;
	otCLI.blob = static_cast<void *>(&fus);
	auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
	if (rc == 0) {
		fus.summarize(otCLI);
		return fus.numErrors;
	}
	return rc;
}
