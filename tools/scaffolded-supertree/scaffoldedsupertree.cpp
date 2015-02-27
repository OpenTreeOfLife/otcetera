#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<typename Node_t::data_type, RTreeOttIDMapping<typename Node_t::data_type>> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
bool handleTabooOTTIdListFile(OTCLI & otCLI, const std::string &nextArg);

struct NodePairing {
	Node_t * scaffoldNode;
	Node_t * phyloNode;
	NodePairing(Node_t *taxo, Node_t *phylo)
		:scaffoldNode(taxo),
		phyloNode(phylo) {
		assert(taxo != nullptr);
		assert(phylo != nullptr);
	}
};
struct PathPairing {
	const Node_t * const phyloChild;
	const Node_t * const phyloParent;
	const Node_t * const scaffoldDes;
	const Node_t * const scaffoldAnc;

	PathPairing(const NodePairing & parent, const NodePairing & child)
		:phyloChild(child.phyloNode),
		phyloParent(parent.phyloNode),
		scaffoldDes(child.scaffoldNode),
		scaffoldAnc(parent.scaffoldNode) {
	}
};
struct AlignmentThreading {
	std::map<Tree_t *, std::set<NodePairing *> > nodeAlignments;
	std::map<Tree_t *, std::set<PathPairing *> > edgeBelowAlignments;
	std::map<Tree_t *, std::set<PathPairing *> > loopAlignments;
	unsigned long getTotalNumNodeMappings() const {
		unsigned long t = 0U;
		for (auto i : nodeAlignments) {
			t += i.second.size();
		}
		return t;
	}
	unsigned long getTotalNumLoops() const {
		unsigned long t = 0U;
		for (auto i : loopAlignments) {
			t += i.second.size();
		}
		return t;
	}
	unsigned long getTotalNumEdgeBelowTraversals() const {
		unsigned long t = 0U;
		for (auto i : edgeBelowAlignments) {
			t += i.second.size();
		}
		return t;
	}
	bool isContested() const {
		for (auto i : edgeBelowAlignments) {
			if (i.second.size() > 1) {
				return true;
			}
		}
		return false;
	}
	std::list<Tree_t *> getContestingTrees() const {
		std::list<Tree_t *> r;
		for (auto i : edgeBelowAlignments) {
			if (i.second.size() > 1) {
				r.push_back(i.first);
			}
		}
		return r;
	}
};

struct RemapToDeepestUnlistedState {
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;
	std::set<long> ottIds;
	std::set<long> tabooIds;
	std::list<std::unique_ptr<Tree_t> > inputTrees;
	std::list<NodePairing> nodePairings;
	std::list<PathPairing> pathPairings;
	std::map<const Node_t*, AlignmentThreading> taxoToAlignment;

	RemapToDeepestUnlistedState()
		:taxonomy(nullptr),
		 numErrors(0) {
	}

	void summarize(const OTCLI &otCLI) {
		assert (taxonomy != nullptr);
		std::map<std::size_t, unsigned long> nodeMappingDegree;
		std::map<std::size_t, unsigned long> passThroughDegree;
		std::map<std::size_t, unsigned long> loopDegree;
		unsigned long totalContested = 0;
		unsigned long redundContested = 0;
		unsigned long totalNumNodes = 0;
		for (auto nd : InternalNodeIter<Tree_t>(*taxonomy)) {
			const auto & thr = taxoToAlignment[nd];
			nodeMappingDegree[thr.getTotalNumNodeMappings()] += 1;
			passThroughDegree[thr.getTotalNumEdgeBelowTraversals()] += 1;
			loopDegree[thr.getTotalNumLoops()] += 1;
			totalNumNodes += 1;
			if (thr.isContested()) {
				otCLI.err << nd->getOttId() << " \"" << nd->getName() << "\" contested by";
				auto c = thr.getContestingTrees();
				for (auto ct : c) {
					otCLI.err << " \"" << ct->getName() << "\"";
				}
				otCLI.err << "\n";
				totalContested += 1;
				if (nd->getOutDegree() == 1) {
					redundContested += 1;
				}
			}
		}
		unsigned long m = std::max(loopDegree.rbegin()->first, passThroughDegree.rbegin()->first);
		m = std::max(m, nodeMappingDegree.rbegin()->first);
		otCLI.out << "Degree\tNodeMaps\tEdgeMaps\tLoops\n";
		for (unsigned long i = 0 ; i <= m; ++i) {
			otCLI.out << i << '\t' << nodeMappingDegree[i]<< '\t' << passThroughDegree[i] << '\t' << loopDegree[i]<< '\n';
		}
		otCLI.out << totalNumNodes << " internals\n" << totalContested << " contested\n" << (totalNumNodes - totalContested) << " uncontested\n";
		otCLI.out << redundContested << " monotypic contested\n";
	}

	bool processTaxonomyTree(OTCLI & otCLI) {
		ottIds = taxonomy->getRoot()->getData().desIds;
		suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
		for (auto nd : NodeIter<Tree_t>(*taxonomy)) {
			taxoToAlignment.emplace(nd, AlignmentThreading{});
		}
		otCLI.getParsingRules().ottIdValidator = &ottIds;
		otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
		return true;
	}

	NodePairing * _addNodeMapping(Node_t *taxo, Node_t *nd, Tree_t *tree) {
		assert(taxo != nullptr);
		assert(nd != nullptr);
		nodePairings.emplace_back(NodePairing(taxo, nd));
		auto ndPairPtr = &(*nodePairings.rbegin());
		auto & athreading = taxoToAlignment[taxo];
		athreading.nodeAlignments[tree].insert(ndPairPtr);
		return ndPairPtr;
	}
	PathPairing * _addPathMapping(NodePairing * parentPairing, NodePairing * childPairing, Tree_t *tree) {
		pathPairings.emplace_back(*parentPairing, *childPairing);
		auto pathPairPtr = &(*pathPairings.rbegin());
		// register a pointer to the path at each traversed...
		auto currTaxo = pathPairPtr->scaffoldDes;
		auto ancTaxo = pathPairPtr->scaffoldAnc;
		if (currTaxo != ancTaxo) {
			while (currTaxo != ancTaxo) {
				taxoToAlignment[currTaxo].edgeBelowAlignments[tree].insert(pathPairPtr);
				currTaxo = currTaxo->getParent();
				if (currTaxo == nullptr) {
					break;
				}
			}
		} else {
			taxoToAlignment[currTaxo].loopAlignments[tree].insert(pathPairPtr);
		}
		return pathPairPtr;
	}
	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<Tree_t> treeup) {
		assert(treeup != nullptr);
		assert(taxonomy != nullptr);
		inputTrees.emplace_back(std::move(treeup));
		Tree_t & tree{**inputTrees.rbegin()};
		tree.setName(otCLI.currentFilename);
		std::map<Node_t *, NodePairing *> currTreeNodePairings;
		std::set<NodePairing *> tipPairings;
		for (auto nd : PostorderIter<Tree_t>(tree)) {
			auto par = nd->getParent();
			if (par == nullptr) {
				continue;
			}
			NodePairing * ndPairPtr = nullptr;
			Node_t * taxoDes = nullptr;
			if (nd->isTip()) {
				assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
				assert(nd->hasOttId());
				auto ottId = nd->getOttId();
				taxoDes = taxonomy->getData().getNodeForOttId(ottId);
				assert(taxoDes != nullptr);
				ndPairPtr = _addNodeMapping(taxoDes, nd, &tree);
				for (auto former : tipPairings) {
					if (areLinearlyRelated(taxoDes, former->scaffoldNode)) {
						std::string m = "Repeated or nested OTT ID in tip mapping of an input tree: \"";
						m += nd->getName();
						m += "\" and \"";
						m += former->phyloNode->getName();
						m += "\" found.";
						throw OTCError(m);
					}
				}
				tipPairings.insert(ndPairPtr);
				currTreeNodePairings[nd] = ndPairPtr;
			} else {
				auto reuseNodePairingIt = currTreeNodePairings.find(nd);
				assert(reuseNodePairingIt != currTreeNodePairings.end());
				ndPairPtr = reuseNodePairingIt->second;
				taxoDes = ndPairPtr->scaffoldNode;
				assert(taxoDes != nullptr);
			}
			NodePairing * parPairPtr = nullptr;
			auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
			if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
				const auto & parDesIds = par->getData().desIds;
				auto taxoAnc = searchAncForMRCAOfDesIds(taxoDes, parDesIds);
				assert(taxoAnc != nullptr);
				parPairPtr = _addNodeMapping(taxoAnc, par, &tree);
				currTreeNodePairings[par] = parPairPtr;
			} else {
				parPairPtr = prevAddedNodePairingIt->second;
			}
			_addPathMapping(parPairPtr, ndPairPtr, &tree);
		}
		otCLI.out << "# pathPairings = " << pathPairings.size() << '\n';
		return true;
	}


	bool parseAndTabooOTTIdListFile(const std::string &fp) {
		auto t = parseListOfOttIds(fp);
		tabooIds.insert(begin(t), end(t));
		return true;
	}
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	RemapToDeepestUnlistedState * ctsp = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->taxonomy == nullptr) {
		ctsp->taxonomy = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	return ctsp->processSourceTree(otCLI, std::move(tree));
}

bool handleTabooOTTIdListFile(OTCLI & otCLI, const std::string &nextArg) {
	RemapToDeepestUnlistedState * fusp = static_cast<RemapToDeepestUnlistedState *>(otCLI.blob);
	assert(fusp != nullptr);
	assert(!nextArg.empty());
	return fusp->parseAndTabooOTTIdListFile(nextArg);
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcdetectcontested",
				"takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Writes the OTT IDs of clades in the taxonomy whose monophyly is questioned by at least one input",
				"taxonomy.tre inp1.tre inp2.tre");
	RemapToDeepestUnlistedState fus;
	otCLI.blob = static_cast<void *>(&fus);
	otCLI.getParsingRules().includeInternalNodesInDesIdSets = true;
	otCLI.addFlag('m',
				  "ARG=a file containing a list of taboo OTT ids.",
				  handleTabooOTTIdListFile,
				  true);
	auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 2);
	if (rc == 0) {
		fus.summarize(otCLI);
		return fus.numErrors;
	}
	return rc;
}
