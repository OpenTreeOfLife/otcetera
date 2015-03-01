#include "otc/otcli.h"
using namespace otc;

template<typename T, typename U>
struct NodePairing {
	T * scaffoldNode;
	U * phyloNode;
	NodePairing(T *taxo, U *phylo)
		:scaffoldNode(taxo),
		phyloNode(phylo) {
		assert(taxo != nullptr);
		assert(phylo != nullptr);
	}
};

template<typename T, typename U>
struct PathPairing {
	const T * const phyloChild;
	const T * const phyloParent;
	const U * const scaffoldDes;
	const U * const scaffoldAnc;

	PathPairing(const NodePairing<T,U> & parent, const NodePairing<T,U> & child)
		:phyloChild(child.phyloNode),
		phyloParent(parent.phyloNode),
		scaffoldDes(child.scaffoldNode),
		scaffoldAnc(parent.scaffoldNode) {
	}
};

template<typename T, typename U>
struct NodeThreading {
	using NodePairSet = std::set<NodePairing<T, U> *>;
	using PathPairSet = std::set<PathPairing<T, U> *>;
	
	std::map<std::size_t, NodePairSet> nodeAlignments;
	std::map<std::size_t, PathPairSet> edgeBelowAlignments;
	std::map<std::size_t, PathPairSet > loopAlignments;
	std::size_t getTotalNumNodeMappings() const {
		unsigned long t = 0U;
		for (auto i : nodeAlignments) {
			t += i.second.size();
		}
		return t;
	}
	std::size_t getTotalNumLoops() const {
		unsigned long t = 0U;
		for (auto i : loopAlignments) {
			t += i.second.size();
		}
		return t;
	}
	std::size_t getTotalNumEdgeBelowTraversals() const {
		unsigned long t = 0U;
		for (auto i : edgeBelowAlignments) {
			t += i.second.size();
		}
		return t;
	}
	bool isContested() const {
		for (auto i : edgeBelowAlignments) {
			if (i.second.size() > 1) {
				const T * firstSrcPar = nullptr;
				for (auto pp : i.second) {
					auto sp = pp->phyloParent;
					if (sp != firstSrcPar) {
						if (firstSrcPar == nullptr) {
							firstSrcPar = sp;
						} else {
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	std::list<std::size_t> getContestingTrees() const {
		std::list<std::size_t> r;
		for (auto i : edgeBelowAlignments) {
			if (i.second.size() > 1) {
				r.push_back(i.first);
			}
		}
		return r;
	}
};

using NodePairingWithSplits = NodePairing<NodeWithSplits, NodeWithSplits>;
using PathPairingWithSplits = PathPairing<NodeWithSplits, NodeWithSplits>;
using NodeThreadingWithSplits = NodeThreading<NodeWithSplits, NodeWithSplits>;

class ThreadedTree {
	protected:
	std::list<NodePairingWithSplits> nodePairings;
	std::list<PathPairingWithSplits> pathPairings;
	std::map<const NodeWithSplits*, NodeThreadingWithSplits> taxoToAlignment;
	public:
	NodePairingWithSplits * _addNodeMapping(NodeWithSplits *taxo, NodeWithSplits *nd, std::size_t treeIndex) {
		assert(taxo != nullptr);
		assert(nd != nullptr);
		nodePairings.emplace_back(NodePairingWithSplits(taxo, nd));
		auto ndPairPtr = &(*nodePairings.rbegin());
		auto & athreading = taxoToAlignment[taxo];
		athreading.nodeAlignments[treeIndex].insert(ndPairPtr);
		return ndPairPtr;
	}
	PathPairingWithSplits * _addPathMapping(NodePairingWithSplits * parentPairing,
											NodePairingWithSplits * childPairing,
											std::size_t treeIndex) {
		pathPairings.emplace_back(*parentPairing, *childPairing);
		auto pathPairPtr = &(*pathPairings.rbegin());
		// register a pointer to the path at each traversed...
		auto currTaxo = pathPairPtr->scaffoldDes;
		auto ancTaxo = pathPairPtr->scaffoldAnc;
		if (currTaxo != ancTaxo) {
			while (currTaxo != ancTaxo) {
				taxoToAlignment[currTaxo].edgeBelowAlignments[treeIndex].insert(pathPairPtr);
				currTaxo = currTaxo->getParent();
				if (currTaxo == nullptr) {
					break;
				}
			}
		} else {
			taxoToAlignment[currTaxo].loopAlignments[treeIndex].insert(pathPairPtr);
		}
		return pathPairPtr;
	}
	void threadNewTree(TreeMappedWithSplits & scaffoldTree, TreeMappedWithSplits & tree, std::size_t treeIndex) {
		// do threading
		std::map<NodeWithSplits *, NodePairingWithSplits *> currTreeNodePairings;
		std::set<NodePairingWithSplits *> tipPairings;
		for (auto nd : PostorderIter<TreeMappedWithSplits>(tree)) {
			auto par = nd->getParent();
			if (par == nullptr) {
				continue;
			}
			NodePairingWithSplits * ndPairPtr = nullptr;
			NodeWithSplits * taxoDes = nullptr;
			if (nd->isTip()) {
				assert(currTreeNodePairings.find(nd) == currTreeNodePairings.end()); // TMP, Remove this to save time?
				assert(nd->hasOttId());
				auto ottId = nd->getOttId();
				taxoDes = scaffoldTree.getData().getNodeForOttId(ottId);
				assert(taxoDes != nullptr);
				ndPairPtr = _addNodeMapping(taxoDes, nd, treeIndex);
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
			NodePairingWithSplits * parPairPtr = nullptr;
			auto prevAddedNodePairingIt = currTreeNodePairings.find(par);
			if (prevAddedNodePairingIt == currTreeNodePairings.end()) {
				const auto & parDesIds = par->getData().desIds;
				auto taxoAnc = searchAncForMRCAOfDesIds(taxoDes, parDesIds);
				assert(taxoAnc != nullptr);
				parPairPtr = _addNodeMapping(taxoAnc, par, treeIndex);
				currTreeNodePairings[par] = parPairPtr;
			} else {
				parPairPtr = prevAddedNodePairingIt->second;
			}
			_addPathMapping(parPairPtr, ndPairPtr, treeIndex);
		}
	}

};

struct RemapToDeepestUnlistedState
	: public TaxonomyDependentTreeProcessor<TreeMappedWithSplits>,
	public ThreadedTree {
	int numErrors;
	std::set<long> tabooIds;
	std::map<std::unique_ptr<TreeMappedWithSplits>, std::size_t> inputTreesToIndex;
	std::vector<TreeMappedWithSplits *> treePtrByIndex;

	virtual ~RemapToDeepestUnlistedState(){}
	RemapToDeepestUnlistedState()
		:TaxonomyDependentTreeProcessor<TreeMappedWithSplits>(),
		 numErrors(0) {
	}

	bool summarize(const OTCLI &otCLI) override {
		assert (taxonomy != nullptr);
		std::map<std::size_t, unsigned long> nodeMappingDegree;
		std::map<std::size_t, unsigned long> passThroughDegree;
		std::map<std::size_t, unsigned long> loopDegree;
		unsigned long totalContested = 0;
		unsigned long redundContested = 0;
		unsigned long totalNumNodes = 0;
		for (auto nd : InternalNodeIter<TreeMappedWithSplits>(*taxonomy)) {
			const auto & thr = taxoToAlignment[nd];
			nodeMappingDegree[thr.getTotalNumNodeMappings()] += 1;
			passThroughDegree[thr.getTotalNumEdgeBelowTraversals()] += 1;
			loopDegree[thr.getTotalNumLoops()] += 1;
			totalNumNodes += 1;
			if (thr.isContested()) {
				auto c = thr.getContestingTrees();
				for (auto cti : c) {
					auto ctree = treePtrByIndex.at(cti);
					otCLI.err << nd->getOttId() << " \"" << nd->getName() << "\" contested by \"" << ctree->getName() << "\"\n";
				}
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
		return true;
	}

	bool processTaxonomyTree(OTCLI & otCLI) override {
		TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
		suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
		for (auto nd : NodeIter<TreeMappedWithSplits>(*taxonomy)) {
			taxoToAlignment.emplace(nd, NodeThreadingWithSplits{});
		}
		return true;
	}

	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
		assert(treeup != nullptr);
		assert(taxonomy != nullptr);
		// Store the tree pointer with a map to its index, and an alias for fast index->tree.
		std::size_t treeIndex = inputTreesToIndex.size();
		assert(treeIndex == treePtrByIndex.size());
		TreeMappedWithSplits * raw = treeup.get();
		inputTreesToIndex[std::move(treeup)] = treeIndex;
		treePtrByIndex.push_back(raw);
		// Store the tree's filename
		raw->setName(otCLI.currentFilename);
		threadNewTree(*taxonomy, *raw, treeIndex);
		otCLI.out << "# pathPairings = " << pathPairings.size() << '\n';
		return true;
	}
};

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcscaffoldedsupertree",
				"takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Crashes or emits bogus output.",
				"taxonomy.tre inp1.tre inp2.tre");
	RemapToDeepestUnlistedState proc;
	return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
