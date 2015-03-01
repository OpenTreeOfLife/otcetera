#include "otc/otcli.h"
using namespace otc;

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg);
void extendSupportedToRedundantNodes(const TreeMappedWithSplits & tree, std::set<const NodeWithSplits *> & supportedNodes);
bool singleDesSupportedOrNamed(const NodeWithSplits *nd, const std::set<const NodeWithSplits *> & supportedNodes);

void extendSupportedToRedundantNodes(const TreeMappedWithSplits & tree, std::set<const NodeWithSplits *> & supportedNodes) {
	for (auto nd : ConstPostorderInternalIter<TreeMappedWithSplits>(tree)) {
		if (nd->isOutDegreeOneNode()) {
			auto c = nd->getFirstChild();
			if (c->hasOttId() || supportedNodes.find(c) != supportedNodes.end()) {
				supportedNodes.insert(nd);
			}
		}
	}
}

bool singleDesSupportedOrNamed(const NodeWithSplits *nd, const std::set<const NodeWithSplits *> & supportedNodes) {
	if (supportedNodes.find(nd) != supportedNodes.end()) {
		return true;
	}
	if (nd->isOutDegreeOneNode()) {
		return nd->hasOttId() || singleDesSupportedOrNamed(nd->getFirstChild(), supportedNodes);
	}
	return false;
}

struct FindUnsupportedState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
	std::unique_ptr<TreeMappedWithSplits> toCheck;
	int numErrors;
	std::map<const NodeWithSplits *, std::set<long> > aPrioriProblemNodes;
	std::set<const NodeWithSplits *> supportedNodes;
	virtual ~FindUnsupportedState(){}
	FindUnsupportedState()
		:toCheck(nullptr),
		 numErrors(0) {
		}

	int describeUnnamedUnsupported(std::ostream &out, const TreeMappedWithSplits & tree,
								   const std::set<const NodeWithSplits *> & supported) const {
		auto ig = ConstPreorderInternalIter<TreeMappedWithSplits>(tree);
		auto nIt = ig.begin();
		const auto eIt = ig.end();
		int numUnsupported = 0;
		++nIt; //skip the root
		for (; nIt != eIt; ++nIt) {
			auto nd = *nIt;
			if (supported.find(nd) != supported.end()) {
				continue;
			}
			if (nd->includesOnlyOneLeaf()) {
				continue;
			}
			auto outDegree = nd->getOutDegree();
			if (outDegree == 1 && singleDesSupportedOrNamed(nd, supported)) {
				continue;
			}
			if (!nd->hasOttId()) { //we don't check named nodes.
				if (aPrioriProblemNodes.empty()) {
					out << "Unsupported node ";
				} else {
					auto gaIt = aPrioriProblemNodes.find(nd);
					if (gaIt == aPrioriProblemNodes.end()) {
						out << "Novel unsupported node ";
					} else {
						out << "Confirmation of unsupported node (designators =";
						writeOttSet(out, "", gaIt->second, " ");
						out << ") ";
					}
				}
				describeUnnamedNode(*nd, out, 0, false);
				numUnsupported += 1;
			}
		}
		return numUnsupported;
	}

	bool summarize(const OTCLI &otCLI) override {
		extendSupportedToRedundantNodes(*toCheck, supportedNodes);
		auto & out = otCLI.out;
		int numUnsupported = describeUnnamedUnsupported(otCLI.out, *toCheck, supportedNodes);
		for (auto gaIt : aPrioriProblemNodes) {
			if (supportedNodes.find(gaIt.first) != supportedNodes.end()) {
				out << "Claim of unsupported apparently refuted for designators: ";
				writeOttSet(out, "", gaIt.second, " ");
				out << ". See standard error stream for details.\n";
			}
		} 
		int supNonNamed = 0;
		int numSupportedInternals = 0;
		for (auto snd : supportedNodes) {
			if (!snd->isTip()) {
				numSupportedInternals += 1;
				if (!snd->hasOttId()) {
					supNonNamed += 1;
				}
			}
		}
		out << "Final summary:\n";
		//out << gRefTreeNumNamedInternalsNodes << " internal nodes were named in the reference tree. These were not rigorously checked against the taxonomy. They may not be detected as errors.\n";
		out << numSupportedInternals << " internal nodes where flagged as being supported by an input (including taxonomy).\n";
		int supNamed = numSupportedInternals - supNonNamed;
		out << "    " << supNamed << " of these were named (some of the support could just be the taxonomic expansion of tips).\n";
		out << "    " << supNonNamed << " of these were unnamed.\n";
		out << numUnsupported << " unsupported nodes.";
		out << std::endl;
		if (numErrors < 0) {
			numErrors -= numUnsupported;
		} else {
			numErrors = numUnsupported;
		}
		return numErrors == 0;
	}

	void parseAndProcessMRCADesignatorsFile(const std::string &fp) {
		if (toCheck == nullptr) {
			throw OTCError("Designator files (if used) must be passed in after the tree to check");
		}
		std::list<std::set<long> > dl = parseDesignatorsFile(fp);
		for (auto d : dl) {
			markSuspectNode(d);
		}
	}

	void markSuspectNode(const std::set<long> & designators) {
		const NodeWithSplits * mrca = findMRCAFromIDSet(*toCheck, designators, -1);
		aPrioriProblemNodes[mrca] = designators;
	}

	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
		assert(taxonomy != nullptr);
		if (toCheck == nullptr) {
			toCheck = std::move(tree);
			return true;
		}
		expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
		return processExpandedTree(otCLI, *tree);
	}

	bool processExpandedTree(OTCLI & otCLI, const TreeMappedWithSplits & tree) {
		assert(toCheck != nullptr);
		std::map<const NodeWithSplits *, std::set<long> > prunedDesId;
		for (auto nd : ConstLeafIter<TreeMappedWithSplits>(tree)) {
			auto ottId = nd->getOttId();
			markPathToRoot(*toCheck, ottId, prunedDesId);
		}
		std::map<std::set<long>, const NodeWithSplits *> sourceClades;
		for (auto nd : ConstPostorderInternalIter<TreeMappedWithSplits>(tree)) {
			if (nd->getParent() != nullptr && !nd->isTip()) {
				sourceClades[nd->getData().desIds] = nd;
			}
		}
		recordSupportedNodes(otCLI, prunedDesId, sourceClades, supportedNodes);
		return true;
	}

	void recordSupportedNodes(OTCLI & otCLI,
							  const std::map<const NodeWithSplits *, std::set<long> > & prunedDesId,
							  const std::map<std::set<long>, const NodeWithSplits *> & sourceClades,
							  std::set<const NodeWithSplits *> & supported) {
		for (auto pd : prunedDesId) {
			auto nd = pd.first;
			auto par = nd->getParent();
			if (par == nullptr) {
				//otCLI.out << "  par null\n";
				continue;
			}
			auto firstBranchingAnc = findFirstBranchingAnc<const NodeWithSplits>(nd);
			if (firstBranchingAnc == nullptr) {
				//otCLI.out << "  firstBranchingAnc null\n";
				continue;
			}
			auto nm = pd.second;
			auto ancIt = prunedDesId.find(firstBranchingAnc);
			assert(ancIt != prunedDesId.end());
			auto anm = ancIt->second;
			const NodeWithSplits * firstNdPtr; // just used to match call
			if (!multipleChildrenInMap(*nd, prunedDesId, &firstNdPtr)) {
				//otCLI.out << "  multipleChildrenInMap false\n";
				continue;
			}
			if (anm == nm) {
				//otCLI.out << "  anc set == node set\n";
				continue;
			}
			auto scIt = sourceClades.find(nm);
			if (scIt != sourceClades.end()) {
				//otCLI.out << "  supported\n";
				if (aPrioriProblemNodes.find(nd) != aPrioriProblemNodes.end()) {
					auto apIt = aPrioriProblemNodes.find(nd);
					otCLI.out << "ERROR!: a priori unsupported node found. Designators were ";
					writeOttSet(otCLI.out, "", apIt->second, " ");
					otCLI.out << ". A node was found, which (when pruned to the leaf set of an input tree) contained:\n";
					writeOttSet(otCLI.out, "    ", nm, " ");
					otCLI.out << "\nThe subtree from the source was: ";
					auto srcNd = scIt->second;
					writePrunedSubtreeNewickForMarkedNodes(otCLI.out, *srcNd, prunedDesId);
					numErrors += 1;
				}
				supported.insert(nd);
			} else {
				//otCLI.out << "  unsupported\n";
			}
		}
	}
};

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg) {
	FindUnsupportedState * fusp = static_cast<FindUnsupportedState *>(otCLI.blob);
	assert(fusp != nullptr);
	assert(!nextArg.empty());
	fusp->parseAndProcessMRCADesignatorsFile(nextArg);
	return true;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcfindunsupportednodes",
				"takes at least 2 newick file paths: a full taxonomy tree, a full supertree, and some number of input trees",
				"taxonomy.tre synth.tre inp1.tre inp2.tre");
	FindUnsupportedState proc;
	otCLI.addFlag('m',
				  "ARG=a designators file. Each line is a list of (white-space separated) OTT ids used to designate the node that is the MRCA of them.",
				  handleDesignator,
				  true);
	return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, false);
}
