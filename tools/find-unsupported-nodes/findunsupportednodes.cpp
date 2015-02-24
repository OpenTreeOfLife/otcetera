#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> MyNodeType;
typedef RTreeOttIDMapping<RTSplits> RootedTreeForNodeType;
typedef otc::RootedTree<RTSplits, RootedTreeForNodeType> Tree_t;

bool handleDesignator(OTCLI & otCLI, const std::string &nextArg);
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
extern const char * badNTreesMessage;
const char * badNTreesMessage = "Expecting a full tree, a full taxonomy, and some number of input trees";
void extendSupportedToRedundantNodes(const Tree_t & tree, std::set<const MyNodeType *> & supportedNodes);
bool singleDesSupportedOrNamed(const MyNodeType *nd, const std::set<const MyNodeType *> & supportedNodes);

void extendSupportedToRedundantNodes(const Tree_t & tree, std::set<const MyNodeType *> & supportedNodes) {
	for (auto nd : ConstPostorderInternalNode<Tree_t>(tree)) {
		if (nd->isOutDegreeOneNode()) {
			auto c = nd->getFirstChild();
			if (c->hasOttId() || supportedNodes.find(c) != supportedNodes.end()) {
				supportedNodes.insert(nd);
			}
		}
	}
}

bool singleDesSupportedOrNamed(const MyNodeType *nd, const std::set<const MyNodeType *> & supportedNodes) {
	if (supportedNodes.find(nd) != supportedNodes.end()) {
		return true;
	}
	if (nd->isOutDegreeOneNode()) {
		return nd->hasOttId() || singleDesSupportedOrNamed(nd->getFirstChild(), supportedNodes);
	}
	return false;
}

struct FindUnsupportedState {
	std::unique_ptr<Tree_t> toCheck;
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;
	std::map<const MyNodeType *, std::set<long> > aPrioriProblemNodes;
	std::set<long> ottIds;
	std::set<const MyNodeType *> supportedNodes;

	FindUnsupportedState()
		:toCheck(nullptr),
		 taxonomy(nullptr),
		 numErrors(0) {
		}

	int describeUnnamedUnsupported(std::ostream &out, const Tree_t & tree,
								   const std::set<const MyNodeType *> & supported) const {
		auto ig = ConstPreorderInternalNode<Tree_t>(tree);
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

	void summarize(const OTCLI &otCLI) {
		if (toCheck == nullptr) {
			otCLI.err << badNTreesMessage << '\n';
			numErrors = 1;
			return;
		}
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
		const MyNodeType * mrca = findMRCAFromIDSet(*toCheck, designators, -1);
		aPrioriProblemNodes[mrca] = designators;
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
	bool processExpandedTree(OTCLI & otCLI, const Tree_t & tree) {
		assert(toCheck != nullptr);
		std::map<const MyNodeType *, std::set<long> > prunedDesId;
		for (auto nd : ConstLeafIter<Tree_t>(tree)) {
			auto ottId = nd->getOttId();
			markPathToRoot(*toCheck, ottId, prunedDesId);
		}
		std::map<std::set<long>, const MyNodeType *> sourceClades;
		for (auto nd : ConstPostorderInternalNode<Tree_t>(tree)) {
			if (nd->getParent() != nullptr && !nd->isTip()) {
				sourceClades[nd->getData().desIds] = nd;
			}
		}
		recordSupportedNodes(otCLI, prunedDesId, sourceClades, supportedNodes);
		return true;
	}
	void recordSupportedNodes(OTCLI & otCLI,
							  const std::map<const MyNodeType *, std::set<long> > & prunedDesId,
							  const std::map<std::set<long>, const MyNodeType *> & sourceClades,
							  std::set<const MyNodeType *> & supported) {
		//otCLI.out << "sourceClades\n";
		//for (auto sc : sourceClades) {
		//	writeOttSet(otCLI.out, " ", sc.first, " ");
		//	otCLI.out << '\n';
		//}
		for (auto pd : prunedDesId) {
			//otCLI.out << "pruned el:";
			//writeOttSet(otCLI.out, " ", pd.second, " ");
			//otCLI.out << "\n";
			auto nd = pd.first;
			auto par = nd->getParent();
			if (par == nullptr) {
				//otCLI.out << "  par null\n";
				continue;
			}
			auto ls = pd.second;
			auto firstBranchingAnc = findFirstBranchingAnc<const MyNodeType>(nd);
			if (firstBranchingAnc == nullptr) {
				//otCLI.out << "  firstBranchingAnc null\n";
				continue;
			}
			auto nm = pd.second;
			auto ancIt = prunedDesId.find(firstBranchingAnc);
			assert(ancIt != prunedDesId.end());
			auto anm = ancIt->second;
			const MyNodeType * firstNdPtr; // just used to match call
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

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	FindUnsupportedState * ctsp = static_cast<FindUnsupportedState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->taxonomy == nullptr) {
		ctsp->taxonomy = std::move(tree);
		return ctsp->processTaxonomyTree(otCLI);
	}
	if (ctsp->toCheck == nullptr) {
		ctsp->toCheck = std::move(tree);
		return true;
	}
	return ctsp->processSourceTree(otCLI, std::move(tree));
}

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
	FindUnsupportedState fus;
	otCLI.blob = static_cast<void *>(&fus);
	otCLI.addFlag('m',
				  "ARG=a designators file. Each line is a list of (white-space separated) OTT ids used to designate the node that is the MRCA of them.",
				  handleDesignator,
				  true);
	auto rc = treeProcessingMain<RTSplits, RootedTreeForNodeType>(otCLI,
																	  argc,
																	  argv,
																	  processNextTree,
																	  nullptr);
	if (rc == 0) {
		fus.summarize(otCLI);
		return fus.numErrors;
	}
	return rc;
}
