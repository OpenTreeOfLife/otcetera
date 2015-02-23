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

	void summarize(const OTCLI &otCLI) {}

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
		for (auto nd : ConstLeafNodeIter<RTSplits, RootedTreeForNodeType>(tree)) {
			auto ottId = nd->getOttId();
			markPathToRoot(*toCheck, ottId, prunedDesId);
		}
		std::map<std::set<long>, const MyNodeType *> sourceClades;
		for (auto nd : ConstPostorderInternalNode<RTSplits, RootedTreeForNodeType>(tree)) {
			if (nd->getParent() != nullptr && !nd->isTip()) {
				sourceClades[nd->getData().desIds] = nd;
			}
		}
		recordSupportedNodes(otCLI, prunedDesId, sourceClades, supportedNodes);
		return true;
	}
	void recordSupportedNodes(OTCLI & otCLI,
							  const std::map<const MyNodeType *, std::set<long> > &prunedDesId,
							  const std::map<std::set<long>, const MyNodeType *> & sourceClades,
							  std::set<const MyNodeType *> & supported) {

		for (auto pd : prunedDesId) {
			auto nd = pd.first;
			auto par = nd->getParent();
			if (par == nullptr) {
				continue;
			}
			auto ls = pd.second;
			auto firstBranchingAnc = findFirstBranchingAnc<const MyNodeType>(nd);
			if (firstBranchingAnc == nullptr) {
				continue;
			}
			auto nm = pd.second;
			auto ancIt = prunedDesId.find(firstBranchingAnc);
			assert(ancIt != prunedDesId.end());
			auto anm = ancIt->second;
			const MyNodeType * firstNdPtr; // just used to match call
			if ((!multipleChildrenInMap(*nd, prunedDesId, &firstNdPtr)) || anm == nm) {
				continue;
			}
			auto scIt = sourceClades.find(nm);
			if (scIt != sourceClades.end()) {
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

#if 0


#include <fstream>
#include "ncl/othelpers.h"
OTCLI gOTCli("findunsupportededgs: takes a tree, an optional file (.txt extension) of MRCA designators, a taxonomy, and the some # of newick files. It reports any edges in the first tree that have no support in the last trees. If the designators file is sent, then the nodes that are the MRCAs of the numbers on each line of that file will be treated as a priori problems to be checked. So you will get a report for those nodes whether or not they are unsupported.",
			 "findunsupportededgs synth.tre probs.txt taxonomy.tre one.tre two.tre three.tre");
using namespace std;
void processContent(PublicNexusReader & nexusReader, ostream *out);

bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);

/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
NxsSimpleTree * gTaxonTree = 0;
map<long, const NxsSimpleNode *> gOttID2RefNode;
map<const NxsSimpleNode *, string> gRefTipToName;
map<long, const NxsSimpleNode *> gOttID2TaxNode;
map<const NxsSimpleNode *, long> gTaxNode2ottID;
set<const NxsSimpleNode *> gSupportedNodes;
string gCurrentFilename;
string gCurrTmpFilepath;
ostream * gCurrTmpOstream = 0L;
bool gReadingTxtFile = false;
map<long, set<long> > gNonMono;
const bool gTrustNamedNodes = true;
map<const NxsSimpleNode *, long> gExpanded;
map<long, const NxsSimpleNode *> gTabooLeaf;
set<long> gTaxLeafOTTIDs;
map<const NxsSimpleNode *, set<long> > gAPrioriProblemNode;
bool gNoAprioriTests = true;
int gRefTreeNumNamedInternalsNodes = 0;
int gExitCode = 0;


void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, set<const NxsSimpleNode *> & gSupportedNodes) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1) {
			if (gSupportedNodes.find(children[0]) != gSupportedNodes.end() || children[0]->GetName().length() > 0) {
				if (gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end()) {
					assert(false); // shouldn't get out-degree one nodes w/ our designators
				}
				gSupportedNodes.insert(nd);
			}
		}
	}
}

bool singleDesSupportedOrNamed(const NxsSimpleNode *nd) {
	if (gSupportedNodes.find(nd) != gSupportedNodes.end()) {
		return true;
	}
	if (nd->GetOutDegree() == 1) {
		if (!nd->GetName().empty()) {
			return true;
		} else {
			return singleDesSupportedOrNamed(nd->GetFirstChild());
		}
	}
	return false;
}

int describeUnnamedUnsupported(ostream &out, const NxsSimpleTree * tree, const set<const NxsSimpleNode *> & supported) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
	int numUnsupported = 0;
	++nIt; //skip the root
	for (;nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0 && supported.find(nd) == supported.end()) {
			if (IsRedundantNodeAroundTip(nd)) {
				//pass
			} else if (outDegree == 1 && gTrustNamedNodes && singleDesSupportedOrNamed(nd)) {
				//pass
			} else if (nd->GetName().length() == 0) { //assume that it is from the taxonomy
				if (gNoAprioriTests) {
					out << "Unsupported node ";
				} else {
					map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.find(nd);
					if (gaIt == gAPrioriProblemNode.end()) {
						out << "Novel unsupported node ";
					} else {
						out << "Confirmation of unsupported node (designators =";
						writeOttSet(out, "", gaIt->second, " ");
						out << ") ";
					}
				}
				describeUnnamedNode(nd, out, 0, gRefTipToName, false);
				numUnsupported += 1;
			}
		}
	}
	return numUnsupported;
}


void summarize(ostream & out) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	extendSupportedToRedundantNodes(gRefTree, gSupportedNodes);
	int numUnsupported = describeUnnamedUnsupported(out, gRefTree, gSupportedNodes);
	for (map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.begin(); gaIt != gAPrioriProblemNode.end(); ++gaIt) {
		if (gSupportedNodes.find(gaIt->first) != gSupportedNodes.end()) {
			out << "Claim of unsupported apparently refuted for designators: ";
			writeOttSet(out, "", gaIt->second, " ");
			out << ". See standard error stream for details.\n";
		}
	} 
	int supNonNamed = 0;
	int numSupportedInternals = 0;
	for (set<const NxsSimpleNode *>::const_iterator rIt = gSupportedNodes.begin(); rIt != gSupportedNodes.end(); ++rIt) {
		if (!(*rIt)->IsTip()) {
			numSupportedInternals += 1;
			if ((*rIt)->GetName().length() == 0) {
				supNonNamed += 1;
			}
		}
	}
	out << "\n\nFinal summary:\n";
	out << gRefTreeNumNamedInternalsNodes << " internal nodes were named in the reference tree. These were not rigorously checked against the taxonomy. They may not be detected as errors.\n";
	out << numSupportedInternals << " internal nodes where flagged as being supported by an input (including taxonomy).\n";
	int supNamed = numSupportedInternals - supNonNamed;
	out << "    " << supNamed << " of these were named (some of the support could just be the taxonomic expansion of tips).\n";
	out << "    " << supNonNamed << " of these were unnamed.\n";
	out << numUnsupported << " unsupported nodes.\n";
	out << endl;
	if (gExitCode < 0) {
		gExitCode -= numUnsupported;
	} else {
		gExitCode = numUnsupported;
	}
}






const NxsSimpleNode * findMRCA(const map<long, const NxsSimpleNode *> & ref,
	                           const map<long, const NxsSimpleNode *> &taxonomy, long ottID) {
	set<long> tipOTTIDs;
	fillTipOTTIDs(taxonomy, ottID, tipOTTIDs, gTaxNode2ottID);
	const unsigned nTips = tipOTTIDs.size();
	if (nTips < 2) {
		cerr << "findMRCA called on " << ottID << '\n';
		assert(false);
	}
	gNonMono[ottID] = tipOTTIDs;
	return findMRCAFromIDSet(ref, tipOTTIDs, ottID);
}






#endif
