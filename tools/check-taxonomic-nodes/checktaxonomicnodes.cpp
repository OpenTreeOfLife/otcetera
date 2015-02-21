#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;
typedef RTreeOttIDMapping<RTNodeNoData> RootedTreeForNodeType;
typedef RootedTree<RTNodeNoData, RootedTreeForNodeType> Tree_t;

bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

struct CheckTaxonState {
	std::unique_ptr<Tree_t> toCheck;
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;

	CheckTaxonState()
		:toCheck(nullptr),
		 taxonomy(nullptr),
		 numErrors(0) {
		}
	void summarize(const OTCLI &otCLI) {

	}
};



inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
	CheckTaxonState * ctsp = static_cast<CheckTaxonState *>(otCLI.blob);
	assert(ctsp != nullptr);
	assert(tree != nullptr);
	if (ctsp->toCheck == nullptr) {
		ctsp->toCheck = std::move(tree);
	} else if (ctsp->taxonomy == nullptr) {
		ctsp->taxonomy = std::move(tree);
	} else {
		otCLI.err << "Expecting only 2 trees a tree to check and the taxonomy!\n";
		return false;
	}
	return true;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcchecktaxonomicnodes",
				"takes 2 newick file paths: to a tree with internal node names and a tree form of a taxonomy. Reports on any taxonomic names that do not map to the correct place in the first tree.",
				"some.tre taxonomy.tre");
	CheckTaxonState cts;
	otCLI.blob = static_cast<void *>(&cts);
	auto rc = treeProcessingMain<RTNodeNoData, RootedTreeForNodeType>(otCLI,
																	  argc,
																	  argv,
																	  processNextTree,
																	  nullptr);
	if (rc == 0) {
		cts.summarize(otCLI);
		return cts.numErrors;
	}
	return rc;
}

#if 0
/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
NxsSimpleTree * gTaxonTree = 0;
std::map<long, const NxsSimpleNode *> gOttID2RefNode;
std::map<long, const NxsSimpleNode *> gOttID2TaxNode;
std::map<const NxsSimpleNode *, long> gTaxNode2ottID;
std::set<const NxsSimpleNode *> gSupportedNodes;

std::map<const NxsSimpleNode *, std::set<long> > gRefNdp2mrca;
std::map<const NxsSimpleNode *, std::set<long> > gTaxNdp2mrca;
set<long> gRefLeafSet;
set<long> gTaxLeafSet;
map<const NxsSimpleNode *, long> gRefNamedNodes;


void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, *nd);
		if (nd->IsTip()) {
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const std::string tn = tb->GetTaxonLabel(ind);
			assert(ottID >= 0);
			gRefNdp2mrca[nd].insert(ottID);
			gRefLeafSet.insert(ottID);
		} else {
			useChildrenToFillMRCASet(nd, gRefNdp2mrca);
			if (ottID > 0) {
				gRefNamedNodes[nd] = ottID;
			}
		}
		if (ottID >= 0) {
			assert(gOttID2RefNode.find(ottID) == gOttID2RefNode.end());
			gOttID2RefNode[ottID] = nd;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafSet.insert(ottID);
			gTaxNdp2mrca[nd].insert(ottID);
		} else {
			useChildrenToFillMRCASet(nd, gTaxNdp2mrca);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (std::map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}

bool isProperSubset(const set<long> & small, const set<long> & big) {
	if (big.size() <= small.size()) {
		return false;
	}
	for (set<long>::const_iterator rIt = small.begin(); rIt != small.end(); ++rIt) {
		if (big.find(*rIt) == big.end()) {
			return false;
		}
	}
	return true;
}

bool doCheckEquivalent(std::ostream &out, long ottID, const NxsSimpleNode * snode, std::map<const NxsSimpleNode *, std::set<long> > & srcLookup,
										  const NxsSimpleNode * tnode, std::map<const NxsSimpleNode *, std::set<long> > & taxLookup,
										  bool topLevel, bool climbSynth, bool climbTax) {
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator streeLSIt = srcLookup.find(snode);
	assert(streeLSIt != srcLookup.end());
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator taxtreeLSIt = taxLookup.find(tnode);
	assert(taxtreeLSIt != taxLookup.end());
	const std::set<long> & streeMRCA = streeLSIt->second;
	const std::set<long> & taxtreeMRCA = taxtreeLSIt->second;
	if (streeMRCA != taxtreeMRCA) {
		if (topLevel) {
			out << "ottID " << ottID << " incorrect:\n";
			writeOttSetDiff(out, "    ", streeMRCA, "synth", taxtreeMRCA, "taxonomy");
		}
		if (climbSynth && isProperSubset(streeMRCA, taxtreeMRCA)) {
			return doCheckEquivalent(out, ottID, snode->GetEdgeToParent().GetParent(), srcLookup, tnode, taxLookup, false, true, false);
		} else if (climbTax && isProperSubset(taxtreeMRCA, streeMRCA)) {
			return doCheckEquivalent(out, ottID, snode, srcLookup, tnode->GetEdgeToParent().GetParent(), taxLookup, false, false, true);
		} else {
			return false;
		}
	} else if (!topLevel) {
		out << "        Found identical leaf sets for the synthetic tree \"" << snode->GetName() << "\" and the taxonomic node \"" << tnode->GetName() << "\".\n";
	}
	return true;
}

void summarize(std::ostream & out) {
	for (map<const NxsSimpleNode *, long>::const_iterator rnit = gRefNamedNodes.begin(); rnit != gRefNamedNodes.end(); ++rnit) {
		const NxsSimpleNode * nd = rnit->first;
		const long ottID = rnit->second;
		std::map<long, const NxsSimpleNode *>::const_iterator tID2nd = gOttID2TaxNode.find(ottID);
		assert(tID2nd != gOttID2TaxNode.end());
		const NxsSimpleNode *taxNd = tID2nd->second;
		if (!doCheckEquivalent(out, ottID, nd, gRefNdp2mrca, taxNd, gTaxNdp2mrca, true, true, true)) {
			out << "        Could not find this set of leaves in the synth \"" << nd->GetName() <<"\" in any taxonomic node.\n";
		}
	}
	if (gTaxLeafSet != gRefLeafSet) {
		writeOttSetDiff(out, "", gRefLeafSet, "synth", gTaxLeafSet, "taxonomy");
	}
}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB) {
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	if (gRefTree == 0) {
		gRefTree = nst;
		processRefTree(taxa, nst);
	} else if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		const char * msg = "Exepting only 2 files: the tree file, and then the taxonomy\n";
		std::cerr << msg;
		throw NxsException(msg);
	}
	if (gRefTree != nst && gTaxonTree != nst) {
		delete nst;
	}
	return false;
}

int main(int argc, char *argv[]) {
	std::vector<std::string> args;
	if (!gOTCli.parseArgs(argc, argv, args)) {
		return 1;
	}
	if (args.size() != 2) {
		cerr << "Expecting a tree file and taxonomy tree.\n";
		return 1;
	}
	try{
		gOTCli.readFilepath(args[0], newTreeHook);
		gOTCli.readFilepath(args[1], newTreeHook);
	} catch (NxsException &x) {
		std::cerr << x.what() << "\n";
		return 1;
	}
	if (gOTCli.exitCode == 0) {
		summarize(std::cout);
	}
	return gOTCli.exitCode;
}

#endif
