#include <set>
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_data.h"
using namespace otc;

typedef otc::RootedTreeNode<RTSplits> TaxoCheckNodeType;
typedef RTreeOttIDMapping<RTSplits> RootedTreeForNodeType;
typedef otc::RootedTree<RTSplits, RootedTreeForNodeType> Tree_t;

bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
extern const char * badNTreesMessage;

const char * badNTreesMessage = "Expecting only 2 trees a tree to check and the taxonomy!\n";

struct CheckTaxonState {
	std::unique_ptr<Tree_t> toCheck;
	std::unique_ptr<Tree_t> taxonomy;
	int numErrors;

	CheckTaxonState()
		:toCheck(nullptr),
		 taxonomy(nullptr),
		 numErrors(0) {
		}
	bool doCheckEquivalent(std::ostream &out,
						   long ottID,
						   const TaxoCheckNodeType * snode,
						   const TaxoCheckNodeType * tnode,
						   bool topLevel,
						   bool climbSynth,
						   bool climbTax) {
		const std::set<long> & streeDes = snode->getData().desIds;
		const std::set<long> & taxTreeDes = tnode->getData().desIds;
		if (streeDes != taxTreeDes) {
			if (topLevel) {
				numErrors += 1;
				out << "ottID " << ottID << " incorrect:\n";
				writeOttSetDiff(out, "    ", streeDes, "toCheck", taxTreeDes, "taxonomy");
			}
			if (climbSynth && isProperSubset(streeDes, taxTreeDes)) {
				return doCheckEquivalent(out, ottID, snode->getParent(), tnode, false, true, false);
			} else if (climbTax && isProperSubset(taxTreeDes, streeDes)) {
				return doCheckEquivalent(out, ottID, snode, tnode->getParent(), false, false, true);
			} else {
				return false;
			}
		} else if (!topLevel) {
			out << "        Found identical leaf sets for the synthetic tree \"" << snode->getName() << "\" and the taxonomic node \"" << tnode->getName() << "\".\n";
		}
		return true;
	}


	void summarize(const OTCLI &otCLI) {
		if (taxonomy == nullptr) {
			numErrors = 1;
			otCLI.err << badNTreesMessage;
			return;
		}
		auto taxOttIdToNode = taxonomy->getData().ottIdToNode;
		const std::set<long> taxOttIds = taxonomy->getRoot()->getData().desIds;
		const std::set<long> toCheckOttIds = toCheck->getRoot()->getData().desIds;
		auto extras = set_difference_as_set(toCheckOttIds, taxOttIds);
		if (!extras.empty()) {
			otCLI.err << "OTT Ids found in the tree to check but not in the taxonomy:\n";
			writeOttSet(otCLI.err, "  ", extras, "\n");
			numErrors = (extras.size() > INT_MAX ? INT_MAX : (int) extras.size());
			return;
		}
		// now check for taxonomic identity
		for (auto toCheckNd: ConstPostorderInternalNode<RTSplits, RootedTreeForNodeType>(*toCheck)) {
			if (!toCheckNd->hasOttId()) {
				continue;
			}
			auto toCheckId = toCheckNd->getOttId();
			assert(contains(taxOttIdToNode, toCheckId));
			auto taxNd = taxOttIdToNode.find(toCheckId)->second;
			if (!doCheckEquivalent(otCLI.out, toCheckId, toCheckNd, taxNd, true, true, true)) {
				otCLI.out << "        Could not find this set of leaves named \"" << toCheckNd->getName() <<"\" in the tree to check when searching through the taxonomic node.\n";
			}
		}
		if (toCheckOttIds != taxOttIds) {
			writeOttSetDiff(otCLI.out, "", toCheckOttIds, "toCheck", taxOttIds, "taxonomy");
			auto omitted = set_difference_as_set(taxOttIds, toCheckOttIds);
			otCLI.err << "OTT Ids found in the taxonomy but not the tree toCheck:\n";
			writeOttSet(otCLI.err, "  ", omitted, "\n");
			for (auto oid: omitted) {
				otCLI.err << "id = " << oid << " name = \"" << taxOttIdToNode.find(oid)->second->getName() << "\"\n";
			}
			return;
		}
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
		otCLI.err << badNTreesMessage;
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
	auto rc = treeProcessingMain<RTSplits, RootedTreeForNodeType>(otCLI,
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



void summarize(std::ostream & out) {
	for (map<const NxsSimpleNode *, long>::const_iterator rnit = gRefNamedNodes.begin(); rnit != gRefNamedNodes.end(); ++rnit) {
		const NxsSimpleNode * nd = rnit->first;
		const long ottID = rnit->second;
		std::map<long, const NxsSimpleNode *>::const_iterator tID2nd = gOttID2TaxNode.find(ottID);
		assert(tID2nd != gOttID2TaxNode.end());
		const NxsSimpleNode *taxNd = tID2nd->second;
		if (!doCheckEquivalent(out, ottID, nd, gRefNdp2mrca, taxNd, gTaxNdp2mrca, true, true, true)) {
			out << "        Could not find this set of leaves in the synth \"" << nd->getName() <<"\" in any taxonomic node.\n";
		}
	}
	if (gTaxLeafSet != gRefLeafSet) {
		writeOttSetDiff(out, "", gRefLeafSet, "synth", gTaxLeafSet, "taxonomy");
	}
}


#endif
