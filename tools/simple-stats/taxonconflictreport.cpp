#include "otc/otcli.h"
using namespace otc;
struct ConflictReporterProc : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
	virtual ~ConflictReporterProc(){}
	bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) {
		assert(tree != nullptr);
		assert(taxonomy != nullptr);
		reportOnInducedConflicts(otCLI.out, *taxonomy, *tree, true);
		return true;
	}
};

int main(int argc, char *argv[]) {
	OTCLI otCLI("otctaxonconflictreport",
				"takes at least 2 newick file paths: a full tree, and some number of input trees. Writes a summary of the difference in taxonomic inclusion for nodes that are in conflict.",
				"taxonomy.tre inp1.tre inp2.tre");
	ConflictReporterProc proc;
	return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);

}
