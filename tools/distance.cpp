#include "otc/otcli.h"
using namespace otc;

// Note that the "taxonomy" data member here will be the first tree (the supertree)
struct DistanceState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    unsigned long totalRF;
    unsigned long totalNumNotDisplayed;
    unsigned long totalNumInternals;
    unsigned long totalNumDisplayed;
    bool showRF;
    bool showNumNotDisplayed;
    bool showNumInternals;
    bool showNumDisplayed;
    std::string prevTreeFilename;
    std::size_t numComparisons;
    std::size_t numTreesInThisTreefile;


    virtual ~DistanceState(){}
    DistanceState()
        :totalRF(0U),
        totalNumNotDisplayed(0U),
        totalNumInternals(0U),
        totalNumDisplayed(0U),
        showRF(false),
        showNumNotDisplayed(false),
        showNumInternals(false),
        showNumDisplayed(false),
        numComparisons(0U), 
        numTreesInThisTreefile(0U) {
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) {
        numComparisons += 1;
        std::string nameToPrint = otCLI.currentFilename;
        if (nameToPrint == prevTreeFilename) {
            numTreesInThisTreefile += 1;
            nameToPrint += "-tree#" + std::to_string(numTreesInThisTreefile);
        } else {
            prevTreeFilename = nameToPrint;
            numTreesInThisTreefile = 1;
        }
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        std::set<std::set<long> > inducedSplits;
        std::set<std::set<long> > tree2Splits;
        inducedCladeSets(*taxonomy, *tree, inducedSplits, tree2Splits, true);
        unsigned long rf = 0;
        unsigned long numNotDisplayed = 0;
        const unsigned long numInternals = tree2Splits.size();
        unsigned long numDisplayed = 0;
        totalNumInternals += numInternals;
        if (showRF) {
            rf = sizeOfSymmetricDifference(tree2Splits, inducedSplits);
            totalRF += rf;
        }
        if (showNumDisplayed || showNumNotDisplayed) {
            for (auto ics : tree2Splits) {
                if (contains(inducedSplits, ics)) {
                    numDisplayed += 1;
                } else {
                    numNotDisplayed += 1;
                }
            }
            totalNumNotDisplayed += numNotDisplayed;
        }
        // display
        // header if first comparison
        if (numComparisons == 1) {
            otCLI.out << "treename";
            if (showRF) {
                otCLI.out << "\tRF";
            }
            if (showNumNotDisplayed) {
                otCLI.out << "\tNumNotDisplayed";
            }
            if (showNumDisplayed) {
                otCLI.out << "\tNumDisplayed";
            }
            if (showNumNotDisplayed) {
                otCLI.out << "\tNumInternals";
            }
            otCLI.out << '\n';
        }
        otCLI.out << nameToPrint;
        if (showRF) {
            otCLI.out << '\t' << rf;
        }
        if (showNumNotDisplayed) {
            otCLI.out << '\t' << numNotDisplayed;
        }
        if (showNumDisplayed) {
            otCLI.out << '\t' << numDisplayed;
        }
        if (showNumNotDisplayed) {
            otCLI.out << '\t' << numInternals;
        }
        otCLI.out << '\n';
        return true;
    }
    bool summarize(OTCLI & otCLI) override {
        otCLI.out << "TOTALS";
        if (showRF) {
            otCLI.out << '\t' << totalRF;
        }
        if (showNumNotDisplayed) {
            otCLI.out << '\t' << totalNumNotDisplayed;
        }
        if (showNumDisplayed) {
            otCLI.out << '\t' << totalNumDisplayed;
        }
        if (showNumNotDisplayed) {
            otCLI.out << '\t' << totalNumInternals;
        }
        otCLI.out << '\n';
        return true;
    }

};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-distance",
                "takes at least 2 newick file paths: a supertree and some number of input trees. Writes one line for each input tree with the statistics requested for the comparison of the supertree to each input tree.",
                "synth.tre inp1.tre inp2.tre");
    DistanceState proc;
    auto rc = taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
    return rc;
}
