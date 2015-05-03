#include "otc/otcli.h"
#include "otc/supertree_util.h"
using namespace otc;

struct FindResolutionState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> toCheck;
    std::size_t numIncludable;
    int numErrors;
    virtual ~FindResolutionState(){}
    FindResolutionState()
        :toCheck(nullptr),
        numIncludable(0U),
        numErrors(0) {
    }

    bool summarize(OTCLI &otCLI) override {
        otCLI.out << numIncludable << " clades found which could be added to the tree.\n";
        numErrors = static_cast<int>(numIncludable);
        return numErrors == 0;
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = true;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set
        clearAndfillDesIdSets(*taxonomy);
        return true;
    }
    
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (toCheck == nullptr) {
            toCheck = std::move(tree);
            otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
            return true;
        }
        const OttIdSet & treeLeafSet = tree->getRoot()->getData().desIds;
        for (auto nd : iter_pre_internal_const(*tree)) {
            const OttIdSet & incGroup = nd->getData().desIds;
            auto mrca = findMRCAUsingDesIds(*toCheck, incGroup);
            OttIdSet md = mrca->getData().desIds;
            OttIdSet rmd;
            for (auto i : md) {
                if (contains(treeLeafSet, i)) {
                    rmd.insert(i);
                }
            }
            if (incGroup == rmd) {
                continue;
            }
            OttIdSet excGroup = treeLeafSet;
            for (auto i :incGroup) {
                excGroup.erase(i);
            }
            if (canBeResolvedToDisplayExcGroup(mrca, incGroup, excGroup)) {
                numIncludable += 1;
                otCLI.out << otCLI.currentFilename << " node " << getDesignator(*nd) << " could be added.\n";
            }
        }
        return true;
    }

};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-find-unsupported-nodes",
                "takes at least 2 newick file paths: a full supertree, and some number of input trees",
                "synth.tre inp1.tre inp2.tre ...");
    FindResolutionState proc;
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
