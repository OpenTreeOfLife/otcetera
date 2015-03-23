#ifndef OTCETERA_EMBEDDING_CLI_H
#define OTCETERA_EMBEDDING_CLI_H
// base class for OTCLI tools that create embeddings of the taxonomy.
#include <tuple>
#include <map>
#include <set>
#include "otc/otc_base_includes.h"
#include "otc/supertree_util.h"
#include "otc/embedded_tree.h"
#include "otc/greedy_forest.h"
#include "otc/node_embedding.h"
#include "otc/tree_iter.h"
#include "otc/tree_data.h"
#include "otc/debug.h"
#include "otc/otcli.h"
#include "otc/util.h"
namespace otc {

class EmbeddingCLI
    : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits>,
    public EmbeddedTree {
    public:
    int numErrors;
    std::map<std::unique_ptr<TreeMappedWithSplits>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedWithSplits *> treePtrByIndex;
    TreeMappedWithSplits * taxonomyAsSource;
    bool debuggingOutput;
    std::map<long, long> monotypicRemapping;

    virtual ~EmbeddingCLI(){}
    EmbeddingCLI()
        :TaxonomyDependentTreeProcessor<TreeMappedWithSplits>(),
         numErrors(0),
         taxonomyAsSource(nullptr),
         debuggingOutput(false) {
    }

    bool processTaxonomyTree(OTCLI & otCLI) override {
        debuggingOutput = otCLI.verbose;
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        checkTreeInvariants(*taxonomy);
        suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy);
        monotypicRemapping = generateIdRemapping(*taxonomy);
        checkTreeInvariants(*taxonomy);
        for (NodeWithSplits * nd : iter_node(*taxonomy)) {
            _getEmdeddingForNode(nd);
        }
        otCLI.getParsingRules().setOttIdForInternals = false;
        otCLI.getParsingRules().idRemapping = &monotypicRemapping;
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
        embedNewTree(*taxonomy, *raw, treeIndex);
        otCLI.err << "# pathPairings = " << pathPairings.size() << '\n';
        return true;
    }

    bool cloneTaxonomyAsASourceTree() {
        assert(taxonomy != nullptr);
        assert(taxonomyAsSource == nullptr);
        std::unique_ptr<TreeMappedWithSplits> tree = std::move(cloneTree(*taxonomy));
        taxonomyAsSource = tree.get();
        std::size_t treeIndex = inputTreesToIndex.size();
        inputTreesToIndex[std::move(tree)] = treeIndex;
        treePtrByIndex.push_back(taxonomyAsSource);
        // suppress the internal node OTT IDs from the des
        OttIdSet internalIDs;
        for (auto nd : iter_post_internal(*taxonomyAsSource)) {
            if (nd->hasOttId()) {
                internalIDs.insert(nd->getOttId());
            }
            auto & d = nd->getData().desIds;
            for (auto o : internalIDs) {
                d.erase(o);
            }
        }
        // Store the tree's filename
        taxonomyAsSource->setName("TAXONOMY");
        embedScaffoldClone(*taxonomy, *taxonomyAsSource, treeIndex);
        return true;
    }
};

} // namespace otc
#endif
