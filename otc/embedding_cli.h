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

/// EmbeddingCLI is intended as a base class for tools that need to 
//      produce an EmbeddedTree with tree 0 serving as the 
//      "outer" tree and all subsequent trees to be embedded in a taxonomic tree.
// It provides overrides for the TaxonomyDependentTreeProcessor hooks process_taxonomy_tree,
//  and process_source_tree. It also proveds a clone_taxonomy_as_a_source_tree method that allows
//  the taxonomy tree to be embedded inside itself (this feature is used by the uncontested-decompose
//  so that it can correctly emit the fracments of the input trees (including the taxonomic tree)
//  when it writes a subproblem).
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

    bool process_taxonomy_tree(OTCLI & otCLI) override {
        debuggingOutput = otCLI.verbose;
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::process_taxonomy_tree(otCLI);
        //check_tree_invariants(*taxonomy);
        suppressMonotypicTaxaPreserveDeepestDangle(*taxonomy, false);
        monotypicRemapping = generateIdRemapping(*taxonomy);
        //check_tree_invariants(*taxonomy);
        for (NodeWithSplits * nd : iter_node(*taxonomy)) {
            _get_embedding_for_node(nd); // side effect is introducint a new, empty embedding
        }
        otCLI.get_parsing_rules().set_ott_idForInternals = false;
        otCLI.get_parsing_rules().id_remapping = &monotypicRemapping;
        return true;
    }

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
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
        suppressMonotypicTaxaPreserveShallowDangle(*raw);
        embed_new_tree(*taxonomy, *raw, treeIndex);
        otCLI.err << "# pathPairings = " << pathPairings.size() << '\n';
        return true;
    }

    bool clone_taxonomy_as_a_source_tree() {
        assert(taxonomy != nullptr);
        assert(taxonomyAsSource == nullptr);
        std::unique_ptr<TreeMappedWithSplits> tree = cloneTree(*taxonomy);
        taxonomyAsSource = tree.get();
        std::size_t treeIndex = inputTreesToIndex.size();
        inputTreesToIndex[std::move(tree)] = treeIndex;
        treePtrByIndex.push_back(taxonomyAsSource);
        // suppress the internal node OTT IDs from the des
        OttIdSet internalIDs;
        for (auto nd : iter_post_internal(*taxonomyAsSource)) {
            if (nd->has_ott_id()) {
                internalIDs.insert(nd->get_ott_id());
            }
            auto & d = nd->get_data().desIds;
            for (auto o : internalIDs) {
                d.erase(o);
            }
        }
        // Store the tree's filename
        taxonomyAsSource->setName("TAXONOMY");
        embed_scaffold_clone(*taxonomy, *taxonomyAsSource, treeIndex);
        return true;
    }
};

} // namespace otc
#endif
