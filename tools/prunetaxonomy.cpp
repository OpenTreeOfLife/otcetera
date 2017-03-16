#include "otc/otcli.h"
using namespace otc;

struct PruneTaxonomyState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    bool reportStats;
    int numErrors;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    std::set<const RootedTreeNodeNoData *> directlyIncludedNodes;
    virtual ~PruneTaxonomyState(){}

    PruneTaxonomyState()
        :reportStats(false),
        numErrors(0) {
    }

    bool summarize(OTCLI &otCLI) override {
        if (reportStats) {
            OttIdSet oids;
            OttIdSet ntoids;

            std::size_t numNonTerminals = 0;
            for (auto tn : directlyIncludedNodes) {
                if (!tn->is_tip()) {
                    numNonTerminals++;
                    ntoids.insert(tn->get_ott_id());
                }
                oids.insert(tn->get_ott_id());
            }
            otCLI.out << numNonTerminals << " non-terminal taxa in OTT that are mapped by at least 1 input.\n";
            otCLI.out << (directlyIncludedNodes.size() - numNonTerminals) << " terminal taxa in OTT that are mapped by at least 1 input\n";
            otCLI.out << directlyIncludedNodes.size() << " total taxa in OTT that are mapped by at least 1 input.\n";
            if (otCLI.verbose) {
                otCLI.out << "total included OTT Ids\n";
                for (const auto & oid : oids) {
                    otCLI.out << oid << '\n';
                }
                otCLI.err << "non-terminal OTT Ids\n";
                for (const auto & oid : ntoids) {
                    otCLI.err << oid << '\n';
                }
            }
            return true;
        }
        assert(taxonomy != nullptr && !includedNodes.empty());
        std::set<RootedTreeNodeNoData *> toPrune;
        std::size_t numLeavesPruned = 0;
        std::size_t numInternalsPruned = 0;
        for (auto nd : iter_node(*taxonomy)) {
            const RootedTreeNodeNoData *  c = const_cast<const RootedTreeNodeNoData *>(nd);
            if (!contains(includedNodes, c)) {
                if (contains(includedNodes, c->get_parent())) {
                    toPrune.insert(nd);
                }
                if (c->is_tip()) {
                    numLeavesPruned += 1;
                } else {
                    numInternalsPruned += 1;
                }
            }
        }
        for (auto nd : toPrune) {
            prune_and_delete(*taxonomy, nd);
        }
        write_tree_as_newick(otCLI.out, *taxonomy);
        otCLI.out << '\n';
        otCLI.err << numLeavesPruned << " terminal taxa pruned\n";
        otCLI.err << numInternalsPruned << " non-terminal taxa pruned\n";
        return true;
    }

    bool process_source_tree(OTCLI & , const std::unique_ptr<TreeMappedEmptyNodes> treePtr) override {
        assert(taxonomy != nullptr);
        std::map<const RootedTreeNodeNoData *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(*treePtr)) {
            auto ottId = nd->get_ott_id();
            auto taxoNode = taxonomy->get_data().get_node_by_ott_id(ottId);
            assert(taxoNode != nullptr);
            if (reportStats) {
                directlyIncludedNodes.insert(taxoNode);
            } else {
                if (!contains(includedNodes, taxoNode)) {
                    includedNodes.insert(taxoNode);
                    insert_ancestors_to_paraphyletic_set(taxoNode, includedNodes);
                }
                insert_descendants_of_unincluded_subtrees(taxoNode, includedNodes);
            }
        }
        return true;
    }
};

bool handleReport(OTCLI & otCLI, const std::string &);

bool handleReport(OTCLI & otCLI, const std::string &) {
    PruneTaxonomyState * proc = static_cast<PruneTaxonomyState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->reportStats  = true;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-prune-taxonomy",
                "takes at least 2 newick file paths: a full taxonomy tree some number of input trees. Prune subtrees from the taxonomy if they are not represented in the inputs",
                "taxonomy.tre inp1.tre inp2.tre");
    PruneTaxonomyState proc;
    otCLI.add_flag('r',
                  "Just report stats on how many tips are included in the inputs.",
                  handleReport,
                  false);
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, false);
}
