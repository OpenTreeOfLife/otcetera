#include "otc/otcli.h"
using namespace otc;

struct DetectContestedState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    int numErrors;
    std::set<const NodeWithSplits *> contestedNodes;
    bool doShortcircuit;
    virtual ~DetectContestedState(){}

    DetectContestedState()
        :numErrors(0),
        doShortcircuit(false) {
        }

    bool summarize(OTCLI &otCLI) override {
        assert (taxonomy != nullptr);
        for (auto nd : contestedNodes) {
            otCLI.out << nd->get_name() << '\n';
        }
        return true;
    }

    bool process_taxonomy_tree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::process_taxonomy_tree(otCLI);
        return true;
    }

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        assert(tree != nullptr);
        expand_ott_internals_which_are_leaves(*tree, *taxonomy);
        return processExpandedTree(otCLI, *tree);
    }

    bool processExpandedTree(OTCLI &otCLI, TreeMappedWithSplits & tree) {
        std::map<const NodeWithSplits *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->get_ott_id();
            mark_path_to_root(*taxonomy, ottId, prunedDesId);
        }
        std::map<std::set<long>, std::list<const NodeWithSplits *> > taxCladesToTaxNdList;
        for (auto & nodeSetPair : prunedDesId) {
            auto nd = nodeSetPair.first;
            auto & ds = nodeSetPair.second;
            if (ds.size() < 2) {
                continue;
            }
            auto tctnlIt = taxCladesToTaxNdList.find(ds);
            if (tctnlIt == taxCladesToTaxNdList.end()) {
                std::list<const NodeWithSplits *> sel{1, nd};
                taxCladesToTaxNdList.emplace(ds, sel);
            } else {
                tctnlIt->second.push_back(nd);
            }
        }
        std::set<std::set<long> > sourceClades;
        for (auto nd : iter_post_internal(tree)) {
            if (nd->get_parent() != nullptr && !nd->is_tip()) {
                sourceClades.insert(std::move(nd->get_data().des_ids));
            }
        }
        auto numLeaves = tree.get_root()->get_data().des_ids.size();
        recordContested(taxCladesToTaxNdList, sourceClades, contestedNodes, numLeaves, otCLI.currentFilename);
        return true;
    }

    void recordContested(const std::map<std::set<long>, std::list<const NodeWithSplits *> > & prunedDesId,
                         const std::set<std::set<long> > & sourceClades,
                         std::set<const NodeWithSplits *> & contestedSet,
                         std::size_t numLeaves,
                         const std::string &treeName) {
        for (const auto & pd : prunedDesId) {
            // shortcircuit taxon nodes that are already marked as contested
            const auto & ndlist = pd.second;
            bool allContested = true;
            for (auto nd : ndlist) {
                if (!contains(contestedSet, nd)) {
                    allContested = false;
                    break;
                }
            }
            if (doShortcircuit && allContested) {
                continue;
            }
            const auto & taxNodesDesSets = pd.first;
            const auto nss = taxNodesDesSets.size();
            if (nss == 1 || nss == numLeaves) {
                continue;
            }
            for (const auto & sc : sourceClades) {
                if (!are_compatible_des_id_sets(taxNodesDesSets, sc)) {
                    for (auto nd : ndlist) {
                        contestedSet.insert(nd);
                        std::cerr << get_contested_preamble_from_name(*nd, treeName) << '\n';
                    }
                    break;
                }
            }
        }
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-detect-contested",
                "takes at least 2 newick file paths: a full taxonomy tree, and some number of input trees. Writes the OTT IDs of clades in the taxonomy whose monophyly is questioned by at least one input",
                "taxonomy.tre inp1.tre inp2.tre");
    DetectContestedState proc;
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, true);
}
