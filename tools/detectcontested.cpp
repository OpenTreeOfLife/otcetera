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

    bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        return true;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        assert(tree != nullptr);
        expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
        return processExpandedTree(otCLI, *tree);
    }

    bool processExpandedTree(OTCLI &otCLI, TreeMappedWithSplits & tree) {
        std::map<const NodeWithSplits *, std::set<long> > prunedDesId;
        for (auto nd : iter_leaf_const(tree)) {
            auto ottId = nd->getOttId();
            markPathToRoot(*taxonomy, ottId, prunedDesId);
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
            if (nd->getParent() != nullptr && !nd->isTip()) {
                sourceClades.insert(std::move(nd->get_data().desIds));
            }
        }
        auto numLeaves = tree.getRoot()->get_data().desIds.size();
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
                if (!areCompatibleDesIdSets(taxNodesDesSets, sc)) {
                    for (auto nd : ndlist) {
                        contestedSet.insert(nd);
                        std::cerr << getContestedPreambleFromName(*nd, treeName) << '\n';
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
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, true);
}
