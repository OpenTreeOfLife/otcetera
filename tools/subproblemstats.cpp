#include "otc/otcli.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
struct PerTreeStats {
    OttIdSet leaf_set; // full leaf set
    OttIdSet leavesInIngroupSet; // leaves whose parents are not the root
    std::size_t numInformativeGroups;
};
struct PerSubproblemStats {
    std::list<PerTreeStats> pts;
};
struct SubproblemStatsState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    int numErrors;
    std::map<std::string, PerSubproblemStats> name2Stats;

};

// global - will inhibit parallelization, but this is not time-critical...
static SubproblemStatsState * gSubproblemStatsState = nullptr;
int summarize(OTCLI & otCLI);
template<typename T>
bool recordSubproblemStats(OTCLI & otCLI, std::unique_ptr<T> tree) {
    auto & pss = gSubproblemStatsState->name2Stats[otCLI.currentFilename];
    pss.pts.push_back(PerTreeStats());
    auto & pts = *(pss.pts.rbegin());
    auto root = tree->get_root();
    for (auto nd : iter_node_const(*tree)) {
        if (nd == root) {
            continue;
        }
        if (nd->is_tip()) {
            assert(nd->has_ott_id());
            const auto ottId = nd->get_ott_id();
            pts.leaf_set.insert(ottId);
            if (nd->get_parent() != root) {
                pts.leavesInIngroupSet.insert(ottId);
            }
        } else {
            pts.numInformativeGroups += 1;
        }
    }
    return true;
}

int summarize(OTCLI & otCLI) {
    auto & nts = gSubproblemStatsState->name2Stats;
    otCLI.err << nts.size() << " subproblems\n";
    otCLI.out << "Subproblem" ;
    otCLI.out << '\t' << "InSp";
    otCLI.out << '\t' << "LSS";
    otCLI.out << '\t' << "ILSS";
    otCLI.out << '\t' << "NT";
    otCLI.out << '\t' << "TreeSummaryName";
    otCLI.out << '\n';
    
    for (auto ntsP : nts) {
        const auto & name = ntsP.first;
        const auto & pss = ntsP.second;
        std::size_t i = 0;
        OttIdSet leaf_set;
        OttIdSet leavesInIngroupSet;
        std::size_t numInformativeGroups = 0;
        OttIdSet allButTaxLeafSet;
        OttIdSet allButTaxLeavesInIngroupSet;
        std::size_t allButTaxNumInformativeGroups = 0;

        for (const auto & pts : pss.pts) {
            allButTaxLeafSet = leaf_set;
            allButTaxLeavesInIngroupSet = leavesInIngroupSet;
            allButTaxNumInformativeGroups = numInformativeGroups;
            otCLI.out << name ;
            otCLI.out << '\t' << pts.numInformativeGroups;
            otCLI.out << '\t' << pts.leaf_set.size();
            otCLI.out << '\t' << pts.leavesInIngroupSet.size();
            otCLI.out << '\t'; // number of trees
            otCLI.out << '\t' << i++;
            otCLI.out << '\n';
            leaf_set.insert(begin(pts.leaf_set), end(pts.leaf_set));
            leavesInIngroupSet.insert(begin(pts.leavesInIngroupSet), end(pts.leavesInIngroupSet));
            numInformativeGroups += pts.numInformativeGroups;
        }
        otCLI.out << name ;
        otCLI.out << '\t' << allButTaxNumInformativeGroups;
        otCLI.out << '\t' << allButTaxLeafSet.size();
        otCLI.out << '\t' << allButTaxLeavesInIngroupSet.size();
        otCLI.out << '\t' << pss.pts.size() - 1;
        otCLI.out << '\t' << "Phylo-only";
        otCLI.out << '\n';
        otCLI.out << name ;
        otCLI.out << '\t' << numInformativeGroups;
        otCLI.out << '\t' << leaf_set.size();
        otCLI.out << '\t' << leavesInIngroupSet.size();
        otCLI.out << '\t' << pss.pts.size();
        otCLI.out << '\t' << "Total";
        otCLI.out << '\n';
   }
    return 0;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-subproblem-stats",
                "takes at series of tree files. Each is treated as a subproblem.\n" \
                "For each subproblem there is tab-separated report for:\n" \
                "    Subproblem name\n" \
                "    InSp = # of informative (nontrivial) splits\n" \
                "    LSS = size of the leaf label set\n" \
                "    ILSS = size of the set of labels included in at least one \"ingroup\"\n" \
                "    NT = The number of trees.\n" \
                "    TreeSummaryName = tree index or summary name\n" \
                "where the summary name can be Phylo-only or Total. \n" \
                "  \"Total\" summarizes info all trees in the file (including the taxonomy).\n" \
                "  \"Phylo-only\" former summarizes all of the phylogenetic inputs",
                "taxonomy.tre supertree/step_7_scratch/export-sub-temp/ott1000236.tre supertree/step_7_scratch/export-sub-temp/ott1000250.tre");
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> rss = recordSubproblemStats<Tree_t>;
    SubproblemStatsState sss;
    gSubproblemStatsState = &sss;
    SubproblemStatsState proc;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, rss, summarize, 1);

}
