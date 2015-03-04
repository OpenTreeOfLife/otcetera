#include "otc/otcli.h"
using namespace otc;

template<typename T, typename U>
void getInducedInformativeGroupingMaps(const T & tree1,
                                       std::map<std::set<long>, std::list<const typename T::node_type *> > & inducedSplitMaps,
                                       const U & tree2) {
    const auto inducingLabels = getOttIdSetForLeaves(tree2);
    auto mrca = findMRCAUsingDesIds(tree1, inducingLabels);
    std::function<bool(const typename T::node_type &)> sf = [inducingLabels](const typename T::node_type &nd){
        return haveIntersection(inducingLabels, nd.getData().desIds);
    };
    for (auto n : iter_pre_filter_n_const(mrca, sf)) {
        //std::cout << n->getOttId() << '\n';
        if (n == mrca) {
            continue;
        }
        auto inducedDesIds = n->getData().desIds;
        const auto x = intersectionOfSets(inducingLabels, inducedDesIds);
        if (x.size() > 1) {
            inducedSplitMaps[std::move(x)].push_back(n);
        }
    }
}


template<typename T>
inline void getInformativeGroupingMaps(const T & tree2,
                                      std::map<std::set<long>, const typename T::node_type *> & tree2Splits) {
    auto t2r = tree2.getRoot();
    for (auto n : iter_pre_internal_const(tree2)) {
        if (n == t2r) {
            continue;
        }
        const std::set<long> & x = n->getData().desIds;
        if (x.size() > 1) {
            tree2Splits[x] = n;
        }
    }
}


template<typename T, typename U>
unsigned long reportOnInducedConflicts(std::ostream & out,
                                       const T & tree1,
                                       const U & tree2,
                                       bool firstIsSuperset) {
    assert(firstIsSuperset);
    std::map<std::set<long>, std::list<const typename T::node_type *> > inducedSplitMap;
    std::map<std::set<long>, const typename U::node_type *> tree2Splits;
    getInducedInformativeGroupingMaps(tree1, inducedSplitMap, tree2);
    getInformativeGroupingMaps(tree2, tree2Splits);
    unsigned long nm = 0;
    for (const auto & icsm : inducedSplitMap) {
        const auto & ics = icsm.first;
        bool found = false;
        std::list<std::set<long> > extraIds;
        std::list<std::set<long> > missingIds;
        std::list<const typename U::node_type *> nodeList;
        for (const auto & t2sP : tree2Splits) {
            const auto & t2s = t2sP.first;
            if (t2s == ics) {
                found = true;
            } else {
                if (!areCompatibleDesIdSets(t2s, ics)) {
                    std::set<long> e = set_difference_as_set(t2s, ics);
                    std::set<long> m = set_difference_as_set(ics, t2s);
                    assert(!e.empty() || !m.empty());
                    extraIds.push_back(e);
                    missingIds.push_back(m);
                    nodeList.push_back(t2sP.second);
                }
            }
        }
        if (!extraIds.empty() || !missingIds.empty()) {
            for (auto taxonNode : icsm.second) {
                assert(extraIds.size() == missingIds.size());
                assert(extraIds.size() == nodeList.size());
                auto eIt = begin(extraIds);
                auto mIt = begin(missingIds);
                auto nIt = begin(nodeList);
                for (; nIt != end(nodeList); ++nIt, ++eIt, ++mIt) {
                    out << getContestedPreamble(*taxonNode, tree2);
                    emitConflictDetails(out, **nIt, *eIt, *mIt);
                }
                nm += 1;
            }
        }
    }
    return nm;
}


struct ConflictReporterProc : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    virtual ~ConflictReporterProc(){}
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        tree->setName(otCLI.currentFilename);
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
