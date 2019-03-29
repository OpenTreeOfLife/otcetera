#include "otc/otcli.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool writeDegreeDistribution(OTCLI & , std::unique_ptr<T> tree);

template<typename T>
inline bool writeDegreeDistribution(OTCLI & otCLI, std::unique_ptr<T> tree) {
    std::map<unsigned long, unsigned long> degreeDistribution;
    for (auto nd : iter_pre_const(*tree)) {
        auto od = nd->get_out_degree();
        degreeDistribution[od] += 1;
    }
    otCLI.out << "Out-degree\tCount\n";
    std::size_t numNodes = 0U;
    for (auto p: degreeDistribution) {
        otCLI.out << p.first << "\t" << p.second << "\n";
        numNodes += p.second;
    }
    const auto numLeaves = degreeDistribution[0UL];
    const auto numElbows = (contains(degreeDistribution, 1UL) ? degreeDistribution[1UL]: 0UL);
    otCLI.err << numLeaves << " leaves\n";
    otCLI.err << numNodes - numLeaves << " non-leaf nodes (internals including the root)\n";
    otCLI.err << numNodes - numLeaves - numElbows<< " non-leaf nodes with out-degree > 1 (internals including the root)\n";
    otCLI.err << numNodes << " total nodes\n";
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-degree-distribution",
                 "takes a filepath to a newick file and reports the number of nodes of each out-degree",
                 "some.tre");
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wdd = writeDegreeDistribution<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wdd, nullptr, nullptr, 1);
}

