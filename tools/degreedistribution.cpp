#include "otc/otcli.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool writeDegreeDistribution(OTCLI & , std::unique_ptr<T> tree);

template<typename T>
inline bool writeDegreeDistribution(OTCLI & otCLI, std::unique_ptr<T> tree) {
    writeDegDist(otCLI.out, &(otCLI.err), *tree);
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-degree-distribution",
                 "takes a filepath to a newick file and reports the number of nodes of each out-degree",
                 "some.tre");
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wdd = writeDegreeDistribution<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wdd, nullptr, nullptr, 1);
}

