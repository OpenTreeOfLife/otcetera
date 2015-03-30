#include "otc/otcli.h"
#pragma clang diagnostic ignored "-Wunused-variable"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool writeNumLeaves(OTCLI & , std::unique_ptr<T> tree);

template<typename T>
bool writeNumLeaves(OTCLI & , std::unique_ptr<T> tree) {
    auto c = 0U;
    for (auto nd : iter_leaf_const(*tree)) {
        c += 1;
    }
    std::cout << c << '\n';
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-count-leaves",
                 "takes a filepath to a newick file and reports the number of leaves",
                 {"some.tre"});
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wnl = writeNumLeaves<Tree_t>;
    return treeProcessingMain<Tree_t>(otCLI, argc, argv, wnl, nullptr, 1);
}

