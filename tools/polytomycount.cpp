#include "otc/otcli.h"
#include "otc/tree_operations.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
template<typename T>
bool writeNumPolytomies(OTCLI & otCLI, std::unique_ptr<T> tree);

template<typename T>
inline bool writeNumPolytomies(OTCLI & otCLI, std::unique_ptr<T> tree) {
    otCLI.out << count_polytomies(*tree) << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-polytomy-count",
                 "takes a filepath to a newick file and reports the number of polytomies in each tree (one line per tree)",
                 {"some.tre"});
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wnp = writeNumPolytomies<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wnp, nullptr, 1);
}

