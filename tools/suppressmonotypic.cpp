#include "otc/otcli.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_operations.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool suppressMonotypicAndWrite(OTCLI & , std::unique_ptr<T> tree);

template<typename T>
inline bool suppressMonotypicAndWrite(OTCLI & otCLI, std::unique_ptr<T> tree) {
    T * p = tree.get();
    suppressMonotypicTaxaPreserveShallowDangle(*p);
    writeTreeAsNewick(otCLI.out, *p);
    otCLI.out << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otcsuppressmonotypic",
                 "takes a filepath to a newick file and writes a newick without any nodes that have just one child",
                 {"some.tre"});
    otCLI.getParsingRules().setOttIdForInternals = false;
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wdd = suppressMonotypicAndWrite<Tree_t>;
    return treeProcessingMain<Tree_t>(otCLI, argc, argv, wdd, nullptr, 1);
}

