#include "otc/otcli.h"
#pragma clang diagnostic ignored "-Wunused-variable"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
bool handleListTips(OTCLI & , const std::string &);
template<typename T>
bool writeNumLeaves(OTCLI & , std::unique_ptr<T> tree);
template<typename T>
bool listTipOttIds(OTCLI & , std::unique_ptr<T> tree);

static bool listTips = false;

template<typename T>
bool listTipOttIds(OTCLI & , T * tree) {
    for (auto nd : iter_leaf_const(*tree)) {
        assert(nd->hasOttId());
        std::cout << nd->getOttId() << '\n';
    }
    return true;
}

template<typename T>
bool writeNumLeaves(OTCLI & otCLI, std::unique_ptr<T> tree) {
    if (listTips) {
        return listTipOttIds(otCLI, tree.get());
    }
    auto c = 0U;
    // the nd loop var is unused
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    for (auto nd : iter_leaf_const(*tree)) {
        c += 1;
    }
#pragma GCC diagnostic pop
    std::cout << c << '\n';
    return true;
}


bool handleListTips(OTCLI & , const std::string &) {
    listTips = true;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-count-leaves",
                 "takes a filepath to a newick file and reports the number of leaves",
                 {"some.tre"});
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wnl = writeNumLeaves<Tree_t>;
    otCLI.addFlag('l',
                  "If present, list the tip OTT IDs rather than counting them",
                  handleListTips,
                  false);
    return treeProcessingMain<Tree_t>(otCLI, argc, argv, wnl, nullptr, 1);
}

