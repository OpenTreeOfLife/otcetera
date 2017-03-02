#include "otc/otcli.h"
#include "otc/tree.h"
#include "otc/tree_data.h"
#include "otc/tree_operations.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool suppressMonotypicAndWrite(OTCLI & , std::unique_ptr<T> tree);

// BDR - I copied this (with some modifications) from name-unnamed-nodes.cpp
//       For trees with no attributes this seems to be much faster than
//         suppressMonotypicTaxaPreserveShallowDangle(*p), which I commented out.
//       We do not need to construct a set of removed nodes here, and we probably
//         do not want to check invariants.
//       I believe this is a "shallow dangle" version, because we delete monotypic nodes
//         instead of preserving the one closest to the root.
void suppressMonotypicFast(Tree_t& tree)
{
    std::vector<Tree_t::node_type*> remove;
    for (auto nd:iter_pre(tree)) {
        if (nd->isOutDegreeOneNode()) {
            remove.push_back(nd);
        }
    }
    for (auto nd: remove) {
        auto child = nd->getFirstChild();
        child->detachThisNode();
        nd->addSibOnRight(child);
        nd->detachThisNode();
        delete nd;
    }
}

template<typename T>
inline bool suppressMonotypicAndWrite(OTCLI & otCLI, std::unique_ptr<T> tree) {
    T * p = tree.get();
//    suppressMonotypicTaxaPreserveShallowDangle(*p);
    suppressMonotypicFast(*p);
    writeTreeAsNewick(otCLI.out, *p);
    otCLI.out << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-suppress-monotypic",
                 "takes a filepath to a newick file and writes a newick without any nodes that have just one child",
                 {"some.tre"});
    otCLI.get_parsing_rules().set_ott_idForInternals = false;
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wdd = suppressMonotypicAndWrite<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wdd, nullptr, 1);
}

