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
//         suppress_monotypic_taxa_preserve_shallow_dangle(*p), which I commented out.
//       We do not need to construct a set of removed nodes here, and we probably
//         do not want to check invariants.
//       I believe this is a "shallow dangle" version, because we delete monotypic nodes
//         instead of preserving the one closest to the root.
void suppress_monotypic_fast(Tree_t& tree) {
    std::vector<Tree_t::node_type*> remove;
    for (auto nd:iter_pre(tree)) {
        if (nd->is_outdegree_one_node()) {
            remove.push_back(nd);
        }
    }
    for (auto nd: remove) {
        auto child = nd->get_first_child();
        child->detach_this_node();
        nd->add_sib_on_right(child);
        nd->detach_this_node();
        delete nd;
    }
}

template<typename T>
inline bool suppressMonotypicAndWrite(OTCLI & otCLI, std::unique_ptr<T> tree) {
    T * p = tree.get();
//    suppress_monotypic_taxa_preserve_shallow_dangle(*p);
    suppress_monotypic_fast(*p);
    write_tree_as_newick(otCLI.out, *p);
    otCLI.out << std::endl;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-suppress-monotypic",
                 "takes a filepath to a newick file and writes a newick without any nodes that have just one child",
                 "some.tre");
    otCLI.get_parsing_rules().set_ott_idForInternals = false;
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wdd = suppressMonotypicAndWrite<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wdd, nullptr, 1);
}

