#include "otc/otcli.h"
#include <unordered_set>
#include "otc/induced_tree.h"
#include "otc/tree_operations.h"

using namespace otc;
using std::vector;


struct RTNodeDepth {
    int depth = 0; // depth = number of nodes to the root of the tree including the  endpoints (so depth of root = 1)
    RootedTreeNode<RTNodeDepth>* summary_node;
};

using Tree_t = RootedTree<RTNodeDepth, RTreeNoData>;

struct InducedSubtreeState
  : public TaxonomyDependentTreeProcessor<Tree_t> {
    std::unordered_set<long> inducingIds;
    virtual ~InducedSubtreeState(){}

    // write the induced tree to the output stream
    virtual bool summarize(OTCLI &otCLI) override {
        auto tax_node_map = get_ottid_to_const_node_map(*taxonomy);
        vector<const Tree_t::node_type*> leaves;
        for(auto id: inducingIds)
            leaves.push_back(tax_node_map.at(id));
        compute_depth(*taxonomy);
        auto induced_tree = get_induced_tree<Tree_t>(leaves);
        write_tree_as_newick(otCLI.out, *induced_tree);
        otCLI.out << std::endl;
        return true;
    }

    // accumulate the set of leaves to include in inducingIds
    virtual bool process_source_tree(OTCLI &, std::unique_ptr<Tree_t> tree) override {
        assert(tree != nullptr);
        assert(taxonomy != nullptr);
        for(auto leaf: iter_leaf_const(*tree))
            inducingIds.insert(leaf->get_ott_id());
        return true;
    }
};

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-induced-subtree",
                "takes at least 2 newick file paths: a full tree, and some number of input trees. Writes the topology of the first tree if it is pruned down to the leafset of the inputs (without removing internal nodes)",
                "taxonomy.tre inp1.tre inp2.tre");
    InducedSubtreeState proc;
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, true);
}


