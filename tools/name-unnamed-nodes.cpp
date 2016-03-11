#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <sstream>

#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/tree_iter.h"
#include "otc/node_naming.h"
using namespace otc;
using std::vector;
using std::unique_ptr;
using std::string;

struct RTNodeSmallestChild {
    int smallestChild = 0;
};

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

inline int smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline int& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-name-unnamed-nodes",
                "Takes a series of tree files and writes each tree with names computed for unnamed nodes.\n"
                "Nodes with OTT ids are written as ott#####, with longer names suppressed.\n",
                "tree1.tre [tree2.tre ... treeN.tree]");
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1)) {
        return 1;
    }
    if (trees.empty()) {
        throw OTCError("No trees loaded!");
    }
    // Add names to unnamed nodes
    for(const auto& tree: trees) {
        nameUnamedNodes(*tree);
    }
    for(const auto& tree: trees) {
        writeTreeAsNewick(std::cout, *tree);
        std::cout<<"\n";
    }
    return 0;
}
