#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"

using namespace otc;
using std::vector;
using std::unique_ptr;

struct RTNodePartialDesSet {
    long smallestChild  = 0;
};

using Tree_t = RootedTree<RTNodePartialDesSet, RTreeNoData>;

inline long smallestChild(const Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

inline long& smallestChild(Tree_t::node_type* node) {
    return node->getData().smallestChild;
}

template<typename T>
void writeTreesFromSampledTipsToEncodeStructure(T & tree, std::ostream & outStream);
template<typename T>
void writeTreesFromSampledTipsToEncodeStructureDecoratedTree(const T & tree, std::ostream & outStream);
template<typename N>
void writeSmallestIDForEachChild(const N &ndRef, std::ostream & outStream);


template<typename T>
void writeTreesFromSampledTipsToEncodeStructure(T & tree, std::ostream & outStream) {
    calculateSmallestChild(tree);
    writeTreesFromSampledTipsToEncodeStructureDecoratedTree(tree, outStream);
}


template<typename N>
void writeSmallestIDForEachChild(const N &ndRef, std::ostream & outStream) {
    outStream << '(';
    bool first = true;
    for (auto c: iter_child_const(ndRef)) {
        if (first) {
            first = false;
        } else {
            outStream << ',';
        }
        outStream << "ott" << smallestChild(c);
    }
    outStream << ')';
}

template<typename T>
void writeTreesFromSampledTipsToEncodeStructureDecoratedTree(const T & tree, std::ostream & outStream) {
    for (auto nd: iter_post_const(tree)){
        if (nd->isTip()) {
            continue;
        }
        auto anc = findFirstForkingAnc(nd);
        if (anc == nullptr) {
            // write a trivial statement for the root, just to make sure that all of the tips show up in one
            // tree. Without this leaves attached to the root are omitted.
            writeSmallestIDForEachChild(*nd, outStream);
            outStream << ";\n";
            continue;
        }
        outStream << '(';
        writeSmallestIDForEachChild(*nd, outStream);
        outStream << ',';
        bool found_sib = false;
        for (auto possible_sib: iter_child_const(*anc)) {
            if ((possible_sib != nd) and (!isAncestorDesNoIter(possible_sib, nd))) {
                outStream << "ott" << smallestChild(possible_sib);
                found_sib = true;
                break;
            }
        }
        assert(found_sib);
        outStream << ");\n";
    }
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-sample-tips",
                "Takes a rooted phylogenetic estimate (the first tree file) and writes a series of trees that are induced\n" \
                "trees generated from the tree by sampling subsets of leaves.\n" \
                "The default behavior is to write one tree per internal branch in the original tree.\n" \
                "For each branch the induced tree will encode the descendant node.\n" \
                "This encoding is done by contains a descendant of each child of the node\n" \
                "    and descendant of one of the node\'s siblings.",
                "tree-to-sample.tre");
    vector<unique_ptr<Tree_t>> trees;
    auto get = [&trees](OTCLI &, unique_ptr<Tree_t> nt) {trees.push_back(std::move(nt)); return true;};
    if (argc < 2) {
        throw OTCError("No trees provided!");
    }
    if (treeProcessingMain<Tree_t>(otCLI, argc, argv, get, nullptr, 1)) {
        return 1;
    }
    if (trees.size() != 1) {
        std::cerr << "Supplied " << trees.size() << " trees for sampling, should be 1 tree!";
        return 2;
    }
    auto & solution = *(trees.at(0));
    writeTreesFromSampledTipsToEncodeStructure(solution, std::cout);
    return 0;
}
