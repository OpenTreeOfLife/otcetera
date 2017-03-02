#include "otc/otcli.h"
#include "otc/node_naming.h"
#include "otc/tree_iter.h"

using namespace otc;
using std::vector;
using std::unique_ptr;

using Tree_t = RootedTree<RTNodeSmallestChild, RTreeNoData>;

template<typename T>
void writeTreesFromSampledTipsToEncodeStructure(T & tree, std::ostream & outStream);
template<typename T>
void writeTreesFromSampledTipsToEncodeStructureDecoratedTree(const T & tree, std::ostream & outStream);
template<typename N>
void writeSmallestIDForEachChild(const N &ndRef, std::ostream & outStream);
template<typename T>
void writeTreesFromSampledTipsOnePerInternalAllChildrenOneSib(const T & tree, std::ostream & outStream);
template<typename T>
void writeTreesFromSampledTipsNamedInternalsOneDesPerOutgoing(const T & tree, std::ostream & outStream);

enum SampleMode {
    TREE_PER_INTERNAL_ALL_CHILDREN_ONE_SIB = 0,
    NAMED_INTERNALS_ONE_DES_PER_OUTGOING = 1
};

SampleMode gSampleMode = SampleMode::TREE_PER_INTERNAL_ALL_CHILDREN_ONE_SIB;

template<typename T>
void writeTreesFromSampledTipsToEncodeStructure(T & tree, std::ostream & outStream, SampleMode sm) {
    calculate_smallest_child(tree);
    writeTreesFromSampledTipsToEncodeStructureDecoratedTree(tree, outStream, sm);
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
        outStream << "ott" << smallest_child(c);
    }
    outStream << ')';
}

template<typename T>
void writeTreesFromSampledTipsToEncodeStructureDecoratedTree(const T & tree, std::ostream & outStream, SampleMode sm) {
    if (sm == SampleMode::TREE_PER_INTERNAL_ALL_CHILDREN_ONE_SIB) {
        writeTreesFromSampledTipsOnePerInternalAllChildrenOneSib(tree, outStream);
    } else {
        assert(sm == SampleMode::NAMED_INTERNALS_ONE_DES_PER_OUTGOING);
        writeTreesFromSampledTipsNamedInternalsOneDesPerOutgoing(tree, outStream);
    }
}

template<typename T>
void writeTreesFromSampledTipsOnePerInternalAllChildrenOneSib(const T & tree, std::ostream & outStream) {
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
                outStream << "ott" << smallest_child(possible_sib);
                found_sib = true;
                break;
            }
        }
        assert(found_sib);
        outStream << ");\n";
    }
}


template<typename N>
void emitNewickForSubtreeIndirectRecursion(const N * node, std::ostream & outStream);
template<typename N>
void emitNewickForNamedInternal(const N * node, std::ostream & outStream, bool isRootOfSubtree);

template<typename N>
void emitNewickForSubtreeIndirectRecursion(const N * node, std::ostream & outStream) {
    outStream << '(';
    bool first_child = true;
    for (auto c: iter_child_const(*node)) {
        if (first_child) {
            first_child = false;
        } else {
            outStream << ',';
        }
        emitNewickForNamedInternal(c, outStream, false);
    }
    outStream << ')';
}

template<typename N>
void emitNewickForNamedInternal(const N * node, std::ostream & outStream, bool isRootOfSubtree) {
    if (isRootOfSubtree) {
        emitNewickForSubtreeIndirectRecursion(node, outStream);
        outStream << ";\n";
    } else {
        if (node->has_ott_id()) {
            if (node->isTip()) {
                outStream << "ott" << node->get_ott_id();
            } else {
                outStream << "ott" << smallest_child(node);
            }
        } else {
            emitNewickForSubtreeIndirectRecursion(node, outStream);
        }
    }
}

template<typename T>
void writeTreesFromSampledTipsNamedInternalsOneDesPerOutgoing(const T & tree, std::ostream & outStream) {
    for (auto nd: iter_post_const(tree)){
        if (nd->isTip()) {
            continue;
        }
        auto anc = findFirstForkingAnc(nd);
        if (nd->has_ott_id() || anc == nullptr) {
            emitNewickForNamedInternal(nd, outStream, true);
        }
    }
}


bool handleSubproblemSampling(OTCLI &, const std::string &) {
    gSampleMode = SampleMode::NAMED_INTERNALS_ONE_DES_PER_OUTGOING;
    return true;
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
    otCLI.addFlag('s',
                  "Subproblem sampling. An output tree for every named input with all outgoing arcs sampled.",
                  handleSubproblemSampling,
                  false);
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
    writeTreesFromSampledTipsToEncodeStructure(solution, std::cout, gSampleMode);
    return 0;
}
