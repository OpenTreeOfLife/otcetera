#include "otc/otcli.h"
#include "otc/util.h"
using namespace otc;
typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef RTreeOttIDMapping<RTSplits> RootedTreeForNodeType;
typedef otc::RootedTree<RTSplits, RootedTreeForNodeType> Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);
extern const char * badNTreesMessage;
const char * badNTreesMessage = "Expecting only 1 tree to prune\n";

struct SubtreePrunerState {
    std::unique_ptr<Tree_t> toPrune;
    std::set<long> designators;
    bool emitParent;
    int numErrors;
    SubtreePrunerState()
        :toPrune(nullptr),
        emitParent(false),
        numErrors(0) {
    }

    void summarize(OTCLI &otCLI) {
        if (toPrune == nullptr) {
            otCLI.err << badNTreesMessage;
            numErrors = 1;
            return;
        }
        if (designators.empty()) {
            otCLI.err << "Expecting a -n or -p arg to specify the IDs whose MRCA defines the subtree";
            numErrors = 1;
            return;
        }
        const Node_t * mrca = findMRCAFromIDSet(*toPrune, designators, -1);
        if (mrca == nullptr) {
            otCLI.err << "MRCA could not be found\n";
            numErrors = 1;
            return;
        }
        const Node_t * toEmit = (emitParent ? mrca->getParent() : mrca);
        if (toEmit == nullptr) {
            otCLI.err << "MRCA had no parent\n";
            numErrors = 1;
            return;
        }
        writeNewick<Node_t>(otCLI.out, toEmit);
        otCLI.out << ";\n";
    }
};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
    SubtreePrunerState * ctsp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(ctsp != nullptr);
    assert(tree != nullptr);
    if (ctsp->toPrune == nullptr) {
        ctsp->toPrune = std::move(tree);
    } else {
        otCLI.err << badNTreesMessage;
        return false;
    }
    return true;
}

bool handleNodeFlag(OTCLI & otCLI, const std::string &nextArg);
bool handleParentFlag(OTCLI & otCLI, const std::string &nextArg);

bool handleNodeFlag(OTCLI & otCLI, const std::string &nextArg) {
    SubtreePrunerState * fusp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(fusp != nullptr);
    assert(!nextArg.empty());
    fusp->designators = parseDelimSeparatedIDs(nextArg, ',');
    return true;
}

bool handleParentFlag(OTCLI & otCLI, const std::string &nextArg) {
    if (!handleNodeFlag(otCLI, nextArg)) {
        return false;
    }
    SubtreePrunerState * fusp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->emitParent = true;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-prune-to-subtree",
                "takes a -p or -n argument followed by one filepath to a newick. Writes the specified subtree to standard output stream.",
                "-n5315,3512 some.tre");
    SubtreePrunerState cts;
    otCLI.blob = static_cast<void *>(&cts);
    otCLI.addFlag('p',
                  "ARG=a comma separated list of OTT numbers. The parent of the MRCA of the IDs will be printed.",
                  handleParentFlag,
                  true);
    otCLI.addFlag('n',
                  "ARG=a comma separated list of OTT numbers. The MRCA of the IDs will be printed.",
                  handleNodeFlag,
                  true);
    auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 1);
    if (rc == 0) {
        cts.summarize(otCLI);
        return cts.numErrors;
    }
    return rc;
}
