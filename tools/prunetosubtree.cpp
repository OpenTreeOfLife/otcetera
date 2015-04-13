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
    bool emitChildren;
    bool emitSiblings;
    int numErrors;
    SubtreePrunerState()
        :toPrune(nullptr),
        emitParent(false),
        emitChildren(false),
        emitSiblings(false),
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
        const Node_t * toEmit = ((emitParent || emitSiblings) ? mrca->getParent() : mrca);
        if (toEmit == nullptr) {
            otCLI.err << "MRCA had no parent\n";
            numErrors = 1;
            return;
        }
        if (emitChildren && toEmit->isTip()) {
            otCLI.err << "MRCA has no children\n";
            numErrors = 1;
            return;
        }
        if ((!emitSiblings) && (!emitChildren)) {
            writeNewick<Node_t>(otCLI.out, toEmit);
            otCLI.out << ";\n";
        } else {
            for (auto c : iter_child_const(*toEmit)) {
                if (emitSiblings && c == mrca) {
                    continue;
                }
                writeNewick<Node_t>(otCLI.out, c);
                otCLI.out << ";\n";
            }
        }
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
bool handleChildrenFlag(OTCLI & otCLI, const std::string &nextArg);
bool handleSiblingsFlag(OTCLI & otCLI, const std::string &nextArg);

bool handleNodeFlag(OTCLI & otCLI, const std::string &nextArg) {
    SubtreePrunerState * fusp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(fusp != nullptr);
    assert(!nextArg.empty());
    if (!fusp->designators.empty()) {
        otCLI.err << "Expecting only 1 flag listing designators.\n";
        return false;
    }
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

bool handleChildrenFlag(OTCLI & otCLI, const std::string &nextArg) {
    if (!handleNodeFlag(otCLI, nextArg)) {
        return false;
    }
    SubtreePrunerState * fusp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->emitChildren = true;
    return true;
}

bool handleSiblingsFlag(OTCLI & otCLI, const std::string &nextArg) {
    if (!handleNodeFlag(otCLI, nextArg)) {
        return false;
    }
    SubtreePrunerState * fusp = static_cast<SubtreePrunerState *>(otCLI.blob);
    assert(fusp != nullptr);
    fusp->emitSiblings = true;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-prune-to-subtree",
                "takes a -c, -n, -s or -p argument followed by one filepath to a newick. Writes the specified subtree to standard output stream.",
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
    otCLI.addFlag('s',
                  "ARG=a comma separated list of OTT numbers. The siblings of the MRCA of the IDs will be printed.",
                  handleSiblingsFlag,
                  true);
    otCLI.addFlag('c',
                  "ARG=a comma separated list of OTT numbers. The children of the MRCA of the IDs will be printed.",
                  handleChildrenFlag,
                  true);
    auto rc = treeProcessingMain<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 1);
    if (rc == 0) {
        cts.summarize(otCLI);
        return cts.numErrors;
    }
    return rc;
}
