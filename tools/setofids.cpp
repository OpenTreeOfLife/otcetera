#include "otc/otcli.h"
#include "otc/util.h"
using namespace otc;
typedef otc::RootedTreeNode<RTSplits> Node_t;
typedef otc::RootedTree<RTSplits, RTreeOttIDMapping<RTSplits> > Tree_t;
bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree);

struct SetOfIdsState {
    OttIdSet idsEncountered;
    bool includeInternals;
    int numErrors;
    bool noTreesYet;
    bool showIntersection;
    bool asNewick;
    SetOfIdsState()
        :includeInternals(true),
        numErrors(0),
        noTreesYet(true),
        showIntersection(false),
        asNewick(false) {
    }
    void summarize(OTCLI &otCLI) {
        if (asNewick) {
            otCLI.out << '(';
            bool f = true;
            for (const auto & oid : idsEncountered) {
                if (!f) {
                    otCLI.out << ',';
                }
                otCLI.out << "ott" << oid;
                f = false;
            }
            otCLI.out << ");\n";
        } else {
            for (const auto & oid : idsEncountered) {
                otCLI.out << oid << '\n';
            }
        }
    }
    bool intersectionProcessNextTree(const Tree_t & tree) {
        if (idsEncountered.empty()) {
            return false;
        }
        OttIdSet ier;
        fillWithIds(tree, ier);
        OttIdSet tmp = set_intersection_as_set(ier, idsEncountered);
        idsEncountered.swap(tmp);
        return !idsEncountered.empty();
    }

    void fillWithIds(const Tree_t & tree, OttIdSet & ier) {
        for (const auto nd : iter_node_const(tree)) {
            if (nd->has_ott_id() && (includeInternals || nd->is_tip())) {
                ier.insert(nd->get_ott_id());
            }
        }
    }

};

inline bool processNextTree(OTCLI & otCLI, std::unique_ptr<Tree_t> tree) {
    SetOfIdsState * ctsp = static_cast<SetOfIdsState *>(otCLI.blob);
    assert(ctsp != nullptr);
    if (ctsp->showIntersection && !ctsp->noTreesYet) {
        return ctsp->intersectionProcessNextTree(*tree);
    }
    ctsp->noTreesYet = false;
    OttIdSet & ier = ctsp->idsEncountered;
    ctsp->fillWithIds(*tree, ier);
    return true;
}

bool handleTipsOnly(OTCLI & otCLI, const std::string &);
bool handleIntersectionOnly(OTCLI & otCLI, const std::string &);
bool handleNewick(OTCLI & otCLI, const std::string &);

bool handleTipsOnly(OTCLI & otCLI, const std::string &) {
    SetOfIdsState * proc = static_cast<SetOfIdsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->includeInternals = false;
    return true;
}
bool handleIntersectionOnly(OTCLI & otCLI, const std::string &) {
    SetOfIdsState * proc = static_cast<SetOfIdsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->showIntersection = true;
    return true;
}
bool handleNewick(OTCLI & otCLI, const std::string &) {
    SetOfIdsState * proc = static_cast<SetOfIdsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->asNewick = true;
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-set-of-ids",
                "takes a series of newick file paths and writes the union of the OTT Ids found in all trees",
                "some.tre");
    SetOfIdsState cts;
    otCLI.blob = static_cast<void *>(&cts);
    otCLI.add_flag('t',
                  "Only consider OTT ids from tips",
                  handleTipsOnly,
                  false);
    otCLI.add_flag('i',
                  "write the intersection of the OTT Ids rather than the union.",
                  handleIntersectionOnly,
                  false);
    otCLI.add_flag('n',
                  "write output as a newick polytomy.",
                  handleNewick,
                  false);
    auto rc = tree_processing_main<Tree_t>(otCLI, argc, argv, processNextTree, nullptr, 1);
    if (rc == 0) {
        cts.summarize(otCLI);
        return cts.numErrors;
    }
    return rc;
}
