#include "otc/otcli.h"
#include "otc/supertree_util.h"
using namespace otc;

template <typename T, typename U>
void resolveNode(T & tree, U & parent, const OttIdSet & newInc) {
    std::set<U *> childrenToMove;
    for (auto nd : iter_child(parent)) {
        if (!areDisjoint(nd->getData().desIds, newInc)) {
            childrenToMove.insert(nd);
        }
    }
    assert(!childrenToMove.empty());
    U * phPar = &parent;
   assert(phPar != nullptr);
    U * insertedNodePtr = tree.createNode(phPar);
    for (auto phChild : childrenToMove) {
        phChild->_detachThisNode();
        insertedNodePtr->addChild(phChild);
        const OttIdSet & cd = phChild->getData().desIds;
        insertedNodePtr->getData().desIds.insert(cd.begin(), cd.end());
    }
    assert(phPar->getOutDegree() > 1);
}

struct FindResolutionState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    std::unique_ptr<TreeMappedWithSplits> toCheck;
    std::size_t numIncludable;
    bool addGroups;
    int numErrors;
    virtual ~FindResolutionState(){}
    FindResolutionState()
        :toCheck(nullptr),
        numIncludable(0U),
        addGroups(false),
        numErrors(0) {
    }

    bool summarize(OTCLI &otCLI) override {
        if (addGroups) {
            otCLI.err << numIncludable << " nodes added to the tree.\n";
            if (numIncludable == 0) {
                return false;
            }
            writeTreeAsNewick(otCLI.out, *toCheck);
            otCLI.out << '\n';
            return true;
        }
        otCLI.out << numIncludable << " clades found which could be added to the tree.\n";
        numErrors = static_cast<int>(numIncludable);
        return numErrors == 0;
    }

    virtual bool processTaxonomyTree(OTCLI & otCLI) override {
        TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::processTaxonomyTree(otCLI);
        otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
        // now we get a little cute and reprocess the taxonomy desIds so that they 
        // exclude internals. So that when we expand source trees, we expand just
        // to the taxonomy's leaf set
        clearAndfillDesIdSets(*taxonomy);
        return true;
    }
    
    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> tree) override {
        assert(taxonomy != nullptr);
        if (toCheck == nullptr) {
            toCheck = std::move(tree);
            otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
            return true;
        }
        expandOTTInternalsWhichAreLeaves(*tree, *taxonomy);
        clearAndfillDesIdSets(*tree);
        const OttIdSet & treeLeafSet = tree->getRoot()->getData().desIds;
        std::set<const NodeWithSplits *> nodesFromInp;
        for (auto nd : iter_pre_internal_const(*tree)) {
            if (nd->getParent() != nullptr) {
                nodesFromInp.insert(nd);
            }
        }
        while (!nodesFromInp.empty()) {
            std::map<NodeWithSplits *, std::list<const NodeWithSplits *> > toCheckNodeToResolves;
            for (const auto nd : nodesFromInp) {
                const OttIdSet & incGroup = nd->getData().desIds;
                auto r = findMRCAUsingDesIds(*toCheck, incGroup);
                NodeWithSplits * mrca = const_cast<NodeWithSplits *>(r);
                assert(mrca);
                OttIdSet md = mrca->getData().desIds;
                OttIdSet rmd;
                for (auto i : md) {
                    if (contains(treeLeafSet, i)) {
                        rmd.insert(i);
                    }
                }
                if (incGroup == rmd) {
                    continue;
                }
                OttIdSet excGroup = treeLeafSet;
                for (auto i :incGroup) {
                    excGroup.erase(i);
                }
                if (canBeResolvedToDisplayExcGroup(mrca, incGroup, excGroup)) {
                    if (addGroups) {
                        toCheckNodeToResolves[mrca].push_back(nd);
                    } else {
                        numIncludable += 1;
                        otCLI.out << otCLI.currentFilename << " node " << getDesignator(*nd) << " could be added.\n";
                    }
                }
            }
            // rather than deal with the logic for adding to the correct node in the resolved tree, we'll just iterate
            //  over the whole loop again. NOT EFFICIENT. MTH is just being lazy...
            std::set<const NodeWithSplits *> toDealWith;
            for (auto mrcaListIncPair : toCheckNodeToResolves) {
                auto mrcaPtr = mrcaListIncPair.first;
                const NodeWithSplits * toAdd = nullptr;
                for (auto tac : mrcaListIncPair.second) {
                    if (toAdd == nullptr) {
                        toAdd = tac;
                    } else {
                        toDealWith.insert(tac);
                    }
                }
                if (toAdd != nullptr) {
                    numIncludable += 1;
                    const OttIdSet & incGroup = toAdd->getData().desIds;
                    resolveNode(*toCheck, *mrcaPtr, incGroup);
                }
            }
            nodesFromInp = toDealWith;
        }
        return true;
    }

};
bool handleResolve(OTCLI & otCLI, const std::string &);

bool handleResolve(OTCLI & otCLI, const std::string &) {
    FindResolutionState * proc = static_cast<FindResolutionState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->addGroups = true;
    return true;
}


int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-find-unsupported-nodes",
                "takes at least 2 newick file paths: a full supertree, and some number of input trees",
                "synth.tre inp1.tre inp2.tre ...");
    FindResolutionState proc;
    otCLI.addFlag('r',
                  "Resolve the supertree rather than counting number of groups that could be added.",
                  handleResolve,
                  false);
    
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 3, true);
}
