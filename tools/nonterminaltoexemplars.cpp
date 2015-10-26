#include "otc/otcli.h"
using namespace otc;

template<typename T>
bool writeTreeOrDie(OTCLI & otCLI, const std::string & fp, const T & tree, bool useStdOut) {
    std::ostream *outPtr = &std::cout;
    std::ofstream outp;
    if (!useStdOut) {
        LOG(INFO) << "writing \"" << fp << "\"";
        outp.open(fp);
        if (!outp.good()) {
            otCLI.err << "Could not open \"" << fp << "\" to write output.\n";
            return false;
        }
        outPtr = &outp;
    } else {
        *outPtr << fp << '\n';
    }
    writeTreeAsNewick(*outPtr, tree);
    *outPtr << "\n";
    if (!useStdOut) {
        outp.close();
    }
    return true;
}

template<typename T, typename Y>
inline OttIdSet findIncludedTipIds(const T & nd, const Y & container) {
    OttIdSet r;
    for (auto t : iter_leaf_n_const(nd)) {
        if (contains(container, t)) {
            assert(t->hasOttId());
            r.insert(t->getOttId());
        }
    }
    return r;
}



template<typename T, typename Y>
inline void replaceTipWithSet(T & tree, Y * nd, const OttIdSet & oids) {
    Y * p = nd->getParent();
    Y * nps = nd->getPrevSib();
    Y * nns = nd->getNextSib();
    assert(p != nullptr);
    assert(!oids.empty());
    Y * firstC = nullptr;
    Y * lastC = nullptr;
    for (auto oid : oids) {
        lastC = tree.createChild(p);
        if (firstC == nullptr) {
            firstC = lastC;
        }
        lastC->setOttId(oid);
    }
    assert(firstC != nullptr);
    assert(lastC != nullptr);
    assert(lastC->getNextSib() == nullptr);
    Y * fcps = firstC->getPrevSib();
    assert(fcps != nullptr);
    // remove firstChild from the sib array
    fcps->_setNextSib(nullptr);
    // Add it in place of nd
    if (nps == nullptr) {
        p->_setFirstChild(firstC);
    } else {
        assert(nps->getNextSib() == nd);
        nps->_setNextSib(firstC);
    }
    // make sure that the next sib of the last child is the same as the incoming exit node.
    lastC->_setNextSib(nns);
}

struct NonTerminalsToExemplarsState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    int numErrors;
    bool useStdOut;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    std::map<std::unique_ptr<TreeMappedEmptyNodes>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedEmptyNodes *> treePtrByIndex;
    using TreeNdPair = std::pair<TreeMappedEmptyNodes *, RootedTreeNodeNoData *>;
    using ListTreeNdPair = std::list<TreeNdPair>;
    std::map<RootedTreeNodeNoData *, ListTreeNdPair> nonTermToMappedPhylo;
    std::string exportDir;
    std::ofstream nonEmptyFileStream;
    std::string outputNonEmptyTreeOutput;

    virtual ~NonTerminalsToExemplarsState(){}
    NonTerminalsToExemplarsState()
        :numErrors(0),
        useStdOut(false) {
    }

    // here we walk through the taxonomy in postorder in case we have tips like:
    //  tree 1: Homo
    //  tree 2: Hominidae
    // and no other tree with a tip that is under either Homo or Hominidae.
    // In such cases, we want to expand both 'Homo' and 'Hominidae' to the same terminal taxon.
    // this is not hard to do if we are walking in post order
    void exemplifyNonterminals() {
        for (auto nd : iter_post(*taxonomy)) {
            auto mappedPhyloListIt = nonTermToMappedPhylo.find(nd);
            if (mappedPhyloListIt == nonTermToMappedPhylo.end()) {
                continue;
            }
            const ListTreeNdPair & treeNdPairList{mappedPhyloListIt->second};
            assert(!nd->isTip());
            bool hasIncludedDes = false;
            for (auto c : iter_child(*nd)) {
                if (contains(includedNodes, c)) {
                    hasIncludedDes = true;
                    break;
                }
            }
            const auto nid = nd->getOttId();
            
            OttIdSet exemplarIDs;
            if (hasIncludedDes) {
                exemplarIDs = findIncludedTipIds(*nd, includedNodes);
            } else {
                const RootedTreeNodeNoData * n = findLeftmostInSubtree(nd);
                includedNodes.insert(n);
                insertAncestorsToParaphyleticSet(n, includedNodes);
                exemplarIDs.insert(n->getOttId());
            }
            LOG(INFO) << "Exemplifying OTT-ID" << nid << " with:";
            for (auto rid : exemplarIDs) {
                LOG(INFO) << "    " << rid;
            }
            for (auto treeNdPair : treeNdPairList) {
                TreeMappedEmptyNodes * treeP = treeNdPair.first;
                RootedTreeNodeNoData * nodeP = treeNdPair.second;
                replaceTipWithSet(*treeP, nodeP, exemplarIDs);
            }
        }
    }

    void pruneTaxonomyToIncludedLeaves() {
        assert(taxonomy != nullptr && !includedNodes.empty());
        std::set<RootedTreeNodeNoData *> toPrune;
        for (auto nd : iter_node(*taxonomy)) {
            const RootedTreeNodeNoData *  c = const_cast<const RootedTreeNodeNoData *>(nd);
            if ((!contains(includedNodes, c)) && contains(includedNodes, c->getParent())) {
                toPrune.insert(nd);
            }
        }
        for (auto nd : toPrune) {
            pruneAndDelete(*taxonomy, nd);
        }
    }

    bool summarize(OTCLI &otCLI) override {
        if ((!useStdOut) && exportDir.empty()) {
            otCLI.err << "Either the -e flag to specify and export directory or the -o flag mandatory\n";
            return false;
        }
        if ((!useStdOut) && exportDir[exportDir.length() - 1] != '/') {
            exportDir += '/';
        }
        // replace non-terminal tips with their expansion
        exemplifyNonterminals();
        // prune down the taxonomy to the set of used leaves
        pruneTaxonomyToIncludedLeaves();
        // write the output
        const std::string tn = "taxonomy.tre";
        auto tp = exportDir + tn;
        if (!writeTreeOrDie(otCLI, tp, *taxonomy, useStdOut)) {
            return false;
        }
        for (auto treePtr : treePtrByIndex) {
            auto pp = exportDir + treePtr->getName();
            if (!writeTreeOrDie(otCLI, pp, *treePtr, useStdOut)) {
                return false;
            }
        }
        nonEmptyFileStream.close();
        return true;
    }
    
    bool processTaxonomyTree(OTCLI & otCLI) override {
        bool r = TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes>::processTaxonomyTree(otCLI);
        // we can ignore the internal node labels for the non-taxonomic trees
        otCLI.getParsingRules().setOttIdForInternals = false;
        if (!outputNonEmptyTreeOutput.empty()) {
            nonEmptyFileStream.open(outputNonEmptyTreeOutput.c_str());
        }
        return r;
    }

    bool processSourceTree(OTCLI & otCLI, std::unique_ptr<TreeMappedEmptyNodes> treeup) override {
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedEmptyNodes * raw = treeup.get();
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        treePtrByIndex.push_back(raw);
        // Store the tree's filename
        raw->setName(otCLI.currentFilename);
        std::map<const RootedTreeNodeNoData *, std::set<long> > prunedDesId;
        auto nleaves = 0;
        for (auto nd : iter_leaf(*raw)) {
            nleaves += 1;
            auto ottId = nd->getOttId();
            auto taxoNode = taxonomy->getData().getNodeForOttId(ottId);
            assert(taxoNode != nullptr);
            if (!contains(includedNodes, taxoNode)) {
                includedNodes.insert(taxoNode);
                insertAncestorsToParaphyleticSet(taxoNode, includedNodes);
            }
            if (!taxoNode->isTip()) {
                TreeNdPair tnp{raw, nd};
                nonTermToMappedPhylo[taxoNode].push_back(tnp);
            }
        }
        if (!outputNonEmptyTreeOutput.empty()) {
            nonEmptyFileStream << otCLI.currentFilename << '\n';
            nonEmptyFileStream.flush();
        }
        return true;
    }
};

bool handleNonemptyTreeOutput(OTCLI & otCLI, const std::string &narg);
bool handleExportModified(OTCLI & otCLI, const std::string &narg);
bool handleStdout(OTCLI & otCLI, const std::string &narg);

bool handleStdout(OTCLI & otCLI, const std::string &) {
    NonTerminalsToExemplarsState * proc = static_cast<NonTerminalsToExemplarsState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->useStdOut = true;
    return true;
}

bool handleExportModified(OTCLI & otCLI, const std::string &narg) {
    NonTerminalsToExemplarsState * proc = static_cast<NonTerminalsToExemplarsState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -e argument.");
    }
    proc->exportDir = narg;
    return true;
}

bool handleNonemptyTreeOutput(OTCLI & otCLI, const std::string &narg) {
    NonTerminalsToExemplarsState * proc = static_cast<NonTerminalsToExemplarsState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting a list of IDs after the -n argument.");
    }
    proc->outputNonEmptyTreeOutput = narg;
    return true;
}


int main(int argc, char *argv[]) {
    const char * helpMsg = "takes an -e flag specifying an export diretory and at least 2 newick file paths: " \
        "a full taxonomy tree some number of input trees. Any tip in non-taxonomic input that is mapped to " \
        "non-terminal taoxn will be remapped such that the parent of the non-terminal tip will hold all of " \
        "the expanded exemplars. The exemplars will be the union of tips that (a) occur below this non-terminal " \
        "taxon in the taxonomy and (b) occur, or are used as an exemplar, in another input tree. The modified " \
        "version of each input will be written in the export directory. Trees with no non-terminal tips should " \
        "be unaltered. The taxonomy written out will be the taxonomy restricted to the set of leaves that are " \
        "leaves of the exported trees";
    OTCLI otCLI("otc-nonterminals-to-exemplars",
                helpMsg,
                "-estep_5 taxonomy.tre inp1.tre inp2.tre");
    NonTerminalsToExemplarsState proc;
    otCLI.addFlag('e',
                  "ARG should be the name of a directory. A .tre file will be written to this directory for each input tree",
                  handleExportModified,
                  true);
    otCLI.addFlag('o',
                  " requests that standard output stream, rather than the export directory, be used for all output.",
                  handleStdout,
                  false);
    otCLI.addFlag('n',
                  "ARG is and output file that will list the filename of input tree that was not empty",
                  handleNonemptyTreeOutput,
                  true);
    return taxDependentTreeProcessingMain(otCLI, argc, argv, proc, 2, false);
}
