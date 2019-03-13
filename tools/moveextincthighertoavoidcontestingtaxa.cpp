#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "json.hpp"
using namespace otc;
using json = nlohmann::json;

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
    write_tree_as_newick(*outPtr, tree);
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
            assert(t->has_ott_id());
            r.insert(t->get_ott_id());
        }
    }
    return r;
}


static bool transferNodeNameToExemplars = true; //use -t command line flag to set this to false

template<typename T, typename Y>
inline void replaceTipWithSet(T & tree, Y * nd, const OttIdSet & oids) {
    bool hasNodeName = false;
    std::string noden;
    if (transferNodeNameToExemplars) {
        std::string onn = nd->get_name();
        auto opts = get_source_node_name(onn);
        if (opts) {
            noden = *opts;
            hasNodeName = true;
        }
    }
    for (auto oid : oids) {
        auto x = tree.create_node(nullptr);
        x->set_ott_id(oid);
        if (hasNodeName) {
            std::string n = noden + " ott" + std::to_string(oid);
            x->set_name(n);
        }
        nd->add_sib_on_left(x);
    }
    nd->detach_this_node();
}

void writeJSONLogForExemplifications(const std::string & outfilename,
                                     const std::map<OttId, OttIdSet> & exemplificationsForJSONLog,
                                     const std::map<OttId, std::list<std::string> > & exemplifedTaxonToTreeNamesForJSONLog) {
    std::ofstream outstream(outfilename);
    if (!outstream.good()) {
        throw OTCError() << "Could not open the output filepath \"" << outfilename << "\"\n";
    }
    json document;
    json exemp;
    for (auto p: exemplificationsForJSONLog) {
        const auto & ottId = p.first;
        const auto & expandedSet = p.second;
        json thisTaxon;
        json exemplarsJSON = json::array();
        for (auto eid : expandedSet) {
            std::string eOttIdStr = "ott" + std::to_string(eid);
            exemplarsJSON.push_back(eOttIdStr);
        }
        thisTaxon["exemplars_used"] = exemplarsJSON;
        const auto & treeList = exemplifedTaxonToTreeNamesForJSONLog.at(ottId);
        json treeNamesJSON = json::array();
        for (auto treeName : treeList) {
            treeNamesJSON.push_back(treeName);
        }
        thisTaxon["trees_modified"] = treeNamesJSON;
        std::string ottIdStr = "ott" + std::to_string(ottId);
        exemp[ottIdStr] = thisTaxon;
    }
    document["taxa_exemplified"] = exemp;
    outstream << document.dump(1) << std::endl;
}

struct MoveExtinctHigherState : public TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes> {
    int numErrors;
    bool useStdOut;
    std::set<const RootedTreeNodeNoData *> includedNodes;
    std::map<std::unique_ptr<TreeMappedEmptyNodes>, std::size_t> inputTreesToIndex;
    std::vector<TreeMappedEmptyNodes *> treePtrByIndex;
    using TreeNdPair = std::pair<TreeMappedEmptyNodes *, RootedTreeNodeNoData *>;
    using ListTreeNdPair = std::list<TreeNdPair>;
    std::map<RootedTreeNodeNoData *, ListTreeNdPair> nonTermToMappedPhylo;
    std::string exportTreeFile;
    std::ofstream nonEmptyFileStream;
    std::string outputNonEmptyTreeOutput;
    std::string outputJSONLogFile;
    std::string extinctInputJSON;
    std::string extinctInputJSONInline;
    bool storeLogInfo = false;
    std::map<OttId, OttIdSet> exemplificationsForJSONLog;
    std::map<OttId, std::list<std::string> > exemplifedTaxonToTreeNamesForJSONLog;

    virtual ~MoveExtinctHigherState(){}
    MoveExtinctHigherState()
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
            assert(!nd->is_tip());
            bool hasIncludedDes = false;
            for (auto c : iter_child(*nd)) {
                if (contains(includedNodes, c)) {
                    hasIncludedDes = true;
                    break;
                }
            }
            const auto nid = nd->get_ott_id();
            
            OttIdSet exemplarIDs;
            if (hasIncludedDes) {
                exemplarIDs = findIncludedTipIds(*nd, includedNodes);
            } else {
                const RootedTreeNodeNoData * n = find_leftmost_in_subtree(nd);
                includedNodes.insert(n);
                insert_ancestors_to_paraphyletic_set(n, includedNodes);
                exemplarIDs.insert(n->get_ott_id());
            }
            LOG(INFO) << "Exemplifying OTT-ID" << nid << " with:";
            for (auto rid : exemplarIDs) {
                LOG(INFO) << "    " << rid;
            }
            if (storeLogInfo) {
                exemplificationsForJSONLog[nid] = exemplarIDs;
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
            if ((!contains(includedNodes, c)) && contains(includedNodes, c->get_parent())) {
                toPrune.insert(nd);
            }
        }
        for (auto nd : toPrune) {
            prune_and_delete(*taxonomy, nd);
        }
    }

    bool summarize(OTCLI &otCLI) override {
        if ((!useStdOut) && exportTreeFile.empty()) {
            otCLI.err << "Either the -e flag to specify and export directory or the -o flag mandatory\n";
            return false;
        }
        if ((!useStdOut) && exportTreeFile[exportTreeFile.length() - 1] != '/') {
            exportTreeFile += '/';
        }
        // replace non-terminal tips with their expansion
        exemplifyNonterminals();
        // prune down the taxonomy to the set of used leaves
        pruneTaxonomyToIncludedLeaves();
        // write the output
        const std::string tn = "taxonomy.tre";
        auto tp = exportTreeFile + tn;
        if (!writeTreeOrDie(otCLI, tp, *taxonomy, useStdOut)) {
            return false;
        }
        for (auto treePtr : treePtrByIndex) {
            auto pp = exportTreeFile + treePtr->get_name();
            if (!writeTreeOrDie(otCLI, pp, *treePtr, useStdOut)) {
                return false;
            }
        }
        nonEmptyFileStream.close();
        if (storeLogInfo) {
            writeJSONLogForExemplifications(outputJSONLogFile, exemplificationsForJSONLog, exemplifedTaxonToTreeNamesForJSONLog);
        }
        return true;
    }
    
    bool process_taxonomy_tree(OTCLI & otCLI) override {
        bool r = TaxonomyDependentTreeProcessor<TreeMappedEmptyNodes>::process_taxonomy_tree(otCLI);
        // we can ignore the internal node labels for the non-taxonomic trees
        otCLI.get_parsing_rules().set_ott_idForInternals = false;
        if (!outputNonEmptyTreeOutput.empty()) {
            nonEmptyFileStream.open(outputNonEmptyTreeOutput.c_str());
        }
        return r;
    }

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedEmptyNodes> treeup) override {
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedEmptyNodes * raw = treeup.get();
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        treePtrByIndex.push_back(raw);
        // Store the tree's filename
        raw->set_name(otCLI.currentFilename);
        std::map<const RootedTreeNodeNoData *, OttIdSet > prunedDesId;
        auto nleaves = 0;
        for (auto nd : iter_leaf(*raw)) {
            nleaves += 1;
            auto ottId = nd->get_ott_id();
            auto taxoNode = taxonomy->get_data().get_node_by_ott_id(ottId);
            assert(taxoNode != nullptr);
            if (!contains(includedNodes, taxoNode)) {
                includedNodes.insert(taxoNode);
                insert_ancestors_to_paraphyletic_set(taxoNode, includedNodes);
            }
            if (!taxoNode->is_tip()) {
                if (storeLogInfo) {
                    exemplifedTaxonToTreeNamesForJSONLog[ottId].push_back(raw->get_name());
                }
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


bool handleExport(OTCLI & otCLI, const std::string &narg);
bool handleStdout(OTCLI & otCLI, const std::string &narg);
bool handleJSONOutput(OTCLI & otCLI, const std::string &narg);
bool handleExtinctJSON(OTCLI & , const std::string &);
bool handleExtinctJSONInline(OTCLI & , const std::string &) ;

bool handleStdout(OTCLI & otCLI, const std::string &) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    proc->useStdOut = true;
    return true;
}

bool handleExport(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("output filepath after the -e argument.");
    }
    proc->exportTreeFile = narg;
    return true;
}


bool handleJSONOutput(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath of output JSON after the -j argument.");
    }
    proc->outputJSONLogFile = narg;
    proc->storeLogInfo = true;
    return true;
}


bool handleExtinctJSON(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath to the exinction ID JSON file  the -t argument.");
    }
    proc->extinctInputJSON = narg;
    return true;
}

bool handleExtinctJSONInline(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting JSON content inline after the -i argument.");
    }
    proc->extinctInputJSONInline = narg;
    return true;
}


int main(int argc, char *argv[]) {
    const char * helpMsg =  "This tool treats the fossil taxa in the input taxonomic tree as poorly placed. " \
        "It will produce a set of edits (encoded as JSON) that can be applied to the taxonomy. If these " \
        "edits were applied to the taxonomy, the none of the fossil taxa included in the phylogenetic trees " \
        "would increase the set of taxa that are contested by input trees.\n" \
        "Input is a full taxonomy tree some number of input trees, and a (JSON-formatted) list of taxon IDs that are extinct.\n" \
        "Extinct taxon info can be given as a file name (-t flag) or on the command-line (-i flag)\n" \
        "The output JSON can be specified as stdout (-o flag) or as filename (-j flag)\n";
        
    OTCLI otCLI("otc-move-extinct-higher-to-avoid-contesting-taxa",
                helpMsg,
                "-jedits.json -textinct.json taxonomy.tre inp1.tre inp2.tre");
    MoveExtinctHigherState proc;
    otCLI.add_flag('o',
                  " requests that standard output stream, rather than the export directory, be used for all output.",
                  handleStdout,
                  false);
    otCLI.add_flag('j',
                  "Name of an output JSON file that summarizes the set of IDs used to exemplify each taxon.",
                  handleJSONOutput,
                  true);
    otCLI.add_flag('t',
                  "Name of input JSON file listing the exinct OTT IDs",
                  handleExtinctJSON,
                  true);
    otCLI.add_flag('i',
                  "exinct OTT IDs as JSON argument, not the filepath to JSON",
                  handleExtinctJSONInline,
                  true);
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, false);
}
