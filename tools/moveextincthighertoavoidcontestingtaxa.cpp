#include "otc/otcli.h"
#include "otc/supertree_util.h"
#include "json.hpp"
#include <sstream>
using namespace otc;
using json = nlohmann::json;



struct MoveExtinctHigherState : public TaxonomyDependentTreeProcessor<TreeMappedWithSplits> {
    int numErrors;
    bool useStdOut;
    bool useCmdLine;
    std::string extinctInputJSONFilepath;
    std::string extinctInputJSONInline;
    std::string outputJSONLogFile;
    std::istream * extinctInStream;
    std::ostream * jsonOutStream;
    std::ofstream jsonOutFileIfUsed;
    std::set<OttId> extinctIDSet;
    std::map<std::unique_ptr<TreeMappedWithSplits>, std::size_t> inputTreesToIndex;
    //std::vector<TreeMappedWithSplits *> treePtrByIndex;
    std::size_t phylogeniesProcessed = 0;
    std::vector<std::string> namesOfTreesWithoutExtinct;
    /*
    std::set<const RootedTreeNodeNoData *> includedNodes;
    using TreeNdPair = std::pair<TreeMappedWithSplits *, RootedTreeNodeNoData *>;
    using ListTreeNdPair = std::list<TreeNdPair>;
    std::map<RootedTreeNodeNoData *, ListTreeNdPair> nonTermToMappedPhylo;
    std::string exportTreeFile;
    std::ofstream nonEmptyFileStream;
    std::string outputNonEmptyTreeOutput;
    std::string outputJSONLogFile;
    std::string extinctInputJSONInline;
    bool storeLogInfo = false;
    std::map<OttId, OttIdSet> exemplificationsForJSONLog;
    std::map<OttId, std::list<std::string> > exemplifedTaxonToTreeNamesForJSONLog;
    */
    virtual ~MoveExtinctHigherState(){}
    MoveExtinctHigherState()
        :numErrors(0),
        useStdOut(true),
        useCmdLine(false) {
    }
    bool readExtinctIDsFromStream(OTCLI & otCLI, std::istream &inp) {
        json j;
        try {
            inp >> j;
        } catch (json::parse_error & x) {
            LOG(ERROR) << "Could not parse input JSON list of IDs.:\n" << x.what() << " at byte " << x.byte ;
            return false;
        }
        if (j.type() == json::value_t::null) {
            //pass
        } else if (j.type() == json::value_t::array) {
            unsigned index = 0;
            for (auto i : j) {
                auto itype = i.type();
                if (itype != json::value_t::number_integer && itype != json::value_t::number_unsigned) {
                    LOG(ERROR) << "Could not parse all elements of input list of as integer:\n" << i ;
                    return false;
                }
                unsigned val = i.get<OttId>();
                std::cerr << "Element[" << index++ << "] = " << val << '\n';
                extinctIDSet.insert(val);
            }

        } else {
            LOG(ERROR) << "Expecting an array of integers as the JSON content.";
        }
        return true;
    }

    bool readExtinctIDs(OTCLI & otCLI) {
        std::ifstream extinctJSONStreamLocal;
        std::istringstream extinctJSONStringWrapper;
        if (useCmdLine) {
            if (!extinctInputJSONFilepath.empty()) {
                LOG(ERROR) << "Cannot use command line JSON and input JSON file options at the same time!";
                return false;
            }
            extinctJSONStringWrapper.str(extinctInputJSONInline);
            extinctInStream = &extinctJSONStringWrapper;
        } else {
            if (extinctInputJSONFilepath.empty()) {
                LOG(ERROR) << "Must use either the command line JSON or input JSON file options!";
                return false;
            }
            extinctJSONStreamLocal.open(extinctInputJSONFilepath);
            if (!extinctJSONStreamLocal.good()) {
                LOG(ERROR) << "Could not open input extinct JSON file \"" << extinctInputJSONFilepath << "\"";
                return false;    
            }
            extinctInStream = &extinctJSONStreamLocal;
        }
        if (useStdOut) {
            jsonOutStream = &std::cout;
        } else {
            jsonOutFileIfUsed.open(outputJSONLogFile);
            if (!jsonOutFileIfUsed.good()) {
                jsonOutFileIfUsed.close();
                LOG(ERROR) << "Could not open JSON output path at \"" << outputJSONLogFile << "\"";
                return false;
            }
            jsonOutStream = &jsonOutFileIfUsed;
        }
        return readExtinctIDsFromStream(otCLI, *extinctInStream);
    }

    bool pretree_read_hook(OTCLI & otCLI) override {
        return readExtinctIDs(otCLI);
    }
    
    bool process_taxonomy_tree(OTCLI & otCLI) override {
        bool r = TaxonomyDependentTreeProcessor<TreeMappedWithSplits>::process_taxonomy_tree(otCLI);
        // we can ignore the internal node labels for the non-taxonomic trees
        otCLI.get_parsing_rules().set_ott_idForInternals = false;

        for (auto nd : iter_pre(*taxonomy)) {
            std::cerr << "nd id = " << nd->get_ott_id() << '\n';
        }
        return r;
    }
    /*
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
                TreeMappedWithSplits * treeP = treeNdPair.first;
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
    */
    
    bool summarize(OTCLI &otCLI) override {
        json document;
        json jTreesWoExt = json::array(); 
        for (auto i : namesOfTreesWithoutExtinct) {
            jTreesWoExt.push_back(i);
        }
        document["trees_with_no_extinct_tips"] = jTreesWoExt;
        assert(jsonOutStream != nullptr);
        *jsonOutStream << document << std::endl;
        if (jsonOutStream == & jsonOutFileIfUsed) {
            jsonOutFileIfUsed.close();
            jsonOutStream = nullptr;  
        }
        return true;
    }
 

    bool process_source_tree(OTCLI & otCLI, std::unique_ptr<TreeMappedWithSplits> treeup) override {
        std::size_t treeIndex = phylogeniesProcessed++;
        TreeMappedWithSplits * raw = treeup.get();
        raw->set_name(otCLI.currentFilename);
        std::set<const TreeMappedWithSplits::node_type *> extinctTips; 
        for (auto nd : iter_leaf(*raw)) {
            auto ottId = nd->get_ott_id();
            if (contains(extinctIDSet, ottId)) {
                extinctTips.insert(nd);
                LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " has extinct tip with ID=" << ottId << ".\n";
            }
        }
        if (extinctTips.empty()) {
            LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " contains no extinct tips. flushing...\n";
            namesOfTreesWithoutExtinct.push_back(raw->get_name());
            return true;
        }
        inputTreesToIndex[std::move(treeup)] = treeIndex;
        

        /*
        assert(treeup != nullptr);
        assert(taxonomy != nullptr);
        // Store the tree pointer with a map to its index, and an alias for fast index->tree.
        std::size_t treeIndex = inputTreesToIndex.size();
        assert(treeIndex == treePtrByIndex.size());
        TreeMappedWithSplits * raw = treeup.get();
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
        */
        return true;
    }
};


bool handleJSONOutput(OTCLI & otCLI, const std::string &narg);
bool handleExtinctJSON(OTCLI & , const std::string &);
bool handleExtinctJSONInline(OTCLI & , const std::string &) ;


bool handleJSONOutput(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath of output JSON after the -j argument.");
    }
    proc->outputJSONLogFile = narg;
    proc->useStdOut = false;
    return true;
}


bool handleExtinctJSON(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting filepath to the exinction ID JSON file  the -t argument.");
    }
    proc->extinctInputJSONFilepath = narg;
    return true;
}

bool handleExtinctJSONInline(OTCLI & otCLI, const std::string & narg) {
    MoveExtinctHigherState * proc = static_cast<MoveExtinctHigherState *>(otCLI.blob);
    assert(proc != nullptr);
    if (narg.empty()) {
        throw OTCError("Expecting JSON content inline after the -i argument.");
    }
    proc->extinctInputJSONInline = narg;
    proc->useCmdLine = true;
    return true;
}


int main(int argc, char *argv[]) {
    const char * helpMsg =  "This tool treats the fossil taxa in the input taxonomic tree as poorly placed. " \
        "It will produce a set of edits (encoded as JSON) that can be applied to the taxonomy. If these " \
        "edits were applied to the taxonomy, the none of the fossil taxa included in the phylogenetic trees " \
        "would increase the set of taxa that are contested by input trees.\n" \
        "Input is a full taxonomy tree some number of input trees, and a (JSON-formatted) list of taxon IDs that are extinct.\n" \
        "Extinct taxon info can be given as a file name (-t flag) or on the command-line (-i flag)\n" \
        "The output JSON can be specified as stdout (default) or as filename (-j flag)\n";
        
    OTCLI otCLI("otc-move-extinct-higher-to-avoid-contesting-taxa",
                helpMsg,
                "-jedits.json -textinct.json taxonomy.tre inp1.tre inp2.tre");
    MoveExtinctHigherState proc;
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
    return tax_dependent_tree_processing_main(otCLI, argc, argv, proc, 2, true);
}
