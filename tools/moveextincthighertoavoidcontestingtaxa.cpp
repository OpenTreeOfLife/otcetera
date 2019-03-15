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
    std::size_t phylogeniesProcessed = 0;
    std::vector<std::string> namesOfTreesWithoutExtinct;
    
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
            RTSplits & d = nd->get_data();
            auto p = nd->get_parent();
            if (p == nullptr) {
                d.depth = 0;
            } else {
                d.depth = 1 + p->get_data().depth;
            }
            //db_write_ott_id_set(nd->get_name().c_str(), d.des_ids);
            std::cerr << "nd id = " << nd->get_ott_id() << " depth = " << d.depth << '\n';
        }
        return r;
    }
    
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
                const RTSplits & d = nd->get_data();
                LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " has extinct tip with ID=" << ottId << ".\n";
            }
        }
        if (extinctTips.empty()) {
            LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " contains no extinct tips. flushing...\n";
            namesOfTreesWithoutExtinct.push_back(raw->get_name());
            return true;
        }
        RTreeOttIDMapping<RTSplits> & taxdata = taxonomy->get_data();
        auto phylo_root = raw->get_root();
        const RTSplits & phylo_root_data = phylo_root->get_data();
        const auto & phylo_tip_ids = phylo_root_data.des_ids;
        db_write_ott_id_set("phylo root des_ids", phylo_tip_ids);
        for (auto nd : iter_leaf(*raw)) {
            auto ottId = nd->get_ott_id();
            auto & nd_des_ids = nd->get_data().des_ids;
            auto tax_nd_for_tip = taxdata.ott_id_to_node.at(ottId);
            auto & tax_nd_des_ids = tax_nd_for_tip->get_data().des_ids;
            auto tax_tip_des_overlap = set_intersection_as_set(tax_nd_des_ids, phylo_tip_ids);
            if (tax_tip_des_overlap.size() != 1) {
                db_write_ott_id_set("evidence of overlapping Ids: ", tax_tip_des_overlap);
                LOG(ERROR) << "Node \"" << nd->get_name() << "\" has taxonomic descendants that overlap with other tips in the tree.\n";
                throw OTCError("Non disjunct taxa as tips of tree");
            }
            if (contains(extinctIDSet, ottId)) {
                extinctTips.insert(nd);
                const RTSplits & d = nd->get_data();
                LOG(INFO) << "Tree " << treeIndex << " name=" <<  raw->get_name() << " has extinct tip with ID=" << ottId << ".\n";
            }
        }
        
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
