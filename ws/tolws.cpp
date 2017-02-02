#include <memory>
#include <thread>
#include <cstdlib>
#include <restbed>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem/operations.hpp>

#include "otc/util.h"
#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

#include "json.hpp"

using namespace std;
using namespace restbed;
using namespace otc;
namespace fs = boost::filesystem;
using json = nlohmann::json;

typedef std::set<fs::path> fp_set;
typedef std::pair<bool, fp_set > bool_fp_set; 

using TaxTree_t = RootedTree<RTTaxNodeData, RTreeNoData>;
using SummaryTree_t = RootedTree<RTNodeNoData, RTreeNoData>;
namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;


variables_map parse_cmd_line(int argc, char* argv[]) { 
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description output("Server options");
    output.add_options()
        ("tree-dir,D",value<string>(),"Filepath to directory that will hold synthetic tree output")
        ("port,P",value<long>(),"Port to bind to.")
        ("num-threads,t",value<long>(),"number of threads")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-tol-ws <taxonomy-dir> -D<dir> [OPTIONS]\n"
                                                    "Load taxonomy, check dir for synth tree outputs, and serve",
                                                    visible, invisible, p);

    return vm;
}

bool_fp_set get_subdirs(const fs::path & dirname) {
    if (!fs::is_directory(dirname)) {
        LOG(ERROR) << "\"" << dirname << "\" is not a directory.\n";
        return bool_fp_set(false, fp_set());
    }
    fp_set fps;
    for (auto sub : fs::directory_iterator(dirname)) {
        if (fs::is_directory(sub)) {
            fps.insert(sub);
        }
        
    }
    return bool_fp_set(true, fps);
}


struct SourceTreeId {
    std::string tree_id;
    std::string study_id;
    std::string git_sha;
};

struct SummaryTreeAnnotation {
    public:
        bool initialized = false;
        std::string date_completed;
        std::string filtered_flags;
        std::vector<std::string> filtered_flags_vec;
        std::string generated_by;
        OttId num_leaves_in_exemplified_taxonomy;
        OttId num_source_studies;
        OttId num_source_trees;
        OttId num_tips;
        OttId root_ott_id;
        std::string root_taxon_name;
        //source_id_map
        //sources
        std::string synth_id;
        std::string taxonomy_version;
        std::string tree_id;
        std::map<std::string, SourceTreeId> source_id_map;
        std::vector<std::string> sources;
        json full_source_id_map_json;

};

void from_json(const json &j, SummaryTreeAnnotation & sta);
void from_json(const json &j, SourceTreeId & sti);
void to_json(const json &j, SourceTreeId & sti);

class TreesToServe {
        list< SummaryTreeAnnotation> annotation_list;
        list< unique_ptr<SummaryTree_t> > tree_list;
        map<string, const SummaryTree_t *> id_to_tree;
        map<string, const SummaryTreeAnnotation *> id_to_annotations;
        string default_synth_id;
        const Taxonomy * taxonomy_ptr = nullptr;
        OttIdSet ott_id_set;
        unique_ptr<TaxTree_t> taxonomy_tree;

    public:
        void setTaxonomy(const Taxonomy &taxonomy) {
            assert(taxonomy_ptr == nullptr);
            taxonomy_ptr = &taxonomy;
            ott_id_set.clear();
            auto nodeNamer = [](const auto&){return string();};
            taxonomy_tree = taxonomy.getTree<TaxTree_t>(nodeNamer);
        }
        const Taxonomy & getTaxonomy() const {
            assert(taxonomy_ptr != nullptr);
            return *taxonomy_ptr;
        }
        void fillOttIdSet(const std::bitset<32> & flags) {
            ott_id_set.clear();
            for (const auto nd : iter_node_const(*taxonomy_tree)) {
                const auto & tax_record = nd->getData().taxonomy_line;
                auto intersection = flags & tax_record->flags;
                if (!intersection.any()) {
                    ott_id_set.insert(tax_record->id);
                }
            }
        }
        pair<SummaryTree_t &, SummaryTreeAnnotation &> getNewTreeAndAnnotations(const string & configfilename,
                                                                                const string & filename) {
            
            auto cleaning_flags = cleaning_flags_from_config_file(configfilename);
            fillOttIdSet(cleaning_flags);

            assert(taxonomy_ptr != nullptr);
            ParsingRules parsingRules;
            parsingRules.ottIdValidator = &ott_id_set;
            parsingRules.includeInternalNodesInDesIdSets = true;
            parsingRules.setOttIdForInternals = true;
            parsingRules.requireOttIds = true;
            parsingRules.setOttIds = true;
            std::ifstream inp;
            if (!openUTF8File(filename, inp)) {
                throw OTCError("Could not open \"" + filename + "\"");
            }
            LOG(INFO) << "reading \"" << filename << "\"...";
            ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
            FilePosStruct pos(filenamePtr);
            std::unique_ptr<SummaryTree_t> nt = readNextNewick<SummaryTree_t>(inp, pos, parsingRules);
            tree_list.push_back(move(nt));
            annotation_list.push_back(SummaryTreeAnnotation());
            return {*(tree_list.back()), annotation_list.back()};
        }
        void registerLastTreeAndAnnotations() {
            const SummaryTreeAnnotation & sta = annotation_list.back();
            const SummaryTree_t & tree = *(tree_list.back());
            default_synth_id = sta.synth_id; // @ TODO need a better system for deciding the default synth ID.
            id_to_tree[sta.synth_id] = &tree;
            id_to_annotations[sta.synth_id] = &sta;
        }
        void freeLastTreeAndAnnotations() {
            tree_list.back()->clear();
            tree_list.pop_back();
            annotation_list.pop_back();
        }

        const SummaryTreeAnnotation * getAnnotations(string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_annotations.find(key);
            return mit == id_to_annotations.end() ? nullptr : mit->second;
        }

        const SummaryTree_t * getSummaryTree(string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_tree.find(key);
            return mit == id_to_tree.end() ? nullptr : mit->second;
        }
        size_t getNumTrees() const {
            return id_to_tree.size();
        }
};

bool read_tree_and_annotations(const fs::path & configpath,
                               const fs::path & treepath,
                               const fs::path & annotationspath,
                               TreesToServe & tts);

// Globals. TODO: lock if we read twice
fp_set checked_dirs;
fp_set known_tree_dirs;

bool read_trees(const fs::path & dirname, TreesToServe & tts) {
    auto sdr = get_subdirs(dirname);
    if (!sdr.first) {
        return false;
    }
    const auto & subdir_set = sdr.second;
    for (auto p : subdir_set) {
        if (!contains(checked_dirs, p)) {
            checked_dirs.insert(p);
            fs::path configpath = p;
            configpath /= "config";
            fs::path treepath = p;
            treepath /= "labelled_supertree";
            treepath /= "labelled_supertree.tre";
            fs::path annotationspath = p;
            annotationspath /= "annotated_supertree";
            annotationspath /= "annotations.json";
            bool was_tree_par = false;
            try {
                if (fs::is_regular_file(treepath)
                    && fs::is_regular_file(annotationspath)
                    && fs::is_regular_file(configpath)) {
                    if (read_tree_and_annotations(configpath, treepath, annotationspath, tts)) {
                        known_tree_dirs.insert(p);
                        was_tree_par = true;
                    }
                }
            } catch (const std::exception & x) {
                LOG(WARNING) << "Exception while reading summary tree directory:\n   ";
                LOG(WARNING) << x.what() << '\n';
            }
            if (!was_tree_par) {
                LOG(WARNING) << "Rejected \"" << p << "\" due to lack of " << treepath << " or lack of " << annotationspath << " or parsing error.\n";
            }
        }
    }
    return true;
}


bool read_tree_and_annotations(const fs::path & config_path,
                               const fs::path & tree_path,
                               const fs::path & annotations_path,
                               TreesToServe & tts) {
    std::string annot_str = annotations_path.native();
    std::ifstream annotations_stream(annot_str.c_str());
    json annotations_obj;
    try {
        annotations_stream >> annotations_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << annotations_path << "\" as JSON.\n";
        throw;
    }
    auto tree_and_ann = tts.getNewTreeAndAnnotations(config_path.native(), tree_path.native());
    try {
        SummaryTree_t & tree = tree_and_ann.first;
        SummaryTreeAnnotation & sta = tree_and_ann.second;
        sta = annotations_obj;
        sta.initialized = true;
        tts.registerLastTreeAndAnnotations();
    } catch (...) {
        tts.freeLastTreeAndAnnotations();
        throw;
    }
    return true;
}

std::string extract_string(const json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_string()) {
        return dc_el->get<std::string>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}

OttId extract_unsigned_long(const json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_number()) {
        return dc_el->get<OttId>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a non-negative integer.\n";
}

void from_json(const json &j, SummaryTreeAnnotation & sta) {
    sta.date_completed = extract_string(j, "date_completed");
    sta.filtered_flags = extract_string(j, "filtered_flags");
    auto splitff = split_string(sta.filtered_flags, ',');
    sta.filtered_flags_vec.assign(splitff.begin(), splitff.end());
    // generated_by gets converted back to a string
    auto gb_el = j.find("generated_by");
    if (gb_el == j.end()) {
        throw OTCError() << "Missing generated_by field.\n";
    }
    sta.generated_by = gb_el->dump();
    sta.num_leaves_in_exemplified_taxonomy = extract_unsigned_long(j, "num_leaves_in_exemplified_taxonomy");
    sta.num_source_studies = extract_unsigned_long(j, "num_source_studies");
    sta.num_source_trees = extract_unsigned_long(j, "num_source_trees");
    sta.num_tips = extract_unsigned_long(j, "num_tips");
    sta.root_ott_id = extract_unsigned_long(j, "root_ott_id");
    sta.root_taxon_name = extract_string(j, "root_taxon_name");
    sta.synth_id = extract_string(j, "synth_id");
    sta.taxonomy_version = extract_string(j, "taxonomy_version");
    sta.tree_id = extract_string(j, "tree_id");
    auto sim_el = j.find("source_id_map");
    if (sim_el == j.end()) {
        throw OTCError() << "Missing source_id_map field.\n";
    }
    if (!sim_el->is_object()) {
        throw OTCError() << "Expected \"source_id_map\" field to be an object.\n";
    }
    try {
        for (json::const_iterator sim_it = sim_el->begin(); sim_it != sim_el->end(); ++sim_it) {
            sta.source_id_map[sim_it.key()] = sim_it.value(); 
        }
    } catch (OTCError & x) {
        throw OTCError() << "Error reading source_id_map field: " << x.what();
    }
    sta.full_source_id_map_json = *sim_el;
    auto s_el = j.find("sources");
    if (s_el == j.end()) {
        throw OTCError() << "Missing sources field.\n";
    }
    if (!s_el->is_array()) {
        throw OTCError() << "Expected \"sources\" field to be an array.\n";
    }
    sta.sources.resize(s_el->size());
    for (auto i = 0U; i < s_el->size(); ++i) {
        try {
            sta.sources[i] = s_el->at(i).get<std::string>();
        } catch (OTCError & x) {
            throw OTCError() << "Error expected each element of the sources array to be a string: " << x.what();
        }
    }
}

void to_json(json &j, SourceTreeId & sti) {
    j = json{{"git_sha", sti.git_sha}, {"study_id", sti.study_id}, {"tree_id", sti.tree_id}};
}

void from_json(const json &j, SourceTreeId & sti) {
    sti.git_sha = extract_string(j, "git_sha");
    sti.study_id = extract_string(j, "study_id");
    sti.tree_id = extract_string(j, "tree_id");
}


TreesToServe tts;
template<typename T>
bool extract_from_request(json & j, string opt_name, T & setting, string & response, int & status_code);

template<>
bool extract_from_request(json & j, string opt_name, bool & setting, string & response, int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_boolean()) {
            setting = opt->get<bool>();
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be a boolean.\n";
        status_code = 400;
    }
    return false;
}

template<>
bool extract_from_request(json & j, string opt_name, string & setting, string & response, int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_boolean()) {
            setting = opt->get<string>();
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be a string.\n";
        status_code = 400;
    }
    return false;
}


void about_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     bool include_sources,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["date_created"] = sta->date_completed;
    response["num_source_trees"] = sta->num_source_trees;
    response["taxonomy_version"] = sta->taxonomy_version;
    response["filtered_flags"] = sta->filtered_flags_vec;
    response["synth_id"] = sta->synth_id;
    if (include_sources) {
        response["source_id_map"] = sta->full_source_id_map_json;
        response["sources"] = sta->sources;
    }
    json root;
    auto root_node = tree_ptr->getRoot();
    auto root_id = root_node->getOttId();
    const auto & root_taxon = taxonomy.record_from_id(root_id);
    root["node_id"] = root_node->getName();
    json taxon;
    taxon["tax_sources"] = root_taxon.sourceinfoAsVec();
    auto un = string(root_taxon.uniqname);
    auto n = string(root_taxon.name);
    taxon["name"] = n;
    taxon["uniqname"] = un;
    auto r = string(root_taxon.rank);
    taxon["rank"] = r;
    taxon["ott_id"] = root_taxon.id;
    
    root["taxon"] = taxon;
    response["root"] = root;
    response_str = response.dump(1);
}

void about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );

    size_t content_length = request->get_header( "Content-Length", 0 );

    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        stringstream id;
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        try {
            if (!body.empty()) {
                parsedargs = json::parse(body);
            }
        } catch (...) {
            rbody = "Could not parse body of call as JSON.\n";
            status_code = 400;
        }
        bool include_sources = false;
        string synth_id;
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_source_list", include_sources, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "synth_id", synth_id, rbody, status_code);
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.getAnnotations(synth_id);
            const SummaryTree_t * treeptr = tts.getSummaryTree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                about_ws_method(tts, treeptr, sta, include_sources, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}




int main( const int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
        variables_map args = parse_cmd_line(argc,argv);

        int num_threads = 4;
        int port_number = 1984;
        if (args.count("num-threads")) {
            port_number = args["num-threads"].as<int>();
        }
        if (args.count("port")) {
            port_number = args["port"].as<int>();
        }
        if (!args.count("tree-dir")) {
            cerr << "Expecting a tree-dir argument for a path to a directory of synth outputs.\n";
            return 1;
        }
        
        const fs::path topdir{args["tree-dir"].as<string>()};
        // Must load taxonomy before trees
        auto taxonomy = load_taxonomy(args);
        tts.setTaxonomy(taxonomy);
        if (!read_trees(topdir, tts)) {
            return 2;
        }
        //
        if (tts.getNumTrees() == 0) {
            cerr << "No tree to serve. Exiting...\n";
            return 3;
        }
        ////// ROUTES
        auto resource = make_shared< Resource >( );
        resource->set_path( "/tree_of_life/about" );
        resource->set_method_handler( "POST", about_method_handler );
        /////  SETTINGS
        auto settings = make_shared< Settings >( );
        settings->set_port( port_number );
        settings->set_worker_limit( num_threads );
        settings->set_default_header( "Connection", "close" );
        
        Service service;
        service.publish( resource );
        LOG(INFO) << "starting service with " << num_threads << " on port " << port_number << "...\n";
        service.start( settings );
        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        LOG(ERROR) <<"otc-tol-ws: Error! " << e.what() << std::endl;
        exit(1);
    }

}