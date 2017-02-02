#include <memory>
#include <thread>
#include <cstdlib>
#include <restbed>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
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

using TaxTree_t = RootedTree<RTNodeNoData, RTreeNoData>;
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

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<long>(), "OTT id of root node of subtree to keep")
        ;

    options_description output("Server options");
    output.add_options()
        ("tree-dir,D",value<string>(),"Filepath to directory that will hold synthetic tree output")
        ("port,P",value<long>(),"Port to bind to.")
        ("num-threads,t",value<long>(),"number of threads")
        ;

    options_description visible;
    visible.add(taxonomy).add(output).add(otc::standard_options());

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
        cerr << "\"" << dirname << "\" is not a directory.\n";
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
    public:
        SummaryTreeAnnotation sta;
};

bool read_tree_and_annotations(const fs::path & treepath, const fs::path & annotationspath, TreesToServe & tts);

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
            fs::path treepath = p;
            treepath /= "labelled_supertree";
            treepath /= "labelled_supertree.tre";
            fs::path annotationspath = p;
            annotationspath /= "annotated_supertree";
            annotationspath /= "annotations.json";
            bool was_tree_par = false;
            try {
                if (fs::is_regular_file(treepath) && fs::is_regular_file(annotationspath)) {
                    if (read_tree_and_annotations(treepath, annotationspath, tts)) {
                        known_tree_dirs.insert(p);
                        was_tree_par = true;
                    }
                }
            } catch (const std::exception & x) {
                std::cerr << "Exception while reading summary tree directory:\n   ";
                std::cerr << x.what() << '\n';
            }
            if (!was_tree_par) {
                std::cerr << "Rejected \"" << p << "\" due to lack of " << treepath << " or lack of " << annotationspath << " or parsing error.\n";
            }
        }
    }
    return true;
}


bool read_tree_and_annotations(const fs::path & tree_path, const fs::path & annotations_path, TreesToServe & tts) {
    std::string annot_str = annotations_path.native();
    std::ifstream annotations_stream(annot_str.c_str());
    json annotations_obj;
    try {
        annotations_stream >> annotations_obj;
    } catch (...) {
        cerr << "Could not read \"" << annotations_path << "\" as JSON.\n";
        throw;
    }
    tts.sta = annotations_obj;
    tts.sta.initialized = true;
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

void about_method_handler( const shared_ptr< Session > session )
{
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
        if (status_code == OK) {
            auto opt = parsedargs.find("include_source_list");
            if (opt != parsedargs.end()) {
                if (opt->is_boolean()) {
                    include_sources = opt->get<bool>();
                } else {
                    rbody = "Expecting include_source_list to be a boolean.\n";
                    status_code = 400;
                }
            }
        }
        if (status_code == OK) {
            SummaryTreeAnnotation & sta = tts.sta;
            json response;
            response["date_created"] = sta.date_completed;
            response["num_source_trees"] = sta.num_source_trees;
            response["taxonomy_version"] = sta.taxonomy_version;
            response["filtered_flags"] = sta.filtered_flags_vec;
            response["synth_id"] = sta.synth_id;
            if (include_sources) {
                response["source_id_map"] = sta.full_source_id_map_json;
                response["sources"] = sta.sources;
            }
            rbody = response.dump(1);
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
        auto taxonomy = load_taxonomy(args);

        if (!read_trees(topdir, tts)) {
            return 2;
        }
        if (!tts.sta.initialized) {
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
        std::cerr << "starting service with " << num_threads << " on port " << port_number << "...\n";
        service.start( settings );
        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        cerr<<"otc-tol-ws: Error! " << e.what() << std::endl;
        exit(1);
    }

}