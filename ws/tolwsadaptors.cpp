#include <thread>
// PID reporting from https://github.com/Corvusoft/restbed/blob/master/example/signal_handling/source/example.cpp
#ifdef _WIN32
    #include <process.h>
#else
    #include <unistd.h>
#endif
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <restbed>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "ws/tolwsadaptors.h"
#include "otc/otcli.h"

// unlike most headers, we'll go ahead an use namespaces
//    because this is an implementation file
using namespace std;
namespace po = boost::program_options;
using namespace otc;
namespace fs = boost::filesystem;
using json = nlohmann::json;
using namespace restbed;

namespace otc {
    // global
    TreesToServe tts;
}

bool parse_body_or_err(const Bytes & body, json & parsedargs, std::string & rbody, int & status_code) {
    try {
        if (!body.empty()) {
            parsedargs = json::parse(body);
            return true;
        }
    } catch (...) {
        rbody = "Could not parse body of call as JSON.\n";
        status_code = 400;
    }
    return false;
}

///////////////////////
// handlers that are registered as callback

void about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        bool include_sources = false;
        string synth_id;
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_source_list", include_sources, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "synth_id", synth_id, rbody, status_code);
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.get_annotations(synth_id);
            const SummaryTree_t * treeptr = tts.get_summary_tree(synth_id);
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


void get_synth_and_node_id(const json &j, string & synth_id, string & node_id, string & rbody, int & status_code) {
    if (status_code == OK) {
        extract_from_request(j, "synth_id", synth_id, rbody, status_code);
    }
    if (status_code == OK) {
        if (!extract_from_request(j, "node_id", node_id, rbody, status_code)) {
            rbody = "node_id is required.\n";
            status_code = 400;
        }
    }
}

void get_synth_and_node_id_vec(const json &j, string & synth_id, vector<string> & vec_node_ids, string & rbody, int & status_code) {
    if (status_code == OK) {
        extract_from_request(j, "synth_id", synth_id, rbody, status_code);
    }
    if (status_code == OK) {
        if (!extract_from_request(j, "node_ids", vec_node_ids, rbody, status_code)) {
            rbody = "node_ids is required.\n";
            status_code = 400;
        }
    }
}


void get_label_format(const json &j, NodeNameStyle & nns, string & rbody, int & status_code) {
    string label_format = "name_and_id";
    if (status_code == OK) {
        if (extract_from_request(j, "label_format", label_format, rbody, status_code)) {
            if (label_format != "name_and_id"
                && label_format != "id"
                && label_format != "name") {
                rbody = "label_format must be \"name_and_id\", \"name\" or \"id\".\n";
                status_code = 400;
            }
            if (label_format == "id") {
                nns = NodeNameStyle::NNS_ID_ONLY;
            } else if (label_format == "name") {
                nns = NodeNameStyle::NNS_NAME_ONLY;
            }
        }
    }
}

void node_info_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        bool include_lineage = false;
        string synth_id;
        string node_id;
        get_synth_and_node_id(parsedargs, synth_id, node_id, rbody, status_code);
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_lineage", include_lineage, rbody, status_code);
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.get_annotations(synth_id);
            const SummaryTree_t * treeptr = tts.get_summary_tree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                node_info_ws_method(tts, treeptr, sta, node_id, include_lineage, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

void mrca_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        string synth_id;
        vector<string> node_id_vec;
        get_synth_and_node_id_vec(parsedargs, synth_id, node_id_vec, rbody, status_code);
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.get_annotations(synth_id);
            const SummaryTree_t * treeptr = tts.get_summary_tree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                mrca_ws_method(tts, treeptr, sta, node_id_vec, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

void subtree_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        string synth_id;
        string node_id;
        string format = "newick"; // : (string) Defines the tree format; either "newick" or "arguson"; default="newick"
        int height_limit = -1; // :
        get_synth_and_node_id(parsedargs, synth_id, node_id, rbody, status_code);
        if (status_code == OK) {
            if (extract_from_request(parsedargs, "format", format, rbody, status_code)) {
                if (format != "newick" && format != "arguson") {
                    rbody = "format must be \"newick\" or \"arguson\".\n";
                    status_code = 400;
                } else if (format == "arguson") {
                    height_limit = 3;
                }
            }
        }
        NodeNameStyle nns = NodeNameStyle::NNS_NAME_AND_ID;
        if (status_code == OK && format == "newick") {
            get_label_format(parsedargs, nns, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "height_limit", height_limit, rbody, status_code);
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.get_annotations(synth_id);
            const SummaryTree_t * treeptr = tts.get_summary_tree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else if (format == "newick") {
                newick_subtree_ws_method(tts, treeptr, sta,
                                         node_id, nns, height_limit,
                                         rbody, status_code);
            } else {
                arguson_subtree_ws_method(tts, treeptr, sta,
                                         node_id, height_limit,
                                         rbody, status_code);    
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}


void induced_subtree_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        string synth_id;
        vector<string> node_id_vec;
        get_synth_and_node_id_vec(parsedargs, synth_id, node_id_vec, rbody, status_code);
        NodeNameStyle nns = NodeNameStyle::NNS_NAME_AND_ID;
        get_label_format(parsedargs, nns, rbody, status_code);
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.get_annotations(synth_id);
            const SummaryTree_t * treeptr = tts.get_summary_tree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                induced_subtree_ws_method(tts, treeptr, sta,
                                         node_id_vec, nns,
                                         rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

void tax_about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        if (status_code == OK) {
            tax_about_ws_method(tts, rbody, status_code);
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

// looks for ott_id or source_id args to find a node
const RTRichTaxNode * extract_taxon_node_from_args(const json & parsedargs,
                                                   std::string & rbody,
                                                   int & status_code) {
    if (status_code != OK) {
        return nullptr;
    }
    OttId ott_id;
    string source_id;
    bool supplied_ott_id = false;
    bool supplied_source_id = false;
    const auto & taxonomy = tts.get_taxonomy();
    const auto & taxonomy_tree = taxonomy.getTaxTree();
    const auto & taxonomy_tree_data = taxonomy_tree.get_data();
    supplied_ott_id = extract_from_request(parsedargs, "ott_id", ott_id, rbody, status_code);
    if (status_code != OK) {
        return nullptr;
    }
    supplied_source_id = extract_from_request(parsedargs, "source_id", source_id, rbody, status_code);
    if (status_code != OK) {
        return nullptr;
    }
    if (supplied_ott_id && supplied_source_id) {
        rbody = "ott_id or source_id arguments cannot both be supplied.";
        status_code = 400;
        return nullptr;
    }
    if ((!supplied_source_id) && (!supplied_ott_id)) {
        rbody = "An ott_id or source_id argument is required.";
        status_code = 400;
        return nullptr;
    }
    if (supplied_source_id) {
        auto pref_id = split_string(source_id, ':');
        if (pref_id.size() != 2) {
            rbody = "Expecting exactly 1 colon in a source ID string. Found: \"";
            rbody += source_id;
            rbody += "\".";
            status_code = 400;
            return nullptr;
        }
        string source_prefix = *pref_id.begin();
        if (indexed_source_prefixes.count(source_prefix) == 0) {
            rbody = "IDs from source ";
            rbody += source_prefix ;
            rbody += " are not known or not indexed for searching.";
            status_code = 400;
            return nullptr;
        }
        string id_str = *pref_id.rbegin();
        try {
            std::size_t pos;
            long raw_foreign_id  = std::stol(id_str.c_str(), &pos);
            if (pos < id_str.length()) {
                rbody = "Expecting the ID portion of the source_id to be numeric. Found: ";
                rbody += id_str;
                status_code = 400;
                return nullptr;
            }
            OttId foreign_id;
            try {
                foreign_id = check_ott_id_size(raw_foreign_id);
            } catch (OTCError &) {
                rbody = "The ID portion of the source_id was too large. Found: ";
                rbody += id_str;
                status_code = 400;
                return nullptr;
            }
            try {
                if (source_prefix == "ncbi") {
                    return taxonomy_tree_data.ncbi_id_map.at(foreign_id);
                } else if (source_prefix == "gbif") {
                    return taxonomy_tree_data.gbif_id_map.at(foreign_id);
                } else if (source_prefix == "worms") {
                    return taxonomy_tree_data.worms_id_map.at(foreign_id);
                } else if (source_prefix == "if") {
                    return taxonomy_tree_data.if_id_map.at(foreign_id);
                } else if (source_prefix == "irmng") {
                    return taxonomy_tree_data.irmng_id_map.at(foreign_id);
                } else {
                    assert(false);
                }
            } catch (std::out_of_range & x) {
                rbody = "No taxon in the taxonomy is associated with source_id of ";
                rbody += source_id;
                status_code = 400;
            }
        } catch (...) {
            rbody = "Expecting the ID portion of the source_id to be numeric. Found: ";
            rbody += id_str;
            status_code = 400;
        }   
        return nullptr;
    }
    const RTRichTaxNode * taxon_node = nullptr;
    if (status_code == OK && supplied_ott_id) {
        taxon_node = taxonomy.taxon_from_id(ott_id);
        if (taxon_node == nullptr) {
            rbody = "Unrecognized OTT ID: ";
            rbody += to_string(ott_id);
            status_code = 400;
        }
    }
    return taxon_node;
}

void taxon_info_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        bool include_lineage = false;
        bool include_children = false;
        bool include_terminal_descendants = false;
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_lineage", include_lineage, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_children", include_children, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_terminal_descendants", include_terminal_descendants, rbody, status_code);
        }
        const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, rbody, status_code);
        if (status_code == OK) {
            assert(taxon_node != nullptr);
            taxon_info_ws_method(tts, taxon_node, include_lineage, include_children, include_terminal_descendants, rbody, status_code);
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}
        
void taxon_mrca_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        string ott_version;
        OttIdSet ott_id_set;
        if (!extract_from_request(parsedargs, "ott_ids", ott_id_set, rbody, status_code)) {
            if (status_code == OK) {
                rbody = "The ott_ids argument is required.";
                status_code = 400;
            }
        }
        if (status_code == OK) {
            taxonomy_mrca_ws_method(tts, ott_id_set, rbody, status_code);
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

void taxon_subtree_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        parse_body_or_err(body, parsedargs, rbody, status_code);
        NodeNameStyle nns = NodeNameStyle::NNS_NAME_AND_ID;
        if (status_code == OK) {
            get_label_format(parsedargs, nns, rbody, status_code);
        }
        const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, rbody, status_code);
        if (status_code == OK) {
            assert(taxon_node != nullptr);
            taxon_subtree_ws_method(tts, taxon_node, nns, rbody, status_code);
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

/// End of method_handler. Start of global service related code
Service * global_service_ptr = nullptr;

void kill_service_worker() {
    Service * gp = global_service_ptr;
    global_service_ptr = nullptr;
    if (gp == nullptr) {
        return;
    }
    LOG(WARNING) <<  "Stopping service...";
    gp->stop();
    LOG(WARNING) <<  "Service stopped...";
}

void sigterm_handler( const int signal_number ) {
    LOG(WARNING) <<  "Received signal number " << signal_number;
    if (global_service_ptr != nullptr) {
        std::thread killer_thread(kill_service_worker);
        killer_thread.detach();
    }
}

static string pidfile;

void ready_handler( Service& ) {
#ifdef _WIN32
    auto pid = _getpid();
#else
    auto pid = getpid();
#endif
    LOG(INFO) << "Service is ready. PID is " << pid;
    if (!pidfile.empty()) {
        ofstream pstream(pidfile);
        if (pstream.good()) {
            pstream << pid << '\n';
            pstream.close();
            LOG(INFO) << "PID written to " << pidfile;
        } else {
            LOG(WARNING) << "PID file " << pidfile << " could not be opened!";    
        }
    } else {
        LOG(INFO) << "No PID file was requested.";
    }
}

int run_server(const po::variables_map & args) {
    time_t start_time;
    time(&start_time);
    int num_threads = 4;
    int port_number = 1984;
    if (args.count("num-threads")) {
        num_threads = args["num-threads"].as<int>();
    }
    if (args.count("port")) {
        port_number = args["port"].as<int>();
    }
    if (args.count("pidfile")) {
        pidfile = args["pidfile"].as<string>();
    }
    if (!args.count("tree-dir")) {
        cerr << "Expecting a tree-dir argument for a path to a directory of synth outputs.\n";
        return 1;
    }
    const fs::path topdir{args["tree-dir"].as<string>()};
    // Must load taxonomy before trees
    LOG(INFO) << "reading taxonomy...";
    RichTaxonomy taxonomy = std::move(load_rich_taxonomy(args));
    time_t post_tax_time;
    time(&post_tax_time);
    tts.set_taxonomy(taxonomy);
    if (!read_trees(topdir, tts)) {
        return 2;
    }
    time_t post_trees_time;
    time(&post_trees_time);
    //
    if (tts.get_num_trees() == 0) {
        cerr << "No tree to serve. Exiting...\n";
        return 3;
    }
    ////// ROUTES
    // tree web services
    auto r_about = make_shared< Resource >( );
    r_about->set_path( "/tree_of_life/about" );
    r_about->set_method_handler( "POST", about_method_handler );
    auto r_node_info = make_shared< Resource >( );
    r_node_info->set_path( "/tree_of_life/node_info" );
    r_node_info->set_method_handler( "POST", node_info_method_handler );
    auto r_mrca = make_shared< Resource >( );
    r_mrca->set_path( "/tree_of_life/mrca" );
    r_mrca->set_method_handler( "POST", mrca_method_handler );
    auto r_subtree = make_shared< Resource >( );
    r_subtree->set_path( "/tree_of_life/subtree" );
    r_subtree->set_method_handler( "POST", subtree_method_handler );
    auto r_induced_subtree = make_shared< Resource >( );
    r_induced_subtree->set_path( "/tree_of_life/induced_subtree" );
    r_induced_subtree->set_method_handler( "POST", induced_subtree_method_handler );
    // taxonomy web services
    auto r_tax_about = make_shared< Resource >( );
    r_tax_about->set_path( "/taxonomy/about" );
    r_tax_about->set_method_handler( "POST", tax_about_method_handler );
    auto r_taxon_info = make_shared< Resource >( );
    r_taxon_info->set_path( "/taxonomy/taxon_info" );
    r_taxon_info->set_method_handler( "POST", taxon_info_method_handler );
    auto r_taxon_mrca = make_shared< Resource >( );
    r_taxon_mrca->set_path( "/taxonomy/mrca" );
    r_taxon_mrca->set_method_handler( "POST", taxon_mrca_method_handler );
    auto r_taxon_subtree = make_shared< Resource >( );
    r_taxon_subtree->set_path( "/taxonomy/subtree" );
    r_taxon_subtree->set_method_handler( "POST", taxon_subtree_method_handler );
    /////  SETTINGS
    auto settings = make_shared< Settings >( );
    settings->set_port( port_number );
    settings->set_worker_limit( num_threads );
    settings->set_default_header( "Connection", "close" );
    
    Service service;
    global_service_ptr = &service;
    service.set_ready_handler( ready_handler );
    service.publish( r_about );
    service.publish( r_node_info );
    service.publish( r_mrca );
    service.publish( r_subtree );
    service.publish( r_induced_subtree );
    service.publish( r_tax_about );
    service.publish( r_taxon_info );
    service.publish( r_taxon_mrca );
    service.publish( r_taxon_subtree );
    service.set_signal_handler( SIGINT, sigterm_handler );
    service.set_signal_handler( SIGTERM, sigterm_handler );
    LOG(INFO) << "starting service with " << num_threads << " threads on port " << port_number << "...";
    time_t service_prep_time;
    time(&service_prep_time);
    LOG(INFO) << "Taxonomy reading took " << difftime(post_tax_time, start_time) << " seconds.";
    LOG(INFO) << "Tree reading took " << difftime(post_trees_time, post_tax_time) << " seconds.";
    LOG(INFO) << "Service prep took " << difftime(service_prep_time, post_trees_time) << " seconds.";
    LOG(INFO) << "Total boot took " << difftime(service_prep_time, start_time) << " seconds.";
    
    try {
        service.start( settings );
    } catch (std::exception & x) {
        LOG(ERROR) << "Exiting due to an exception after service.start: " << x.what();
        return 1;
    }
    return EXIT_SUCCESS;
}



po::variables_map parse_cmd_line(int argc, char* argv[]) { 
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description output("Server options");
    output.add_options()
        ("tree-dir,D",value<string>(),"Filepath to directory that will hold synthetic tree output")
        ("port,P",value<int>(),"Port to bind to.")
        ("pidfile,p",value<string>(),"filepath for PID")
        ("num-threads,n",value<int>(),"number of threads")
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