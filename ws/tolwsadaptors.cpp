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
using boost::optional;
using std::pair;


namespace otc {
// global
TreesToServe tts;

#if defined(DEBUGGING_THREAD_LOCKING)
    std::mutex otc::ParallelReadSerialWrite::cout_mutex;
#endif

}// namespace otc

json parse_body_or_throw(const Bytes & body) {
    // This line is necessary for the otcetera test for tree_of_life/about to succeed, apparently.
    if (body.empty()) {
        return json();
    }
    try {
        return json::parse(body);
    }
    catch (...) {
        throw OTCBadRequest("Could not parse body of call as JSON.\n");
    }
}

template<typename T>
optional<T> convert_to(const json & j);

template <>
optional<bool> convert_to(const json & j) {
    if (j.is_boolean()) {
        return j.get<bool>();
    }
    return boost::none;
}

template <>
optional<string> convert_to(const json & j) {
    if (j.is_string()) {
        return j.get<string>();
    }
    return boost::none;
}

template <>
optional<int> convert_to(const json & j) {
    if (j.is_number()) {
        return j.get<int>();
    }
    return boost::none;
}

template <>
optional<vector<string>> convert_to(const json & j) {
    if (not j.is_array()) {
        return boost::none;
    }
    vector<string> v;
    for(auto& xj: j) {
        auto x = convert_to<string>(xj);
        if (not x) return boost::none;
        v.push_back(*x);
    }
    return v;
}

#if defined(LONG_OTT_ID)
template <>
optional<OttId> convert_to(const json & j) {
    return (j.is_number() ? j.get<OttId>() : boost::none);
}
#endif

template <>
optional<OttIdSet> convert_to(const json & j) {
    if (not j.is_array()) {
        return boost::none;
    }
    OttIdSet ids;
    for(auto& jid: j) {
        auto id = convert_to<OttId>(jid);
        if (not id) return boost::none;
        ids.insert(*id);
    }
    return ids;
}



template <typename T>
constexpr const char* type_name_with_article();

template <> constexpr const char* type_name_with_article<bool>() {
    return "a boolean";
}
template <> constexpr const char* type_name_with_article<int>() {
    return "an integer";
}
template <> constexpr const char* type_name_with_article<string>() {
    return "a string";
}
#if defined(LONG_OTT_ID)
template <> constexpr const char* type_name_with_article<OttId>() {
    return "an OttId";
}
#endif
template <> constexpr const char* type_name_with_article<vector<string>>() {
    return "an array of strings";
}
template <> constexpr const char* type_name_with_article<OttIdSet>() {
    return "an array of integers";
}

template<typename T>
optional<T> extract_argument(const json & j, const std::string& opt_name, bool required=false) {
    auto opt = j.find(opt_name);
    if (opt == j.end()) {
        if (required) {
            throw OTCBadRequest("expecting ") << type_name_with_article<T>() << " argument called \"" << opt_name << "\"\n";
        }
        return boost::none;
    }
    auto arg = convert_to<T>(*opt);
    if (not arg) {
        throw OTCBadRequest("expecting argument '") << opt_name << "' to be " << type_name_with_article<T>() <<"! Found \"" << *opt << "\"\n";
    }
    return arg;
}

template<typename T>
T extract_argument_or_default(const json & j, const std::string& opt_name, const T& _default_) {
    auto arg = extract_argument<T>(j, opt_name);
    return (arg ? *arg : _default_);
}

template<typename T>
T extract_required_argument(const json & j, const std::string& opt_name) {
    auto arg = extract_argument<T>(j, opt_name, true);
    assert(arg);
    return *arg;
}

///////////////////////
// handlers that are registered as callback

const SummaryTreeAnnotation * get_annotations(const TreesToServe& tts, const string& synth_id) {
    auto sta = tts.get_annotations(synth_id);
    if (not sta) {
        throw OTCBadRequest() << "Did not recognize synth_id '" << synth_id << "'";
    }
    return sta;
}

const SummaryTree_t * get_summary_tree(const TreesToServe& tts, const string& synth_id) {
    auto treeptr = tts.get_summary_tree(synth_id);
    if (not treeptr) {
        throw OTCBadRequest() << "Did not recognize synth_id '" << synth_id << "'";
    }
    return treeptr;
}


void about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            auto parsedargs = parse_body_or_throw(body);

            bool include_sources = extract_argument_or_default<bool>  (parsedargs, "include_source_list", false);
            string synth_id      = extract_argument_or_default<string>(parsedargs, "synth_id",            ""   );

            const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
            const SummaryTree_t * treeptr     = get_summary_tree(tts, synth_id);

            string rbody = about_ws_method(tts, treeptr, sta, include_sources);

            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[/tree_of_life/about] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}


pair<string,string> get_synth_and_node_id(const json &j) {
    auto synth_id = extract_argument_or_default<string>(j, "synth_id", "");
    auto node_id = extract_required_argument<string>(j, "node_id");
    return pair<string,string>(synth_id, node_id);
}

pair<string,vector<string>> get_synth_and_node_id_vec(const json &j) {
    auto synth_id =  extract_argument_or_default<string>(j, "synth_id", "");
    auto node_id_vec = extract_required_argument<vector<string>>(j, "node_ids");
    return pair<string, vector<string>>(synth_id, node_id_vec);
}


NodeNameStyle get_label_format(const json &j) {
    auto label_format = extract_argument_or_default<string>(j, "label_format", "name_and_id");
    if (label_format == "id") {
        return NodeNameStyle::NNS_ID_ONLY;
    } else if (label_format == "name") {
        return NodeNameStyle::NNS_NAME_ONLY;
    } else if (label_format == "name_and_id") {
        return NodeNameStyle::NNS_NAME_AND_ID;
    } else {
        throw OTCBadRequest("label_format must be \"name_and_id\", \"name\" or \"id\".\n");
    }
}

void node_info_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            json parsedargs = parse_body_or_throw(body);
            string synth_id;
            string node_id;
            tie(synth_id, node_id) = get_synth_and_node_id(parsedargs);
            bool include_lineage = extract_argument_or_default<bool>(parsedargs, "include_lineage", false);
            const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
            const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
            string rbody = node_info_ws_method(tts, treeptr, sta, node_id, include_lineage);
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[/tree_of_life/node_info] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

void mrca_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {    
            auto parsedargs = parse_body_or_throw(body);
            string synth_id;
            vector<string> node_id_vec;
            tie(synth_id, node_id_vec) = get_synth_and_node_id_vec(parsedargs);
            const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
            const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
            string rbody = mrca_ws_method(tts, treeptr, sta, node_id_vec);
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[/tree_of_life/mrca] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

std::string process_subtree(const json& parsedargs)
{
    string synth_id;
    string node_id;
    tie(synth_id, node_id) = get_synth_and_node_id(parsedargs);
    auto format = extract_argument_or_default<string>(parsedargs, "format", "newick");
    if (format != "newick" && format != "arguson") {
	throw OTCBadRequest("format must be \"newick\" or \"arguson\".\n");
    }
    NodeNameStyle nns = get_label_format(parsedargs);
    int height_limit = extract_argument_or_default<int>(parsedargs, "height_limit", (format == "arguson")? 3 : -1);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);

    if (format == "newick") {
	return newick_subtree_ws_method(tts, treeptr, node_id, nns, height_limit);
    } else {
	return arguson_subtree_ws_method(tts, treeptr, sta, node_id, height_limit);
    }
}

void induced_subtree_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            auto parsedargs = parse_body_or_throw(body);
            string synth_id;
            vector<string> node_id_vec;
            tie(synth_id, node_id_vec) = get_synth_and_node_id_vec(parsedargs);
            NodeNameStyle nns = get_label_format(parsedargs);
            const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
            const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
            auto rbody = induced_subtree_ws_method(tts, treeptr, node_id_vec, nns);
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[tree_of_life/induced_subtree] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

void tax_about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & ) {
        try {
            string rbody;
            {
                auto locked_taxonomy = tts.get_readable_taxonomy();
                const auto & taxonomy = locked_taxonomy.first;
                rbody = tax_about_ws_method(taxonomy);
            }
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[taxonomy/about] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

// looks for ott_id or source_id args to find a node
const RTRichTaxNode * extract_taxon_node_from_args(const json & parsedargs, const RichTaxonomy & taxonomy) {
    const auto & taxonomy_tree = taxonomy.get_tax_tree();
    const auto & taxonomy_tree_data = taxonomy_tree.get_data();

    auto ott_id = extract_argument<OttId>(parsedargs, "ott_id");
    auto source_id = extract_argument<string>(parsedargs, "source_id");

    if (ott_id and source_id) {
        throw OTCBadRequest("'ott_id' and 'source_id' arguments cannot both be supplied.");
    } else if (not ott_id and not source_id) {
        throw OTCBadRequest("An 'ott_id' or 'source_id' argument is required.");
    } else if (source_id) {
        auto pref_id = split_string(*source_id, ':');
        if (pref_id.size() != 2) {
            throw OTCBadRequest() << "Expecting exactly 1 colon in a source ID string. Found: \"" << *source_id << "\".";
        }
        string source_prefix = *pref_id.begin();
        if (indexed_source_prefixes.count(source_prefix) == 0) {
            throw OTCBadRequest() << "IDs from source " << source_prefix << " are not known or not indexed for searching.";
        }
        string id_str = *pref_id.rbegin();
        try {
            std::size_t pos;
            long raw_foreign_id  = std::stol(id_str.c_str(), &pos);
            if (pos < id_str.length()) {
                throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " <<  id_str;
            }
            OttId foreign_id = std::stol(id_str.c_str(), &pos);
            if (pos < id_str.length()) {
                throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " << id_str;
            }
            try {
                foreign_id = check_ott_id_size(raw_foreign_id);
            } catch (OTCError &) {
                throw OTCBadRequest() << "The ID portion of the source_id was too large. Found: " << id_str;
            }
            try {
#               if defined(MAP_FOREIGN_TO_POINTER)
                    const RTRichTaxNode * in_ott = nullptr;
#               else
                    OttId in_ott = 0;
#               endif

                if (source_prefix == "ncbi") {
                    in_ott = taxonomy_tree_data.ncbi_id_map.at(foreign_id);
                } else if (source_prefix == "gbif") {
                    in_ott = taxonomy_tree_data.gbif_id_map.at(foreign_id);
                } else if (source_prefix == "worms") {
                    in_ott =  taxonomy_tree_data.worms_id_map.at(foreign_id);
                } else if (source_prefix == "if") {
                    in_ott =  taxonomy_tree_data.if_id_map.at(foreign_id);
                } else if (source_prefix == "irmng") {
                    in_ott =  taxonomy_tree_data.irmng_id_map.at(foreign_id);
                } else {
                    assert(false); // We should catch an unknown source_prefix above
                    throw OTCBadRequest() << "Don't recognized source_prefix = '" << source_prefix << "' - but we shouldn't get here.";
                }
#               if defined(MAP_FOREIGN_TO_POINTER)
                    return in_ott;
#               else
                    auto taxon_node =  taxonomy.included_taxon_from_id(in_ott);
                    if (taxon_node == nullptr) {
                        rbody = "Foreign ID " ;
                        rbody += supplied_source_id;
                        rbody += " mapped to unknown OTT ID: ";
                        rbody += to_string(ott_id);
                        status_code = 400;
                    }
                    return taxon_node;
#               endif
            } catch (std::out_of_range & x) {
                throw OTCBadRequest() << "No taxon in the taxonomy is associated with source_id of " << source_id;
            }
        } catch (...) {
            throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " << id_str;
        }
    } else {
        assert(ott_id);
        auto taxon_node = taxonomy.included_taxon_from_id(*ott_id);
        if (taxon_node == nullptr) {
            throw OTCBadRequest() << "Unrecognized OTT ID: " << *ott_id;
        }
        return taxon_node;
    }
}

void taxon_info_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            auto parsedargs = parse_body_or_throw(body);
            auto include_lineage = extract_argument_or_default<bool>(parsedargs, "include_lineage", false);
            auto include_children = extract_argument_or_default<bool>(parsedargs, "include_children", false);
            auto include_terminal_descendants = extract_argument_or_default<bool>(parsedargs, "include_terminal_descendants", false);       
            string rbody;
            {
                auto locked_taxonomy = tts.get_readable_taxonomy();
                const auto & taxonomy = locked_taxonomy.first;
                const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
                rbody = taxon_info_ws_method(taxonomy, taxon_node, include_lineage, include_children, include_terminal_descendants);
            }
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[taxonomy/taxon_info] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}
        
void taxon_mrca_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            auto parsedargs = parse_body_or_throw(body);
            OttIdSet ott_id_set = extract_required_argument<OttIdSet>(parsedargs, "ott_ids");
            string rbody;
            {
                auto locked_taxonomy = tts.get_readable_taxonomy();
                const auto & taxonomy = locked_taxonomy.first;
                rbody = taxonomy_mrca_ws_method(taxonomy, ott_id_set);
            }
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[taxonomy/mrca] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

void taxon_subtree_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            string rbody;
            int status_code = OK;
            json parsedargs = parse_body_or_throw(body);
            NodeNameStyle nns = get_label_format(parsedargs);
            {
                auto locked_taxonomy = tts.get_readable_taxonomy();
                const auto & taxonomy = locked_taxonomy.first;
                const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
                rbody = taxon_subtree_ws_method(taxonomy, taxon_node, nns);
            }
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[taxonomy/subtree] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
    });
}

void conflict_conflict_status_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        try {
            auto parsed_args = parse_body_or_throw(body);
            string tree1 = extract_required_argument<string>(parsed_args, "tree1");
            string tree2 = extract_required_argument<string>(parsed_args, "tree2");
            const auto& summary = *tts.get_summary_tree("");
            string rbody;
            {
                auto locked_taxonomy = tts.get_readable_taxonomy();
                const auto & taxonomy = locked_taxonomy.first;
                rbody = conflict_ws_method(summary, taxonomy, tree1, tree2);
            }
            session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        } catch (OTCWebError& e) {
            string rbody = string("[conflict-status] Error: ") + e.what();
            session->close( e.status_code(), rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
        }
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

multimap<string,string> request_headers(const string& rbody)
{
    multimap<string,string> headers;
    headers.insert({ "Access-Control-Allow-Credentials", "true" });
    headers.insert({ "Access-Control-Allow-Methods", "POST"});
    headers.insert({ "Access-Control-Allow-Origin", "*" });
    headers.insert({ "Access-Control-Max-Age","86400" });
    headers.insert({ "Connection", "Keep-Alive"});
    headers.insert({ "Content-Length", ::to_string(rbody.length())});
    return headers;
}

multimap<string,string> options_headers()
{
    multimap<string,string> headers;
    headers.insert({ "Access-Control-Allow-Credentials", "true" });
    headers.insert({ "Access-Control-Allow-Methods", "POST"});
    headers.insert({ "Access-Control-Allow-Origin", "*" });
    headers.insert({ "Access-Control-Max-Age","86400" });
    headers.insert({ "Connection", "Keep-Alive"});
    headers.insert({ "Content-Type", "text/html; charset=UTF-8"});
    headers.insert({ "Content-Length", "0"});
    return headers;
}

std::function<void(const shared_ptr< Session > session)>
create_method_handler(const string& path, const std::function<std::string(const json&)> process_request)
{
    return [=](const shared_ptr< Session > session ) {
	const auto request = session->get_request( );
	size_t content_length = request->get_header( "Content-Length", 0 );
	session->fetch( content_length, [ path, process_request, request ]( const shared_ptr< Session > session, const Bytes & body ) {
		try {
		    json parsedargs = parse_body_or_throw(body);
		    auto rbody = process_request(parsedargs);
		    session->close( OK, rbody, request_headers(rbody) );
		} catch (OTCWebError& e) {
		    string rbody = string("[") + path + ("] Error: ") + e.what();
		    session->close( e.status_code(), rbody, request_headers(rbody) );
		}
	    });
    };
}

void options_method_handler( const shared_ptr< Session > session ) {
    session->close( OK, "", options_headers() );
}

shared_ptr< Resource > path_handler(const string& path, std::function<std::string(const json &)> process_request)
{
    auto r_subtree = make_shared< Resource >( );
    r_subtree->set_path( path );
    r_subtree->set_method_handler( "POST", create_method_handler(path,process_request));
    r_subtree->set_method_handler( "OPTIONS", options_method_handler);
    return r_subtree;
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

    auto r_subtree = path_handler("/tree_of_life/subtree", process_subtree);

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
    // conflict
    auto r_conflict_conflict_status = make_shared< Resource >( );
    r_conflict_conflict_status->set_path( "/conflict/conflict-status" );
    r_conflict_conflict_status->set_method_handler( "POST", conflict_conflict_status_method_handler );
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
    service.publish( r_conflict_conflict_status );
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
