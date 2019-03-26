#include <thread>
// PID reporting from https://github.com/Corvusoft/restbed/blob/master/example/signal_handling/source/example.cpp
#ifdef _WIN32
    #include <process.h>
#else
    #include <unistd.h>
#endif
#include <boost/program_options.hpp>
#include <restbed>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "ws/tolwsadaptors.h"
#include "otc/otcli.h"

// unlike most headers, we'll go ahead an use namespaces
//    because this is an implementation file

using std::vector;
using std::map;
using std::multimap;
using std::shared_ptr;
using std::pair;
using std::set;
using std::string;
using std::ostringstream;
using std::optional;
using std::unique_ptr;

using std::make_shared;
using std::to_string;

namespace po = boost::program_options;
using namespace otc;
namespace fs = boost::filesystem;
using json = nlohmann::json;
using namespace restbed;
namespace chrono = std::chrono;
using std::pair;


namespace otc {
// global
TreesToServe tts;

#if defined(DEBUGGING_THREAD_LOCKING)
    std::mutex otc::ParallelReadSerialWrite::cout_mutex;
#endif

}// namespace otc

template<typename T>
optional<T> convert_to(const json & j);

template <>
optional<bool> convert_to(const json & j) {
    if (j.is_boolean()) {
        return j.get<bool>();
    }
    return {};
}

template <>
optional<string> convert_to(const json & j) {
    if (j.is_string()) {
        return j.get<string>();
    }
    return {};
}

template <>
optional<int> convert_to(const json & j) {
    if (j.is_number()) {
        return j.get<int>();
    }
    return {};
}

template <>
optional<vector<string>> convert_to(const json & j) {
    if (not j.is_array()) {
        return {};
    }
    vector<string> v;
    for(auto& xj: j) {
        auto x = convert_to<string>(xj);
        if (not x) return {};
        v.push_back(*x);
    }
    return v;
}

#if defined(LONG_OTT_ID)
template <>
optional<OttId> convert_to(const json & j) {
    return (j.is_number() ? j.get<OttId>() : {});
}
#endif

template <>
optional<OttIdSet> convert_to(const json & j) {
    if (not j.is_array()) {
        return {};
    }
    OttIdSet ids;
    for(auto& jid: j) {
        auto id = convert_to<OttId>(jid);
        if (not id) return {};
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
        return {};
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


string available_trees_method_handler(const json&)
{
    return available_trees_ws_method(tts);
}

string about_method_handler(const json& parsedargs)
{
    bool include_sources = extract_argument_or_default<bool>  (parsedargs, "include_source_list", false);
    string synth_id      = extract_argument_or_default<string>(parsedargs, "synth_id",            ""   );

    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr     = get_summary_tree(tts, synth_id);

    return about_ws_method(tts, treeptr, sta, include_sources);
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

auto lookup_source_id(const string& source_prefix, OttId foreign_id, const RichTaxonomy& taxonomy, const string& source_id)
{
    const auto & taxonomy_tree = taxonomy.get_tax_tree();
    const auto & taxonomy_tree_data = taxonomy_tree.get_data();

    try
    {
	if (source_prefix == "ncbi")
	    return taxonomy_tree_data.ncbi_id_map.at(foreign_id);
	else if (source_prefix == "gbif")
	    return taxonomy_tree_data.gbif_id_map.at(foreign_id);
	else if (source_prefix == "worms")
	    return taxonomy_tree_data.worms_id_map.at(foreign_id);
	else if (source_prefix == "if")
	    return taxonomy_tree_data.if_id_map.at(foreign_id);
	else if (source_prefix == "irmng")
	    return taxonomy_tree_data.irmng_id_map.at(foreign_id);
	else
	    throw OTCBadRequest() << "Don't recognize source_prefix = '" << source_prefix << "' - but we shouldn't get here.";
    }
    catch (std::out_of_range & x) {
	throw OTCBadRequest() << "No taxon in the taxonomy is associated with source_id of '"<<source_id<<"'";
    }
}


// Get an OttId from a string.
//
// The fact that this function is so long (uses exceptions) is ridiculous.
// + We could use std::strtol, which sets errno.
// + c++17 has a function from_chars( ) which seems less ridiculous.
//
OttId id_from_string(const string& id_str)
{
    std::size_t pos;
    long raw_id;
    try
    {
	raw_id  = std::stol(id_str.c_str(), &pos);
    }
    catch(const std::out_of_range&)
    {
	throw OTCBadRequest() << "The ID portion of the source_id was too large. Found: " << id_str;
    }
    catch(std::invalid_argument&)
    {
	throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " <<  id_str;
    }
    if (pos < id_str.length())
	throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " <<  id_str;

    auto id = to_OttId(raw_id);

    if (not id)
	throw OTCBadRequest() << "The ID portion of the source_id was too large. Found: " << id_str;

    return *id;
}

const RTRichTaxNode* taxon_from_source_id(const string& source_id, const RichTaxonomy& taxonomy)
{
    auto pref_id = split_string(source_id, ':');
    if (pref_id.size() != 2) {
	throw OTCBadRequest() << "Expecting exactly 1 colon in a source ID string. Found: \"" << source_id << "\".";
    }

    string source_prefix = *pref_id.begin();
    if (indexed_source_prefixes.count(source_prefix) == 0) {
	throw OTCBadRequest() << "IDs from source " << source_prefix << " are not known or not indexed for searching.";
    }

    string id_str = *pref_id.rbegin();

    auto foreign_id = id_from_string(id_str);

    auto in_ott = lookup_source_id(source_prefix, foreign_id, taxonomy, source_id);

#   if defined(MAP_FOREIGN_TO_POINTER)
        return in_ott;
#   else
	if (auto taxon_node = taxonomy.included_taxon_from_id(in_ott))
	    return taxon_node;
	else
	    throw OTCBadRequest("Foreign ID '"+ source_id+"' mapped to unknown OTT ID: "+ to_string(in_ott));
#   endif
}

string node_info_method_handler( const json& parsed_args)
{
    string synth_id = extract_argument_or_default<string>(parsed_args, "synth_id", "");
    auto node_id = extract_argument<string>(parsed_args,"node_id");
    auto source_id = extract_argument<string>(parsed_args,"source_id");

    if (node_id and source_id)
        throw OTCBadRequest("'node_id' and 'source_id' arguments cannot both be supplied.");
    else if (not node_id and not source_id)
        throw OTCBadRequest("An 'node_id' or 'source_id' argument is required.");

    if (source_id)
    {
	auto locked_taxonomy = tts.get_readable_taxonomy();
	const auto & taxonomy = locked_taxonomy.first;
	auto tax_node = taxon_from_source_id(*source_id, taxonomy);
	node_id = "ott"+std::to_string(tax_node->get_ott_id());
    }

    bool include_lineage = extract_argument_or_default<bool>(parsed_args, "include_lineage", false);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);

    return node_info_ws_method(tts, treeptr, sta, *node_id, include_lineage);
}

string mrca_method_handler( const json& parsedargs)
{
    string synth_id;
    vector<string> node_id_vec;
    tie(synth_id, node_id_vec) = get_synth_and_node_id_vec(parsedargs);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
    return mrca_ws_method(tts, treeptr, sta, node_id_vec);
}

std::string process_subtree(const json& parsedargs)
{
    // FIXME: According to treemachine/ws-tests/tests.subtree, there is an "include_all_node_labels"
    //        argument.  Unless this is explicitly set to true, we are supposed to not write node labels
    //        for non-ottids.  At least in Newick.

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

    bool all_node_labels = extract_argument_or_default<bool>(parsedargs, "include_all_node_labels", false);

    if (format == "newick") {
	return newick_subtree_ws_method(tts, treeptr, node_id, nns, all_node_labels, height_limit);
    } else {
	return arguson_subtree_ws_method(tts, treeptr, sta, node_id, height_limit);
    }
}

string induced_subtree_method_handler( const json& parsedargs )
{
    string synth_id;
    vector<string> node_id_vec;
    tie(synth_id, node_id_vec) = get_synth_and_node_id_vec(parsedargs);
    NodeNameStyle nns = get_label_format(parsedargs);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
    return induced_subtree_ws_method(tts, treeptr, node_id_vec, nns);
}

string tax_about_method_handler( const json& )
{
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tax_about_ws_method(taxonomy);
}

// looks for ott_id or source_id args to find a node
const RTRichTaxNode * extract_taxon_node_from_args(const json & parsedargs, const RichTaxonomy & taxonomy)
{
    auto ott_id = extract_argument<OttId>(parsedargs, "ott_id");
    auto source_id = extract_argument<string>(parsedargs, "source_id");

    if (ott_id and source_id)
        throw OTCBadRequest("'ott_id' and 'source_id' arguments cannot both be supplied.");
    else if (not ott_id and not source_id)
        throw OTCBadRequest("An 'ott_id' or 'source_id' argument is required.");

    if (source_id)
	return taxon_from_source_id(*source_id, taxonomy);
    else {
        assert(ott_id);
        auto taxon_node = taxonomy.included_taxon_from_id(*ott_id);
        if (taxon_node == nullptr) {
            throw OTCBadRequest() << "Unrecognized OTT ID: " << *ott_id;
        }
        return taxon_node;
    }
}

string taxon_info_method_handler( const json& parsedargs )
{
    auto include_lineage = extract_argument_or_default<bool>(parsedargs, "include_lineage", false);
    auto include_children = extract_argument_or_default<bool>(parsedargs, "include_children", false);
    auto include_terminal_descendants = extract_argument_or_default<bool>(parsedargs, "include_terminal_descendants", false);       

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
    return taxon_info_ws_method(taxonomy, taxon_node, include_lineage, include_children, include_terminal_descendants);
}

string taxon_flags_method_handler( const json& )
{
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return taxonomy_flags_ws_method(taxonomy);
}

string taxon_mrca_method_handler( const json& parsedargs )
{
    OttIdSet ott_id_set = extract_required_argument<OttIdSet>(parsedargs, "ott_ids");
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return taxonomy_mrca_ws_method(taxonomy, ott_id_set);
}

string taxon_subtree_method_handler( const json& parsedargs )
{
    NodeNameStyle nns = get_label_format(parsedargs);
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
    return taxon_subtree_ws_method(taxonomy, taxon_node, nns);
}

// See taxomachine/src/main/java/org/opentree/taxonomy/plugins/tnrs_v3.java

// 10,000 queries at .0016 second per query = 16 seconds
const int MAX_NONFUZZY_QUERY_STRINGS = 10000;
// 250 queries at .3 second per query = 75 seconds
const int MAX_FUZZY_QUERY_STRINGS = 250;

static string LIFE_NODE_NAME = "life";

string tnrs_match_names_handler( const json& parsedargs )
{
    // 1. Requred argument: "names"
    vector<string> names = extract_required_argument<vector<string>>(parsedargs, "names");

    // 2. Optional argunments
    optional<string> context_name = extract_argument<string>(parsedargs, "context_name");
    bool do_approximate_matching  = extract_argument_or_default(parsedargs, "do_approximate_matching", false);
    vector<string> ids            = extract_argument_or_default(parsedargs, "ids",                     names);
    bool include_suppressed       = extract_argument_or_default(parsedargs, "include_suppressed",      false);

    // 3. Check that "ids" have the same length as "names", if supplied
    if (ids.size() != names.size())
    {
	throw OTCBadRequest()<<"The number of names and ids does not match. If you provide ids, then you "
			     <<"must provide exactly as many ids as names.";
    }

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;

    return tnrs_match_names_ws_method(names, context_name, do_approximate_matching, ids, include_suppressed, taxonomy);
}

string tnrs_autocomplete_name_handler( const json& parsedargs )
{
    string name              = extract_required_argument<string>(parsedargs, "name");
    string context_name      = extract_argument_or_default(parsedargs, "context_name",            LIFE_NODE_NAME);
    bool include_suppressed  = extract_argument_or_default(parsedargs, "include_suppressed",      false);

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tnrs_autocomplete_name_ws_method(name, context_name, include_suppressed, taxonomy);
}

string tnrs_contexts_handler( const json& )
{
    return tnrs_contexts_ws_method();
}

string tnrs_infer_context_handler( const json& parsedargs )
{
    vector<string> names = extract_required_argument<vector<string>>(parsedargs, "names");

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tnrs_infer_context_ws_method(names, taxonomy);
}

string conflict_status_method_handler( const json& parsed_args )
{
    auto tree1newick = extract_argument<string>(parsed_args, "tree1newick");
    auto tree1 = extract_argument<string>(parsed_args, "tree1");

    string tree2 = extract_required_argument<string>(parsed_args, "tree2");

    const auto& summary = *tts.get_summary_tree("");
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;

    if (tree1newick)
	return newick_conflict_ws_method(summary, taxonomy, *tree1newick, tree2);
    else if (tree1)
	return phylesystem_conflict_ws_method(summary, taxonomy, *tree1, tree2);
    else
	throw OTCBadRequest()<<"Expecting argument 'tree1' or argument 'tree1newick'";
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
	std::ofstream pstream(pidfile);
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

// this is a hack.  Also it doesn't include the time zone
string ctime(const chrono::system_clock::time_point& t)
{
    time_t t2 = chrono::system_clock::to_time_t(t);
    char* c = ctime(&t2);
    string tt = c;
    tt.pop_back(); // remove newline
    return tt;
}

multimap<string,string> request_headers(const string& rbody)
{
    multimap<string,string> headers;

    headers.insert({ "Access-Control-Allow-Credentials", "true" });
    headers.insert({ "Access-Control-Allow-Origin", "*" });
    headers.insert({ "Access-Control-Max-Age","86400" });
    headers.insert({ "Cache-Control", "no-store, no-cache, must-revalidate, post-check=0, pre-check=0"});
// Connection:Keep-Alive  -- I 
//  We're calling 'close' so this doesn't make sense, I think...
    headers.insert({ "Content-Length", ::to_string(rbody.length())});

//  All of our replies should be JSON, so that users can unconditionally parse the response a JSON.
    headers.insert({ "Content-Type", "application/json;"});

    headers.insert({ "Date", ctime(chrono::system_clock::now())});
    headers.insert({ "Expires", ctime(chrono::system_clock::now())});
// Keep-Alive:timeout=5, max=99
    headers.insert({ "Pragma", "no-cache"});
// Server:Apache/2.4.10 (Debian)
// Set-Cookie:session_id_phylesystem=152.3.12.201-82a8df2c-2a3d-4c10-aca3-3c79ca8ebdd1; Path=/
    headers.insert({ "Vary", "Accept-Encoding"} );   // cache separately for each encoding?
    headers.insert({ "X-Powered-By","otc-tol-ws"});  //X-Powered-By:web2py
    return headers;
}

multimap<string,string> options_headers()
{
    multimap<string,string> headers;
    headers.insert({ "Access-Control-Allow-Credentials", "true" });
    headers.insert({ "Access-Control-Allow-Headers", "content-type" });
    headers.insert({ "Access-Control-Allow-Methods", "POST"});
    headers.insert({ "Access-Control-Allow-Origin", "*" });
    headers.insert({ "Access-Control-Max-Age","86400" });
//    headers.insert({ "Connection", "Keep-Alive"});

//    There is no content, to don't include a Content-Type header
//    headers.insert({ "Content-Type", "text/html; charset=UTF-8"});
    headers.insert({ "Content-Length", "0"});
    headers.insert({ "Date", ctime(chrono::system_clock::now())});  //    Date:Mon, 22 May 2017 20:55:05 GMT
    headers.insert({ "X-Powered-By","otc-tol-ws"});  //X-Powered-By:web2py
//Keep-Alive:timeout=5, max=100
//Server:Apache/2.4.10 (Debian)
//Set-Cookie:session_id_phylesystem=152.3.12.201-0006b381-2288-41df-b7f9-d4285aad48a2; Path=/

    return headers;
}

std::string error_response(const string& path, const std::exception& e)
{
    string msg = string("[") + path + ("] Error: ") + e.what();
    LOG(DEBUG)<<msg;
    json j = { {"message", msg} };
    return j.dump(4)+"\n";
}

std::function<void(const shared_ptr< Session > session)>
create_method_handler(const string& path, const std::function<std::string(const json&)> process_request)
{
    return [=](const shared_ptr< Session > session ) {
	const auto request = session->get_request( );
	size_t content_length = request->get_header( "Content-Length", 0 );
	session->fetch( content_length, [ path, process_request, request ]( const shared_ptr< Session > session, const Bytes & body ) {
		try {
		    LOG(DEBUG)<<"request: "<<path;
		    json parsedargs = parse_body_or_throw(body);
		    LOG(DEBUG)<<"   argument "<<parsedargs.dump(1);
		    auto rbody = process_request(parsedargs);
		    LOG(DEBUG)<<"request: DONE";
		    session->close( OK, rbody, request_headers(rbody) );
		} catch (OTCWebError& e) {
		    string rbody = error_response(path,e);
		    session->close( e.status_code(), rbody, request_headers(rbody) );
		} catch (OTCError& e) {
		    string rbody = error_response(path,e);
		    session->close( 500, rbody, request_headers(rbody) );
		}
	    });
    };
}

json request_to_json(const Request& request)
{
    LOG(DEBUG)<<"GET "<<request.get_path();
    json query;
    for(auto& key_value_pair: request.get_query_parameters())
	query[key_value_pair.first] = key_value_pair.second;
    return query;
}

std::function<void(const shared_ptr< Session > session)>
create_GET_method_handler(const string& path, const std::function<std::string(const json&)> process_request)
{
    return [=](const shared_ptr< Session > session )
    {
	try
	{
	    LOG(DEBUG)<<"request: "<<path;
	    const auto& request = session->get_request( );
	    json parsedargs = request_to_json(*request);
	    LOG(DEBUG)<<"   argument "<<parsedargs.dump(1);
	    auto rbody = process_request(parsedargs);
	    LOG(DEBUG)<<"request: DONE";
	    session->close( OK, rbody, request_headers(rbody) );
	}
	catch (OTCWebError& e)
	{
	    string rbody = error_response(path, e);
	    session->close( e.status_code(), rbody, request_headers(rbody) );
	} catch (OTCError& e) {
	    string rbody = error_response(path, e);
	    session->close( 500, rbody, request_headers(rbody) );
	}
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
    string v3_prefix = "/v3";
    string v4_prefix = "/v4-beta";

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
	std::cerr << "Expecting a tree-dir argument for a path to a directory of synth outputs.\n";
        return 1;
    }
    const fs::path topdir{args["tree-dir"].as<string>()};

    // Must load taxonomy before trees
    LOG(INFO) << "reading taxonomy...";
    RichTaxonomy taxonomy = load_rich_taxonomy(args);
    time_t post_tax_time;
    time(&post_tax_time);
    tts.set_taxonomy(taxonomy);

    // Now load trees
    if (!read_trees(topdir, tts)) {
        return 2;
    }
    time_t post_trees_time;
    time(&post_trees_time);
    if (tts.get_num_trees() == 0) {
	std::cerr << "No tree to serve. Exiting...\n";
        return 3;
    }

    ////// v3 ROUTES
    // tree web services
    auto v3_r_about            = path_handler(v3_prefix + "/tree_of_life/about", about_method_handler);
    auto v3_r_node_info        = path_handler(v3_prefix + "/tree_of_life/node_info", node_info_method_handler );
    auto v3_r_mrca             = path_handler(v3_prefix + "/tree_of_life/mrca", mrca_method_handler );
    auto v3_r_subtree          = path_handler(v3_prefix + "/tree_of_life/subtree", process_subtree);
    auto v3_r_induced_subtree  = path_handler(v3_prefix + "/tree_of_life/induced_subtree", induced_subtree_method_handler );

    // taxonomy web services
    auto v3_r_tax_about        = path_handler(v3_prefix + "/taxonomy/about", tax_about_method_handler );
    auto v3_r_taxon_info       = path_handler(v3_prefix + "/taxonomy/taxon_info", taxon_info_method_handler );
    auto v3_r_taxon_flags      = path_handler(v3_prefix + "/taxonomy/flags", taxon_flags_method_handler );
    auto v3_r_taxon_mrca       = path_handler(v3_prefix + "/taxonomy/mrca", taxon_mrca_method_handler );
    auto v3_r_taxon_subtree    = path_handler(v3_prefix + "/taxonomy/subtree", taxon_subtree_method_handler );

    // tnrs
    auto v3_r_tnrs_match_names       = path_handler(v3_prefix + "/tnrs/match_names", tnrs_match_names_handler );
    auto v3_r_tnrs_autocomplete_name = path_handler(v3_prefix + "/tnrs/autocomplete_name", tnrs_autocomplete_name_handler );
    auto v3_r_tnrs_contexts          = path_handler(v3_prefix + "/tnrs/contexts", tnrs_contexts_handler );
    auto v3_r_tnrs_infer_context     = path_handler(v3_prefix + "/tnrs/infer_context", tnrs_infer_context_handler );

    // conflict
    auto v3_r_conflict_status  = path_handler(v3_prefix + "/conflict/conflict-status", conflict_status_method_handler );

    // v2 conflict --
    auto v3_r_old_conflict_status = make_shared< Resource >( );
    {
	string path = v3_prefix + "/conflict/old-conflict-status";
	v3_r_old_conflict_status->set_path( path );
	v3_r_old_conflict_status->set_method_handler( "GET", create_GET_method_handler(path, conflict_status_method_handler) );
	v3_r_old_conflict_status->set_method_handler( "OPTIONS", options_method_handler);
    }

    ////// v4 ROUTES

    // tree web services
    auto v4_r_available_trees  = path_handler(v4_prefix + "/tree_of_life/available_trees", available_trees_method_handler);
    auto v4_r_about            = path_handler(v4_prefix + "/tree_of_life/about", about_method_handler);
    auto v4_r_node_info        = path_handler(v4_prefix + "/tree_of_life/node_info", node_info_method_handler );
    auto v4_r_mrca             = path_handler(v4_prefix + "/tree_of_life/mrca", mrca_method_handler );
    auto v4_r_subtree          = path_handler(v4_prefix + "/tree_of_life/subtree", process_subtree);
    auto v4_r_induced_subtree  = path_handler(v4_prefix + "/tree_of_life/induced_subtree", induced_subtree_method_handler );

    // taxonomy web services
    auto v4_r_tax_about        = path_handler(v4_prefix + "/taxonomy/about", tax_about_method_handler );
    auto v4_r_taxon_info       = path_handler(v4_prefix + "/taxonomy/taxon_info", taxon_info_method_handler );
    auto v4_r_taxon_flags      = path_handler(v4_prefix + "/taxonomy/flags", taxon_flags_method_handler );
    auto v4_r_taxon_mrca       = path_handler(v4_prefix + "/taxonomy/mrca", taxon_mrca_method_handler );
    auto v4_r_taxon_subtree    = path_handler(v4_prefix + "/taxonomy/subtree", taxon_subtree_method_handler );

    // tnrs
    auto v4_r_tnrs_match_names       = path_handler(v4_prefix + "/tnrs/match_names", tnrs_match_names_handler );
    auto v4_r_tnrs_autocomplete_name = path_handler(v4_prefix + "/tnrs/autocomplete_name", tnrs_autocomplete_name_handler );
    auto v4_r_tnrs_contexts          = path_handler(v4_prefix + "/tnrs/contexts", tnrs_contexts_handler );
    auto v4_r_tnrs_infer_context     = path_handler(v4_prefix + "/tnrs/infer_context", tnrs_infer_context_handler );

    // conflict
    auto v4_r_conflict_status  = path_handler(v4_prefix + "/conflict/conflict-status", conflict_status_method_handler );

    /////  SETTINGS
    auto settings = make_shared< Settings >( );
    settings->set_port( port_number );
    settings->set_worker_limit( num_threads );
    settings->set_default_header( "Connection", "close" );
    
    Service service;
    global_service_ptr = &service;
    service.set_ready_handler( ready_handler );
    service.publish( v3_r_about );
    service.publish( v3_r_node_info );
    service.publish( v3_r_mrca );
    service.publish( v3_r_subtree );
    service.publish( v3_r_induced_subtree );
    service.publish( v3_r_tax_about );
    service.publish( v3_r_taxon_info );
    service.publish( v3_r_taxon_flags );
    service.publish( v3_r_taxon_mrca );
    service.publish( v3_r_taxon_subtree );
    service.publish( v3_r_tnrs_match_names );
    service.publish( v3_r_tnrs_autocomplete_name );
    service.publish( v3_r_tnrs_contexts );
    service.publish( v3_r_tnrs_infer_context );
    service.publish( v3_r_conflict_status );
    service.publish( v3_r_old_conflict_status );

    service.publish( v4_r_available_trees );
    service.publish( v4_r_about );
    service.publish( v4_r_node_info );
    service.publish( v4_r_mrca );
    service.publish( v4_r_subtree );
    service.publish( v4_r_induced_subtree );
    service.publish( v4_r_tax_about );
    service.publish( v4_r_taxon_info );
    service.publish( v4_r_taxon_flags );
    service.publish( v4_r_taxon_mrca );
    service.publish( v4_r_taxon_subtree );
    service.publish( v4_r_tnrs_match_names );
    service.publish( v4_r_tnrs_autocomplete_name );
    service.publish( v4_r_tnrs_contexts );
    service.publish( v4_r_tnrs_infer_context );
    service.publish( v4_r_conflict_status );

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
