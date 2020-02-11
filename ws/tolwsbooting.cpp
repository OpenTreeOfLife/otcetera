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
#include "otc/ws/tolwsadaptors.h"
#include "otc/otcli.h"
#include "otc/ctrie/context_ctrie_db.h"
#include "otc/tnrs/context.h"
#include "otc/supertree_util.h"

INITIALIZE_EASYLOGGINGPP

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
            throw OTCBadRequest("expecting ") << type_name_with_article<T>() << " argument called '" << opt_name << "'\n";
        }
        return {};
    }
    auto arg = convert_to<T>(*opt);
    if (not arg) {
        throw OTCBadRequest("expecting argument '") << opt_name << "' to be " << type_name_with_article<T>() <<"! Found '" << opt->dump() << "'\n";
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


namespace otc {
// global
TreesToServe tts;


}// namespace otc

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


string available_trees_method_handler(const json&) {
    return available_trees_ws_method(tts);
}

string about_method_handler(const json& parsedargs) {
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

auto lookup_source_id(const string& source_prefix, OttId foreign_id, const RichTaxonomy& taxonomy, const string& source_id) {
    const auto & taxonomy_tree = taxonomy.get_tax_tree();
    const auto & taxonomy_tree_data = taxonomy_tree.get_data();
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
            throw OTCBadRequest() << "Don't recognize source_prefix = '" << source_prefix << "' - but we shouldn't get here.";
        }
    } catch (std::out_of_range & x) {
        throw OTCBadRequest() << "No taxon in the taxonomy is associated with source_id of '"<< source_id <<"'";
    }
}


// Get an OttId from a string.
//
// The fact that this function is so long (uses exceptions) is ridiculous.
// + We could use std::strtol, which sets errno.
// + c++17 has a function from_chars( ) which seems less ridiculous.
//
OttId id_from_string(const string& id_str) {
    std::size_t pos;
    long raw_id;
    try {
        raw_id  = std::stol(id_str.c_str(), &pos);
    } catch (const std::out_of_range&) {
        throw OTCBadRequest() << "The ID portion of the source_id was too large. Found: " << id_str;
    } catch(std::invalid_argument&) {
        throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " <<  id_str;
    }
    if (pos < id_str.length()) {
        throw OTCBadRequest() << "Expecting the ID portion of the source_id to be numeric. Found: " <<  id_str;
    }
    auto id = to_OttId(raw_id);
    if (not id) {
        throw OTCBadRequest() << "The ID portion of the source_id was too large. Found: " << id_str;
    }
    return *id;
}

const RTRichTaxNode* taxon_from_source_id(const string& source_id, const RichTaxonomy& taxonomy) {
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
        if (auto taxon_node = taxonomy.included_taxon_from_id(in_ott)) {
            return taxon_node;
        } else {
            throw OTCBadRequest("Foreign ID '"+ source_id+"' mapped to unknown OTT ID: "+ to_string(in_ott));
        }
#   endif
}

string node_info_method_handler(const json& parsed_args) {
    string synth_id = extract_argument_or_default<string>(parsed_args, "synth_id", "");
    auto node_id = extract_argument<string>(parsed_args,"node_id");
    auto source_id = extract_argument<string>(parsed_args,"source_id");
    auto node_ids = extract_argument<vector<string>>(parsed_args,"node_ids");
    int count =0;
    count += bool(node_id)?1:0;
    count += bool(node_ids)?1:0;
    count += bool(source_id)?1:0;
    if (count != 1) {
        throw OTCBadRequest("Must supply exactly one of 'node_id', 'node_ids', or 'source_id'.");
    }
    if (source_id) {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        auto tax_node = taxon_from_source_id(*source_id, taxonomy);
        node_id = "ott"+std::to_string(tax_node->get_ott_id());
    }
    bool include_lineage = extract_argument_or_default<bool>(parsed_args, "include_lineage", false);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
    if (node_id) {
        return node_info_ws_method(tts, treeptr, sta, *node_id, include_lineage);
    } else {
        return nodes_info_ws_method(tts, treeptr, sta, *node_ids, include_lineage);
    }
}

string mrca_method_handler( const json& parsedargs)
{
    auto [synth_id, node_id_vec] = get_synth_and_node_id_vec(parsedargs);
    auto excluded_node_ids =  extract_argument_or_default<vector<string>>(parsedargs, "excluded_node_ids", {});
    auto soft_exclude =  extract_argument_or_default<bool>(parsedargs, "soft_exclude", false);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
    return mrca_ws_method(tts, treeptr, sta, node_id_vec, excluded_node_ids, soft_exclude);
}

std::string process_subtree(const json& parsedargs) {
    // FIXME: According to treemachine/ws-tests/tests.subtree, there is an "include_all_node_labels"
    //        argument.  Unless this is explicitly set to true, we are supposed to not write node labels
    //        for non-ottids.  At least in Newick.

    auto [synth_id, node_id] = get_synth_and_node_id(parsedargs);
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

string induced_subtree_method_handler( const json& parsedargs ) {
    auto [synth_id, node_id_vec] = get_synth_and_node_id_vec(parsedargs);
    NodeNameStyle nns = get_label_format(parsedargs);
    const SummaryTreeAnnotation * sta = get_annotations(tts, synth_id);
    const SummaryTree_t * treeptr = get_summary_tree(tts, synth_id);
    return induced_subtree_ws_method(tts, treeptr, node_id_vec, nns);
}

string tax_about_method_handler( const json& ) {
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tax_about_ws_method(taxonomy);
}

// looks for ott_id or source_id args to find a node
const RTRichTaxNode * extract_taxon_node_from_args(const json & parsedargs, const RichTaxonomy & taxonomy) {
    auto ott_id = extract_argument<OttId>(parsedargs, "ott_id");
    auto source_id = extract_argument<string>(parsedargs, "source_id");
    if (ott_id and source_id) {
        throw OTCBadRequest("'ott_id' and 'source_id' arguments cannot both be supplied.");
    } else if (not ott_id and not source_id) {
        throw OTCBadRequest("An 'ott_id' or 'source_id' argument is required.");
    }
    if (source_id) {
        return taxon_from_source_id(*source_id, taxonomy);
    } else {
        assert(ott_id);
        auto taxon_node = taxonomy.included_taxon_from_id(*ott_id);
        if (taxon_node == nullptr) {
            throw OTCBadRequest() << "Unrecognized OTT ID: " << *ott_id;
        }
        return taxon_node;
    }
}

string taxon_info_method_handler( const json& parsedargs ) {
    auto include_lineage = extract_argument_or_default<bool>(parsedargs, "include_lineage", false);
    auto include_children = extract_argument_or_default<bool>(parsedargs, "include_children", false);
    auto include_terminal_descendants = extract_argument_or_default<bool>(parsedargs, "include_terminal_descendants", false);       
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
    return taxon_info_ws_method(taxonomy, taxon_node, include_lineage, include_children, include_terminal_descendants);
}

string taxon_flags_method_handler( const json& ) {
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return taxonomy_flags_ws_method(taxonomy);
}

string taxon_mrca_method_handler( const json& parsedargs ) {
    OttIdSet ott_id_set = extract_required_argument<OttIdSet>(parsedargs, "ott_ids");
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return taxonomy_mrca_ws_method(taxonomy, ott_id_set);
}

string taxon_subtree_method_handler( const json& parsedargs ) {
    NodeNameStyle nns = get_label_format(parsedargs);
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    const RTRichTaxNode * taxon_node = extract_taxon_node_from_args(parsedargs, taxonomy);
    return taxon_subtree_ws_method(tts, taxonomy, taxon_node, nns);
}

// See taxomachine/src/main/java/org/opentree/taxonomy/plugins/tnrs_v3.java

// 10,000 queries at .0016 second per query = 16 seconds
const int MAX_NONFUZZY_QUERY_STRINGS = 10000;
// 250 queries at .3 second per query = 75 seconds
const int MAX_FUZZY_QUERY_STRINGS = 250;

static string LIFE_NODE_NAME = "life";

string tnrs_match_names_handler( const json& parsedargs ) {
    // 1. Requred argument: "names"
    vector<string> names = extract_required_argument<vector<string>>(parsedargs, "names");
    // 2. Optional argunments
    optional<string> context_name = extract_argument<string>(parsedargs, "context_name");
    bool do_approximate_matching  = extract_argument_or_default(parsedargs, "do_approximate_matching", false);
    vector<string> ids            = extract_argument_or_default(parsedargs, "ids",                     names);
    bool include_suppressed       = extract_argument_or_default(parsedargs, "include_suppressed",      false);

    // 3. Check that "ids" have the same length as "names", if supplied
    if (ids.size() != names.size()) {
        throw OTCBadRequest() << "The number of names and ids does not match. If you provide ids, then you "
                              << "must provide exactly as many ids as names.";
    }
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tnrs_match_names_ws_method(names, context_name, do_approximate_matching, ids, include_suppressed, taxonomy);
}

string tnrs_autocomplete_name_handler( const json& parsedargs ) {
    string name              = extract_required_argument<string>(parsedargs, "name");
    string context_name      = extract_argument_or_default(parsedargs, "context_name",            LIFE_NODE_NAME);
    bool include_suppressed  = extract_argument_or_default(parsedargs, "include_suppressed",      false);
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tnrs_autocomplete_name_ws_method(name, context_name, include_suppressed, taxonomy);
}

string tnrs_contexts_handler( const json& ) {
    return tnrs_contexts_ws_method();
}

string tnrs_infer_context_handler( const json& parsedargs ) {
    vector<string> names = extract_required_argument<vector<string>>(parsedargs, "names");
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    return tnrs_infer_context_ws_method(names, taxonomy);
}

string conflict_status_method_handler( const json& parsed_args ) {
    auto tree1newick = extract_argument<string>(parsed_args, "tree1newick");
    auto tree1 = extract_argument<string>(parsed_args, "tree1");
    string tree2 = extract_required_argument<string>(parsed_args, "tree2");

    const auto& summary = *tts.get_summary_tree("");
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;

    if (tree1newick) {
        return newick_conflict_ws_method(summary, taxonomy, *tree1newick, tree2);
    } else if (tree1) {
        return phylesystem_conflict_ws_method(summary, taxonomy, *tree1, tree2);
    } else {
        throw OTCBadRequest() << "Expecting argument 'tree1' or argument 'tree1newick'";
    }
}


// 
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/ctrie/str_utils.h"

using namespace std;
namespace fs = boost::filesystem;
using json = nlohmann::json;
typedef std::set<fs::path> fp_set;
typedef std::pair<bool, fp_set > bool_fp_set; 

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


/// formerly tolwsadaptors.h
inline otc::vec_src_node_ids extract_node_id_vec(otc::TreesToServe & tts,
                                                 const nlohmann::json & sbv
#                                                if defined(JOINT_MAPPING_VEC)
                                                   , otc::SourceEdgeMappingType semt
#                                                endif
                                                 ) {
#   if defined(JOINT_MAPPING_VEC)
       using lel_t = otc::semt_ind_t;
#   else
       using lel_t = std::uint32_t;
#   endif
    std::list<lel_t> lsni;
    for (nlohmann::json::const_iterator jit = sbv.begin(); jit != sbv.end(); ++jit) {
        const std::string * kp = tts.get_stored_string(jit.key());
        const auto & v = jit.value();
        for (nlohmann::json::const_iterator vit = v.begin(); vit != v.end(); ++vit) {
            const std::string * vp = tts.get_stored_string(*vit);
            const auto sni_ind = tts.get_source_node_id_index(otc::src_node_id(kp, vp));
#           if defined(JOINT_MAPPING_VEC)
                lsni.push_back(lel_t(semt, sni_ind));
#           else
                lsni.push_back(sni_ind);
#           endif
        } 
    }
    return otc::vec_src_node_ids(lsni.begin(), lsni.end());
}


inline const nlohmann::json & extract_obj(const nlohmann::json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw otc::OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_object()) {
        nlohmann::json & k = const_cast<nlohmann::json &>(j);
        return k[field];
    }
    throw otc::OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}

int run_server(const boost::program_options::variables_map & args);
boost::program_options::variables_map parse_cmd_line(int argc, char* argv[]);


/// end formerly tolwsadaptors.h
/// formerly tolwsadaptors.cpp

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
string ctime(const chrono::system_clock::time_point& t) {
    time_t t2 = chrono::system_clock::to_time_t(t);
    char* c = ctime(&t2);
    string tt = c;
    tt.pop_back(); // remove newline
    return tt;
}

multimap<string,string> request_headers(const string& rbody) {
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

multimap<string,string> options_headers() {
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

// OK, so I want to make something like OTCWebError, but add the ability to pass back
// some JSON along with the error.

// The difficult thing is how to generically (polymorphically) use the same interface
// for this error and other stuff.

std::string error_response(const string& path, const std::exception& e) {
    string msg = string("[") + path + ("] Error: ") + e.what();
    LOG(DEBUG)<<msg;
    json j = { {"message", msg} };
    return j.dump(4)+"\n";
}

std::string error_response(const string& path, const OTCWebError& e1) {
    OTCWebError e2 = e1;
    e2.prepend(string("[") + path + ("] Error: "));
    LOG(DEBUG)<<e2.what();
    return e2.json().dump(4)+"\n";
}

std::function<void(const shared_ptr< Session > session)>
create_method_handler(const string& path, const std::function<std::string(const json&)> process_request) {
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

json request_to_json(const Request& request) {
    LOG(DEBUG)<<"GET "<<request.get_path();
    json query;
    for(auto& key_value_pair: request.get_query_parameters()) {
        query[key_value_pair.first] = key_value_pair.second;
    }
    return query;
}

std::function<void(const shared_ptr< Session > session)>
create_GET_method_handler(const string& path, const std::function<std::string(const json&)> process_request) {
    return [=](const shared_ptr< Session > session ) {
        try {
            LOG(DEBUG)<<"request: "<<path;
            const auto& request = session->get_request( );
            json parsedargs = request_to_json(*request);
            LOG(DEBUG)<<"   argument "<<parsedargs.dump(1);
            auto rbody = process_request(parsedargs);
            LOG(DEBUG)<<"request: DONE";
            session->close( OK, rbody, request_headers(rbody) );
        } catch (OTCWebError& e) {
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

shared_ptr< Resource > path_handler(const string& path, std::function<std::string(const json &)> process_request) {
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


    const Context * c = determine_context({});
    if (c == nullptr) {
        throw OTCError() << "no context found for entire taxonomy";
    }
    ContextAwareCTrieBasedDB ct{*c, taxonomy};
    taxonomy.set_fuzzy_matcher(&ct);

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

/// end formerly tolwsadaptors.cpp

namespace otc {
void from_json(const nlohmann::json &j, SummaryTreeAnnotation & sta) {
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
    sta.root_ott_id = extract_ott_id(j, "root_ott_id");
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
        for (nlohmann::json::const_iterator sim_it = sim_el->begin(); sim_it != sim_el->end(); ++sim_it) {
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


bool read_tree_and_annotations(const fs::path & configpath,
                               const fs::path & treepath,
                               const fs::path & annotationspath,
                               const fs::path & brokentaxapath,
                               const fs::path & contestingtrees_path,
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
            fs::path brokentaxapath = p;
            brokentaxapath /= "labelled_supertree";
            brokentaxapath /= "broken_taxa.json";
            fs::path annotationspath = p;
            annotationspath /= "annotated_supertree";
            annotationspath /= "annotations.json";

            fs::path contestingtrees_path = p / "subproblems" / "contesting-trees.json";

            bool was_tree_par = false;
            try {
                if (fs::is_regular_file(treepath)
                    && fs::is_regular_file(annotationspath)
                    && fs::is_regular_file(configpath)) {
                    if (read_tree_and_annotations(configpath, treepath, annotationspath, brokentaxapath, contestingtrees_path, tts)) {
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
    tts.final_tree_added();
    return true;
}

#if defined(REPORT_MEMORY_USAGE)

template<>
inline std::size_t calc_memory_used(const RTRichTaxTreeData &d, MemoryBookkeeper &mb) {
    std::size_t sz_el_size = sizeof(OttId) + sizeof(const RTRichTaxNode *);
    std::size_t nm_sz = calc_memory_used_by_map_eqsize(d.ncbi_id_map, sz_el_size, mb);
    std::size_t gm_sz = calc_memory_used_by_map_eqsize(d.gbif_id_map, sz_el_size, mb);
    std::size_t wm_sz = calc_memory_used_by_map_eqsize(d.worms_id_map, sz_el_size, mb);
    std::size_t fm_sz = calc_memory_used_by_map_eqsize(d.if_id_map, sz_el_size, mb);
    std::size_t im_sz = calc_memory_used_by_map_eqsize(d.irmng_id_map, sz_el_size, mb);
    std::size_t f2j_sz = calc_memory_used_by_map_simple(d.flags2json, mb);
    std::size_t in_sz = calc_memory_used_by_map_simple(d.id_to_node, mb);
    std::size_t nn_sz = calc_memory_used_by_map_simple(d.name_to_node, mb);
    std::size_t nutn_sz = calc_memory_used_by_map_simple(d.non_unique_taxon_names, mb);
    std::size_t htn_sz = 0;
    for (auto el : d.homonym_to_node) {
        htn_sz += sizeof(std::string_view);
        htn_sz += calc_memory_used_by_vector_eqsize(el.second, sizeof(const RTRichTaxNode *), mb);
    }
    mb["taxonomy data ncbi map"] += nm_sz;
    mb["taxonomy data gbif map"] += gm_sz;
    mb["taxonomy data worms map"] += wm_sz;
    mb["taxonomy data indexfungorum map"] += fm_sz;
    mb["taxonomy data irmng map"] += im_sz;
    mb["taxonomy data flags2json"] += f2j_sz;
    mb["taxonomy data id_to_node"] += in_sz;
    mb["taxonomy data name_to_node"] += nn_sz;
    mb["taxonomy data non_unique_taxon_names"] += nutn_sz;
    mb["taxonomy data homonym_to_node"] += htn_sz;
    return nm_sz + gm_sz + wm_sz + fm_sz + im_sz + f2j_sz + in_sz + nn_sz + nutn_sz + htn_sz;
}

template<>
inline std::size_t calc_memory_used(const TaxonomicJuniorSynonym &d, MemoryBookkeeper &mb) {
    return calc_memory_used(d.name, mb) + calc_memory_used(d.source_string, mb) + sizeof(char *);
}

template<>
inline std::size_t calc_memory_used(const std::list<TaxonomicJuniorSynonym> &d, MemoryBookkeeper &mb) {
    std::size_t total = sizeof(size_t) * (2 + 2*d.size()) ; // start and capacity
    for (auto i = d.begin(); i != d.end(); ++i) {
        total += calc_memory_used(*i, mb);
    }
    return total;
}
template<>
inline std::size_t calc_memory_used(const RTRichTaxNodeData &rtn, MemoryBookkeeper &mb) {
    size_t total = 0;
    size_t x = calc_memory_used_by_vector_eqsize(rtn.junior_synonyms, sizeof(char *), mb);
    mb["taxonomy node data junior_synonyms"] += x; total += x;
    x = 2*sizeof(std::uint32_t);
    mb["taxonomy node data traversal"] += x; total += x;
    x = sizeof(TaxonomicRank);
    mb["taxonomy node data rank"] += x; total += x;
    x = sizeof(int32_t) + sizeof(std::bitset<32>);
    mb["taxonomy node data flags"] += x; total += x;
    x = calc_memory_used(rtn.source_info, mb);
    mb["taxonomy node data source_info"] += x; total += x;
    x = sizeof(std::string_view);
    mb["taxonomy node data nonunique name"] += x; total += x;
    return total;
}

template<>
inline std::size_t calc_memory_used(const RichTaxonomy &rt, MemoryBookkeeper &mb) {
    const auto & tax_tree = rt.get_tax_tree();
    std::size_t ttsz = calc_memory_used_by_tree(tax_tree, mb);
    const auto & syn_list = rt.get_synonyms_list();
    std::size_t slsz = calc_memory_used(syn_list, mb);
    const auto & suppress_trns_set = rt.get_ids_to_suppress_from_tnrs();
    std::size_t stssz = calc_memory_used(suppress_trns_set, mb);
    mb["taxonomy tree"] += ttsz;
    mb["taxonomy synonyms list"] += slsz;
    mb["taxonomy set of ids to suppress from tnrs"] += stssz;
    return ttsz + slsz + stssz;
}

#endif

void mark_summary_tree_nodes_extinct(SummaryTree_t& tree, const RichTaxonomy& taxonomy)
{
    // compute extinctness for each node.  Post means that a node is only visited after all its children.
    for (auto node: iter_post(tree))
    {
        auto& node_data = node->get_data();
        if (node->is_tip())
        {
            auto id = node->get_ott_id();
            auto& taxon = taxonomy.included_taxon_from_id(id)->get_data();
            node_data.extinct_mark = taxon.is_extinct();
        }
        else
        {
            // If any child is not extinct, then this node is not extinct either.
            node_data.extinct_mark = true;
            for (auto c : iter_child_const(*node))
                if (not c->get_data().is_extinct())
                    node_data.extinct_mark = false;

            // Complain about higher taxa with extinctness that doesn't match the computed extinctness.
            if (node->has_ott_id())
            {
                auto id = node->get_ott_id();
                auto& taxon = taxonomy.included_taxon_from_id(id)->get_data();
                if (node_data.is_extinct() != taxon.is_extinct())
                {
                    LOG(WARNING)<<"Higher taxon "<<taxon.possibly_nonunique_name<<" is extinct="<<taxon.is_extinct()<<"  but the computed extinctness is extinct="<<node_data.is_extinct();
                    for (auto c : iter_child_const(*node))
                        if (not c->get_data().is_extinct())
                            LOG(WARNING)<<"    Child "<<c->get_name()<<" is NOT extinct!";
                        else
                            LOG(WARNING)<<"    Child "<<c->get_name()<<" is EXTINCT!";
                }
            }
        }
    }
}

bool read_tree_and_annotations(const fs::path & config_path,
                               const fs::path & tree_path,
                               const fs::path & annotations_path,
                               const fs::path & brokentaxa_path,
                               const fs::path & contestingtrees_path,
                               TreesToServe & tts)
{
    std::ifstream contestingtrees_stream(contestingtrees_path.native().c_str());
    json contestingtrees_obj;
    try {
        contestingtrees_stream >> contestingtrees_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << contestingtrees_path << "\" as JSON.\n";
        throw;
    }

    std::string annot_str = annotations_path.native();
    std::ifstream annotations_stream(annot_str.c_str());
    json annotations_obj;
    try {
        annotations_stream >> annotations_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << annotations_path << "\" as JSON.\n";
        throw;
    }
    std::string bt_str = brokentaxa_path.native();
    std::ifstream brokentaxa_stream(bt_str.c_str());
    json brokentaxa_obj;
    try {
        brokentaxa_stream >> brokentaxa_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << brokentaxa_path << "\" as JSON.\n";
        throw;
    }
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
#   if defined(REPORT_MEMORY_USAGE)
        MemoryBookkeeper tax_mem_b;
        std::size_t tree_mem = 0;
        auto tax_mem = calc_memory_used(taxonomy, tax_mem_b);
        write_memory_bookkeeping(LOG(INFO), tax_mem_b, "taxonomy", tax_mem);
#   endif
    auto [tree,sta] = tts.get_new_tree_and_annotations(config_path.native(), tree_path.native());
    try {

        mark_summary_tree_nodes_extinct(tree, taxonomy);

        sta = annotations_obj;
        json tref;
        tref["taxonomy"] = taxonomy.get_version();
        sta.full_source_id_map_json[taxonomy.get_version()] = tref;

        // read the node annotations and add them to the tree.
        auto node_obj = extract_obj(annotations_obj, "nodes");
        auto & sum_tree_data = tree.get_data();
        //const auto & n2n = sum_tree_data.name_to_node;
        for (json::const_iterator nit = node_obj.begin(); nit != node_obj.end(); ++nit) {
            string k = nit.key();
            auto result = find_node_by_id_str(tree, k);
            //auto stnit = n2n.find(k);
            if (result.node == nullptr) {
                throw OTCError() << "Node " << k << " from annotations not found in tree.";
            }
            //const SumTreeNode_t * stn = stnit->second;
            SumTreeNode_t * mstn = const_cast<SumTreeNode_t *>(result.node);
            SumTreeNodeData & mstnd = mstn->get_data();
            const auto & supportj = nit.value();
#           if defined(JOINT_MAPPING_VEC)
                vec_src_node_ids tmpv;
#           endif
            for (json::const_iterator sbit = supportj.begin(); sbit != supportj.end(); ++sbit) {
                const auto & sbk = sbit.key();
                const auto & sbv = sbit.value();
                if (sbk == "supported_by") {            
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::SUPPORTED_BY_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.supported_by = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "terminal") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::TERMINAL_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.terminal = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "conflicts_with") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::CONFLICTS_WITH_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.conflicts_with = extract_node_id_vec(tts, sbv)
#                   endif
                } else if (sbk == "partial_path_of") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::PARTIAL_PATH_OF_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.partial_path_of = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "resolves") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::RESOLVES_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.resolves = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "was_uncontested") {
                    if (sbv.is_boolean()) {
                        mstnd.was_uncontested =  sbv.get<bool>();
                    } else {
                        throw OTCError() << "Expected was_uncontested to be a boolean.";
                    }
                } else {
                    if (sbk != "was_constrained") {
                        throw OTCError() << "Unrecognized annotations key " << sbit.key();
                    }
                }
            }
#           if defined(JOINT_MAPPING_VEC)
                mstnd.source_edge_mappings.clear();
                std::swap(mstnd.source_edge_mappings, tmpv);
                //LOG(INFO) << "mstnd.source_edge_mappings size = " << mstnd.source_edge_mappings.size() << " for k = " << k;
#           endif

        }

        // Read in the 'contesting-trees.json' file.  We aren't using all the info yet.
        auto& contesting_trees_for_taxon = sum_tree_data.contesting_trees_for_taxon;
        for(auto& [taxon,trees]: contestingtrees_obj.items())
        {
            vector<contesting_tree_t> contesting_trees;
            for(auto& [tree,attachment_points_json]: trees.items())
            {
                contesting_tree_t contesting_tree;
                // Remove extension ".tre"
                assert(tree.substr(tree.size()-4) == ".tre");
                contesting_tree.tree = tree.substr(0,tree.size()-4);
                for(auto& attachment_point_json: attachment_points_json)
                {
                    attachment_point_t A;
                    assert(attachment_point_json.count("parent"));
                    assert(attachment_point_json.count("children_from_taxon")>0);

                    A.parent = attachment_point_json["parent"].get<string>();
                    A.parent = strip_surrounding_whitespace(A.parent);
                    A.parent = get_source_node_name_if_available(A.parent);

                    for(auto& child_json: attachment_point_json["children_from_taxon"])
                    {
                        string child = child_json.get<string>();
                        child = strip_surrounding_whitespace(child);
                        child = get_source_node_name_if_available(child);
                        A.children_from_taxon.push_back(child);
                    }
                    contesting_tree.attachment_points.push_back(A);
                }
                contesting_trees.push_back(contesting_tree);
            }
            contesting_trees_for_taxon.insert({taxon, contesting_trees});
        }

        auto & tree_broken_taxa = sum_tree_data.broken_taxa;
        // read the info from the broken taxa file
        if (brokentaxa_obj.count("non_monophyletic_taxa")
            && (!brokentaxa_obj["non_monophyletic_taxa"].is_null())) {
            auto & nmt_obj = extract_obj(brokentaxa_obj, "non_monophyletic_taxa");
            for (json::const_iterator btit = nmt_obj.begin(); btit != nmt_obj.end(); ++btit) {
                string broken_ott = btit.key();
                auto & dest_obj = btit.value();
                string mrca_id = extract_string(dest_obj, "mrca");
                auto & attach_obj = extract_obj(dest_obj, "attachment_points");
                list<string> attach_id_list;
                for (json::const_iterator ai_it = attach_obj.begin(); ai_it != attach_obj.end(); ++ai_it) {
                    attach_id_list.push_back(ai_it.key());
                }
                auto mrca_result = find_node_by_id_str(tree, mrca_id);
                vector<const SumTreeNode_t *> avec;
                avec.reserve(attach_id_list.size());
                for (auto attach_id : attach_id_list) {
                    auto [anptr, _] = find_node_by_id_str(tree, attach_id);
                    assert(anptr != nullptr);
                    avec.push_back(anptr);
                }
                tree_broken_taxa[broken_ott] = BrokenMRCAAttachVec(mrca_result.node, avec);
            }
        }
        sta.initialized = true;
        tts.register_last_tree_and_annotations();
#       if defined(REPORT_MEMORY_USAGE)
            MemoryBookkeeper tree_mem_b;
            tree_mem += calc_memory_used_by_tree(tree, tree_mem_b);
            write_memory_bookkeeping(LOG(INFO), tree_mem_b, "tree", tree_mem);
            LOG(INFO) << "tax + tree memory = " << tax_mem << " + " << tree_mem << " = " << tax_mem + tree_mem;
#       endif
    } catch (...) {
        tts.free_last_tree_and_annotations();
        throw;
    }
    return true;
}

}// namespace otc

int main( const int argc, char** argv) {
    if (otc::set_global_conv_facet() != 0) {
        return 1;
    }
    
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc,argv);
        return run_server(args);
    } catch (std::exception& e) {
        LOG(ERROR) <<"otc-tol-ws: Error! " << e.what() << std::endl;
        return 1;
    }
}
