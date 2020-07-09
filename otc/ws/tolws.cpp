#include <regex>
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/ws/node_namer_supported_by_stasher.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/node_naming.h"       // for make_name(prefix,OttId)
#include <optional>
#include <string_view>
#include "otc/tnrs/context.h"

#include "otc/ctrie/str_utils.h"


using std::vector;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::optional;
using std::ostringstream;
using std::unique_ptr;

using namespace boost::property_tree;
using json=nlohmann::json;


using otc::OttId;

namespace otc {
const int OK = restbed::OK;

string get_synth_node_label(const SumTreeNode_t* node)
{
    if (node->has_ott_id())
        return make_name("ott", node->get_ott_id());
    else
        return node->get_name();
}

string_view taxon_nonuniquename(const RichTaxonomy& taxonomy, const SumTreeNode_t& nd)
{
    if (not nd.has_ott_id())
        throw OTCError()<<"Node "<<nd.get_name()<<" has no OTT id";

    auto id = nd.get_ott_id();
    auto nd_taxon = taxonomy.included_taxon_from_id(id);
    auto& taxon_data = nd_taxon->get_data();
    return taxon_data.get_nonuniqname();
}

// Corresponds to getNamesOfRepresentativeDescendants( ) in treemachine/src/main/java/opentree/GraphExplorer.java
// FIXME - Clean up to eliminate recursive calls.
// FIXME - Looking at more names on each level seems better because it would find higher-ranking descendant names
//       - is there a reason we weren't doing this?  e.g. scanning all children could go too slow?
// FIXME - We could cache representative descendants (which takes memory) if too slow.
json get_descendant_names(const RichTaxonomy& taxonomy, const SumTreeNode_t& nd)
{
    json names = json::array();
    if (nd.has_children())
    {
        auto first = nd.get_first_child();
        auto last  = nd.get_last_child();
        if (first->has_ott_id())
            names.push_back(taxon_nonuniquename(taxonomy, *first));
        else
        {
            auto names2 = get_descendant_names(taxonomy, *first);
            if (not names2.empty())
                names.push_back(names2[0]);
        }

        if (last != first)
        {
            if (last->has_ott_id())
                names.push_back(taxon_nonuniquename(taxonomy, *last));
            else
            {
                auto names3 = get_descendant_names(taxonomy, *last);
                if (not names3.empty())
                    names.push_back(names3.back());
            }
        }
    }
    return names;
}

// Corresponds to getNodeBlob( ) and getNodeBlobArguson( ) in treemachine/src/main/java/opentree/GraphExplorer.java
void add_basic_node_info(const RichTaxonomy & taxonomy, const SumTreeNode_t & nd, json & noderepr, bool is_arguson = false) {
    noderepr["node_id"] = node_id_for_summary_tree_node(nd);

    // The number of descendant tips (.e.g not including this node).
    if (nd.is_tip())
        noderepr["num_tips"] = 0;
    else
        noderepr["num_tips"] = nd.get_data().num_tips;

    if (is_arguson)
        noderepr["extinct"] = nd.get_data().is_extinct();

    if (nd.has_ott_id()) {
        auto nd_id = nd.get_ott_id();
        const auto * nd_taxon = taxonomy.included_taxon_from_id(nd_id);
        if (nd_taxon == nullptr) {
            throw OTCError() << "OTT Id " << nd_id << " not found in taxonomy! Please report this bug";
        }
        json taxon;
        add_taxon_info(taxonomy, *nd_taxon, taxon);
        noderepr["taxon"] = taxon;
    }
    else if (is_arguson)
        noderepr["descendant_name_list"] = get_descendant_names(taxonomy, nd);
}

inline void add_str_to_vec_string(json & o, const string& first, const string& second) {
    const char * studyc = first.c_str();
    if (not o.count(studyc))
        o[studyc] = json::array();
    o[studyc].push_back(second);
}

inline void add_str_to_str(json & o, const string& first, const string& second) {
    const char * studyc = first.c_str();
    assert(not o.count(studyc) or o[studyc].get<string>() == second);
    o[studyc] = second;
}

inline void add_str_to_str_or_vec_string(json & o, const string& first, const string& second) {
    const char * studyc = first.c_str();
    if (o.count(studyc)) {
        if (o[studyc].is_array()) {
            o[studyc].push_back(second);
        } else {
            string prev = o[studyc].get<string>();
            json v = {prev, second};
            o[studyc] = v;
        }
    } else {
        o[studyc] = second;
    }
}

#if !defined(JOINT_MAPPING_VEC)
void add_support_info_vec(const TreesToServe & tts, 
                          const char * tag,
                          const vec_src_node_ids & v,
                          json & noderepr,
                          set<string> & usedSrcIds,
                          const optional<string>& extra_src = {},
                          const optional<string>& extra_node_id = {}) {
    json o;
    for (const auto & sni : v ) {
        const auto study_node_pair = tts.decode_study_node_id_index(sni);
        usedSrcIds.insert(*study_node_pair.first);
        add_str_to_vec_string(o, study_node_pair.first, study_node_pair.second);
    }
    if (extra_src && extra_node_id) {
        add_str_to_vec_string(o, *extra_src, *extra_node_id);
    }
    noderepr[tag] = o;
}

void add_support_info_single_element(const TreesToServe & tts,
                                     const char * tag,
                                     const vec_src_node_ids & v,
                                     json & noderepr,
                                     set<string> & usedSrcIds,
                                     const optional<string>& extra_src = {},
                                     const optional<string>& extra_node_id = {}) {
    json o;
    for (const auto & sni : v ) {
        const auto study_node_pair = tts.decode_study_node_id_index(sni);
        usedSrcIds.insert(*study_node_pair.first);
        add_str_to_str(o, study_node_pair.first, study_node_pair.second);
    }
    if (extra_src && extra_node_id) {
        add_str_to_str(o, *extra_src, *extra_node_id);
    }
    noderepr[tag] = o;
}
#endif

void add_node_support_info(const TreesToServe & tts,
                           const SumTreeNode_t & nd,
                           json & noderepr,
                           set<string> & usedSrcIds) {
    const auto & d = nd.get_data();
    optional<string> extra_src;
    optional<string> extra_node_id;
    string tmp;
    if (nd.has_ott_id()) {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        extra_src = string("ott") + taxonomy.get_version();
        extra_node_id = node_id_for_summary_tree_node(nd);
        usedSrcIds.insert(*extra_src);
    }

#if defined(JOINT_MAPPING_VEC)
    json supported_j; bool had_conflicts = false;
    json conflicts_j; bool had_supported = false;
    json partial_path_j; bool had_partial_path = false;
    json resolves_j; bool had_resolves = false;
    json terminal_j; bool had_terminal = false;
    for (auto el : d.source_edge_mappings) {
        const auto study_node_pair = tts.decode_study_node_id_index(el.second);
        usedSrcIds.insert(*study_node_pair.first);
        switch (el.first) {
            // array
            case SourceEdgeMappingType::CONFLICTS_WITH_MAPPING:
                add_str_to_vec_string(conflicts_j, *study_node_pair.first, *study_node_pair.second);
                had_conflicts = true;
                break;
            // single element
            case SourceEdgeMappingType::PARTIAL_PATH_OF_MAPPING:
                add_str_to_str(partial_path_j, *study_node_pair.first, *study_node_pair.second);
                had_partial_path = true;
                break;
            // single element
            case SourceEdgeMappingType::RESOLVES_MAPPING:
                add_str_to_str(resolves_j, *study_node_pair.first, *study_node_pair.second);
                had_resolves = true;
                break;
            // single element
            case SourceEdgeMappingType::SUPPORTED_BY_MAPPING:
                add_str_to_str(supported_j, *study_node_pair.first, *study_node_pair.second);
                had_supported = true;
                break;
            // single element
            case SourceEdgeMappingType::TERMINAL_MAPPING:
                add_str_to_str(terminal_j, *study_node_pair.first, *study_node_pair.second);
                had_terminal = true;
                break;
            // resolved_by: array
        }
    }
    if (extra_src && extra_node_id) {
        add_str_to_str(supported_j, *extra_src, *extra_node_id);
        had_supported = true;
    }
    if (had_conflicts) {
        noderepr["conflicts_with"] = conflicts_j;
    }
    if (had_partial_path) {
        noderepr["partial_path_of"] = partial_path_j;
    }
    if (had_resolves) {
        noderepr["resolves"] = resolves_j;
    }
    if (had_supported) {
        noderepr["supported_by"] = supported_j;
    }
    if (had_terminal) {
        noderepr["terminal"] = terminal_j;
    }

#else
    if (extra_src || !d.supported_by.empty()) {
        add_support_info_single_element("supported_by", d.supported_by, noderepr, usedSrcIds, extra_src, extra_node_id);
    }
    if (!d.conflicts_with.empty()) {
        add_support_info_vec("conflicts_with", d.conflicts_with, noderepr, usedSrcIds);
    }
    if (!d.resolves.empty()) {
        add_support_info_single_element("resolves", d.resolves, noderepr, usedSrcIds);
    }
    if (!d.partial_path_of.empty()) {
        add_support_info_single_element("partial_path_of", d.partial_path_of, noderepr, usedSrcIds);
    }
    if (!d.terminal.empty()) {
        add_support_info_single_element("terminal", d.terminal, noderepr, usedSrcIds);
    }
#endif
    if (d.was_uncontested) {
        noderepr["was_uncontested"] = true;
        noderepr["was_constrained"] = true;
    }
}


const std::regex ott_id_pattern("^ott(\\d+)$");

optional<OttId> is_ott_id(const string& node_id)
{
    std::smatch matches;
    if (std::regex_match(node_id, matches, ott_id_pattern))
    {
        long raw_ott_id = long_ott_id_from_name(node_id);
        if (raw_ott_id >= 0)
            return to_OttId(raw_ott_id);
    }
    return {};
}

const std::regex mrca_id_pattern("^mrca(ott\\d+)(ott\\d+)$");

node_lookup_t broken_taxon(const SumTreeNode_t* n)
{
    node_lookup_t result(n);
    result.was_broken = true;
    return result;
}

template<typename N>
bool is_ancestor_of(N* f, N* s)
{
    assert(f);
    assert(s);
    const auto * fdata = &(f->get_data());
    const auto sec_ind = s->get_data().trav_enter;
    return (sec_ind >= fdata->trav_enter and sec_ind <= fdata->trav_exit);
}

node_lookup_t find_node_by_id_str(const SummaryTree_t & tree, const string & node_id)
{
    const auto & tree_data = tree.get_data();

    std::smatch matches;
    if (std::regex_match(node_id, matches, ott_id_pattern))
    {
        // Get the OTT ID
        long raw_ott_id = long_ott_id_from_name(node_id);

        // Try to find the OTT ID in the summary tree.
        if (raw_ott_id >= 0) {
            OttId ott_id = check_ott_id_size(raw_ott_id);
            auto i2nit = tree_data.id_to_node.find(ott_id);
            if (i2nit != tree_data.id_to_node.end()) {
                return {i2nit->second};
            }
            LOG(WARNING) << "not finding " << ott_id << " extracted from " << node_id;
        }

        // We didn't find a summary tree node for this OTT ID.  Is this node listed as broken?
        if (auto bt_it = tree_data.broken_taxa.find(node_id); bt_it != tree_data.broken_taxa.end())
        {
            // if found we return the MRCA pointer in the first slot of the pair.
            return broken_taxon(bt_it->second.first);
        }
        else
            return {};
    }

    if (std::regex_match(node_id, matches, mrca_id_pattern))
    {
        // Does this match a canonical mrcaottXottY name?
        auto n2nit = tree_data.broken_name_to_node.find(node_id);
        if (n2nit != tree_data.broken_name_to_node.end()) {
            return {n2nit->second};
        }

        // If it does not match a canonical name, then lookup ottX and ottY
        assert(matches.size() >= 2);
        std::string first_id = matches[1];
        std::string second_id = matches[2];
        auto result1 = find_node_by_id_str(tree, first_id);
        if (result1.node == nullptr) {
            return result1;
        }
        auto result2 = find_node_by_id_str(tree, second_id);
        if (result2.node == nullptr) {
            return result2;
        }
        return {find_mrca_via_traversal_indices(result1.node, result2.node)};
    }
    return {};
}

node_lookup_t find_required_node_by_id_str(const SummaryTree_t & tree, const string & node_id)
{
    auto result = find_node_by_id_str(tree, node_id);
    if (not result.node) {
        throw OTCBadRequest() << "node_id '" << node_id << "' was not found!";
    }
    return result;
}

// See API docs at https://github.com/OpenTreeOfLife/germinator/wiki/Synthetic-tree-API-v3

string available_trees_ws_method(const TreesToServe &tts) {
    json response;
    json trees = json::array();
    for(auto& synth_id: tts.get_available_trees()) {
        trees.push_back(synth_id);
    }
    response["synth_ids"] = trees;
    response["default"] = tts.get_default_tree();
    return response.dump(1);
}

string about_ws_method(const TreesToServe &tts,
                       const SummaryTree_t * tree_ptr,
                       const SummaryTreeAnnotation * sta,
                       bool include_sources) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    json response;
    response["date_created"] = sta->date_completed;
    response["num_source_trees"] = sta->num_source_trees;
    response["num_source_studies"] = sta->num_source_studies;
    response["taxonomy_version"] = sta->taxonomy_version;
    response["filtered_flags"] = sta->filtered_flags_vec;
    response["synth_id"] = sta->synth_id;
    if (include_sources) {
        response["source_id_map"] = sta->full_source_id_map_json;
        response["source_list"] = sta->sources;
    }
    json root;
    auto root_node = tree_ptr->get_root();
    {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        add_basic_node_info(taxonomy, *root_node, root);
    }
    response["root"] = root;
    return response.dump(1);
}

json tax_about_json(const RichTaxonomy & taxonomy) {
    json response;
    string weburl;
    response["author"] = "open tree of life project";
    response["name"] = "ott";
    weburl = "https://tree.opentreeoflife.org/about/taxonomy-version/ott";
    weburl += taxonomy.get_version_number();
    response["source"] = string("ott") + taxonomy.get_version();
    response["version"] = taxonomy.get_version_number();
    response["weburl"] = weburl;
    return response;
}

string tax_about_ws_method(const RichTaxonomy & taxonomy) {
    json response = tax_about_json(taxonomy);
    return response.dump(1);
}


inline void add_lineage(const TreesToServe & tts,
                        json & j,
                        const SumTreeNode_t * focal,
                        const RichTaxonomy & taxonomy,
                        set<string> & usedSrcIds, bool is_arguson = false) {
    json lineage_arr;
    const SumTreeNode_t * anc = focal->get_parent();
    if (!anc) {
        vector<string> c;
        j["lineage"] = c;
        return;
    }
    while (anc) {
        json ancj;
        add_basic_node_info(taxonomy, *anc, ancj, is_arguson);
        add_node_support_info(tts, *anc, ancj, usedSrcIds);
        lineage_arr.push_back(ancj);
        anc = anc->get_parent();
    }
    j["lineage"] = lineage_arr;
}

inline void add_source_id_map(json & j,
                              const set<string> & usedSrcIds,
                              const RichTaxonomy & taxonomy,
                              const SummaryTreeAnnotation * sta ) {
    json sim;
    for (auto srcTag : usedSrcIds) {
        json jt;
        if (srcTag == string("ott")+taxonomy.get_version()) {
            jt["taxonomy"] = string("ott")+taxonomy.get_version();
        } else {
            const auto & simentry = sta->source_id_map.at(srcTag);
            jt = simentry;
        }
        sim[srcTag] = jt;   
    }
    j["source_id_map"] = sim;
}

json node_info_json(const TreesToServe & tts,
                    const SummaryTreeAnnotation * sta,
                    const SumTreeNode_t* focal,
                    bool include_lineage)
{
    assert(focal != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;

    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        add_basic_node_info(taxonomy, *focal, response);
        add_node_support_info(tts, *focal, response, usedSrcIds);
        if (include_lineage) {
            add_lineage(tts, response, focal, taxonomy, usedSrcIds);
        }
        add_source_id_map(response, usedSrcIds, taxonomy, sta);
    }
    if (was_broken) {
        response["response_for_mrca_of_broken_taxon"] = true; //@TODO: discuss and document
    }
    return response;
}

string node_info_ws_method(const TreesToServe & tts,
                           const SummaryTree_t * tree_ptr,
                           const SummaryTreeAnnotation * sta,
                           const string & node_id,
                           bool include_lineage) {
    LOG(DEBUG)<<"Got to node_info_ws_method( )";

    auto result = find_required_node_by_id_str(*tree_ptr, node_id);

    auto response = node_info_json(tts, sta, result.node, include_lineage);
    response["query"] = node_id;
    return response.dump(1);
}

pair<vector<const SumTreeNode_t*>,json> find_nodes_for_id_strings(const RichTaxonomy& taxonomy,
                                                                  const SummaryTree_t* tree_ptr,
                                                                  const vector<string>& node_ids,
                                                                  bool fail_broken = false) {
    vector<const SumTreeNode_t *> nodes;
    json unknown;
    json broken = json::object();
    optional<string> bad_node_id;
    for (auto node_id : node_ids) {
        auto result = find_node_by_id_str(*tree_ptr, node_id);
        if (not result.node or (result.was_broken and fail_broken)) {
            // Possible statuses:
            //  - invalid    (never minted id)
            //  - pruned     (valid but not in synth)
            //  - deprecated (previously valid, not forwarded)
            string reason = "unknown_id";
            if (result.node and result.was_broken) {
                reason = "broken";
            } else if (auto id = is_ott_id(node_id)) {
                auto taxon = taxonomy.included_taxon_from_id(*id);
                // Not currently implemented...
                // if (id == -2)
                //  reason = "deprecated";
                reason = (not taxon ? "invalid_ott_id": "pruned_ott_id");
            }
            unknown[node_id] = reason;
            bad_node_id = node_id;
        }
        if (result.was_broken) {
            broken[node_id] = node_id_for_summary_tree_node(*result.node);
        }
        // Current default strategy means that we include MRCAs for broken taxa.
        nodes.push_back(result.node);
    }
    if (unknown.size()) {
        throw OTCBadRequest()<<"node_id '"<< *bad_node_id << "' was not found!"<<json{ {"unknown", unknown} };
    }
    return {nodes, broken};
}

string nodes_info_ws_method(const TreesToServe & tts,
                            const SummaryTree_t * tree_ptr,
                            const SummaryTreeAnnotation * sta,
                            const vector<string> & node_ids,
                            bool include_lineage) {
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    auto [nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_ids);
    json response;
    for(auto i = 0U; i < nodes.size(); i++) {
        auto j = node_info_json(tts, sta, nodes[i], include_lineage);
        j["query"] = node_ids[i];
        response.push_back(j);
    }
    return response.dump(1);
}

void add_nearest_taxon(const RichTaxonomy& taxonomy, const SumTreeNode_t& node, json& j) {
    if (node.has_ott_id()) {
        return;
    }
    auto anc = node.get_parent();
    assert(anc != nullptr);
    while (!anc->has_ott_id()) {
        anc = anc->get_parent();
        if (anc == nullptr) {
            throw OTCWebError("No ancestors were taxa. That is odd.\n");
        }
    }
    const RTRichTaxNode * anc_taxon = taxonomy.included_taxon_from_id(anc->get_ott_id());
    if (anc_taxon == nullptr) {
        throw OTCWebError() << "Anc OTT Id " << anc->get_ott_id() << " not found in taxonomy! Please report this bug";
    }
    json nt;
    add_taxon_info(taxonomy, *anc_taxon, nt);
    j["nearest_taxon"] = nt;
}

template <typename T>
const SumTreeNode_t * mrca(const T & nodes) {
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    for (auto n : nodes) {
        if (first) {
            first = false;
            focal = n;
        } else {
            focal = find_mrca_via_traversal_indices(focal, n);
            if (focal == nullptr) {
                break;
            }
        }
    }
    return focal;
}


// Check if the excluded node is inside the include group, and update the MRCA of the exclude group if not.
bool check_node_and_update_excluded_ancestor(const SumTreeNode_t* mrca_included, const SumTreeNode_t* excluded_node, const SumTreeNode_t* &closest_excluded_ancestor)
{
    auto mrca = find_mrca_via_traversal_indices(mrca_included, excluded_node);

    // Return false if the excluded node is actually inside the include group.
    if (mrca == mrca_included)
        return false;

    assert(not closest_excluded_ancestor or is_ancestor_of(closest_excluded_ancestor, mrca_included) or is_ancestor_of(mrca_included, closest_excluded_ancestor));

    // If this is the first excluded taxon OR the mrca is closer, the update the closest excluded ancestor
    if (not closest_excluded_ancestor or is_ancestor_of(closest_excluded_ancestor, mrca))
        closest_excluded_ancestor = mrca;

    return true;
}


string mrca_ws_method(const TreesToServe & tts,
                      const SummaryTree_t * tree_ptr,
                      const SummaryTreeAnnotation * sta,
                      const vector<string> & node_id_vec,
                      const vector<string> & excluded_node_ids,
                      bool soft_exclude)
{
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);

    auto [taxonomy,_] = tts.get_readable_taxonomy();

    // 1. Find MRCA of include group
    auto [tip_nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_id_vec);

    auto mrca_included = mrca(tip_nodes);

    if (not mrca_included) throw OTCBadRequest("MRCA of taxa was not found.\n");

    const SumTreeNode_t* closest_excluded_ancestor = nullptr;

    // 2. Find MRCA of exclude group
    auto [excluded_nodes, broken_excluded] = find_nodes_for_id_strings(taxonomy, tree_ptr, excluded_node_ids, false);
    json reversals;
    for(int i=0;i<excluded_nodes.size();i++)
    {
        if (not check_node_and_update_excluded_ancestor(mrca_included, excluded_nodes[i], closest_excluded_ancestor))
            reversals.push_back(excluded_node_ids[i]);
    }

    json response;
    if (not reversals.empty())
        response["reversals"] = reversals;
    response["synth_id"] = sta->synth_id;

    // 3a. Error: check for hard reversals.
    if (not reversals.empty() and not soft_exclude)
        throw OTCWebError(404)<<"Excluded taxa were nested within the include group!"<<response;

    // 3b. Error: check for empty exclude group after soft reversals.
    if (not excluded_node_ids.empty() and closest_excluded_ancestor == nullptr)
        throw OTCWebError(404)<<"All excluded taxa were nested within the include group!";
    
    // 4. Do the phylo-ref response here, if we had any excluded ids.
    if (not excluded_node_ids.empty())
    {
        assert(is_ancestor_of(closest_excluded_ancestor, mrca_included));
        json nodes;
        for(auto node = mrca_included; node != closest_excluded_ancestor; node = node->get_parent())
            nodes.push_back(node_id_for_summary_tree_node(*node));
        response["node_ids"] = nodes;
        return response.dump(1);
    }

    // 5. Do the standard mrca response.
    json mrcaj;
    set<string> usedSrcIds;
    add_node_support_info(tts, *mrca_included, mrcaj, usedSrcIds);
    add_basic_node_info(taxonomy, *mrca_included, mrcaj);
    add_nearest_taxon(taxonomy, *mrca_included, response);
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    response["mrca"] = mrcaj;
    return response.dump(1);
}

void to_json(json& j, const attachment_point_t& attachment_point)
{
    j["parent"] = attachment_point.parent;
    j["children_from_taxon"] = attachment_point.children_from_taxon;
}


template <typename Map, typename Key>
auto lookup(const Map& m, const Key& k)
{
    auto it = m.find(k);
    return it != m.end()
               ? std::make_optional(it->second)
               : std::nullopt;
}


const SumTreeNode_t * get_node_for_subtree(const SummaryTree_t * tree_ptr,
                                           const string & node_id, 
                                           int height_limit,
                                           uint32_t tip_limit)
{
    assert(tree_ptr != nullptr);
    auto result = find_required_node_by_id_str(*tree_ptr, node_id);
    if (result.was_broken) {
        json broken;
        broken["mrca"] = get_synth_node_label(result.node);

        auto& contesting_trees_for_taxon = tree_ptr->get_data().contesting_trees_for_taxon;

        if (auto contesting_trees = lookup(contesting_trees_for_taxon, node_id))
        {
            json contesting;
            for(auto& [tree, attachment_points]: *contesting_trees)
            {
                contesting[tree] = {{"attachment_points", attachment_points}};
            }
            broken["contesting_trees"] = contesting;
        }

        json j;
        j["broken"] = broken;
        throw OTCBadRequest("node_id was not found (broken taxon).\n")<<j;
    }
    if (result.node->get_data().num_tips > tip_limit && height_limit < 0) {
        throw OTCBadRequest() << "The requested subtree is too large to be returned via the API. (Tip limit = " << tip_limit << ".) Download the entire tree.\n";
    }
    return result.node;
}


template<typename C, typename T, typename Y>
inline void write_visited_newick_no_semi(std::ostream & out,
                                         const C & visited,
                                         T nd,
                                         Y & nodeNamer) {
    assert(nd != nullptr);
    if (!(nd->is_tip())) {
        bool first = true;
        for (auto c : iter_child_const(*nd)) {
            if (visited.count(c)) {
                if (first) {
                    out << '(';
                    first = false;
                } else {
                    out << ',';
                }
                write_visited_newick_no_semi<C, T, Y>(out, visited, c, nodeNamer);
            }
        }
        if (!first) {
            out << ')';
        }
    }
    write_escaped_for_newick(out, nodeNamer(nd));
}

template<typename C, typename T, typename Y>
inline void write_visited_newick(std::ostream & out,
                                 const C & visited,
                                 T nd,

                                 Y & nodeNamer) {
    write_visited_newick_no_semi<C,T,Y>(out, visited, nd, nodeNamer);
    out << ';';
}

json get_supporting_studies(const set<const string*>& study_id_set) {
    json ss_arr = json::array();
    for (auto study_it_ptr : study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    return ss_arr;
}

string induced_subtree_ws_method(const TreesToServe & tts,
                                 const SummaryTree_t * tree_ptr,
                                 const vector<string> & node_id_vec,
                                 NodeNameStyle label_format) {
    assert(tree_ptr != nullptr);
    const SumTreeNode_t * focal = nullptr;
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    // Check if any of the tip nodes are either (i) broken or (ii) not found.
    auto [tip_nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_id_vec);
    // Find the mrca
    bool first = true;
    for (auto n: tip_nodes) {
        if (first) {
            first = false;
            focal = n;
        } else {
            focal = find_mrca_via_traversal_indices(focal, n);
            if (focal == nullptr) {
                break;
            }
        }
    }
    if (focal == nullptr) {
        throw OTCBadRequest() << "MRCA of taxa was not found.\n";
    }
    // Visit some nodes beween tip_nodes and the mrca
    set<const SumTreeNode_t *> visited;
    visited.insert(focal);
    for (auto tni : tip_nodes) {
        auto cnd = tni;
        while (visited.count(cnd) == 0) {
            visited.insert(cnd);
            cnd = cnd->get_parent(); 
        } 
    }
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy, tts);
    ostringstream out;
    write_visited_newick(out, visited, focal, nnsbs);
    json response;
    response["newick"] = out.str();
    response["supporting_studies"] = get_supporting_studies(nnsbs.study_id_set);
    response["broken"] = broken;
    return response.dump(1);
}

string newick_subtree_ws_method(const TreesToServe & tts,
                                const SummaryTree_t * tree_ptr,
                                const string & node_id,
                                NodeNameStyle label_format, 
                                bool include_all_node_labels,
                                int height_limit)
{
    const uint32_t NEWICK_TIP_LIMIT = 100000;
    auto focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT);
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy, tts);
    ostringstream out;
    write_newick_generic<const SumTreeNode_t *, NodeNamerSupportedByStasher>(out, focal, nnsbs, include_all_node_labels, height_limit);
    json response;
    response["newick"] = out.str();
    response["supporting_studies"] = get_supporting_studies(nnsbs.study_id_set);
    return response.dump(1);
}


template<typename T>
inline void write_arguson(json & j,
                          const TreesToServe & tts,
                          const SummaryTreeAnnotation * sta,
                          const RichTaxonomy & taxonomy,
                          T nd,
                          long height_limit,
                          set<string> & usedSrcIds) {
    assert(nd != nullptr);
    if (!(nd->is_tip()) && height_limit != 0) {
        json c_array;
        const long nhl = height_limit - 1;
        for (auto c : iter_child_const(*nd)) {
            json cj;
            write_arguson<T>(cj, tts, sta, taxonomy, c, nhl, usedSrcIds);
            c_array.push_back(cj);
        }
        j["children"] = c_array;
    }
    add_basic_node_info(taxonomy, *nd, j, true);
    add_node_support_info(tts, *nd, j, usedSrcIds);
}

string arguson_subtree_ws_method(const TreesToServe & tts,
                                 const SummaryTree_t * tree_ptr,
                                 const SummaryTreeAnnotation * sta,
                                 const string & node_id,
                                 int height_limit) {
    const uint32_t NEWICK_TIP_LIMIT = 25000;
    auto focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT);
    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    json a;
    {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        write_arguson(a, tts, sta, taxonomy, focal, height_limit, usedSrcIds);
        add_lineage(tts, a, focal, taxonomy, usedSrcIds, true);
        add_source_id_map(a, usedSrcIds, taxonomy, sta);
    }
    response["arguson"] = a;
    return response.dump(1);
}

template <typename Tree>
void delete_subtree_and_monotypic_ancestors(Tree& tree, typename Tree::node_type* node)
{
    auto parent = node->get_parent();
    tree.prune_and_delete(node);
    if (parent and parent->is_tip()) {
        delete_tip_and_monotypic_ancestors(tree, parent);
    }
}

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



} //namespace otc

