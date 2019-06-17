#include <regex>
#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"
#include "ws/trees_to_serve.h"
#include "otc/conflict.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "nexson/nexson.h"
#include <optional>
#include <string_view>
#include "ws/tnrs/context.h"
INITIALIZE_EASYLOGGINGPP


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

    extern TreesToServe tts;
    typedef RTRichTaxNode Taxon;

    void OTCWebError::prepend(const std::string& s)
    {
        if (not data.count("message"))
        {
                 data["message"] = s;
        }
        else
        {
            string tmp = data["message"];
            data["message"] = s + tmp;
        }
    }

    void OTCWebError::append(const std::string& s)
    {
        if (not data.count("message"))
        {
                 data["message"] = s;
        }
        else
        {
            string tmp = data["message"];
            data["message"] = tmp + s;
        }
    }

    OTCWebError::OTCWebError() noexcept
        : OTCWebError(500,"")
    {}

    OTCWebError::OTCWebError(int c) noexcept
        :OTCWebError(c, "")
    {}

    OTCWebError::OTCWebError(const std::string & msg) noexcept
        :OTCWebError(500,msg)
    { }

    OTCWebError::OTCWebError(int c, const std::string & msg) noexcept
        :status_code_(c)
    {
        prepend(msg);
    }

    template<>
    OTCWebError& OTCWebError::operator<<(const nlohmann::json& j)
    {
        for(auto& [key,value]: j.items())
            data[key] = value;
        return *this;
    }

    template <>
    OTCWebError& OTCWebError::operator<<(const std::string& j)
    {
        append(j);
        return *this;
    }

    inline const string & get_taxon_unique_name(const RTRichTaxNode & nd_taxon) {
        return nd_taxon.get_name();
    }

void add_taxon_info(const RichTaxonomy & , const RTRichTaxNode & nd_taxon, json & taxonrepr) {
    const auto & taxon_data = nd_taxon.get_data();
    taxonrepr["tax_sources"] = taxon_data.get_sources_json();
    taxonrepr["name"] = taxon_data.get_nonuniqname();
    taxonrepr["unique_name"] = get_taxon_unique_name(nd_taxon);
    taxonrepr["rank"] = taxon_data.get_rank();
    taxonrepr["ott_id"] = nd_taxon.get_ott_id();    
}

inline string ott_id_to_idstr(OttId ott_id) {
    string ret;
    ret.reserve(12);
    ret = "ott";
    ret += std::to_string(ott_id);
    return ret;
}

inline string node_id_for_summary_tree_node(const SumTreeNode_t & nd) {
    //return nd.get_name();
    // code below from when we were trying calling clear_names_for_nodes_with_ids to clear
    //    the name strings, but I that saved trivial memory
    return (nd.has_ott_id() ? ott_id_to_idstr(nd.get_ott_id()) : nd.get_name());
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
void add_support_info_vec(const char * tag,
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

void add_support_info_single_element(const char * tag,
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

template<typename N>
N * find_mrca_via_traversal_indices(N *f, N *s) {
    const auto * fdata = &(f->get_data());
    const auto sec_ind = s->get_data().trav_enter;
    while (sec_ind < fdata->trav_enter || sec_ind > fdata->trav_exit) {
        f = f->get_parent();
        if (f == nullptr) {
            assert(false); 
            return nullptr;
        }
        fdata = &(f->get_data());
    }
    return f;
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

const SumTreeNode_t * find_node_by_id_str(const SummaryTree_t & tree,
                                          const string & node_id,
                                          bool & was_broken) {
    was_broken = false;
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
                return i2nit->second;
            }
            LOG(WARNING) << "not finding " << ott_id << " extracted from " << node_id;
        }

        // We didn't find a summary tree node for this OTT ID.  Is this node listed as broken?
        if (auto bt_it = tree_data.broken_taxa.find(node_id); bt_it != tree_data.broken_taxa.end())
        {
            // if found we return the MRCA pointer in the first slot of the pair.
            was_broken = true;
            return bt_it->second.first;
        }
        else
            return nullptr;
    }

    if (std::regex_match(node_id, matches, mrca_id_pattern))
    {
        auto n2nit = tree_data.broken_name_to_node.find(node_id);
        if (n2nit != tree_data.broken_name_to_node.end()) {
            return n2nit->second;
        }
        assert(matches.size() >= 2);
        std::string first_id = matches[1];
        std::string second_id = matches[2];
        bool bogus = false;
        auto fir_nd = find_node_by_id_str(tree, first_id, bogus);
        if (fir_nd == nullptr) {
            return nullptr;
        }
        auto sec_nd = find_node_by_id_str(tree, second_id, bogus);
        if (sec_nd == nullptr) {
            return nullptr;
        }
        return find_mrca_via_traversal_indices(fir_nd, sec_nd);
    }
    return nullptr;
}

const SumTreeNode_t * find_required_node_by_id_str(const SummaryTree_t & tree,
                                                   const string & node_id,
                                                   bool & was_broken) {
    auto node  = find_node_by_id_str(tree, node_id, was_broken);
    if (not node) {
        throw OTCBadRequest() << "node_id '" << node_id << "' was not found!";
    }
    return node;
}

// See API docs at https://github.com/OpenTreeOfLife/germinator/wiki/Synthetic-tree-API-v3

string available_trees_ws_method(const TreesToServe &tts)
{
    json response;
    json trees = json::array();
    for(auto& synth_id: tts.get_available_trees())
        trees.push_back(synth_id);
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

json tax_about_json(const RichTaxonomy & taxonomy)
{
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


inline void add_lineage(json & j, const SumTreeNode_t * focal, const RichTaxonomy & taxonomy, set<string> & usedSrcIds, bool is_arguson = false) {
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
            add_lineage(response, focal, taxonomy, usedSrcIds);
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
                           bool include_lineage)
{
    LOG(DEBUG)<<"Got to node_info_ws_method( )";
    bool was_broken = false;

    auto node = find_required_node_by_id_str(*tree_ptr, node_id, was_broken);

    auto response = node_info_json(tts, sta, node, include_lineage);
    response["query"] = node_id;
    return response.dump(1);
}

pair<vector<const SumTreeNode_t*>,json> find_nodes_for_id_strings(const RichTaxonomy& taxonomy,
                                                                  const SummaryTree_t* tree_ptr,
                                                                  const vector<string>& node_ids,
                                                                  bool fail_broken = false)
{
    vector<const SumTreeNode_t *> nodes;
    json unknown;
    json broken = json::object();

    optional<string> bad_node_id;

    for (auto node_id : node_ids)
    {
        bool was_broken = false;
        const SumTreeNode_t * n = find_node_by_id_str(*tree_ptr, node_id, was_broken);

        if (not n or (was_broken and fail_broken))
        {
            // Possible statuses:
            //  - invalid    (never minted id)
            //  - pruned     (valid but not in synth)
            //  - deprecated (previously valid, not forwarded)
            string reason = "unknown_id";

            if (n and was_broken)
                reason = "broken";
            else if (auto id = is_ott_id(node_id))
            {
                auto taxon = taxonomy.included_taxon_from_id(*id);
//                // Not currently implemented...
//                if (id == -2)
//                    reason = "deprecated";
                if (not taxon)
                    reason = "invalid_ott_id";
                else
                    reason = "pruned_ott_id";
            }

            unknown[node_id] = reason;
            bad_node_id = node_id;
        }

        if (was_broken)
            broken[node_id] = node_id_for_summary_tree_node(*n);

        // Current default strategy means that we include MRCAs for broken taxa.
        nodes.push_back(n);
    }

    if (unknown.size())
        throw OTCBadRequest()<<"node_id '"<< *bad_node_id << "' was not found!"<<json{ {"unknown", unknown} };

    return {nodes, broken};
}

string nodes_info_ws_method(const TreesToServe & tts,
                            const SummaryTree_t * tree_ptr,
                            const SummaryTreeAnnotation * sta,
                            const vector<string> & node_ids,
                            bool include_lineage)
{
    auto [taxonomy,_] = tts.get_readable_taxonomy();

    auto [nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_ids);

    json response;
    for(int i=0;i<nodes.size();i++)
    {
        auto j = node_info_json(tts, sta, nodes[i], include_lineage);
        j["query"] = node_ids[i];
        response.push_back(j);
    }

    return response.dump(1);
}

void add_nearest_taxon(const RichTaxonomy& taxonomy, const SumTreeNode_t& node, json& j)
{
    if (node.has_ott_id()) return;

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
const SumTreeNode_t* mrca(const T& nodes)
{
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    for (auto n : nodes)
    {
        if (first)
        {
            first = false;
            focal = n;
        }
        else
        {
            focal = find_mrca_via_traversal_indices(focal, n);
            if (focal == nullptr) {
                break;
            }
        }
    }
    return focal;
}

string mrca_ws_method(const TreesToServe & tts,
                      const SummaryTree_t * tree_ptr,
                      const SummaryTreeAnnotation * sta,
                      const vector<string> & node_id_vec)
{
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);

    auto [taxonomy,_] = tts.get_readable_taxonomy();

    auto [tip_nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_id_vec);

    auto focal = mrca(tip_nodes);

    if (not focal) throw OTCBadRequest("MRCA of taxa was not found.\n");

    json response;
    response["synth_id"] = sta->synth_id;
    json mrcaj;
    set<string> usedSrcIds;
    add_node_support_info(tts, *focal, mrcaj, usedSrcIds);
    add_basic_node_info(taxonomy, *focal, mrcaj);
    add_nearest_taxon(taxonomy, *focal, response);
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    response["mrca"] = mrcaj;
    return response.dump(1);
}


const SumTreeNode_t * get_node_for_subtree(const SummaryTree_t * tree_ptr,
                                           const string & node_id, 
                                           int height_limit,
                                           uint32_t tip_limit) {
    assert(tree_ptr != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_required_node_by_id_str(*tree_ptr, node_id, was_broken);
    if (was_broken) {
        throw OTCBadRequest("node_id was not found (broken taxon).\n");
    }
    if (focal->get_data().num_tips > tip_limit && height_limit < 0) {
        throw OTCBadRequest() << "The requested subtree is too large to be returned via the API. (Tip limit = " << tip_limit << ".) Download the entire tree.\n";
    }
    return focal;
}

class NodeNamerSupportedByStasher {
    public:
        mutable set<const string *> study_id_set;
        NodeNameStyle nns;
        const RichTaxonomy & taxonomy;
        NodeNamerSupportedByStasher(NodeNameStyle in_nns, const RichTaxonomy &tax)
            :nns(in_nns),
            taxonomy(tax) {
        }

        string operator()(const SumTreeNode_t *nd) const {
            const SumTreeNodeData & d = nd->get_data();
#           if defined(JOINT_MAPPING_VEC)
                for (auto & el : d.source_edge_mappings) {
                    if (el.first == SourceEdgeMappingType::SUPPORTED_BY_MAPPING) {
                        study_id_set.insert(tts.decode_study_node_id_index(el.second).first);
                    }
                }
#           else
                for (auto & sni : d.supported_by) {
                    study_id_set.insert(tts.decode_study_node_id_index(sni).first);
                }
#           endif
            if (nns != NodeNameStyle::NNS_ID_ONLY && nd->has_ott_id()) {
                const auto * tr = taxonomy.included_taxon_from_id(nd->get_ott_id());
                if (tr == nullptr) {
                    throw OTCError() << "OTT Id " << nd->get_ott_id() << " in namer not found in taxonomy! Please report this bug";
                }
                string taxon_name = get_taxon_unique_name(*tr);
                if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
                    string ret;
                    const string id_str = node_id_for_summary_tree_node(*nd);
                    ret.reserve(taxon_name.length() + 1 + id_str.length());
                    ret = taxon_name;
                    ret += ' ';
                    ret += id_str;
                    return ret;    
                } else {
                    return taxon_name;
                }
            }
            return node_id_for_summary_tree_node(*nd);
        }
        string operator()(const RTRichTaxNode *nd) const {
            assert(nd != nullptr);
            if (nns == NodeNameStyle::NNS_ID_ONLY) {
                return ott_id_to_idstr(nd->get_ott_id());
            }
            const string & taxon_name = get_taxon_unique_name(*nd);
            if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
                string ret;
                string id_str = ott_id_to_idstr(nd->get_ott_id());
                ret.reserve(taxon_name.length() + 1 + id_str.length());
                ret = taxon_name;
                ret += ' ';
                ret += id_str;
                return ret;    
            }
            return taxon_name;
        }
};

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
    out<<';';
}

json get_supporting_studies(const set<const string*>& study_id_set)
{
    json ss_arr = json::array();
    for (auto study_it_ptr : study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    return ss_arr;
}

string induced_subtree_ws_method(const TreesToServe & tts,
                                 const SummaryTree_t * tree_ptr,
                                 const vector<string> & node_id_vec,
                                 NodeNameStyle label_format)
{
    assert(tree_ptr != nullptr);
    const SumTreeNode_t * focal = nullptr;

    auto [taxonomy,_] = tts.get_readable_taxonomy();

    // Check if any of the tip nodes are either (i) broken or (ii) not found.
    auto [tip_nodes, broken] = find_nodes_for_id_strings(taxonomy, tree_ptr, node_id_vec);

    // Find the mrca
    bool first = true;
    for (auto n: tip_nodes)
    {
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

    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
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
                                int height_limit) {
    const uint32_t NEWICK_TIP_LIMIT = 25000;
    const SumTreeNode_t * focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT);

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
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
        add_lineage(a, focal, taxonomy, usedSrcIds, true);
        add_source_id_map(a, usedSrcIds, taxonomy, sta);
    }
    response["arguson"] = a;
    return response.dump(1);
}

// taxon_info
void tax_service_add_taxon_info(const RichTaxonomy & taxonomy,
                                const RTRichTaxNode & nd_taxon,
                                json & taxonrepr) {
    add_taxon_info(taxonomy, nd_taxon, taxonrepr);
    taxonrepr["source"] = string("ott") + taxonomy.get_version(); //TBD "source" ?
    const auto & taxon_data = nd_taxon.get_data();
    taxonrepr["flags"] = flags_to_string_vec(taxon_data.get_flags());
    json syn_list = json::array();
    for (auto tjs : taxon_data.junior_synonyms) {
        syn_list.push_back(tjs->get_name());
    }
    taxonrepr["synonyms"] = syn_list;
    const auto & ots = taxonomy.get_ids_to_suppress_from_tnrs();
    const auto ott_id = nd_taxon.get_ott_id();
    const bool is_suppressed = (0 < ots.count(ott_id));
    taxonrepr["is_suppressed"] = is_suppressed;
    auto isfs = taxonomy.get_ids_suppressed_from_summary_tree_alias();
    if (isfs) {
        taxonrepr["is_suppressed_from_synth"] = is_suppressed || (0 < isfs->count(ott_id));
    }
}

string taxon_info_ws_method(const RichTaxonomy & taxonomy,
                            const RTRichTaxNode * taxon_node,
                            bool include_lineage,
                            bool include_children,
                            bool include_terminal_descendants) {
    assert(taxon_node != nullptr);
    json response;
    tax_service_add_taxon_info(taxonomy, *taxon_node, response);
    if (include_lineage) {
        json anc_array = json::array();
        for (auto a : iter_anc_const(*taxon_node)) {
            json a_el;
            tax_service_add_taxon_info(taxonomy, *a, a_el);
            anc_array.push_back(a_el);
        }
        response["lineage"] = anc_array;
    }
    if (include_children) {
        json c_array = json::array();
        for (auto c : iter_child_const(*taxon_node)) {
            json c_el;
            tax_service_add_taxon_info(taxonomy, *c, c_el);
            c_array.push_back(c_el);
        }
        response["children"] = c_array;
    }
    if (include_terminal_descendants) {
        json td_array = json::array();
        for (auto nd : iter_leaf_n_const(*taxon_node)) {
            td_array.push_back(nd->get_ott_id());
        }
        response["terminal_descendants"] = td_array;
    }
    return response.dump(1);
}

// flags: not_otu, environmental, environmental_inherited, viral, hidden, hidden_inherited, was_container
//    are excluded from being returned in TNRS results.

std::string taxonomy_flags_ws_method(const RichTaxonomy & taxonomy)
{
    vector<int> flag_counts(32);
    for(auto  node: iter_node_const(taxonomy.get_tax_tree()))
    {
        for(int i=0;i<32;i++)
            if (node->get_data().flags.test(i))
                flag_counts[i]++;
    }
    json flags;

//    `curl -X POST https://api.opentreeoflife.org/v3/taxonomy/flags` yields something like:

    vector<string> required_flags = {  "environmental_inherited",
                                       "unclassified", //0
                                       "unplaced",
                                       "hidden",
                                       "extinct_inherited",
                                       "unclassified_direct", //0
                                       "extinct",
                                       "unclassified_inherited", //0
                                       "incertae_sedis",
                                       "was_container",
                                       "forced_visible", //0
                                       "hidden_inherited",
                                       "not_otu",
                                       "environmental",
                                       "major_rank_conflict",
                                       "edited", //0
                                       "incertae_sedis_direct", //0
                                       "extinct_direct", //0
                                       "unplaced_inherited",
                                       "merged",
                                       "major_rank_conflict_inherited",
                                       "major_rank_conflict_direct", //0
                                       "inconsistent",
                                       "barren",
                                       "sibling_higher",
                                       "tattered", //0
                                       "infraspecific",
                                       "hybrid",
                                       "incertae_sedis_inherited",
                                       "viral",
                                       "tattered_inherited", //0
                                       "sibling_lower", // 0
    };

    for(auto& flag_name: required_flags)
    {
        flags[flag_name] = 0;
    }

    // Add the actual counts.
    for(int i=0;i<32;i++)
        if (flag_counts[i])
            flags[*string_for_flag(i)] = flag_counts[i];

    return flags.dump(1);
}

// BDR: factored this code out of taxonomy_mrca_ws_method below for use in tnrs
const RTRichTaxNode* taxonomy_mrca(const vector<const RTRichTaxNode*>& nodes)
{
    if (nodes.empty()) return nullptr;

    auto focal = nodes[0];
    for(auto& node: nodes)
    {
	focal = find_mrca_via_traversal_indices(focal, node);
	if (not focal) 
	    throw OTCWebError(400, "MRCA of taxa was not found. Please report this bug!\n");
    }

    return focal;
}

string taxonomy_mrca_ws_method(const RichTaxonomy & taxonomy,
                               const OttIdSet & ott_id_set) {
    const RTRichTaxNode * focal = nullptr;
    bool first = true;
    for (auto ott_id : ott_id_set) {
        const RTRichTaxNode * n = taxonomy.included_taxon_from_id(ott_id);
        if (n == nullptr) {
            throw OTCBadRequest() << "ott_id \"" << ott_id <<  "\" was not recognized.\n";
        }
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
        throw OTCWebError(400, "MRCA of taxa was not found. Please report this bug!\n");
    }
    json mrcaj;
    tax_service_add_taxon_info(taxonomy, *focal, mrcaj);
    json response;
    response["mrca"] = mrcaj;
    return response.dump(1);
}

string taxon_subtree_ws_method(const RichTaxonomy & taxonomy,
                             const RTRichTaxNode * taxon_node,
                             NodeNameStyle label_format) {
    assert(taxon_node != nullptr);
    const auto & taxonomy_tree = taxonomy.get_tax_tree();
    json response;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    int height_limit = -1;
    bool include_all_node_labels = true; // ??
    write_newick_generic<const RTRichTaxNode *, NodeNamerSupportedByStasher>(out, taxon_node, nnsbs, include_all_node_labels, height_limit);
    response["newick"] = out.str();
    return response.dump(1);
}

// It seems that infer_context is supposed to ignore synonyms:
// * Bacteria is NOT ambiguous: it has 1 name match, and 1 synonym match.
// * Firmiscala IS ambiguous: it has 0 name matches, and 1 synonym match.
// * Random gibberish is reported as "ambiguous".
// curl -X POST https://api.opentreeoflife.org/v2/tnrs/infer_context -H "content-type:application/json" -d  '{"names":["Bacteria","Firmiscala"]}'
// 

template <typename T>
bool lcase_string_equals(const string_view& s1, const T& s2)
{
    if (s1.size() != s2.size()) return false;

    for(int i=0;i<s1.size();i++)
	if (std::tolower(s1[i]) != std::tolower(s2[i]))
	    return false;

    return true;
}

bool lcase_match_prefix(const string_view& s, const string_view& prefix)
{
    if (prefix.size() > s.size()) return false;

    for(int i=0;i<prefix.size();i++)
	if (std::tolower(s[i]) != std::tolower(prefix[i]))
	    return false;

    return true;
}

bool taxon_is_specific(const Taxon* taxon)
{
    auto rank = taxon->get_data().rank;

    return rank_is_specific(rank);
}

bool taxon_is_genus(const Taxon* taxon)
{
    return taxon->get_data().rank == TaxonomicRank::RANK_GENUS;
}

bool taxon_is_higher(const Taxon* taxon)
{
    return taxon->get_data().rank < TaxonomicRank::RANK_SPECIES;
}

vector<pair<const Taxon*,const string&>> exact_synonym_search(const Taxon* context_root, string query, std::function<bool(const Taxon*)> ok = [](const Taxon*){return true;})
{
    vector<pair<const Taxon*,const string&>> hits;
    for(auto taxon: iter_post_n_const(*context_root))
    {
	if (not ok(taxon)) continue;

	for(auto& tjs: taxon->get_data().junior_synonyms)
	    if (lcase_string_equals(query, tjs->get_name()))
		hits.push_back({taxon,tjs->get_name()});
    }

    return hits;
}

vector<pair<const Taxon*,const string&>> exact_synonym_search(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					       {
						   if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
						   return true;
					       };

    return exact_synonym_search(context_root, query, ok);
}

vector<pair<const Taxon*,const string&>> exact_synonym_search_higher(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					       {
						   if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
						   if (not taxon_is_higher(taxon)) return false;
						   return true;
					       };

    return exact_synonym_search(context_root, query, ok);
}



vector<const Taxon*> exact_name_search(const Taxon* context_root, string query, std::function<bool(const Taxon*)> ok = [](const Taxon*){return true;})
{
    for(auto& c: query)
	c = std::tolower(c);

    vector<const Taxon*> hits;
    for(auto taxon: iter_post_n_const(*context_root))
    {
	if (not ok(taxon)) continue;

	if (lcase_string_equals(query, taxon->get_data().get_nonuniqname()))
	    hits.push_back(taxon);
    }

    return hits;
}

vector<const Taxon*> exact_name_search(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    return true;
					};

    return exact_name_search(context_root, query, ok);
}

vector<const Taxon*> exact_name_search_species(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    if (not taxon_is_specific(taxon)) return false;
					    return true;
					};

    return exact_name_search(context_root, query, ok);
    
}

vector<const Taxon*> exact_name_search_genus(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    if (not taxon_is_genus(taxon)) return false;
					    return true;
					};

    return exact_name_search(context_root, query, ok);
    
}

vector<const Taxon*> exact_name_search_higher(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    if (not taxon_is_higher(taxon)) return false;
					    return true;
					};

    return exact_name_search(context_root, query, ok);
}

vector<const Taxon*> prefix_name_search(const Taxon* context_root, const string& query, std::function<bool(const Taxon*)> ok = [](const Taxon*){return true;})
{
    vector<const Taxon*> hits;
    for(auto taxon: iter_post_n_const(*context_root))
    {
	if (not ok(taxon)) continue;

	if (lcase_match_prefix(taxon->get_data().get_nonuniqname(), query))
	    hits.push_back(taxon);
    }

    return hits;
}

vector<const Taxon*> prefix_name_search(const RichTaxonomy& taxonomy, const Taxon* context_root, const string& query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    return true;
					};

    return prefix_name_search(context_root, query, ok);
}

vector<const Taxon*> prefix_name_search_higher(const RichTaxonomy& taxonomy, const Taxon* context_root, const string& query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					{
					    if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
					    if (not taxon_is_higher(taxon)) return false;
					    return true;
					};

    return prefix_name_search(context_root, query, ok);
}

vector<pair<const Taxon*,const string&>> prefix_synonym_search(const Taxon* context_root, string query, std::function<bool(const Taxon*)> ok = [](const Taxon*){return true;})
{
    vector<pair<const Taxon*,const string&>> hits;
    for(auto taxon: iter_post_n_const(*context_root))
    {
	if (not ok(taxon)) continue;

	for(auto& tjs: taxon->get_data().junior_synonyms)
	    if (lcase_match_prefix(tjs->get_name(), query))
		hits.push_back({taxon,tjs->get_name()});
    }

    return hits;
}

vector<pair<const Taxon*,const string&>> prefix_synonym_search(const RichTaxonomy& taxonomy, const Taxon* context_root, string query, bool include_suppressed)
{
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon)
					       {
						   if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) return false;
						   return true;
					       };

    return prefix_synonym_search(context_root, query, ok);
}

vector<const Taxon*> exact_name_search(const RichTaxonomy& taxonomy, const string& query, bool include_suppressed)
{
    auto context_root = taxonomy.get_tax_tree().get_root();

    return exact_name_search(taxonomy, context_root, query, include_suppressed);
}

enum match_status {unmatched=0, ambiguous_match=1, unambiguous_match=2};

// Should probably rename ContextSearcher => Context, and Context=> ContextDescription
struct ContextSearcher
{
    const RichTaxonomy& taxonomy;
    const Context& context;
    const Taxon* context_root;

    pair<json,match_status> match_name(string query, bool do_approximate_matching, bool include_suppressed);

    ContextSearcher(const RichTaxonomy& t, const Context& c): taxonomy(t), context(c)
    {
	context_root = taxonomy.included_taxon_from_id(c.ott_id);
    }
};

string escape_query_string(const string& name)
{
    return name;
}

json get_taxon_json(const RichTaxonomy& taxonomy, const Taxon& taxon)
{
    json j;
    tax_service_add_taxon_info(taxonomy, taxon, j);
    // What about the "is_suppressed_from_synth" flag?  Do we want that?
    return j;
}

json exact_name_match_json(const string& query, const Taxon* taxon, const RichTaxonomy& taxonomy)
{
    json result;
    result["taxon"] = get_taxon_json(taxonomy, *taxon);
    result["search_string"] = query;
    result["score"] = 1.0;
    result["is_approximate_match"] = false;
    result["nomenclature_code"] = "code";   // FIXME!
    result["is_synonym"] = false;
    result["matched_name"] = taxon->get_data().get_nonuniqname();
    return result;
}

json exact_synonym_match_json(const string& query, const Taxon* taxon, const string& synonym_name, const RichTaxonomy& taxonomy)
{
    json result;
    result["taxon"] = get_taxon_json(taxonomy, *taxon);
    result["search_string"] = query;
    result["score"] = 1.0;
    result["is_approximate_match"] = false;
    result["nomenclature_code"] = "code";   // FIXME!
    result["is_synonym"] = true;
    result["matched_name"] = synonym_name;
    return result;
}

pair<json,match_status> ContextSearcher::match_name(string query, bool do_approximate_matching, bool include_suppressed)
{
    for(auto& c: query)
	c = std::tolower(c);

    // What does an ambiguous match mean?
    json results;

    match_status status = unmatched;

    // 1. See if we can find an exact name match
    auto exact_name_matches = exact_name_search(taxonomy, context_root, query, include_suppressed);

    for(auto taxon: exact_name_matches)
	results.push_back(exact_name_match_json(query, taxon, taxonomy));

    if (exact_name_matches.size() == 1)
	status = unambiguous_match;

    // 2. See if we can find an exact name match for synonyms
    auto exact_synonym_matches = exact_synonym_search(taxonomy, context_root, query, include_suppressed);

    for(auto& [ taxon, synonym_name ]: exact_synonym_matches)
	results.push_back(exact_synonym_match_json(query, taxon, synonym_name, taxonomy));

    if (status == unmatched and results.size())
	status = ambiguous_match;

    // 3. Do fuzzy matching ONLY for names that we couldn't match
    if (do_approximate_matching and status == unmatched)
    {
	// do fuzzy matching.
    }

    json match_results;
    match_results["name"] = query;
    match_results["matches"] = results;

    return {match_results, status};
}

const Context* determine_context_for_names(const vector<string>& names, const optional<string>& context_name, const RichTaxonomy& taxonomy)
{
    // 1. Determine context
    const Context* context;
    if (not context_name)
	context = infer_context_and_ambiguous_names(taxonomy, names).first;
    else
    {
	if (not name_to_context.count(*context_name))
	    throw OTCError()<<"The context '"<<*context_name<<"' could not be found.";
	context = name_to_context.at(*context_name);
    }
    return context;
}

const Context* determine_context(const optional<string>& context_name)
{
    // 1. Determine context
    const Context* context;
    if (not context_name)
	context = name_to_context.at("All life");
    else
    {
	if (not name_to_context.count(*context_name))
	    throw OTCError()<<"The context '"<<*context_name<<"' could not be found.";
	context = name_to_context.at(*context_name);
    }
    return context;
}

//FIXME: how is "suppressed_names" different from "deprecated_taxa"?

// $ curl -X POST https://api.opentreeoflife.org/v3/tnrs/match_names  -H "content-type:application/json" -d '{"names":["Aster","Symphyotrichum","Barnadesia"]}'
// $ curl -X POST http://localhost:1984/v3/tnrs/match_names  -H "content-type:application/json" -d '{"names":["Aster","Symphyotrichum","Barnadesia"]}'
std::string tnrs_match_names_ws_method(const vector<string>& names,
                                       const optional<string>& context_name,
                                       bool do_approximate_matching,
                                       const optional<vector<string>>& /* ids */, // FIXME: Do we need to implement this?
                                       bool include_suppressed,
                                       const RichTaxonomy& taxonomy)
{

    // This corresponds to a MultiNameContextQuery in taxomachine.
    // * See org/opentree/taxonomy/plugins/tnrs_v3.java
    // * See org/opentree/tnrs/queries/MultiNameContextQuery.java:
        
    // 1. Determine context
    auto context = determine_context_for_names(names, context_name, taxonomy);

    ContextSearcher searcher(taxonomy, *context);

    // 2. Iterate over names and fill arrays `results`, `unmatched_names`, `matched_names`, and `unambiguous_names`.
    json results = json::array();
    json unambiguous_names = json::array();
    json unmatched_names = json::array();
    json matched_names = json::array();

    for(auto& name: names)
    {
	// Do the search
	auto [result,status] = searcher.match_name(name, do_approximate_matching, include_suppressed);

	// Store the result
	results.push_back(result);

	// Classify name as unmatched / matched / unambiguous
	if (status == unmatched)
	    unmatched_names.push_back(name);
	else
	{
	    matched_names.push_back(name);

	    if (status == unambiguous_match)
		unambiguous_names.push_back(name);
	}
    }

    // 3. Construct JSON response.
    json response;
    response["governing_code"] = context->code.name;
    response["context"] = context->name;
    response["includes_approximage_matches"] = do_approximate_matching;
    response["includes_deprecated_taxa"] = false; // ?? How is this different from suppressed_names?
    response["includes_suppressed_names"] = include_suppressed;
    response["taxonomy"] = tax_about_json(taxonomy);
    response["unambiguous_names"] = unambiguous_names;
    response["unmatched_names"] = unmatched_names;
    response["matched_names"] = matched_names;
    response["results"] = results;

    return response.dump(1);
}


// return name.split("\s+",2) if the string has a space or an optional with no value otherwise.
optional<pair<string,string>> split_genus_species(const string& name)
{
    auto first_space = name.find(' ');

    // Quit if there's no space
    if (first_space == string::npos) return {};

    auto genus = name.substr(0,first_space);

    auto non_space = name.find_first_not_of(' ', first_space+1);
    if (non_space == string::npos)
	non_space = name.size();

    auto species = name.substr(non_space, name.size() - non_space);

    return {{genus, species}};
}

/*
Note from taxomachine: tnrs_v3.java:

Assumes the input is a taxon name that may be incomplete (i.e. the beginning of a taxon name such as 'Ast', 
which would match 'Astelia', 'Astilbe', 'Aster', 'Asteroidea', 'Asteraceae', 'Astrantia', etc.). If the input 
string is an exact string match to an existing taxon name, then only the exact match will be returned, (i.e. the
input 'Aster' will produce a single result 'Aster').

If name expansion identifies a valid genus name, the results will
not include species names from within that genus, but if a trailing space exists in the input following a valid
genus name, then species names will be returned. For example, both 'Garcin' and 'Garcinia' will match the genus name
'Garcinia' itself but will not match any species names within the genus, but 'Garcinia ' (note the trailing space)
will match all the species in the genus, and 'Garcinia m' with match all species names in Garcinia with a specific
epithet that starts with 'm'.

**IMPORTANT NOTE: This service should not be used for general purpose TNRS queries.** It is optimized for and
(obviously) intended for use *only* with autocomplete boxes on web forms. For all name matching purposes other than
autocompleting name fields on forms, use the `match_names` service.
*/

json autocomplete_json(const RichTaxonomy& taxonomy, const Taxon* taxon)
{
    json match;
    match["ott_id"] = taxon->get_ott_id();
    match["unique_name"] = get_taxon_unique_name(*taxon);
    match["is_suppressed"] = taxonomy.node_is_suppressed_from_tnrs(taxon);
    match["is_higher"] = taxon_is_higher(taxon);
    return match;
}    

json autocomplete_json(const RichTaxonomy& taxonomy, const pair<const Taxon*,const string&>& p)
{
    auto [taxon,syonym] = p;
    return autocomplete_json(taxonomy,taxon);
}

void add_hits(json& j, const RichTaxonomy& taxonomy, const vector<const Taxon*> taxa)
{
    for(auto taxon:taxa)
	j.push_back(autocomplete_json(taxonomy,taxon));
}
		    
void add_hits(json& j, const RichTaxonomy& taxonomy, const vector<pair<const Taxon*,const string&>> taxa)
{
    for(auto [taxon,synonym]:taxa)
	j.push_back(autocomplete_json(taxonomy,taxon));
}

// Find all species in the genus that have the given prefix
vector<const Taxon*> prefix_search_species_in_genus(const Taxon* genus, const string_view& species_prefix)
{
    vector<const Taxon*> match_species;

    auto genus_name = genus->get_data().get_nonuniqname();
    for(auto species: iter_post_n_const(*genus))
    {
	if (not taxon_is_specific(species)) continue;

	auto species_name = species->get_data().get_nonuniqname().substr(genus_name.size()+1);

	if (lcase_match_prefix(species_name, species_prefix))
	    match_species.push_back(species);
    }

    return match_species;
}

// curl -X POST https://api.opentreeoflife.org/v3/tnrs/autocomplete_name -H "content-type:application/json" -d '{"name":"Endoxyla","context_name":"All life"}'
string tnrs_autocomplete_name_ws_method(const string& name, const string& context_name, bool include_suppressed, const RichTaxonomy& taxonomy)
{
    json response = json::array();

    // We need to escape the query string.
    auto escaped_query = escape_query_string(name);

    // This corresponds to a SingleNamePrefixQuery in taxomachine.
    // * See org/opentree/taxonomy/plugins/tnrs_v3.java
    // * See org/opentree/tnrs/queries/SingleNamePrefixQuery.java


    // 0. Escape the query??
    // lower-case the name?

    LOG(DEBUG)<<"tnrs_autocomplete_name_ws_method: name = '"<<name<<"'  context = '"<<context_name<<"' include_suppressed = "<<include_suppressed;
    // 1. Determine context
    auto context = determine_context(context_name);
    auto context_root = taxonomy.included_taxon_from_id(context->ott_id);

    // 2. If we have a space, then assume the first part is a genus and match species names within the genus
    if (auto query_genus_species = split_genus_species(name))
    {
	// Search against species and synonyms
	add_hits(response, taxonomy, exact_name_search_species(taxonomy, context_root, escaped_query, include_suppressed));
	add_hits(response, taxonomy, exact_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);
	
	// no exact hit against the species index
	auto genus_hits = exact_name_search_genus(taxonomy, context_root, escaped_query, include_suppressed);

	if (not genus_hits.empty()) // the first word was an exact match against the genus index
	{
	    auto [query_genus,query_species] = *query_genus_species;

	    for(auto genus: genus_hits)
		add_hits(response, taxonomy, prefix_search_species_in_genus(genus, query_species));
	}
	if (not response.empty()) return response.dump(1);

	// no exact hit for first word against the genus index

        // Hit query string against the higher taxon index... not sure if this is useful, since it has a space
	add_hits(response, taxonomy, exact_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);

	// Prefix query against the synonyms and higher taxa
	add_hits(response, taxonomy, prefix_name_search(taxonomy, context_root, escaped_query, include_suppressed));
	add_hits(response, taxonomy, prefix_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);
	
	// fuzzy search on names and synonyms
    }
    else // does not contain a space at all
    {
	add_hits(response, taxonomy, exact_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
	add_hits(response, taxonomy, exact_synonym_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);

	// Do a prefix query against the higher taxon index
	add_hits(response, taxonomy, prefix_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);

	// Do a prefix query against the all taxa synonym index
	add_hits(response, taxonomy, prefix_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
	if (not response.empty()) return response.dump(1);
	
	// fuzzy search on higher names and synonyms
    }
    
    return response.dump(1);
}

// curl -X POST https://api.opentreeoflife.org/v3/tnrs/contexts
// curl -X POST http://localhost:1984/v3/tnrs/contexts
std::string tnrs_contexts_ws_method()
{
    json response;
    for(auto& context: all_contexts)
	response[context.group].push_back(context.name);
    return response.dump(1);
}

// curl -X POST https://api.opentreeoflife.org/v3/tnrs/infer_context -H "content-type:application/json" -d '{"names":["Pan","Homo","Mus","Bufo","Drosophila"]}'
// curl -X POST http://localhost:1984/v3/tnrs/infer_context -H "content-type:application/json" -d '{"names":["Pan","Homo","Mus","Bufo","Drosophila"]}'
string tnrs_infer_context_ws_method(const vector<string>& names, const RichTaxonomy& taxonomy)
{
    auto results = infer_context_and_ambiguous_names(taxonomy, names);
    auto& context = results.first;
    auto& ambiguous_names = results.second;
    
    json response;

    response["context_name"] = results.first->name;
    response["context_ott_id"] = results.first->ott_id;
    response["ambiguous_names"] = ambiguous_names;
    
    return response.dump(1);
}


using cnode_type = ConflictTree::node_type;

// input tree node names should all be of the form node###
// ott        node names should all be of the form ott###
// synth tree node names should all be of the form ott###  or  mrcaott###ott###

template <typename N>
string node_name(const N* node)
{
    assert(not node->get_name().empty());
    return node->get_name();
}

int min_depth(const cnode_type* node)
{
    while (node and node->get_parent() and node->get_parent()->is_outdegree_one_node())
        node = node->get_parent();
    return node->get_data().depth;
}

struct conflict_stats
{
    protected:
    template <typename T>
    void throw_otcerror_if_second_not_true(const T & result, const cnode_type* node1) {
        if (not result.second) {
            throw OTCError() << "key "<<node_name(node1) << " occurs twice in a conflict_stats map!";
        }    
    }
    public:
    map<string, string> supported_by;
    map<string, string> partial_path_of;
    map<string, std::set<pair<string, int>>> conflicts_with;
    map<string, string> resolved_by;
    map<string, string> resolves;
    map<string, string> terminal;

    void add_supported_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = supported_by.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `supported_by` node1 for only one node2.
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_partial_path_of(const cnode_type* node2, const cnode_type* node1) {
        /* auto result = */ partial_path_of.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `partial_path_of` node1 for multiple node2s.
        // throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolved_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = resolved_by.insert({node_name(node1), node_name(node2)});
        // We can have the relationship node2 `resolved_by` node1 for only one node2.
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolves(const cnode_type* node2, const cnode_type* node1) {
        /* auto result = */ resolves.insert({node_name(node1), node_name(node2)});
        // We can only have the relationship (node2 `resolves` node1) one multiple node2s.
        // throw_otcerror_if_second_not_true(result, node1);
    }

    void add_terminal(const cnode_type* node2, const cnode_type* node1) {
        /* auto result = */ terminal.insert({node_name(node1), node_name(node2)});
//      We can have the relationship (node2 `terminal` node1) for multiple node2s!
//      Let's keep the FIRST one, since it will be the most tipward-one.
//        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_conflicts_with(const cnode_type* node2, const cnode_type* node1)
    {
        auto name1 = node_name(node1);
        auto name2 = node_name(node2);

        // 1. Get the (possibly empty) set of deepest nodes in tree2 conflicting with node1 in tree1.
        set<pair<string, int>>& nodes = conflicts_with[name1];

        int min_depth_node2 = min_depth(node2);

        // 2. If the newest node (node2) is deeper, then clear the non-minimal-depth nodes.
        if (not nodes.empty() and (min_depth_node2 < nodes.begin()->second)) nodes.clear();

        // 3. If the newest node (node2) is shallower then do nothing
        if (not nodes.empty() and (min_depth_node2 > nodes.begin()->second)) return;

        // 4. If the newest node is not shallower then add it.
        nodes.insert(pair<string, int>({name2, min_depth_node2}));
    }
    json get_json(const ConflictTree&, const RichTaxonomy&) const;
};

using tnode_type = RTRichTaxNode;

int depth(const tnode_type* node) {
    return node->get_data().depth;
}

string extract_node_name_if_present(const string& newick_name)
{
    auto node_name = get_source_node_name(newick_name);
    if (node_name)
        return *node_name;
    else
        return newick_name;
}

json get_node_status(const string& witness, string status, const RichTaxonomy& Tax) {
    json j;
    j["witness"] = extract_node_name_if_present(witness);
    j["status"] = std::move(status);
    long raw_ott_id = long_ott_id_from_name(witness);
    if (raw_ott_id >= 0) {
        OttId id = check_ott_id_size(raw_ott_id);
        auto nd = Tax.included_taxon_from_id(id);
        if (nd != nullptr) {
            j["witness_name"] = nd->get_name();
        }
    }
    return j;
}

json get_conflict_node_status(const set<pair<string,int>>& witnesses, string status, const RichTaxonomy& Tax) {
    json j;
    string first_witness = witnesses.begin()->first;
    j["witness"] = extract_node_name_if_present(witnesses.begin()->first);
    j["status"] = std::move(status);

    // Compute first 3 witnesses as name + name + name + ...
    string witness_names;
    int total = 0;
    for(auto& witness: witnesses)
    {
        long raw_ott_id = long_ott_id_from_name(witness.first);

        if (raw_ott_id < 0) continue;

        OttId id = check_ott_id_size(raw_ott_id);

        auto nd = Tax.included_taxon_from_id(id);

        if (nd == nullptr) continue;

        string witness_name  = nd->get_name();
        
        if (total > 2)
        {
            witness_names += " + ...";
            break;
        }
        else if (total == 0)
        {
            witness_names = witness_name;
        }
        else
        {
            witness_names = witness_names + " + " + witness_name;
        }
        total++;
    }

    if (not witness_names.empty())
    {
        //      witness_names  = witness_names + " (" + std::to_string(witnesses.size()) + ")"; // record number of witnesses
        j["witness_name"] = witness_names;
    }

    return j;
}

// FIXME: Conflict relations currently refer to nodes by a name string.
//        We assume that this name string:
//         1. Is accessed as node->get_name()
//         2. Contains the ottid, if there is one.
//        The curator application also assumes that
//         3. JSON results reference the nexml node name (e.g. nodeYYY) for phylesystem trees.
//
//        Actual node->get_name() strings look like this:
//         ott:   ottXXX
//         synth: ottXXX or mrcaottWWWottZZZ
//         (internally) phylesystem: nodeYYY
//         newick: stuff_ottXXX for leaves, and stuff otherwise.
//         newick from phylesystem (externally): nodeYYY (for internal nodes) or nodeYYY_ottXXX (for leaves).
//
//        I assume that phylesystem node names are transformed from name -> name_ottXXX for leaf nodes, when translating phylesystem trees to newick.
//
//        Assumptions 1-3 are met for ott, synth, and (internal) phylesystem trees.
//
//        For newick trees from phylesystem, we handle this situation by translating these names from
//        stuff_nodeYYY_ottXXX -> nodeYYY before constructing the json, as below. (10/03/2017 - BDR).
//
//        *ALTERNATIVELY*, we could chop any _ottXXX suffix from the leaf names to get the name used in the JSON result, but this could affect non-phylesystem inputs.
//


// PROBLEM: It is possible to have both y1 conflicts with x (x:conflicts_with y1) and y2 resolves x (x:resolves y2)
// PROBLEM: It is possible to have both y1 supported_by x (x:supported_by y1) and y2 resolves x (x:resolves y2)
// PROBLEM: It is possible to have both x resolves y (x:resolved_by y1) and y2 resolves x (x:resolves y2)
// Let's solve this situation by NOT reporting when tree2 resolves tree1.

json conflict_stats::get_json(const ConflictTree& tree, const RichTaxonomy& Tax) const {
    json nodes;
//    for(auto& x: resolves) {
//        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "resolves", Tax);
//    }
    for(auto& x: resolved_by) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "resolved_by", Tax);
    }
    for(auto& x: supported_by) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "supported_by", Tax);
    }
    for(auto& x: partial_path_of) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "partial_path_of", Tax);
    }
    for(auto& x: terminal) {
        nodes[extract_node_name_if_present(x.first)] = get_node_status(x.second, "terminal", Tax);
    }
    for(auto& x: conflicts_with) {
        nodes[extract_node_name_if_present(x.first)] = get_conflict_node_status(x.second, "conflicts_with", Tax);
    }
    // For monotypic nodes in the query, copy annotation from child.
    for(auto it: iter_post_const(tree))
        if (it->is_outdegree_one_node())
        {
            auto name = extract_node_name_if_present(it->get_name());
            auto child_name = extract_node_name_if_present(it->get_first_child()->get_name());
            nodes[name] = nodes.at(child_name);
        }
    return nodes;
}

// Get a list of nodes in T2 that are leaves in T1.
// The nodes in T2 do NOT need to be leaves of T2.

template <typename Tree1_t, typename Tree2_t>
auto get_induced_nodes(const Tree1_t& T1, const Tree2_t& T2)
{
    auto& ott_to_nodes2 = T2.get_data().id_to_node;

    std::vector<const typename Tree2_t::node_type*> nodes;
    for(auto leaf: iter_leaf_const(T1))
    {
        auto id = leaf->get_ott_id();
        auto it = ott_to_nodes2.find(id);

        if (it != ott_to_nodes2.end()) nodes.push_back(it->second);
    }
    return nodes;
}

template <typename Tree_t>
std::size_t n_leaves(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    std::size_t count = 0;
    for(auto nd: iter_leaf_const(T)){
        count++;
    }
    return count;
}

// Get the subtree of T1 connecting the leaves of T1 that are also in T2.
template <typename Tree1, typename Tree2, typename Tree_Out_t>
pair<unique_ptr<Tree_Out_t>,unique_ptr<Tree_Out_t>>
get_induced_trees2(const Tree1& T1,
                  std::function<const typename Tree1::node_type*(const typename Tree1::node_type*,const typename Tree1::node_type*)> MRCA_of_pair1,
                  const Tree2& T2,
                  std::function<const typename Tree2::node_type*(const typename Tree2::node_type*,const typename Tree2::node_type*)> MRCA_of_pair2)
{
    LOG(WARNING)<<"T1 = "<<newick_string(T1);
    LOG(WARNING)<<"n_leaves(T1) = "<<n_leaves(T1);
    LOG(WARNING)<<"n_leaves(T2) = "<<n_leaves(T2);

    // 1. First construct the induced tree for T2.

    // 1a. Find the nodes of T2 that corresponds to leaves of T1.
    //     Note that some of these nodes could be ancestral to other ones in T2.
    auto T2_nodes_from_T1_leaves = get_induced_nodes(T1,T2);
    LOG(WARNING)<<T2_nodes_from_T1_leaves.size()<<" leaves of T1 are in T2";

    // 1b. Actually construct the induced tree for T2.
    //     It might have fewer leaves than in T2_nodes_from_T1_leaves, if some of the nodes are ancestral to others.
    auto induced_tree2 = get_induced_tree<Tree2, Tree_Out_t>(T2_nodes_from_T1_leaves, MRCA_of_pair2);
    LOG(WARNING)<<"n_leaves(induced_tree2) = "<<n_leaves(*induced_tree2);
    LOG(WARNING)<<"induced-tree2a = "<<newick_string(*induced_tree2);

    // 1c. Rename internal nodes of synth/taxonomy to ottXXX instead of taxon name
    for(auto node: iter_post(*induced_tree2))
        if (node->has_ott_id())
            node->set_name("ott"+std::to_string(node->get_ott_id()));
    LOG(WARNING)<<"induced_tree2 = "<<newick_string(*induced_tree2);

    //FIXME - handle cases like (Homo sapiens, Homo) by deleting monotypic nodes at the root.

    // 2. Keep only leaves from T1 that (a) have an ottid in T2 and (b) map to a leaf of induced_tree2
    auto ottid_to_induced_tree2_node = otc::get_ottid_to_const_node_map(*induced_tree2);
    std::vector<const typename Tree1::node_type*> induced_leaves1;
    for(auto leaf: iter_leaf_const(T1))
    {
        if (not leaf->has_ott_id())
        {
            LOG(WARNING)<<"Dropping tip: no ott id";
            continue;
        }

        auto it = ottid_to_induced_tree2_node.find(leaf->get_ott_id());
        if (it == ottid_to_induced_tree2_node.end())
            LOG(WARNING)<<"Dropping tip "<<leaf->get_ott_id()<<": not found in induced taxonomy-or-synth tree.";
        else if (not it->second->is_tip())
            LOG(WARNING)<<"Dropping higher taxon tip "<<leaf->get_ott_id();
        else
            induced_leaves1.push_back(leaf);
    }
    LOG(WARNING)<<"keeping "<<induced_leaves1.size()<<" leaves from T1";

    // 3. Construct the induced tree for T1
    auto induced_tree1 = get_induced_tree<Tree1, Tree_Out_t>(induced_leaves1, MRCA_of_pair1);
    LOG(WARNING)<<"n_leaves(induced_tree1) = "<<n_leaves(*induced_tree1);
    LOG(WARNING)<<"induced_tree1 = "<<newick_string(*induced_tree1);

    assert(n_leaves(*induced_tree1) == n_leaves(*induced_tree2));

    return {std::move(induced_tree1), std::move(induced_tree2)};
}


template <typename N>
bool is_fake_tip(const N* n)
{
    return (n->get_out_degree() > 0) and (n->get_first_child()->get_name().empty());
}

/*
 * other_tree is (currently) either the taxonomy tree, or the synth tree
 */

template<typename QT, typename TT, typename QM, typename TM>
json conflict_with_tree_impl(const QT & query_tree,
                             const TT & other_tree,
                             std::function<const QM*(const QM*,const QM*)> & query_mrca,
                             std::function<const TM*(const TM*,const TM*)> & other_mrca,
                             const RichTaxonomy& Tax)
{
    conflict_stats stats;
    auto log_supported_by = [&stats](const QM* node2, const QM* node1) {
        if (is_fake_tip(node1))
            stats.add_terminal(node2,node1);
        else
            stats.add_supported_by(node2, node1);
    };
    auto log_partial_path_of = [&stats](const QM* node2, const QM* node1) {
        if (is_fake_tip(node1))
            stats.add_terminal(node2,node1);
        else
            stats.add_partial_path_of(node2, node1);
    };
    auto log_conflicts_with = [&stats](const QM* node2, const QM* node1) {
        stats.add_conflicts_with(node2, node1);
    };
    auto log_resolved_by = [&stats](const QM* node2, const QM* node1) {
        stats.add_resolved_by(node2, node1);
    };
    auto log_terminal = [&stats](const QM* node2, const QM* node1) {
        // Node1 might have an empty name if is a fake tip.
        if (not node1->get_name().empty())
            stats.add_terminal(node2, node1);
    };

    {
        auto induced_trees = get_induced_trees2<QT,TT,ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);

        perform_conflict_analysis(*induced_trees.first,
                                  *induced_trees.second,
                                  log_supported_by,
                                  log_partial_path_of,
                                  log_conflicts_with,
                                  log_resolved_by,
                                  log_terminal);
    }
/*
//  See PROBLEM notes above, on why we don't record when node2 resolves node1.
//    auto log_resolves = [&stats](const QM* node1, const QM* node2) {
//      if (not is_fake_tip(node1))
//          stats.add_resolves(node2, node1);
//    };
//    auto do_nothing = [](const QM*, const QM*) {};

    //The induced trees are destructively modified, so we can't reuse them.
    {
        auto induced_trees = get_induced_trees2<QT,TT,ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);

        perform_conflict_analysis(*induced_trees.second,
                                  *induced_trees.first,
                                  do_nothing,
                                  do_nothing,
                                  do_nothing,
                                  log_resolves,
                                  do_nothing);
    }
*/

    {
        auto induced_trees = get_induced_trees2<QT,TT,ConflictTree>(query_tree, query_mrca, other_tree, other_mrca);
        return stats.get_json(*induced_trees.first, Tax);
    }
}

json conflict_with_taxonomy(const ConflictTree& query_tree, const RichTaxonomy& Tax) {
    auto & taxonomy = Tax.get_tax_tree();

    using cfunc = std::function<const cnode_type*(const cnode_type*,const cnode_type*)>;
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;
    cfunc query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    return conflict_with_tree_impl(query_tree, taxonomy, query_mrca, taxonomy_mrca, Tax);
}

json conflict_with_summary(const ConflictTree& query_tree,
                           const SummaryTree_t& summary,
                           const RichTaxonomy& Tax) {
    using snode_type = SummaryTree_t::node_type;
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);};
    std::function<const snode_type*(const snode_type*,const snode_type*)> summary_mrca = [](const snode_type* n1, const snode_type* n2) {
        return find_mrca_via_traversal_indices(n1,n2);
    };
    return conflict_with_tree_impl(query_tree, summary, query_mrca, summary_mrca, Tax);
}

template<typename T>
bool has_internal_node_names(const T& t) {
    for(const auto nd: iter_post_const(t))
        if (not nd->get_name().empty()) return false;
    return true;
}

template <typename Tree>
pair<int,int> prune_unmapped_leaves(Tree& tree, const RichTaxonomy& tax)
{
    int mapped_leaves = 0;
    int unmapped_leaves = 0;
    vector<typename Tree::node_type*> leaves;
    for(auto leaf: iter_leaf(tree))
    {
        if (leaf->has_ott_id())
        {
            auto tax_node = tax.included_taxon_from_id(leaf->get_ott_id());
            if (tax_node)
            {
                // Handle forwards
                auto id = tax_node->get_ott_id();
                if (id != leaf->get_ott_id()) leaf->set_ott_id(id);

                // Count as mapped
                mapped_leaves++;
                continue;
            }
        }

        // Mark leaf for deletion
        leaves.push_back(leaf);
        unmapped_leaves++;
    }

    for(auto leaf: leaves)
    {
        while (leaf and leaf->is_tip())
        {
            auto parent = leaf->get_parent();
            if (parent)
            {
                leaf->detach_this_node();
                delete leaf;
                leaf = parent;
            }
            else
            {
                delete leaf;
                tree._set_root(nullptr);
            }
        }
        assert(tree.get_root());
    }
    return {mapped_leaves, unmapped_leaves};
}

template <typename Tree>
void delete_tip_and_monotypic_ancestors(Tree& tree, typename Tree::node_type* node)
{
    assert(node->is_tip());
    while (node and node->is_tip())
    {
        auto parent = node->get_parent();
        if (not parent)
            tree._set_root(nullptr);
        else
            node->detach_this_node();
        delete node;
        node = parent;
    }
}

template <typename Tree>
void delete_subtree_and_monotypic_ancestors(Tree& tree, typename Tree::node_type* node)
{
    auto parent = node->get_parent();
    tree.prune_and_delete(node);
    if (parent and parent->is_tip())
        delete_tip_and_monotypic_ancestors(tree, parent);
}

template <typename Tree>
void prune_duplicate_ottids(Tree& tree)
{
    vector<typename Tree::node_type*> leaves;
    for(auto leaf: iter_leaf(tree))
        leaves.push_back(leaf);

    map<OttId, typename Tree::node_type*> node_ptrs;
    for(auto leaf: leaves)
    {
        if (not leaf->has_ott_id()) continue;

        auto id = leaf->get_ott_id();

        // If the OTT id is new, then add the node as canonical representative of the OTT id
        if (not node_ptrs.count(id))
            node_ptrs.insert({id, leaf});
        // Otherwise delete the non-canonical OTT id and its ancestors
        else
            delete_tip_and_monotypic_ancestors(tree, leaf);
    }
}

void check_all_leaves_have_ott_ids(const ConflictTree& query_tree)
{
    for(auto leaf: iter_leaf_const(query_tree))
        if (not leaf->has_ott_id())
        {
            if (leaf->get_name().empty())
                throw OTCBadRequest()<<"Un-named leaf has no OTT id!";
            else
                throw OTCBadRequest()<<"Leaf '"<<leaf->get_name()<<"' has no OTT id!";
        }

}

void check_all_nodes_have_node_names(const ConflictTree& query_tree)
{
    for(const auto nd: iter_post_const(query_tree))
    {
        if (nd->get_name().empty())
        {
            auto E = OTCBadRequest();
            E << "Query tree has unnamed node";
            if (nd->has_ott_id()) {
                E << " with OTT Id=" << nd->get_ott_id();
            }
            E << "!";
            throw E;
        }
    }
}

// Find the smallest set C of taxa (leaf or internal) that we need to add as children of `id` so that
// i) all descendents of id are descendants of one of these children
// ii) all taxa in C are present in the summary tree
vector<OttId> extra_children_for_node(OttId id, const SummaryTree_t& summary, const RichTaxonomy& taxonomy)
{
    // We could maybe improve speed by:
    // - following only nodes which are ancestors of synthesis leaves (i.e. unpruned leaves)
    // - following only nodes which are ancestors of input phylogeny leaves, exemplified to turn higher taxon leaves -> leaf leaves.

    auto& id_to_node = summary.get_data().id_to_node;
    auto& tax_id_to_node = taxonomy.get_tax_tree().get_data().id_to_node;

    // If the node is in synth already, then we don't need to add any children.
    if (id_to_node.count(id)) return {};

    // Growing list of descendants of `id` that are not in synth.
    vector<OttId> children;
    vector<OttId> bad_parents({id});

    // Walk a frontier leafward from id:
    for(std::size_t i=0;i<bad_parents.size();i++)
    {
        auto parent_id = bad_parents[i];
        auto tax_node = tax_id_to_node.at(bad_parents[i]);

        auto parent_node = taxonomy.included_taxon_from_id(parent_id);

        for(auto c: iter_child_const(*tax_node))
        {
            auto child_id = c->get_ott_id();
            auto child_node = taxonomy.included_taxon_from_id(child_id);

            std::cerr << "Lost taxon " << parent_node->get_name() << " ("<<parent_id<<") -> ";
            std::cerr << "taxon " << child_node->get_name() << " (" << child_id << ") of degree " << child_node->get_out_degree() << ":";

            if (id_to_node.count(child_id))
            {
                // In synth, we can stop expanding the frontier at this node.
                std::cerr << "GOOD\n"; 
                children.push_back(child_id);
            }
            else
            {
                // Not in synth, we need to consider the children of this node.
                std::cerr << "BAD\n";
                bad_parents.push_back(child_id);
            }
        }
    }
    return children;
}

// Remove leaves (and their monotypic ancestors) in the query tree
// that are the ancestors of other leaves in the query tree.
void prune_ancestral_leaves(ConflictTree& query_tree, const RichTaxonomy& taxonomy)
{
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;
    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };

    auto taxonomy_nodes_from_query_leaves = get_induced_nodes(query_tree, taxonomy.get_tax_tree());
    auto induced_taxonomy = get_induced_tree<RichTaxTree, ConflictTree>(taxonomy_nodes_from_query_leaves,
                                                                        taxonomy_mrca);
    auto ottid_to_induced_tax_node = get_ottid_to_const_node_map(*induced_taxonomy);

    LOG(WARNING)<<"induced taxonomy has "<<n_leaves(*induced_taxonomy)<<" leaves.";

    vector<cnode_type*> nodes_to_prune;
    for(auto leaf: iter_leaf(query_tree))
    {
        auto ottid = leaf->get_ott_id();
        auto tax_node = ottid_to_induced_tax_node.at(ottid);
        if (not tax_node->is_tip())
            nodes_to_prune.push_back(leaf);
    }

    for(auto node: nodes_to_prune)
        delete_tip_and_monotypic_ancestors(query_tree, node);

    LOG(WARNING)<<"query tree pruned down to "<<n_leaves(*induced_taxonomy)<<" leaves.";
}


/*
 * OK, a general approach to normalizing two trees for comparison:
 * 0. All tips must have OTT IDs.
 *    The taxonomy and ids are used to assess if a tip is identical, ancestral, or descended from another tip.
 * 1. Remove duplicate tips.
 * 2. Remove tips that are ancestral to other tips.
 * At this point, we COULD just add as children all OTT leaves that are descendants of each higher taxon.
 *
 * A possibly simpler approach would be to.
 * 3. Remove tips that have no shared descendants with tips in the other tree.
 * 4. For each tip, add as children tips from the other tree that are descendants.
 */

string conflict_ws_method(const SummaryTree_t& summary,
                          const RichTaxonomy & taxonomy,
                          std::unique_ptr<ConflictTree>& query_tree,
                          const string& tree2s)
{
    // 0. Prune unmapped leaves.
    auto leaf_counts = prune_unmapped_leaves(*query_tree, taxonomy);
    if (leaf_counts.first < 3)
        throw OTCBadRequest()<<"Query tree has only "<<leaf_counts.first<<" leaves with an OTT id!";

    // 1. Check that all leaves in input tree have OTT ids
    check_all_leaves_have_ott_ids(*query_tree);

    // 2. Check that all leaves in input tree have node names
    check_all_nodes_have_node_names(*query_tree);

    // 3. Prune leaves with duplicate ott ids
    prune_duplicate_ottids(*query_tree);

    // 4. Prune leaves of the query that are ancestral to other query leaves
    prune_ancestral_leaves(*query_tree, taxonomy);

    // 5. Add children to higher-taxon leaves of query_tree that are not in the synth tree.
    if (tree2s == "synth")
    {
        map<cnode_type*,vector<OttId>> children_to_add;
        for(auto leaf: iter_leaf(*query_tree))
        {
//          LOG(WARNING)<<"Considering leaf "<<leaf->get_ott_id()<<":";
            auto leaf_id = leaf->get_ott_id();
            auto c = extra_children_for_node(leaf_id, summary, taxonomy);
            if (not c.empty())
            {
                children_to_add.insert({leaf,c});
                LOG(WARNING)<<"ottid "<<leaf->get_ott_id()<<" has frontier ";
                for(auto id: c)
                    LOG(WARNING)<<"   "<<id;
            }
        }

        // Add nodes with the specified ottids
        for(const auto& job: children_to_add)
        {
            auto leaf = job.first;
            auto& child_ids = job.second;
            for(auto id: child_ids)
                query_tree->create_child(leaf)->set_ott_id(id);
        }
    }

    // 5. Compute depth and number of tips for query_tree.
    compute_depth(*query_tree);
    compute_tips(*query_tree);

    if (tree2s == "ott") {
        return conflict_with_taxonomy(*query_tree, taxonomy).dump(1);
    } else if (tree2s == "synth") {
        return conflict_with_summary(*query_tree, summary, taxonomy).dump(1);
    }
    throw OTCBadRequest() << "tree2 = '" << tree2s << "' not recognized!";
}

string newick_conflict_ws_method(const SummaryTree_t& summary,
                                 const RichTaxonomy & taxonomy,
                                 const string& tree1s,
                                 const string& tree2s)
{
    try
    {
        LOG(WARNING)<<"newick conflict: tree1s = '"<<tree1s<<"'   tree1s = '"<<tree2s<<"'";
        auto query_tree = tree_from_newick_string<ConflictTree>(tree1s);
        return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
    }
    catch (otc::OTCParsingError& e)
    {
        LOG(WARNING)<<"newick conflict: error parsing newick string '"<<tree1s<<"'";
        throw OTCBadRequest() << "Error parsing newick:  \n:"<<e.what();
    }
}

string phylesystem_conflict_ws_method(const SummaryTree_t& summary,
                                      const RichTaxonomy & taxonomy,
                                      const string& tree1s,
                                      const string& tree2s)
{
    LOG(WARNING)<<"phylesystem conflict: tree1s = '"<<tree1s<<"'   tree1s = '"<<tree2s<<"'";
    auto query_tree = get_phylesystem_tree<ConflictTree>(tree1s);
    return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
}




} //namespace otc

int main( const int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc,argv);
        return run_server(args);
    } catch (std::exception& e) {
        LOG(ERROR) <<"otc-tol-ws: Error! " << e.what() << std::endl;
        return 1;
    }
}
