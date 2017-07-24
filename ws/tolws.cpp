#include <regex>
#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"
#include "otc/conflict.h"
#include "otc/tree_operations.h"
#include "nexson/nexson.h"
#include <boost/optional.hpp>
INITIALIZE_EASYLOGGINGPP


using namespace std;
using namespace boost::property_tree;
using json=nlohmann::json;

using boost::optional;

namespace otc {
const int OK = restbed::OK;
extern TreesToServe tts;

inline const string & get_taxon_unique_name(const RTRichTaxNode & nd_taxon) {
    return nd_taxon.get_name();
}

void add_taxon_info(const RichTaxonomy & , const RTRichTaxNode & nd_taxon, json & taxonrepr) {
    const auto & taxon_data = nd_taxon.get_data();
    taxonrepr["tax_sources"] = taxon_data.get_sources_json();
    taxonrepr["name"] = string(taxon_data.get_nonuniqname());
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

string taxon_nonuniquename(const RichTaxonomy& taxonomy, const SumTreeNode_t& nd)
{
    if (not nd.has_ott_id())
	throw OTCError()<<"Node "<<nd.get_name()<<" has no OTT id";

    auto id = nd.get_ott_id();
    auto nd_taxon = taxonomy.included_taxon_from_id(id);
    auto& taxon_data = nd_taxon->get_data();
    return string(taxon_data.get_nonuniqname());
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
                          const optional<string>& extra_src = boost::none,
                          const optional<string>& extra_node_id = boost::none) {
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
				     const optional<string>& extra_src = boost::none,
				     const optional<string>& extra_node_id = boost::none) {
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

const std::regex mrca_id_pattern("^mrca(ott\\d+)(ott\\d+)$");
const std::regex ott_id_pattern("^ott(\\d+)$");

const SumTreeNode_t * find_node_by_id_str(const SummaryTree_t & tree,
                                          const string & node_id,
                                          bool & was_broken) {
    was_broken = false;
    const auto & tree_data = tree.get_data();
    //auto n2nit = tree_data.name_to_node.find(node_id);
    //if (n2nit != tree_data.name_to_node.end()) {
    //    return n2nit->second;
    //}
    std::smatch matches;
    if (std::regex_match(node_id, matches, ott_id_pattern)) {
        long raw_ott_id = long_ott_id_from_name(node_id);
        if (raw_ott_id >= 0) {
            OttId ott_id = check_ott_id_size(raw_ott_id);
            auto i2nit = tree_data.id_to_node.find(ott_id);
            if (i2nit != tree_data.id_to_node.end()) {
                return i2nit->second;
            }
            LOG(WARNING) << "not finding " << ott_id << " extracted from " << node_id;
        }
        auto bt_it = tree_data.broken_taxa.find(node_id);
        if (bt_it == tree_data.broken_taxa.end()) {
            return nullptr;
        }
        // if found we return the MRCA pointer in the first slot of the pair.
        was_broken = true;
        return bt_it->second.first;
    }
    if (std::regex_match(node_id, matches, mrca_id_pattern)) {
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

string tax_about_ws_method(const RichTaxonomy & taxonomy) {
    json response;
    string weburl;
    response["author"] = "open tree of life project";
    response["name"] = "ott";
    weburl = "https://tree.opentreeoflife.org/about/taxonomy-version/ott";
    weburl += taxonomy.get_version_number();
    response["source"] = string("ott") + taxonomy.get_version();
    response["version"] = taxonomy.get_version_number();
    response["weburl"] = weburl;
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

string node_info_ws_method(const TreesToServe & tts,
                           const SummaryTree_t * tree_ptr,
                           const SummaryTreeAnnotation * sta,
                           const string & node_id,
                           bool include_lineage) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_required_node_by_id_str(*tree_ptr, node_id, was_broken);

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
    return response.dump(1);
}

string mrca_ws_method(const TreesToServe & tts,
                      const SummaryTree_t * tree_ptr,
                      const SummaryTreeAnnotation * sta,
                      const vector<string> & node_id_vec) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_required_node_by_id_str(*tree_ptr, node_id, was_broken);

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
    if (not focal) {
        throw OTCBadRequest("MRCA of taxa was not found.\n");
    }
    json response;
    response["synth_id"] = sta->synth_id;
    json mrcaj;
    set<string> usedSrcIds;
    add_node_support_info(tts, *focal, mrcaj, usedSrcIds);
    {
        auto locked_taxonomy = tts.get_readable_taxonomy();
        const auto & taxonomy = locked_taxonomy.first;
        add_basic_node_info(taxonomy, *focal, mrcaj);
        if (!focal->has_ott_id()) {
            auto anc = focal->get_parent();
            assert(anc != nullptr);
            while (!anc->has_ott_id()) {
                anc = anc->get_parent();
                if (anc == nullptr) {
                    throw OTCWebError("No ancestors were taxa. That is odd.\n");
                }
            }
            json nt;
            const RTRichTaxNode * anc_taxon = taxonomy.included_taxon_from_id(anc->get_ott_id());
            if (anc_taxon == nullptr) {
                throw OTCWebError() << "Ancd OTT Id " << anc->get_ott_id() << " not found in taxonomy! Please report this bug";
            }
            add_taxon_info(taxonomy, *anc_taxon, nt);
            response["nearest_taxon"] = nt;
        }
        add_source_id_map(response, usedSrcIds, taxonomy, sta);
    }
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
inline void writeVisitedNewick(std::ostream & out,
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
                writeVisitedNewick<C, T, Y>(out, visited, c, nodeNamer);
            }
        }
        if (!first) {
            out << ')';
        }
    }
    write_escaped_for_newick(out, nodeNamer(nd));
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
                 NodeNameStyle label_format) {
    assert(tree_ptr != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    set<const SumTreeNode_t *> tip_nodes;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_required_node_by_id_str(*tree_ptr, node_id, was_broken);
        tip_nodes.insert(n);
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
    bool is_broken = false;
    if (focal == nullptr) {
        throw OTCBadRequest() << "MRCA of taxa was not found.\n";
    }
    set<const SumTreeNode_t *> visited;
    visited.insert(focal);
    for (auto tni : tip_nodes) {
        auto cnd = tni;
        while (visited.count(cnd) == 0) {
            visited.insert(cnd);
            cnd = cnd->get_parent(); 
        } 
    }

    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    writeVisitedNewick(out, visited, focal, nnsbs);

    json response;
    response["newick"] = out.str();
    response["supporting_studies"] = get_supporting_studies(nnsbs.study_id_set);
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
    const uint NEWICK_TIP_LIMIT = 25000;
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


string taxonomy_mrca_ws_method(const RichTaxonomy & taxonomy,
                               const OttIdSet & ott_id_set) {
    json response;
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
    tax_service_add_taxon_info(taxonomy, *focal, response);
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

using cnode_type = ConflictTree::node_type;

template <typename N>
boost::optional<string> node_name_or_ottid(const N* node) {
    auto result = get_source_node_name(node->get_name());
    if (result) {
        return result;
    } else if (node->has_ott_id()) {
        return string("ott")+std::to_string(node->get_ott_id());
    } else {
        return boost::none;
    }
}

struct conflict_stats
{
    protected:
    template <typename T>
    void throw_otcerror_if_second_not_true(const T & result, const cnode_type* node1) {
        if (not result.second) {
            throw OTCError() << "key "<<*node_name_or_ottid(node1) << " occurs twice in a conflict_stats map!";
        }    
    }
    public:
    map<string, string> supported_by;
    map<string, string> partial_path_of;
    map<string, pair<string, const cnode_type*>> conflicts_with;
    map<string, string> resolved_by;
    map<string, string> resolves;
    map<string, string> terminal;

    void add_supported_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = supported_by.insert({*node_name_or_ottid(node1), *node_name_or_ottid(node2)});
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_partial_path_of(const cnode_type* node2, const cnode_type* node1) {
        auto result = partial_path_of.insert({*node_name_or_ottid(node1), *node_name_or_ottid(node2)});
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolved_by(const cnode_type* node2, const cnode_type* node1) {
        auto result = resolved_by.insert({*node_name_or_ottid(node1), *node_name_or_ottid(node2)});
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_resolves(const cnode_type* node2, const cnode_type* node1) {
        auto result = resolves.insert({*node_name_or_ottid(node1), *node_name_or_ottid(node2)});
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_terminal(const cnode_type* node2, const cnode_type* node1) {
        auto result = terminal.insert({*node_name_or_ottid(node1), *node_name_or_ottid(node2)});
        throw_otcerror_if_second_not_true(result, node1);
    }

    void add_conflicts_with(const cnode_type* node2, const cnode_type* node1) {
        auto name1 = *node_name_or_ottid(node1);
        auto name2 = *node_name_or_ottid(node2);
        auto present = conflicts_with.insert({name1,{name2,node2}});
        if (present.second) {
            const cnode_type* node2b = present.first->second.second;
            if (depth(node2) < depth(node2b)) {
                conflicts_with.erase(present.first);
                conflicts_with.insert({name1,{name2,node2}});
            }
        }
    }
    json get_json(const RichTaxonomy&) const;
};

using tnode_type = RTRichTaxNode;

int depth(const tnode_type* node) {
    return node->get_data().depth;
}

json get_node_status(const string& witness, string status, const RichTaxonomy& Tax) {
    json j;
    j["witness"] = witness;
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

json conflict_stats::get_json(const RichTaxonomy& Tax) const {
    json nodes;
    for(auto& x: supported_by) {
        nodes[x.first] = get_node_status(x.second, "supported_by", Tax);
    }
    for(auto& x: partial_path_of) {
        nodes[x.first] = get_node_status(x.second, "partial_path_of", Tax);
    }
    for(auto& x: terminal) {
        nodes[x.first] = get_node_status(x.second, "terminal", Tax);
    }
    for(auto& x: conflicts_with) {
        nodes[x.first] = get_node_status(x.second.first, "conflicts_with", Tax);
    }
    for(auto& x: resolves) {
        nodes[x.first] = get_node_status(x.second, "resolves", Tax);
    }
    for(auto& x: resolved_by) {
        nodes[x.first] = get_node_status(x.second, "resolved_by", Tax);
    }
    return nodes;
}

template<typename QT, typename TT, typename QM, typename TM>
void conflict_with_tree_impl(conflict_stats & stats,
                             const QT & query_tree,
                             const TT & other_tree,
                             std::function<const QM*(const QM*,const QM*)> & query_mrca,
                             std::function<const TM*(const TM*,const TM*)> & other_mrca) {
    auto log_supported_by = [&stats](const QM* node2, const QM* node1) {
        stats.add_supported_by(node2, node1);
    };
    auto log_partial_path_of = [&stats](const QM* node2, const QM* node1) {
        stats.add_partial_path_of(node2, node1);
    };
    auto log_conflicts_with = [&stats](const QM* node2, const QM* node1) {
        stats.add_conflicts_with(node2, node1);
    };
    auto log_resolved_by = [&stats](const QM* node2, const QM* node1) {
        stats.add_resolved_by(node2, node1);
    };
    auto log_terminal = [&stats](const QM* node2, const QM* node1) {
        stats.add_terminal(node2, node1);
    };
    auto ottid_to_query_node = otc::get_ottid_to_const_node_map(query_tree);
    auto ottid_to_taxonomy_node = otc::get_ottid_to_const_node_map(other_tree);
    perform_conflict_analysis(query_tree, ottid_to_query_node, query_mrca,
                              other_tree, ottid_to_taxonomy_node, other_mrca,
                              log_supported_by,
                              log_partial_path_of,
                              log_conflicts_with,
                              log_resolved_by,
                              log_terminal);
    auto log_resolves = [&stats](const QM* node1, const QM* node2) {
        stats.add_resolves(node2, node1);
    };
    auto do_nothing = [](const QM*, const QM*) {};
    perform_conflict_analysis(other_tree, ottid_to_taxonomy_node, other_mrca,
                              query_tree, ottid_to_query_node, query_mrca,
                              do_nothing,
                              do_nothing,
                              do_nothing,
                              log_resolves,
                              do_nothing);
}

json conflict_with_taxonomy(const ConflictTree& query_tree,
                            const RichTaxonomy& Tax) {
    auto & taxonomy = Tax.get_tax_tree();
    conflict_stats stats;
    using cfunc = std::function<const cnode_type*(const cnode_type*,const cnode_type*)>;
    using tfunc = std::function<const tnode_type*(const tnode_type*,const tnode_type*)>;
    cfunc query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    tfunc taxonomy_mrca = [](const tnode_type* n1, const tnode_type* n2) {
        return mrca_from_depth(n1,n2);
    };
    conflict_with_tree_impl(stats, query_tree, taxonomy, query_mrca, taxonomy_mrca);
    return stats.get_json(Tax);
}

json conflict_with_summary(const ConflictTree& query_tree,
                           const SummaryTree_t& summary,
                           const RichTaxonomy& Tax) {
    using snode_type = SummaryTree_t::node_type;
    conflict_stats stats;
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {
        return mrca_from_depth(n1,n2);};
    std::function<const snode_type*(const snode_type*,const snode_type*)> summary_mrca = [](const snode_type* n1, const snode_type* n2) {
        return find_mrca_via_traversal_indices(n1,n2);
    };
    conflict_with_tree_impl(stats, query_tree, summary, query_mrca, summary_mrca);
    return stats.get_json(Tax);
}
                           
template<typename T>
bool has_internal_node_names(const T& t) {
    for(const auto nd: iter_post_const(t)) {
        if (not node_name_or_ottid(nd->get_name())) {
            return false;
        }
    }
    return true;
}

template <typename Tree>
pair<int,int> prune_unmapped_leaves(Tree& tree)
{
    int mapped_leaves = 0;
    int unmapped_leaves = 0;
    vector<typename Tree::node_type*> leaves;
    for(auto leaf: iter_leaf(tree))
	if (not leaf->has_ott_id())
	{
	    leaves.push_back(leaf);
	    unmapped_leaves++;
	}
	else
	    mapped_leaves++;

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



string conflict_ws_method(const SummaryTree_t& summary,
			  const RichTaxonomy & taxonomy,
			  std::unique_ptr<ConflictTree>& query_tree,
			  const string& tree2s)
{
    // 0. Remove unmapped leaves.
    auto leaf_counts = prune_unmapped_leaves(*query_tree);
    if (leaf_counts.first < 3)
	throw OTCBadRequest()<<"Query tree has only "<<leaf_counts.first<<" leaves with an OTT id!";

    // 1. Check that all leaves have OTT ids
    for(auto leaf: iter_leaf(*query_tree))
	if (not leaf->has_ott_id())
	{
	    if (leaf->get_name().empty())
		throw OTCBadRequest()<<"Un-named leaf has no OTT id!";
	    else
		throw OTCBadRequest()<<"Leaf '"<<leaf->get_name()<<"' has no OTT id!";
	}

    // 2. Check that all leaves have node names
    for(const auto nd: iter_post_const(*query_tree)) {
        string name = nd->get_name();
        auto node_name = node_name_or_ottid(nd);
        if (not node_name) {
            auto E = OTCBadRequest();
            E << "Newick tree node with name='" << name << "'";
            if (nd->has_ott_id()) {
                E << " and OTT Id=" << nd->get_ott_id();
            }
            E << " does not have node name annotation!";
            throw E;
        }
    }

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
    auto query_tree = tree_from_newick_string<ConflictTree>(tree1s);
    return conflict_ws_method(summary, taxonomy, query_tree, tree2s);
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
