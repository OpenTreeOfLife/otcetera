#include <regex>
#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"
#include "otc/conflict.h"
#include "otc/tree_operations.h"
INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace boost::property_tree;
using json=nlohmann::json;

namespace otc {
const int OK = restbed::OK;
extern TreesToServe tts;

void add_taxon_info(const RichTaxonomy & , const RTRichTaxNode & nd_taxon, json & taxonrepr) {
    const auto & taxon_data = nd_taxon.get_data();
    taxonrepr["tax_sources"] = taxon_data.get_sources_json();
    taxonrepr["name"] = string(taxon_data.get_name());
    taxonrepr["uniqname"] = string(taxon_data.get_uniqname());
    taxonrepr["rank"] = string(taxon_data.get_rank());
    taxonrepr["ott_id"] = nd_taxon.get_ott_id();    
}

void add_basic_node_info(const RichTaxonomy & taxonomy, const SumTreeNode_t & nd, json & noderepr) {
    noderepr["node_id"] = nd.get_name();
    noderepr["num_tips"] = nd.get_data().num_tips;
    if (nd.has_ott_id()) {
        auto nd_id = nd.get_ott_id();
        const auto * nd_taxon = taxonomy.taxon_from_id(nd_id);
        if (nd_taxon == nullptr) {
            throw OTCError() << "OTT Id " << nd_id << " not found in taxonomy! Please report this bug";
        }
        json taxon;
        add_taxon_info(taxonomy, *nd_taxon, taxon);
        noderepr["taxon"] = taxon;
    }
}

inline void add_str_to_str_or_vec_string(json & o, const string * first, const string * second) {
    const char * studyc = first->c_str();
    if (o.count(studyc)) {
        if (o[studyc].is_array()) {
            o[studyc].push_back(*second);
        } else {
            string prev = o[studyc].get<string>();
            json v = {prev, *second};
            o[studyc] = v;
        }
    } else {
        o[studyc] = *second;
    }
}

void add_support_info_vec(const char * tag,
                          const vec_src_node_ids & v,
                          json & noderepr,
                          set<string> & usedSrcIds,
                          const string * extra_src=nullptr,
                          const string * extra_node_id=nullptr) {
    json o;
    for (const auto & study_node_pair : v ) {
        usedSrcIds.insert(*study_node_pair.first);
        add_str_to_str_or_vec_string(o, study_node_pair.first, study_node_pair.second);
    }
    if (extra_src && extra_node_id) {
        add_str_to_str_or_vec_string(o, extra_src, extra_node_id);
    }
    noderepr[tag] = o;
}

void add_node_support_info(const TreesToServe & tts,
                           const SumTreeNode_t & nd,
                           json & noderepr,
                           set<string> & usedSrcIds) {
    const auto & taxonomy = tts.get_taxonomy();
    const auto & d = nd.get_data();
    const string * extra_src = nullptr;
    const string * extra_node_id = nullptr;
    if (nd.has_ott_id()) {
        extra_src = &(taxonomy.get_version());
        extra_node_id = &(nd.get_name());
        usedSrcIds.insert(taxonomy.get_version());
    }
    if (extra_src != nullptr || !d.supported_by.empty()) {
        add_support_info_vec("supported_by", d.supported_by, noderepr, usedSrcIds, extra_src, extra_node_id);
    }
    if (!d.conflicts_with.empty()) {
        add_support_info_vec("conflicts_with", d.conflicts_with, noderepr, usedSrcIds);
    }
    if (!d.resolves.empty()) {
        add_support_info_vec("resolves", d.resolves, noderepr, usedSrcIds);
    }
    if (!d.partial_path_of.empty()) {
        add_support_info_vec("partial_path_of", d.partial_path_of, noderepr, usedSrcIds);
    }
    if (!d.terminal.empty()) {
        add_support_info_vec("terminal", d.terminal, noderepr, usedSrcIds);
    }
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
    auto n2nit = tree_data.name_to_node.find(node_id);
    if (n2nit != tree_data.name_to_node.end()) {
        return n2nit->second;
    }
    std::smatch matches;
    if (std::regex_match(node_id, matches, mrca_id_pattern)) {
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
    if (std::regex_match(node_id, matches, ott_id_pattern)) {
        auto bt_it = tree_data.broken_taxa.find(node_id);
        if (bt_it == tree_data.broken_taxa.end()) {
            return nullptr;
        }
        // if found we return the MRCA pointer in the first slot of the pair.
        was_broken = true;
        return bt_it->second.first;
    }
    return nullptr;
}


void about_ws_method(const TreesToServe &tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     bool include_sources,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    const auto & taxonomy = tts.get_taxonomy();
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
    auto root_node = tree_ptr->get_root();
    add_basic_node_info(taxonomy, *root_node, root);
    response["root"] = root;
    response_str = response.dump(1);
}

void tax_about_ws_method(const TreesToServe &tts,
                         string & response_str,
                         int & status_code) {
    const auto & taxonomy = tts.get_taxonomy();
    status_code = OK;
    json response;
    string weburl;
    weburl = "https://tree.opentreeoflife.org/about/taxonomy-version/ott";
    weburl += taxonomy.get_version_number();
    response["weburl"] = weburl;
    response["author"] = "open tree of life project";
    response["name"] = "ott";
    response["source"] = taxonomy.get_version();
    response["version"] = taxonomy.get_version_number();
    response_str = response.dump(1);
}


inline void add_lineage(json & j, const SumTreeNode_t * focal, const RichTaxonomy & taxonomy, set<string> & usedSrcIds) {
    json lineage_arr;
    const SumTreeNode_t * anc = focal->get_parent();
    if (!anc) {
        vector<string> c;
        j["lineage"] = c;
        return;
    }
    while (anc) {
        json ancj;
        add_basic_node_info(taxonomy, *anc, ancj);
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
        if (srcTag == taxonomy.get_version()) {
            jt["taxonomy"] = taxonomy.get_version();
        } else {
            const auto & simentry = sta->source_id_map.at(srcTag);
            jt = simentry;
        }
        sim[srcTag] = jt;   
    }
    j["source_id_map"] = sim;
}

void node_info_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     bool include_lineage,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_node_by_id_str(*tree_ptr, node_id, was_broken);
    if (focal == nullptr) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return;
    }
    const auto & taxonomy = tts.get_taxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    add_basic_node_info(taxonomy, *focal, response);
    add_node_support_info(tts, *focal, response, usedSrcIds);
    if (include_lineage) {
        add_lineage(response, focal, taxonomy, usedSrcIds);
    }
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    if (was_broken) {
        response["response_for_mrca_of_broken_taxon"] = true; //@TODO: discuss and document
    }
    response_str = response.dump(1);
}

void mrca_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const vector<string> & node_id_vec,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_node_by_id_str(*tree_ptr, node_id, was_broken);
        if (n == nullptr) {
            response_str = "node_id \"";
            response_str += node_id;
            response_str += "\" was not recognized.\n";
            status_code = 400;
            return;
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
    bool is_broken = false;
    if (focal == nullptr) {
        response_str = "MRCA of taxa was not found.\n";
        status_code = 400;
        return;
    }
    const auto & taxonomy = tts.get_taxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    json mrcaj;
    add_basic_node_info(taxonomy, *focal, mrcaj);
    set<string> usedSrcIds;
    add_node_support_info(tts, *focal, mrcaj, usedSrcIds);
    if (!focal->has_ott_id()) {
        auto anc = focal->get_parent();
        assert(anc != nullptr);
        while (!anc->has_ott_id()) {
            anc = anc->get_parent();
            if (anc == nullptr) {
                response_str = "No ancestors were taxa. That is odd.\n";
                status_code = 500;
                return;
            }
        }
        json nt;
        const RTRichTaxNode * anc_taxon = taxonomy.taxon_from_id(anc->get_ott_id());
        if (anc_taxon == nullptr) {
            throw OTCError() << "Ancd OTT Id " << anc->get_ott_id() << " not found in taxonomy! Please report this bug";
        }
        add_taxon_info(taxonomy, *anc_taxon, nt);
        response["nearest_taxon"] = nt;
    }
    response["mrca"] = mrcaj;
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    response_str = response.dump(1);
}


const SumTreeNode_t * get_node_for_subtree(const SummaryTree_t * tree_ptr,
                                           const string & node_id, 
                                           int height_limit,
                                           uint32_t tip_limit,
                                           string & response_str,
                                           int & status_code) {
    assert(tree_ptr != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_node_by_id_str(*tree_ptr, node_id, was_broken);
    if (focal == nullptr) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return nullptr;
    }
    if (was_broken) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return nullptr;
    }
    if (focal->get_data().num_tips > tip_limit && height_limit < 0) {
        response_str = "The requested subtree is too large to be returned via the API. Download the entire tree.\n";
        status_code = 400;
        return nullptr;  
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
            for (auto & p : d.supported_by) {
                study_id_set.insert(p.first);
            }
            const string & id_str = nd->get_name(); // in this tree the "name" is really the "ott###" string
            if (nns != NodeNameStyle::NNS_ID_ONLY && nd->has_ott_id()) {
                const auto * tr = taxonomy.taxon_from_id(nd->get_ott_id());
                if (tr == nullptr) {
                    throw OTCError() << "OTT Id " << nd->get_ott_id() << " in namer not found in taxonomy! Please report this bug";
                }
                string taxon_name = string(tr->get_data().get_uniqname());
                if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
                    string ret;
                    ret.reserve(taxon_name.length() + 1 + id_str.length());
                    ret = taxon_name;
                    ret += ' ';
                    ret += id_str;
                    return ret;    
                } else {
                    return taxon_name;
                }
            }
            return id_str;
        }
        static string ott_id_to_idstr(OttId ott_id) {
            string ret;
            ret = "ott";
            ret += std::to_string(ott_id);
            return ret;
        }
        string operator()(const RTRichTaxNode *nd) const {
            if (nns == NodeNameStyle::NNS_ID_ONLY) {
                return ott_id_to_idstr(nd->get_ott_id());
            }
            const string & taxon_name = string(nd->get_data().get_uniqname());
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


void induced_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const vector<string> & node_id_vec,
                     NodeNameStyle label_format, 
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    set<const SumTreeNode_t *> tip_nodes;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_node_by_id_str(*tree_ptr, node_id, was_broken);
        tip_nodes.insert(n);
        if (n == nullptr) {
            response_str = "node_id \"";
            response_str += node_id;
            response_str += "\" was not recognized.\n";
            status_code = 400;
            return;
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
    bool is_broken = false;
    if (focal == nullptr) {
        response_str = "MRCA of taxa was not found.\n";
        status_code = 400;
        return;
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
    const auto & taxonomy = tts.get_taxonomy();
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    writeVisitedNewick(out, visited, focal, nnsbs);
    json response;
    response["newick"] = out.str();
    json ss_arr;
    for (auto study_it_ptr : nnsbs.study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    response["supporting_studies"] = ss_arr;
    response_str = response.dump(1);
}

void newick_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     NodeNameStyle label_format, 
                     int height_limit,
                     string & response_str,
                     int & status_code) {
    const uint32_t NEWICK_TIP_LIMIT = 25000;
    const SumTreeNode_t * focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT, response_str, status_code);
    if (focal == nullptr) {
        return;
    }
    assert(sta != nullptr);
    const auto & taxonomy = tts.get_taxonomy();
    status_code = OK;
    json response;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    write_newick_generic<const SumTreeNode_t *, NodeNamerSupportedByStasher>(out, focal, nnsbs, height_limit);
    response["newick"] = out.str();
    json ss_arr;
    for (auto study_it_ptr : nnsbs.study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    response["supporting_studies"] = ss_arr;
    response_str = response.dump(1);
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
    add_basic_node_info(taxonomy, *nd, j);
    add_node_support_info(tts, *nd, j, usedSrcIds);
}

void arguson_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     int height_limit,
                     string & response_str,
                     int & status_code) {
    const uint NEWICK_TIP_LIMIT = 25000;
    auto focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT, response_str, status_code);
    if (focal == nullptr) {
        return;
    }
    assert(sta != nullptr);
    const auto & taxonomy = tts.get_taxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    json a;
    write_arguson(a, tts, sta, taxonomy, focal, height_limit, usedSrcIds);
    add_lineage(a, focal, taxonomy, usedSrcIds);
    response["arguson"] = a;
    response_str = response.dump(1);
}

// taxon_info
void tax_service_add_taxon_info(const RichTaxonomy & taxonomy, const RTRichTaxNode & nd_taxon, json & taxonrepr) {
    add_taxon_info(taxonomy, nd_taxon, taxonrepr);
    taxonrepr["source"] = taxonomy.get_version(); //TBD "source" ?
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

void taxon_info_ws_method(const TreesToServe & tts,
                          const RTRichTaxNode * taxon_node,
                          bool include_lineage,
                          bool include_children,
                          bool include_terminal_descendants,
                          string & response_str,
                          int & status_code) {
    const auto & taxonomy = tts.get_taxonomy();
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
    response_str = response.dump(1);
    status_code = OK;
}


void taxonomy_mrca_ws_method(const TreesToServe & tts,
                             const OttIdSet & ott_id_set,
                             string & response_str,
                             int & status_code) {
    const auto & taxonomy = tts.get_taxonomy();
    const RTRichTaxNode * focal = nullptr;
    bool first = true;
    for (auto ott_id : ott_id_set) {
        const RTRichTaxNode * n = taxonomy.taxon_from_id(ott_id);
        if (n == nullptr) {
            response_str = "ott_id \"";
            response_str += ott_id;
            response_str += "\" was not recognized.\n";
            status_code = 400;
            return;
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
    bool is_broken = false;
    if (focal == nullptr) {
        response_str = "MRCA of taxa was not found. Please report this bug!\n";
        status_code = 400;
        return;
    }
    status_code = OK;
    json response;
    tax_service_add_taxon_info(taxonomy, *focal, response);
    response_str = response.dump(1);
}

void taxon_subtree_ws_method(const TreesToServe & tts,
                          const RTRichTaxNode * taxon_node,
                          NodeNameStyle label_format, 
                          string & response_str,
                          int & status_code) {
    const auto & taxonomy = tts.get_taxonomy();
    assert(taxon_node != nullptr);
    const auto & taxonomy_tree = taxonomy.getTaxTree();
    json response;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    int height_limit = -1;
    write_newick_generic<const RTRichTaxNode *, NodeNamerSupportedByStasher>(out, taxon_node, nnsbs, height_limit);
    response["newick"] = out.str();
    response_str = response.dump(1);
    status_code = OK;
}

using cnode_type = ConflictTree::node_type;

struct conflict_stats
{
    map<string, set<const cnode_type*>> supported_by;
    map<string, set<const cnode_type*>> partial_path_of;
    map<string, set<const cnode_type*>> conflicts_with;
    map<string, set<const cnode_type*>> resolved_by;
    map<string, set<const cnode_type*>> resolves;
    map<string, set<const cnode_type*>> terminal;

};

void conflict_with_taxonomy(const ConflictTree& query_tree,
			    const RichTaxTree& taxonomy)
{
    using tnode_type = RTRichTaxNode;

    conflict_stats stats;

    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {return mrca_from_depth(n1,n2);};
    std::function<const tnode_type*(const tnode_type*,const tnode_type*)> taxonomy_mrca;

    auto log_supported_by    = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.supported_by[node1->get_name()].insert(node2);};
    auto log_partial_path_of = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.partial_path_of[node1->get_name()].insert(node2);};
    auto log_conflicts_with  = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.conflicts_with[node1->get_name()].insert(node2);};
    auto log_resolved_by     = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.resolved_by[node1->get_name()].insert(node2);};
    auto log_terminal        = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.terminal[node1->get_name()].insert(node2);};

    auto ottid_to_query_node    = otc::get_ottid_to_const_node_map(query_tree);
    auto ottid_to_taxonomy_node = otc::get_ottid_to_const_node_map(taxonomy);

    perform_conflict_analysis(query_tree, ottid_to_query_node, query_mrca,
			      taxonomy, ottid_to_taxonomy_node, taxonomy_mrca,
			      log_supported_by,
			      log_partial_path_of,
			      log_conflicts_with,
			      log_resolved_by,
			      log_terminal);

    auto log_resolves        = [&stats](const cnode_type* node1, const cnode_type* node2) {stats.resolves[node1->get_name()].insert(node2);};
    auto do_nothing          = [](const cnode_type*, const cnode_type*) {};

    perform_conflict_analysis(taxonomy, ottid_to_taxonomy_node, taxonomy_mrca,
			      query_tree, ottid_to_query_node, query_mrca,
			      do_nothing,
			      do_nothing,
			      do_nothing,
			      log_resolves,
			      do_nothing);
    
}
			   

void conflict_with_summary(const ConflictTree& query_tree, const SummaryTree_t& summary)
{
    using snode_type = SummaryTree_t::node_type;

    conflict_stats stats;
    
    std::function<const cnode_type*(const cnode_type*,const cnode_type*)> query_mrca = [](const cnode_type* n1, const cnode_type* n2) {return mrca_from_depth(n1,n2);};
    std::function<const snode_type*(const snode_type*,const snode_type*)> summary_mrca;

    auto log_supported_by    = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.supported_by[node1->get_name()].insert(node2);};
    auto log_partial_path_of = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.partial_path_of[node1->get_name()].insert(node2);};
    auto log_conflicts_with  = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.conflicts_with[node1->get_name()].insert(node2);};
    auto log_resolved_by     = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.resolved_by[node1->get_name()].insert(node2);};
    auto log_terminal        = [&stats](const cnode_type* node2, const cnode_type* node1) {stats.terminal[node1->get_name()].insert(node2);};

    auto ottid_to_query_node = otc::get_ottid_to_const_node_map(query_tree);
    auto ottid_to_summary_node = otc::get_ottid_to_const_node_map(summary);

    perform_conflict_analysis(query_tree, ottid_to_query_node, query_mrca,
			      summary, ottid_to_summary_node, summary_mrca,
			      log_supported_by,
			      log_partial_path_of,
			      log_conflicts_with,
			      log_resolved_by,
			      log_terminal);

    auto log_resolves        = [&stats](const cnode_type* node1, const cnode_type* node2) {stats.resolves[node1->get_name()].insert(node2);};
    auto do_nothing          = [](const cnode_type*, const cnode_type*) {};

    perform_conflict_analysis(summary, ottid_to_summary_node, summary_mrca,
			      query_tree, ottid_to_query_node, query_mrca,
			      do_nothing,
			      do_nothing,
			      do_nothing,
			      log_resolves,
			      do_nothing);
    
}
			   



void conflict_ws_method(const SummaryTree_t& summary,
			const RichTaxonomy & taxonomy,
			const string& tree1s,
			const string& tree2s,
			string & response_str,
			int& status_code )
{
    ConflictTree query_tree;
    if (tree2s == "ott")
    {
	conflict_with_taxonomy(query_tree, taxonomy.getTaxTree());
    }
    else if (tree2s == "synth")
    {
	conflict_with_summary(query_tree, summary);
    }

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
