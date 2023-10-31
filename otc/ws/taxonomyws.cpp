#include <regex>
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/ws/node_namer_supported_by_stasher.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/ws/nexson/nexson.h"
#include <optional>
#include <string_view>
#include "otc/ws/extract.h"


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


namespace otc {

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

void tax_service_add_suppressed_taxon_info(const RichTaxonomy & taxonomy,
                                          const TaxonomyRecord & record,
                                          nlohmann::json & taxonrepr){
    add_taxon_record_info(taxonomy, record, taxonrepr);
    taxonrepr["source"] = string("ott") + taxonomy.get_version(); //TBD "source" ?
    taxonrepr["flags"] = flags_to_string_vec(record.flags);
    json syn_list = json::array();
    taxonrepr["synonyms"] = syn_list;
    taxonrepr["is_suppressed"] = true;
    taxonrepr["is_suppressed_from_synth"] = true;
}


json taxon_info_ws_method_j(const RichTaxonomy & taxonomy,
                            const RTRichTaxNode * taxon_node,
                            bool include_lineage,
                            bool include_children,
                            bool include_terminal_descendants)
{
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
    return response;
}

string taxon_info_ws_method(const RichTaxonomy & taxonomy,
                            const RTRichTaxNode * taxon_node,
                            bool include_lineage,
                            bool include_children,
                            bool include_terminal_descendants)
{
    auto response = taxon_info_ws_method_j(taxonomy, taxon_node, include_lineage, include_children, include_terminal_descendants);
    return response.dump(1);
}

string taxon_infos_ws_method(const RichTaxonomy & taxonomy,
                             const OttIdSet& ott_ids,
                             bool include_lineage,
                             bool include_children,
                             bool include_terminal_descendants)
{
    json response = json::array();
    for(auto ott_id: ott_ids)
    {
        json j;
        if (auto taxon_node = taxonomy.included_taxon_from_id(ott_id))
            j = taxon_info_ws_method_j(taxonomy, taxon_node, include_lineage, include_children, include_terminal_descendants);
        else
        {
            j["error"] = "unrecognized";
        }
        
        j["query"] = ott_id;
        response.push_back(j);
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
            focal = mrca_from_depth(focal, n);
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

string taxon_subtree_ws_method(const TreesToServe & tts, 
                               const RichTaxonomy & taxonomy,
                               const RTRichTaxNode * taxon_node,
                               NodeNameStyle label_format) {
    assert(taxon_node != nullptr);
    const auto & taxonomy_tree = taxonomy.get_tax_tree();
    json response;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy, tts);
    ostringstream out;
    int height_limit = -1;
    bool include_all_node_labels = true; // ??
    write_newick_generic<const RTRichTaxNode *, NodeNamerSupportedByStasher>(out, taxon_node, nnsbs, include_all_node_labels, height_limit);
    response["newick"] = out.str();
    return response.dump(1);
}

// How do we handle a situation where we add 3 taxa, and the 2nd one fails?
// Right now, re-adding the first taxa will fail.

std::string taxon_addition_ws_method(const TreesToServe & tts,
				     PatchableTaxonomy & taxonomy,
                                     const json& taxa)
{
    const auto & taxonomy_tree = taxonomy.get_tax_tree();

    int added = 0;
    for(auto& taxon: taxa)
    {
        auto name = extract_required_argument<string>(taxon, "name");
        auto ott_id = extract_required_argument<OttId>(taxon, "ott_id");
        auto parent_id = extract_required_argument<OttId>(taxon, "parent");
        auto rank = extract_required_argument<string>(taxon, "rank");

        auto [ok,error] = taxonomy.add_new_taxon(ott_id, parent_id, name, rank, "", "", {});

        if (not ok)
        {
            json j;
            j["ott_id"] = ott_id;
            j["error"] = error;
            j["added"] = added;

            throw OTCWebError()<<"Error adding taxon '"<<name<<"' with ott_id "<<ott_id<<j;
        }
        else
            added++;
    }

    json response = {{"added", added}};
    return response.dump(1);
}

} // namespace otc
