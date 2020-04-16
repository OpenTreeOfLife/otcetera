#include <regex>
#include "otc/ws/tolws.h"
#include "otc/ws/tolwsadaptors.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/tree_operations.h"
#include "otc/supertree_util.h"
#include "otc/ctrie/search_data_models.h"
#include "otc/ctrie/context_ctrie_db.h"
#include "otc/ws/nexson/nexson.h"
#include <optional>
#include <string_view>
#include "otc/tnrs/context.h"


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

using tax_pred_t = std::function<bool(const Taxon*)>;
json get_taxon_json(const RichTaxonomy& taxonomy, const Taxon& taxon);

enum match_status {unmatched=0,
                   ambiguous_match=1,
                   unambiguous_match=2};

// Should probably rename ContextSearcher => Context, and Context=> ContextDescription
struct ContextSearcher {
    const RichTaxonomy& taxonomy;
    const Context& context;
    const Taxon* context_root;

    pair<json,match_status> match_name(const string & query, bool do_approximate_matching, bool include_suppressed);

    ContextSearcher(const RichTaxonomy& t, const Context& c): taxonomy(t), context(c) {
        context_root = taxonomy.included_taxon_from_id(c.ott_id);
    }
    private:
    json fuzzy_name_match_json(const string& query, 
                               const FuzzyQueryResultWithTaxon & fqrwt) const;
    json exact_name_match_json(const string& query,
                               const Taxon* taxon) const;
    json exact_synonym_match_json(const string& query,
                                  const Taxon* taxon,
                                    const string& synonym_name) const;
    json _base_name_match_json(const string& query, const Taxon * taxon) const;
    json _base_name_match_json(const string & query, const TaxonomyRecord * record) const; 
};

string escape_query_string(const string& name) {
    return name;
}

// It seems that infer_context is supposed to ignore synonyms:
// * Bacteria is NOT ambiguous: it has 1 name match, and 1 synonym match.
// * Firmiscala IS ambiguous: it has 0 name matches, and 1 synonym match.
// * Random gibberish is reported as "ambiguous".
// curl -X POST https://api.opentreeoflife.org/v2/tnrs/infer_context -H "content-type:application/json" -d  '{"names":["Bacteria","Firmiscala"]}'
// 


bool taxon_is_specific(const Taxon* taxon) {
    auto rank = taxon->get_data().rank;
    return rank_is_specific(rank);
}

bool taxon_is_genus(const Taxon* taxon) {
    return taxon->get_data().rank == TaxonomicRank::RANK_GENUS;
}

bool taxon_is_higher(const Taxon* taxon) {
    return taxon->get_data().rank < TaxonomicRank::RANK_SPECIES;
}

using vec_tax_str_pair_t = vector<pair<const Taxon*, const string&> >;
vec_tax_str_pair_t exact_synonym_search(const RichTaxonomy& taxonomy,
                                        const Taxon* context_root,
                                        string query, 
                                        tax_pred_t ok = [](const Taxon*){return true;})
{
    query = normalize_query(query);
    vec_tax_str_pair_t hits;
    for(auto taxon: iter_post_n_const(*context_root)) {
        if (not ok(taxon)) {
            continue;
        }
        for(auto& tjs: taxon->get_data().junior_synonyms) {
            if (lcase_string_equals(query, tjs->get_name())) {
                hits.push_back({taxon,tjs->get_name()});
            }
        }
    }


    if (auto ctp = taxonomy.get_fuzzy_matcher())
    {
        auto fuzzy_results = ctp->to_taxa(ctp->exact_query(query), context_root, taxonomy, true);
        vec_tax_str_pair_t hits2;
        for(auto& result: fuzzy_results)
        {
            if (result.is_synonym())
            {
                auto t = result.get_taxon();
                if (ok(t))
                    hits2.push_back({t,query});
            }
        }
//        we can't sort references -- use string_view?
        vector<const Taxon*> taxon_hits1;
        for(auto& [taxon,_]: hits)
            taxon_hits1.push_back(taxon);
        vector<const Taxon*> taxon_hits2;
        for(auto& [taxon,_]: hits2)
            taxon_hits2.push_back(taxon);
        std::sort(taxon_hits1.begin(), taxon_hits1.end());
        std::sort(taxon_hits2.begin(), taxon_hits2.end());
        LOG(INFO)<<"exact_synonym_search: query = '"<<query<<"'  context_id = "<<context_root->get_ott_id();
        if (taxon_hits1 != taxon_hits2)
        {
            LOG(INFO)<<"lcase match:";
            for(int i=0;i<taxon_hits1.size();i++)
                LOG(INFO)<<"   "<<taxon_hits1[i]->get_data().get_nonuniqname();
            LOG(INFO)<<"ctrie match:";
            for(int i=0;i<taxon_hits2.size();i++)
                LOG(INFO)<<"   "<<taxon_hits2[i]->get_data().get_nonuniqname();
        }
        else
            LOG(INFO)<<"exact synonym search: "<<hits.size()<<" names agree";
    }

    return hits;
}

vec_tax_str_pair_t exact_synonym_search(const RichTaxonomy& taxonomy,
                                        const Taxon* context_root,
                                        string query,
                                        bool include_suppressed) {
    if (include_suppressed) {
        return exact_synonym_search(taxonomy, context_root, query);
    }
    tax_pred_t ok = [&](const Taxon* taxon) {
        return not taxonomy.node_is_suppressed_from_tnrs(taxon);
    };
    return exact_synonym_search(taxonomy, context_root, query, ok);
}

vec_tax_str_pair_t exact_synonym_search_higher(const RichTaxonomy& taxonomy,
                                               const Taxon* context_root,
                                               string query,
                                               bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) {
            return false;
        }
        return taxon_is_higher(taxon);
    };
    return exact_synonym_search(taxonomy, context_root, query, ok);
}

vector<const Taxon*> exact_name_search_species(const RichTaxonomy& taxonomy,
                                              const Taxon* context_root,
                                              string query,
                                              bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) {
            return false;
        }
        return taxon_is_specific(taxon);
    };
    return exact_name_search(taxonomy, context_root, query, ok);
    
}

vector<const Taxon*> exact_name_search_genus(const RichTaxonomy& taxonomy,
                                             const Taxon* context_root,
                                             string query,
                                             bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) {
            return false;
        }
        return taxon_is_genus(taxon);
    };
    return exact_name_search(taxonomy, context_root, query, ok);
}

vector<const Taxon*> exact_name_search_higher(const RichTaxonomy& taxonomy,
                                              const Taxon* context_root,
                                              string query,
                                              bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) {
            return false;
        }
        return taxon_is_higher(taxon);
    };
    return exact_name_search(taxonomy, context_root, query, ok);
}

vector<const Taxon*> prefix_name_search(const Taxon* context_root,
                                        const string& query,
                                        tax_pred_t ok = [](const Taxon*){return true;}) {
    vector<const Taxon*> hits;
    for(auto taxon: iter_post_n_const(*context_root)) {
        if (not ok(taxon)) {
            continue;
        }
        if (lcase_match_prefix(taxon->get_data().get_nonuniqname(), query)) {
            hits.push_back(taxon);
        }
    }
    return hits;
}

vector<const Taxon*> prefix_name_search(const RichTaxonomy& taxonomy,
                                        const Taxon* context_root,
                                        const string& query,
                                        bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        return include_suppressed or not taxonomy.node_is_suppressed_from_tnrs(taxon);
    };
    return prefix_name_search(context_root, query, ok);
}

vector<const Taxon*> prefix_name_search_higher(const RichTaxonomy& taxonomy,
                                               const Taxon* context_root,
                                               const string& query,
                                               bool include_suppressed) {
    tax_pred_t ok = [&](const Taxon* taxon) {
        if (not include_suppressed and taxonomy.node_is_suppressed_from_tnrs(taxon)) {
            return false;
        }
        return taxon_is_higher(taxon);
    };
    return prefix_name_search(context_root, query, ok);
}

vec_tax_str_pair_t prefix_synonym_search(const Taxon* context_root,
                                         string query,
                                         tax_pred_t ok = [](const Taxon*){return true;})
{
    vec_tax_str_pair_t hits;
    for(auto taxon: iter_post_n_const(*context_root)) {
        if (not ok(taxon)) {
            continue;
        }
        for(auto& tjs: taxon->get_data().junior_synonyms) {
            if (lcase_match_prefix(tjs->get_name(), query)) {
                hits.push_back({taxon,tjs->get_name()});
            }
        }
    }
    return hits;
}

vec_tax_str_pair_t prefix_synonym_search(const RichTaxonomy& taxonomy,
                                         const Taxon* context_root,
                                         string query,
                                         bool include_suppressed) {
    std::function<bool(const Taxon*)> ok = [&](const Taxon* taxon) {
        return include_suppressed or not taxonomy.node_is_suppressed_from_tnrs(taxon);
    };
    return prefix_synonym_search(context_root, query, ok);
}

inline json get_taxon_json(const RichTaxonomy& taxonomy, const Taxon& taxon) {
    json j;
    tax_service_add_taxon_info(taxonomy, taxon, j);
    // What about the "is_suppressed_from_synth" flag?  Do we want that?
    return j;
}
 
inline json get_taxon_record_json(const RichTaxonomy& taxonomy, const TaxonomyRecord & taxon) {
    json j;
    tax_service_add_suppressed_taxon_info(taxonomy, taxon, j);
    return j;
}


inline json ContextSearcher::_base_name_match_json(const string& query,
                                                   const Taxon * taxon) const {
    json result;
    result["taxon"] = get_taxon_json(taxonomy, *taxon);
    result["search_string"] = query;
    result["nomenclature_code"] = Context::get_code_name(taxonomy, taxon);
    return result;
}

inline json ContextSearcher::_base_name_match_json(const string & query, 
                                                   const TaxonomyRecord * record) const {
    json result;
    result["taxon"] = get_taxon_record_json(taxonomy, *record);
    result["search_string"] = query;
    result["nomenclature_code"] = Context::get_code_name(taxonomy, record);
    return result;
}

inline json ContextSearcher::fuzzy_name_match_json(const string& query, 
                           const FuzzyQueryResultWithTaxon & fqrwt) const {
    json result;
    auto taxon = fqrwt.get_taxon();
    if (taxon != nullptr) {
        result = _base_name_match_json(query, taxon);
    } else {
        result = _base_name_match_json(query, taxon);
    }
    result["score"] = fqrwt.get_score();
    result["is_approximate_match"] = true;
    result["is_synonym"] = false;
    result["matched_name"] = fqrwt.get_matched_name();
    return result;                           
}

inline json ContextSearcher::exact_name_match_json(const string& query,
                                            const Taxon* taxon) const {
    auto result = _base_name_match_json(query, taxon);
    result["score"] = 1.0;
    result["is_approximate_match"] = false;
    result["is_synonym"] = false;
    result["matched_name"] = taxon->get_data().get_nonuniqname();
    return result;
}

inline json ContextSearcher::exact_synonym_match_json(const string& query,
                              const Taxon* taxon,
                              const string& synonym_name) const {
    auto result = _base_name_match_json(query, taxon);
    result["score"] = 1.0;
    result["is_approximate_match"] = false;
    result["is_synonym"] = true;
    result["matched_name"] = synonym_name;
    return result;
}


pair<json,match_status> ContextSearcher::match_name(const string & raw_query,
                                                    bool do_approximate_matching,
                                                    bool include_suppressed) {
    auto query = normalize_query(raw_query);
    json results;
    match_status status = unmatched;
    // 1. See if we can find an exact name match
    auto exact_name_matches = exact_name_search(taxonomy, context_root, query, include_suppressed);
    for(auto taxon: exact_name_matches) {
        results.push_back(exact_name_match_json(query, taxon));
    }
    if (exact_name_matches.size() == 1) {
        status = unambiguous_match;
    }
    // 2. See if we can find an exact name match for synonyms
    auto exact_synonym_matches = exact_synonym_search(taxonomy, context_root, query, include_suppressed);
    for(auto& [ taxon, synonym_name ]: exact_synonym_matches) {
        results.push_back(exact_synonym_match_json(query, taxon, synonym_name));
    }
    if (status == unmatched and results.size()) {
        status = ambiguous_match;
    }
    // 3. Do fuzzy matching ONLY for names that we couldn't match
    if (do_approximate_matching and status == unmatched) {
        // do fuzzy matching.
        auto ctp = taxonomy.get_fuzzy_matcher();
        if (ctp == nullptr) {
            throw OTCError() << "Fuzzy matching has not been enabled in the taxonomy, but was requested in match_name.";
        }
        auto fuzzy_results = ctp->fuzzy_query_to_taxa(query, context_root, taxonomy, include_suppressed);
        if (fuzzy_results.size() > 0) {
            if (fuzzy_results.size() == 1) {
                status = unambiguous_match;
            } else {
                status = ambiguous_match;
            }
            for (auto fqr : fuzzy_results) {
                results.push_back(fuzzy_name_match_json(query, fqr));
            }
        }
    }
    json match_results;
    match_results["name"] = raw_query;
    match_results["matches"] = results;
    return {match_results, status};
}


const Context* determine_context_for_names(const vector<string>& names,
                                           const optional<string>& context_name,
                                           const RichTaxonomy& taxonomy) {
    if (not context_name) {
        return infer_context_and_ambiguous_names(taxonomy, names).first;
    }
    return get_context_by_name(*context_name);
}



// return name.split("\s+",2) if the string has a space or an optional with no value otherwise.
optional<pair<string,string>> split_genus_species(const string& name) {
    auto first_space = name.find(' ');
    // Quit if there's no space
    if (first_space == string::npos) {
        return {};
    }
    auto genus = name.substr(0,first_space);
    auto non_space = name.find_first_not_of(' ', first_space+1);
    if (non_space == string::npos) {
        non_space = name.size();
    }
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


//FIXME: how is "suppressed_names" different from "deprecated_taxa"?

// $ curl -X POST https://api.opentreeoflife.org/v3/tnrs/match_names  -H "content-type:application/json" -d '{"names":["Aster","Symphyotrichum","Barnadesia"]}'
// $ curl -X POST http://localhost:1984/v3/tnrs/match_names  -H "content-type:application/json" -d '{"names":["Aster","Symphyotrichum","Barnadesia"]}'
std::string tnrs_match_names_ws_method(const vector<string>& names,
                                       const optional<string>& context_name,
                                       bool do_approximate_matching,
                                       const optional<vector<string>>& /* ids */, // FIXME: Do we need to implement this?
                                       bool include_suppressed,
                                       const RichTaxonomy& taxonomy) {
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
    for(auto& name: names) {
        // Do the search
        auto [result, status] = searcher.match_name(name, do_approximate_matching, include_suppressed);
        // Store the result
        results.push_back(result);
        // Classify name as unmatched / matched / unambiguous
        if (status == unmatched) {
            unmatched_names.push_back(name);
        } else {
            matched_names.push_back(name);
            if (status == unambiguous_match) {
                unambiguous_names.push_back(name);
            }
        }
    }
    // 3. Construct JSON response.
    json response;
    response["governing_code"] = context->code.name;
    response["context"] = context->name;
    response["includes_approximate_matches"] = do_approximate_matching;
    response["includes_deprecated_taxa"] = false; // ?? How is this different from suppressed_names?
    response["includes_suppressed_names"] = include_suppressed;
    response["taxonomy"] = tax_about_json(taxonomy);
    response["unambiguous_names"] = unambiguous_names;
    response["unmatched_names"] = unmatched_names;
    response["matched_names"] = matched_names;
    response["results"] = results;
    return response.dump(1);
}

json autocomplete_json(const RichTaxonomy& taxonomy, const Taxon* taxon) {
    json match;
    match["ott_id"] = taxon->get_ott_id();
    match["unique_name"] = get_taxon_unique_name(*taxon);
    match["is_suppressed"] = taxonomy.node_is_suppressed_from_tnrs(taxon);
    match["is_higher"] = taxon_is_higher(taxon);
    return match;
}    

inline json autocomplete_json(const RichTaxonomy& taxonomy, const pair<const Taxon*,const string&>& p) {
    return autocomplete_json(taxonomy, p.first);
}

inline void add_hits(json& j, const RichTaxonomy& taxonomy, const vector<const Taxon*> taxa) {
    for(auto taxon: taxa) {
        j.push_back(autocomplete_json(taxonomy, taxon));
    }
}

inline void add_hits(json& j, const RichTaxonomy& taxonomy, const vector<FuzzyQueryResultWithTaxon>& fuzzy_query_results)
{
    for(auto fuzzy_query_result: fuzzy_query_results) {
        j.push_back(autocomplete_json(taxonomy, fuzzy_query_result.get_taxon()));
    }
}

inline void add_hits(json& j, const RichTaxonomy& taxonomy, const vec_tax_str_pair_t taxa) {
    for(auto [taxon, synonym]: taxa) {
        j.push_back(autocomplete_json(taxonomy, taxon));
    }
}

// Find all species in the genus that have the given prefix
vector<const Taxon*> prefix_search_species_in_genus(const Taxon* genus,
                                                    const string_view& species_prefix) {
    vector<const Taxon*> match_species;
    auto genus_name = genus->get_data().get_nonuniqname();
    for(auto species: iter_post_n_const(*genus)) {
        if (not taxon_is_specific(species)) {
            continue;
        }
        auto species_name = species->get_data().get_nonuniqname().substr(genus_name.size()+1);
        if (lcase_match_prefix(species_name, species_prefix)) {
            match_species.push_back(species);
        }
    }
    return match_species;
}

/*
 * Fuzzy matching DOES occur in autocomplete_name:
 *
 * curl -X POST https://api.opentreeoflife.org/v3/tnrs/autocomplete_name -H "content-type:application/json" -d '{"name":"Homo salien"}'
[ {
  "is_suppressed" : false,
  "unique_name" : "Homo sapiens",
  "ott_id" : 770315,
  "is_higher" : false
} ]
 *
 * However, fuzzy matching does not find `Homo sapiens neanderthalensis` and `Homo sapiens subsp. 'Denisova'`, whereas direct matching does.
 *
 * curl -X POST https://api.opentreeoflife.org/v3/tnrs/autocomplete_name -H "content-type:application/json" -d '{"name":"Homo sapien"}'
 *
 * However, it appears that a difference of 2 chars is allowed on an exact match, instead of doing a prefix-query on the fuzzy index:
 *
 * "Hono saliens" -> "Homo sapiens" + "Neobodo saliens"
 * "Hono salien"  -> nothing
 * "Homo salens"  -> "Homo sapiens"
 * "Homo salen"   -> nothing
 */



// curl -X POST https://api.opentreeoflife.org/v3/tnrs/autocomplete_name -H "content-type:application/json" -d '{"name":"Endoxyla","context_name":"All life"}'
string tnrs_autocomplete_name_ws_method(const string& name,
                                        const string& context_name,
                                        bool include_suppressed,
                                        const RichTaxonomy& taxonomy) {
    // An empty response should be `[ ]`, not `null`.
    json response = json::array();
    // We need to escape the query string.
    auto escaped_query = normalize_query(name);
    // This corresponds to a SingleNamePrefixQuery in taxomachine.
    // * See org/opentree/taxonomy/plugins/tnrs_v3.java
    // * See org/opentree/tnrs/queries/SingleNamePrefixQuery.java
    // 0. Escape the query??
    // lower-case the name?
    // 1. Determine context
    auto context = determine_context(context_name);
    auto context_root = taxonomy.included_taxon_from_id(context->ott_id);
    if (auto query_genus_species = split_genus_species(name)) {
        // 2. If we have a space, then assume the first part is a genus and match species names within the genus
        // Search against species and synonyms
        add_hits(response, taxonomy, exact_name_search_species(taxonomy, context_root, escaped_query, include_suppressed));
        add_hits(response, taxonomy, exact_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
        if (response.size()) {
            return response.dump(1);
        }
        // no exact hit against the species index
        auto genus_hits = exact_name_search_genus(taxonomy, context_root, escaped_query, include_suppressed);
        if (not genus_hits.empty()) { // the first word was an exact match against the genus index
            auto [query_genus,query_species] = *query_genus_species;
            for(auto genus: genus_hits) {
                add_hits(response, taxonomy, prefix_search_species_in_genus(genus, query_species));
            }
        }
        if (not response.empty()) {
            return response.dump(1);
        }
        // no exact hit for first word against the genus index
        // Hit query string against the higher taxon index... not sure if this is useful, since it has a space
        add_hits(response, taxonomy, exact_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
        if (not response.empty()) {
            return response.dump(1);
        }
        // Prefix query against the synonyms and higher taxa
        add_hits(response, taxonomy, prefix_name_search(taxonomy, context_root, escaped_query, include_suppressed));
        add_hits(response, taxonomy, prefix_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
        if (not response.empty()) {
            return response.dump(1);
        }
        // fuzzy search on names and synonyms (BDR -- not a prefix search?)
        {
            auto ctp = taxonomy.get_fuzzy_matcher();
            if (ctp == nullptr) {
                throw OTCError() << "Fuzzy matching has not been enabled in the taxonomy, but was requested in match_name.";
            }
            auto fuzzy_results = ctp->fuzzy_query_to_taxa(escaped_query, context_root, taxonomy, include_suppressed);
            add_hits(response, taxonomy, fuzzy_results);
        }
    } else { // does not contain a space at all
        add_hits(response, taxonomy, exact_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
        add_hits(response, taxonomy, exact_synonym_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
        if (not response.empty()) {
            return response.dump(1);
        }
        // Do a prefix query against the higher taxon index
        add_hits(response, taxonomy, prefix_name_search_higher(taxonomy, context_root, escaped_query, include_suppressed));
        if (not response.empty()) {
            return response.dump(1);
        }
        // Do a prefix query against the all taxa synonym index
        add_hits(response, taxonomy, prefix_synonym_search(taxonomy, context_root, escaped_query, include_suppressed));
        if (not response.empty()) {
            return response.dump(1);
        }
        // fuzzy search on HIGHER names and synonyms (BDR -- not a prefix search?)
        {
            auto ctp = taxonomy.get_fuzzy_matcher();
            if (ctp == nullptr) {
                throw OTCError() << "Fuzzy matching has not been enabled in the taxonomy, but was requested in match_name.";
            }
            auto fuzzy_results = ctp->fuzzy_query_to_taxa(escaped_query, context_root, taxonomy, include_suppressed);
            add_hits(response, taxonomy, fuzzy_results);
        }
    }
    return response.dump(1);
}

// curl -X POST https://api.opentreeoflife.org/v3/tnrs/contexts
// curl -X POST http://localhost:1984/v3/tnrs/contexts
std::string tnrs_contexts_ws_method() {
    json response;
    for(auto& context: all_contexts) {
        response[context.group].push_back(context.name);
    }
    return response.dump(1);
}

// curl -X POST https://api.opentreeoflife.org/v3/tnrs/infer_context -H "content-type:application/json" -d '{"names":["Pan","Homo","Mus","Bufo","Drosophila"]}'
// curl -X POST http://localhost:1984/v3/tnrs/infer_context -H "content-type:application/json" -d '{"names":["Pan","Homo","Mus","Bufo","Drosophila"]}'
string tnrs_infer_context_ws_method(const vector<string>& names, const RichTaxonomy& taxonomy)
{
    auto [context, ambiguous_names] = infer_context_and_ambiguous_names(taxonomy, names);

    json response;
    response["context_name"] = context->name;
    response["context_ott_id"] = context->ott_id;
    response["ambiguous_names"] = ambiguous_names;

    return response.dump(1);
}

} // namespace otc
