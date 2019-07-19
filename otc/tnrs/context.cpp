#include "context.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/otc_base_includes.h" // for OttID

using std::string;
using std::map;
using std::pair;
using std::vector;

namespace otc
{

int find_context_index_by_name(const vector<Context> & context_vec, const std::string &n) {
    for (auto i = 0; i < (int) context_vec.size(); ++i) {
        if (context_vec[i].name == n) {
            return i;
        }
    }
    return -1;
}

// NOTE: These context definitions were taken from `taxomachine/src/main/java/org/opentree/taxonomy/contexts/ContextDescription.java`
//       on Jan 15, 2019.

vector<Context> generate_all_contexts() {
    vector<Context> ret = {
        // Name             Group        Index suffix        Node name string    ott id      Code
        {"All life",        "LIFE",      "",                 "life",             805080L,    Nomenclature::Undefined},
        // MICROBES group
        {"Bacteria",        "MICROBES",  "Bacteria",         "Bacteria",         844192L,    Nomenclature::ICNP},
        {"SAR group",       "MICROBES",  "SAR",              "SAR",              5246039L,   Nomenclature::Undefined},
        {"Archaea",         "MICROBES",  "Archaea",          "Archaea",          996421L,    Nomenclature::ICNP},
        {"Excavata",        "MICROBES",  "Excavata",         "Excavata",         2927065L,   Nomenclature::Undefined},
        {"Amoebozoa",       "MICROBES",  "Amoebae",          "Amoebozoa",        1064655L,   Nomenclature::ICZN},
        {"Centrohelida",    "MICROBES",  "Centrohelida",     "Centrohelida",     755852L,    Nomenclature::ICZN},
        {"Haptophyta",      "MICROBES",  "Haptophyta",       "Haptophyta",       151014L,    Nomenclature::Undefined},
        {"Apusozoa",        "MICROBES",  "Apusozoa",         "Apusozoa",         671092L,    Nomenclature::ICZN},
        {"Diatoms",         "MICROBES",  "Diatoms",          "Bacillariophyta",  5342311L,   Nomenclature::ICN},
        {"Ciliates",        "MICROBES",  "Ciliates",         "Ciliophora",       302424L,    Nomenclature::Undefined},
        {"Forams",          "MICROBES",  "Forams",           "Foraminifera",     936399L,    Nomenclature::ICZN},
        // ANIMALS group
        {"Animals",         "ANIMALS",   "Animals",          "Metazoa",          691846L,    Nomenclature::ICZN},
        {"Birds",           "ANIMALS",   "Birds",            "Aves",             81461L,     Nomenclature::ICZN},
        {"Tetrapods",       "ANIMALS",   "Tetrapods",        "Tetrapoda",        229562L,    Nomenclature::ICZN},
        {"Mammals",         "ANIMALS",   "Mammals",          "Mammalia",         244265L,    Nomenclature::ICZN},
        {"Amphibians",      "ANIMALS",   "Amphibians",       "Amphibia",         544595L,    Nomenclature::ICZN},
        {"Vertebrates",     "ANIMALS",   "Vertebrates",      "Vertebrata",       801601L,    Nomenclature::ICZN},
        {"Arthropods",      "ANIMALS",   "Arthopods",        "Arthropoda",       632179L,    Nomenclature::ICZN},
        {"Molluscs",        "ANIMALS",   "Molluscs",         "Mollusca",         802117L,    Nomenclature::ICZN},
        {"Nematodes",       "ANIMALS",   "Nematodes",        "Nematoda",         395057L,    Nomenclature::ICZN},
        {"Platyhelminthes", "ANIMALS",   "PlatyContexthelminthes",  "Platyhelminthes",  555379L,    Nomenclature::ICZN},
        {"Annelids",        "ANIMALS",   "Annelids",         "Annelida",         941620L,    Nomenclature::ICZN},
        {"Cnidarians",      "ANIMALS",   "Cnidarians",       "Cnidaria",         641033L,    Nomenclature::ICZN},
        {"Arachnids",       "ANIMALS",   "Arachnids",        "Arachnida",        511967L,    Nomenclature::ICZN},
        {"Insects",         "ANIMALS",   "Insects",          "Insecta",          1062253L,   Nomenclature::ICZN},
        // FUNGI group
        {"Fungi",           "FUNGI",     "Fungi",            "Fungi",            352914L,    Nomenclature::ICN},
        {"Basidiomycetes",  "FUNGI",     "Basidiomycetes",   "Basidiomycota",    634628L,    Nomenclature::ICN},
        {"Ascomycetes",     "FUNGI",     "Ascomycota",       "Ascomycota",       439373L,    Nomenclature::ICN},
        // PLANTS group
        {"Land plants",     "PLANTS",    "Plants",           "Embryophyta",      5342313L,   Nomenclature::ICN},
        {"Hornworts",       "PLANTS",    "Anthocerotophyta", "Anthocerotophyta", 738980L,    Nomenclature::ICN},
        {"Mosses",          "PLANTS",    "Bryophyta",        "Bryophyta",        246594L,    Nomenclature::ICN},
        {"Liverworts",      "PLANTS",    "Marchantiophyta",  "Marchantiophyta",  56601L,     Nomenclature::ICN},
        {"Vascular plants", "PLANTS",    "Tracheophyta",     "Tracheophyta",     10210L,     Nomenclature::ICN},
        {"Club mosses",     "PLANTS",    "Lycopodiophyta",   "Lycopodiophyta",   144803L,    Nomenclature::ICN},
        {"Ferns",           "PLANTS",    "Moniliformopses",  "Moniliformopses",  166292L,    Nomenclature::ICN},
        {"Seed plants",     "PLANTS",    "Spermatophyta",    "Spermatophyta",    10218L,     Nomenclature::ICN},
        {"Flowering plants","PLANTS",    "Magnoliophyta",    "Magnoliophyta",    99252L,     Nomenclature::ICN},
        {"Monocots",        "PLANTS",    "Monocots",         "Liliopsida",       1058517L,   Nomenclature::ICN},
        {"Eudicots",        "PLANTS",    "Eudicots",         "eudicotyledons",   431495L,    Nomenclature::ICN},
        {"Rosids",          "PLANTS",    "Rosids",           "rosids",           1008296L,   Nomenclature::ICN},
        {"Asterids",        "PLANTS",    "Asterids",         "asterids",         1008294L,   Nomenclature::ICN},
        {"Asterales",       "PLANTS",    "Asterales",        "Asterales",        1042120L,   Nomenclature::ICN},
        {"Asteraceae",      "PLANTS",    "Asteraceae",       "Asteraceae",       46248L,     Nomenclature::ICN},    
        {"Aster",           "PLANTS",    "Aster",            "Aster",            409712L,    Nomenclature::ICN},    
        {"Symphyotrichum",  "PLANTS",    "Symphyotrichum",   "Symphyotrichum",   1058735L,   Nomenclature::ICN},
        {"Campanulaceae",   "PLANTS",    "Campanulaceae",    "Campanulaceae",    1086303L,   Nomenclature::ICN},
        {"Lobelia",         "PLANTS",    "Lobelia",          "Lobelia",          1086294L,   Nomenclature::ICN}
    };

    const std::map<std::string, std::vector<std::string> > parname2childnames = {
        {"All life", {"Bacteria", "Archaea", "SAR group", "Excavata", "Amoebozoa", "Centrohelida", "Haptophyta", "Apusozoa", "Animals", "Fungi", "Land plants"}},
        {"SAR group", {"Diatoms", "Ciliates", "Forams"}},
        {"Animals", {"Vertebrates", "Arthropods", "Molluscs", "Nematodes", "Platyhelminthes", "Annelids", "Cnidarians"}},
        {"Vertebrates", {"Tetrapods"}},
        {"Tetrapods", {"Birds", "Mammals", "Amphibians"}},
        {"Arthropods", {"Arachnids","Insects"}},
        {"Fungi", {"Basidiomycetes", "Ascomycetes"}},
        {"Land plants", {"Hornworts", "Mosses", "Liverworts", "Club mosses", "Vascular plants"}},
        {"Vascular plants", {"Ferns", "Seed plants"}},
        {"Seed plants", {"Monocots", "Eudicots"}},
        {"Eudicots", {"Rosids", "Asterids"}},
        {"Asterids", {"Asterales"}},
        {"Asterales", {"Asteraceae", "Campanulaceae"}},
        {"Asteraceae", {"Aster", "Symphyotrichum"}},
        {"Campanulaceae", {"Lobelia"}}
        };

    for (auto p2cit : parname2childnames) {
        const auto &  par_name = p2cit.first;
        auto par_ind = find_context_index_by_name(ret, par_name);
        if (par_ind < 0) {
            throw OTCError() << "did not find context with name \"" << par_name << "\"";
        }
        const auto & cvec = p2cit.second;
        Context & par_context = ret[par_ind];
        par_context.subcontext_indices.reserve(cvec.size());
        for (const auto & cn : cvec) {
            auto child_ind = find_context_index_by_name(ret, cn);
            if (child_ind < 0) {
                throw OTCError() << "did not find context with name \"" << cn << "\"";
            }
            par_context.subcontext_indices.push_back((std::size_t)child_ind);
        }
    }
    return ret;
}

const vector<Context> all_contexts = generate_all_contexts();

map<OttId, const Context*> make_ottid_to_context(const vector<Context> & all_contexts)
{
    map<OttId, const Context*> ottid_to_context;
    for(auto& context: all_contexts) {
        ottid_to_context.insert({context.ott_id, &context});
    }
    return ottid_to_context;
}

map<OttId, const Context*> ottid_to_context = make_ottid_to_context(all_contexts);



map<string, const Context*> make_name_to_context(const vector<Context>& all_contexts)
{
    map<string, const Context*> name_to_context;
    for(auto& context: all_contexts)
	name_to_context.insert({context.name, &context});
    return name_to_context;
}

map<string, const Context*> name_to_context = make_name_to_context(all_contexts);


using trav_range_t = std::pair<std::uint32_t, std::uint32_t>;
using nom_cod_and_ranges_t = std::pair<const Nomenclature::Code *, std::vector<trav_range_t> >;
using vec_nom_cod_and_ranges_t = std::vector<nom_cod_and_ranges_t>;
vec_nom_cod_and_ranges_t global_nom_code_and_ranges;


using id_name_pair_t = std::pair<OttId, std::string>;
using vec_id_name_pair_t = std::vector<id_name_pair_t>;

std::vector<trav_range_t> _fill_and_sort_ranges(const RichTaxonomy & taxonomy, const vec_id_name_pair_t & vec_ids_names) {
    const auto & tree = taxonomy.get_tax_tree();
    const auto & tax_tree_data = tree.get_data();
    std::vector<trav_range_t> ret;
    for (auto & p : vec_ids_names) {
        const RTRichTaxNode * nd;
        try {
             nd = tax_tree_data.id_to_node.at(p.first);
        } catch (...) {
            throw OTCError() << "Barrier taxon with id=" << p.first << " not found";
        }
        if (nd->get_name() != p.second) {
            throw OTCError() << "Barrier taxon with id=" << p.first << " had name \"" << nd->get_name() << "\", but we expected \"" << p.second << "\".";
        }
        const auto & ndd = nd->get_data();
        ret.push_back(trav_range_t{ndd.trav_enter, ndd.trav_exit});
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}


/* From BarrierNodes.java in taxomachine...
ICN
352914 "Fungi"
361838 "Viridiplantae"
266751 "Alveolata"
878953 "Rhodophyta"
664970 "Glaucocystophyceae"
151014 "Haptophyceae"

Nomenclature.ICNP
844192 "Bacteria"
996421 "Archaea"

ICZN
{691846, "Metazoa"}, {202765, "Choanoflagellida"}
 */

void Context::init_nom_codes_to_traversal(const RichTaxonomy & taxonomy) {
    global_nom_code_and_ranges.clear();
    if (taxonomy.get_tax_tree().get_root()->get_ott_id() != 805080) {
        LOG(WARNING) << "taxonomy root did not have ID=805080, assuming that this is a taxonomy for testing purposes only... detection of nomenclatural code will be non-functional.";
        return;
    }
    const vec_id_name_pair_t iczn{{691846, "Metazoa"},
                                  {202765, "Choanoflagellida"}};
    const vec_id_name_pair_t icnp{{844192, "Bacteria"},
                                  {996421, "Archaea (domain silva:D37982/#1)"}};
    const vec_id_name_pair_t icn{{352914, "Fungi"},
                                 {361838, "Chloroplastida"},
                                 {266751, "Alveolata"},
                                 {878953, "Rhodophyta"},
                                 {664970, "Glaucophyta"},
                                 {151014, "Haptophyta"}};
    const auto viczn = _fill_and_sort_ranges(taxonomy, iczn);
    const auto vinp = _fill_and_sort_ranges(taxonomy, icnp);
    const auto vicn = _fill_and_sort_ranges(taxonomy, icn);
    global_nom_code_and_ranges.push_back(nom_cod_and_ranges_t{&Nomenclature::ICZN, viczn});
    global_nom_code_and_ranges.push_back(nom_cod_and_ranges_t{&Nomenclature::ICNP, vinp});
    global_nom_code_and_ranges.push_back(nom_cod_and_ranges_t{&Nomenclature::ICN, vicn});
}

const vec_nom_cod_and_ranges_t & get_nom_code_to_trav_ranges() {
    return global_nom_code_and_ranges;
}

inline bool contains_trav(const std::vector<trav_range_t> & sorted_trav_pairs, std::uint32_t query) {
    for (const auto & p : sorted_trav_pairs) {
        if (p.first > query) {
            return false;
        }
        if (p.second >= query) {
            return true;
        }
    }
    return false;
}

const Nomenclature::Code & get_nom_code_for_trav(std::uint32_t trav) {
    const auto & nom_code_ranges = get_nom_code_to_trav_ranges();
    for (const auto & c_and_r : nom_code_ranges) {
        if (contains_trav(c_and_r.second, trav)) {
            return *(c_and_r.first);
        }
    }
    return Nomenclature::Undefined;
}

const std::string & Context::get_code_name(const RichTaxonomy & , const RTRichTaxNode * taxon) {
    const auto & nom = get_nom_code_for_trav(taxon->get_data().trav_enter);
    return nom.name;
}

const RTRichTaxNode * get_closest_anc_with_node(const RichTaxonomy & taxonomy, const TaxonomyRecord * record) {
    const auto & tax_data = taxonomy.get_tax_tree().get_data();
    OttId par_id = record->parent_id;
    for (;;) {
        auto ndit = tax_data.id_to_node.find(par_id);
        if (ndit == tax_data.id_to_node.end()) {
            auto taxit = tax_data.id_to_record.find(par_id);
            if (taxit == tax_data.id_to_record.end()) {
                return nullptr;
            } else {
                par_id = taxit->second->parent_id;
            }
        } else {
            return ndit->second;
        }
    }
}

const std::string & Context::get_code_name(const RichTaxonomy & taxonomy, const TaxonomyRecord * record) {
    const auto taxon = get_closest_anc_with_node(taxonomy, record);
    if (taxon == nullptr) {
        return Nomenclature::Undefined.name;
    }
    return get_code_name(taxonomy, taxon);
}

// FIXME - We should return Context* or Context&.
const Context* least_inclusive_context(const vector<const RTRichTaxNode*> & taxa)
{
    if (taxa.size() == 0) {
        return name_to_context.at("All life");
    }
    auto mrca = taxonomy_mrca(taxa);
    while (mrca) {
        auto id = mrca->get_ott_id();
        if (ottid_to_context.count(id)) {
            return ottid_to_context.at(id);
        }
        mrca = mrca->get_parent();
    }
    throw OTCError() << "Can't find least inclusive context for " << taxa.size() << " taxa!";
}

// FIXME - technically we could try and save the exact matches to avoid work.
// FIXME - We should return Context* or Context&.
pair<const Context*, vector<string>> infer_context_and_ambiguous_names(const RichTaxonomy& taxonomy,
                                                                      const vector<string>& names) {
    vector<const RTRichTaxNode*> unique_taxa;
    vector<string> ambiguous_names;
    for(auto& name: names) {
        // This search (i) includes suppressed/dubious/deprecated taxa and (ii) DOES (?NOT) include synonyms.
        auto hits = taxonomy.exact_name_search(name, nullptr, true);
        if (hits.size() == 1) {
            unique_taxa.push_back(hits.front());
        } else {
            ambiguous_names.push_back(name);
        }
    }
    return {least_inclusive_context(unique_taxa), ambiguous_names};
}


}
