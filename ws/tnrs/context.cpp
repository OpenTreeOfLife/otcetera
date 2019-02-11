#include "context.h"
#include "../tolws.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/otc_base_includes.h" // for OttID

using std::string;
using std::map;
using std::pair;
using std::vector;

namespace otc
{

// NOTE: These context definitions were taken from `taxomachine/src/main/java/org/opentree/taxonomy/contexts/ContextDescription.java`
//       on Jan 15, 2019.

const vector<Context> all_contexts = {
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
    {"Platyhelminthes", "ANIMALS",   "Platyhelminthes",  "Platyhelminthes",  555379L,    Nomenclature::ICZN},
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

map<OttId, const Context*> make_ottid_to_context(const vector<Context>& all_contexts)
{
    map<OttId, const Context*> ottid_to_context;
    for(auto& context: all_contexts)
	ottid_to_context.insert({context.ott_id, &context});
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

// FIXME - We should return Context* or Context&.
const Context* least_inclusive_context(const vector<const RTRichTaxNode*>& taxa)
{
    if (taxa.size() == 0)
	return name_to_context.at("All life");

    auto mrca = taxonomy_mrca(taxa);

    while(true)
    {
	auto id = mrca->get_ott_id();
	if (ottid_to_context.count(id))
	    return ottid_to_context.at(id);
	mrca = mrca->get_parent();
    }

    throw OTCError()<<"Can't find least inclusive context for "<<taxa.size()<<" taxa!";
}

// FIXME - technically we could try and save the exact matches to avoid work.
// FIXME - We should return Context* or Context&.
pair<const Context*,vector<string>> infer_context_and_ambiguous_names(const RichTaxonomy& taxonomy, const vector<string>& names)
{
    vector<const RTRichTaxNode*> unique_taxa;

    vector<string> ambiguous_names;

    for(auto& name: names)
    {
	// This search (i) includes suppressed/dubious/deprecated taxa and (ii) DOES (?NOT) include synonyms.
	auto hits = exact_name_search(taxonomy, name, true);
	if (hits.size() == 1)
	    unique_taxa.push_back(hits.front());
	else
	    ambiguous_names.push_back(name);
    }

    return {least_inclusive_context(unique_taxa), ambiguous_names};
}

}
