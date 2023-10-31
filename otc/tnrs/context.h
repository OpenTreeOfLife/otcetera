#ifndef CONTEXT_H
#define CONTEXT_H

#include <string>
#include <utility>
#include <optional>
#include <map>
#include "otc/taxonomy/taxonomy.h"
#include "nomenclature.h"

namespace otc
{

class CompressedTrieBasedDB;

struct Context
{
    std::string name;
    std::string group;
    std::string name_suffix;
    std::string lica_node_name;
    OttId ott_id;
    const Nomenclature::Code & code;
    mutable const CompressedTrieBasedDB * name_matcher; // for fuzzy matching

    Context(std::string name_arg,
            std::string group_arg,
            std::string name_suffix_arg,
            std::string lica_node_name_arg,
            OttId ott_arg,
            const Nomenclature::Code &code_arg)
     :name(name_arg),
     group(group_arg),
     name_suffix(name_suffix_arg),
     lica_node_name(lica_node_name_arg),
     ott_id(ott_arg),
     code(code_arg), 
     name_matcher(nullptr)
     {}

    static void init_nom_codes_boundaries(const RichTaxonomy &);
    static const std::string & get_code_name(const RichTaxonomy & taxonomy, const RTRichTaxNode * taxon);
    static const std::string & get_code_name(const RichTaxonomy & taxonomy, const TaxonomyRecord * record);
    // MUST CALL this FIRST. IMPORTANT GLOBAL SIDE EFFECTS
    static int cull_contexts_to_taxonomy(const RichTaxonomy &); 
};

extern std::vector<Context> all_contexts;
extern std::map<std::string, const Context*> name_to_context;


const Context * determine_context(const std::optional<std::string> & context_name);
const Context * get_context_by_name(const std::string & context_name);

inline const Context* get_context_by_name(const std::string & context_name) {
    if (not name_to_context.count(context_name)) {
        throw OTCError() << "The context '" << context_name << "' could not be found.";
    }
    return name_to_context.at(context_name);
}

inline const Context* determine_context(const std::optional<std::string> & context_name) {
    if (not context_name) {
        return name_to_context.at("All life");
    }
    return get_context_by_name(*context_name);
}

std::pair<const Context*,std::vector<std::string>> infer_context_and_ambiguous_names(const RichTaxonomy& taxonomy, const std::vector<std::string>& names);

const Context* least_inclusive_context(const std::vector<RTRichTaxNode*>& taxa);
}
#endif
