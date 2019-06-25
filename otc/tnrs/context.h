#ifndef CONTEXT_H
#define CONTEXT_H

#include <string>
#include <utility>
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
    std::vector<std::size_t> subcontext_indices;
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
};

extern const std::vector<Context> all_contexts;

extern std::map<OttId, const Context*> ottid_to_context;

extern std::map<std::string, const Context*> name_to_context;

std::pair<const Context*,std::vector<std::string>> infer_context_and_ambiguous_names(const RichTaxonomy& taxonomy, const std::vector<std::string>& names);

const Context* least_inclusive_context(const std::vector<RTRichTaxNode*>& taxa);
}
#endif
