#ifndef CONTEXT_H
#define CONTEXT_H

#include <string>
#include <utility>
#include <map>
#include "otc/taxonomy/taxonomy.h"
#include "nomenclature.h"

namespace otc
{

struct Context
{
    std::string name;
    std::string group;
    std::string name_suffix;
    std::string lica_node_name;
    OttId ott_id;
    Nomenclature::Code code;
};

extern const std::vector<Context> all_contexts;

extern std::map<OttId, const Context*> ottid_to_context;

extern std::map<std::string, const Context*> name_to_context;

std::pair<const Context*,std::vector<std::string>> infer_context_and_ambiguous_names(const RichTaxonomy& taxonomy, const std::vector<std::string>& names);

const Context* least_inclusive_context(const std::vector<RTRichTaxNode*>& taxa);
}
#endif
