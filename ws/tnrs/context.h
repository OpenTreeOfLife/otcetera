#ifndef CONTEXT_H
#define CONTEXT_H

#include <string>
#include <map>
#include "otc/otc_base_includes.h" // for OttId
#include "nomenclature.h"

struct Context
{
    std::string name;
    std::string name_suffix;
    std::string lica_node_name;
    otc::OttId ott_id;
    Nomenclature::Code code;
};

struct AllContexts: public std::map<std::string, Context>
{
    AllContexts(const std::initializer_list<Context>&);
};

extern AllContexts all_contexts;
#endif
