#ifndef OTC_WS_NEXSON_NEXSON_H
#define OTC_WS_NEXSON_NEXSON_H
#endif

#include <string>

#include "json.hpp"

namespace otc
{
    nlohmann::json get_phylesystem_study(const std::string& study_id);
    std::pair<nlohmann::json,nlohmann::json> extract_tree_nexson(const nlohmann::json& nexson, const std::string& treeid);

} // namespace otc
