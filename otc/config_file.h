#ifndef OTC_CONFIG_H
#define OTC_CONFIG_H

#include <string>
#include <vector>
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>

namespace otc {
boost::optional<std::string> load_config(const std::string& filename, const std::string& section, const std::string& name);
boost::optional<std::string> load_config(const std::vector<std::string>& filenames, const std::string& section, const std::string& name);
boost::optional<std::string> dot_opentree();
boost::optional<std::string> interpolate(const boost::property_tree::ptree& pt, const std::string& section_name, const std::string& key);
}
#endif
