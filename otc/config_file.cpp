#include "config_file.h"
#include <regex>
#include "error.h"
#include <fstream>
#include <filesystem>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/functional/hash.hpp>

namespace fs = std::filesystem;

using std::string;
using std::size_t;
using std::vector;
using std::optional;
using boost::property_tree::ptree;

namespace otc {

optional<string> interpolate(const ptree& pt, const string& section_name, const string& key) {
    static std::regex KEYCRE ("%\\(([^)]+)\\)s");
    if (not pt.get_child_optional(section_name)) {
        return {};
    }
    const ptree& section = pt.get_child(section_name);
    if (not section.get_optional<string>(key)) {
        return {};
    }
    string value = section.get<string>(key);
    for(int i = 0; i < 20 and value.find('%') != string::npos; i++) {
        string value2;
        size_t p1 = 0U;
        size_t p2 = 0U;
        while (p1 < value.size()) {
            p2 = value.find('%', p1);
            if (p2 == string::npos) {
                p2 = value.size(); 
            }
            value2 += value.substr(p1, p2-p1);
            if (p2 == value.size()) {
                break;
            }
            if (p2+1 >= value.size()) {
                throw OTCError() << "Found '%' at end of string!";
            }
            char c = value[p2 + 1];
            if (c == '%') {
                value2 += "%";
                p1 = p2 + 1;
            } else if (c == '(') {
                std::cmatch m;
                bool matched = std::regex_search(value.c_str()+p2, value.c_str()+value.size(), m, KEYCRE);
                if (not matched) {
                    throw OTCError() << "Bad interpolation variable reference: '" << value.substr(p2) << "'";
                }
                string name = m[1];
                if (not section.get_optional<string>(name)) {
                    throw OTCError() << "Reference to undefined key '" << name << "' in " << key << " = " << value << "\n";
                }
                value2 += section.get<string>(name);
                p1 = p2 + m.length(0);
            } else {
                throw OTCError() << "Bad interpolation variable reference: '" << value.substr(p2) << "'";
            }
        }
        value = value2;
    }
    return value;
}


optional<string> load_config(const string& filename, const string& section, const string& name) {
    ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    return interpolate(pt, section, name);
}

optional<string> load_config(const vector<string>& filenames, const string& section, const string& name) {
    for(const auto& filename: filenames) {
        auto result = load_config(filename, section, name);
        if (result) {
            return result;
        }
    }
    return {};
}

optional<string> dot_opentree() {
    auto envconfig = std::getenv("OTC_CONFIG");
    if (envconfig) {
        std::ifstream test(envconfig);
        if (not test) {
            return {};
        }
        return envconfig;
    }
    auto homedir = std::getenv("HOME");
    if (not homedir) {
        return {};
    }
    string path = homedir;
    path += "/.opentree";
    std::ifstream test(path);
    if (not test) {
        return {};
    }
    test.close();
    return path;
}

} // namespace otc
