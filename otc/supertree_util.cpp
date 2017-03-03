#include "otc/supertree_util.h"
#include <regex>
namespace otc {

std::string string_between_chars(const std::string & s, char beforeC, char endC) {
    auto start = s.find_last_of(beforeC);
    assert(start != std::string::npos);
    ++start;
    assert(start < s.length() - 1);
    const auto end = s.find_last_of(endC);
    assert(end != std::string::npos);
    assert(end > start);
    return s.substr(start, end - start);
}

std::string study_from_tree_name(const std::string& name) {
    return string_between_chars(name, ' ', '_');
}

std::string tree_in_study_from_tree_name(const std::string& name) {
    return string_between_chars(name, '_', '.');
}

std::string source_from_tree_name(const std::string& name) {
    return string_between_chars(name, ' ', '.');
}

boost::optional<std::string> get_source_node_name(const std::string& name) {
    static std::regex e("(.*[ _])?(node\\d+)([ _].*)?");
    std::smatch matches;
    if (std::regex_match(name,matches,e))
    {
        assert(matches.size() >= 2);
        std::string source = matches[2];
        return source;
    }
    else
        return boost::none;
}

bool culled_and_complete_incompat_wrt_leaf_set(const OttIdSet & culled,
                                                const OttIdSet & complete,
                                                const OttIdSet & leaf_set) {
    //TMP this could be more efficient. See are_compatible_des_id_sets
    const OttIdSet inter = set_intersection_as_set(culled, complete);
    if (inter.empty()) {
        return false;
    }
    if (inter == culled) {
        return false;
    }
    const OttIdSet compCulled = set_intersection_as_set(complete, leaf_set);
    return (inter != compCulled);
}

}// namespace otc

