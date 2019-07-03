#include "otc/supertree_util.h"
#include <regex>

using std::optional;
using std::string;

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

optional<string> get_source_node_name(const std::string& name)
{
    // Here the middle group can match ott12345, so we add ottXXX on the end to avoid this.
    static std::regex with_ott(".*[ _]([a-zA-Z]+\\d+)[ _]ott.*");
    // Then we need another regex to handle the case where there's no ottid.
    static std::regex without_ott(".*[ _]([a-zA-Z]+\\d+)[ _]?");
    std::smatch matches;
    if (std::regex_match(name, matches, with_ott))
    {
        assert(matches.size() >= 2);
        std::string source = matches[1];
        return source;
    }

    if (std::regex_match(name, matches, without_ott)) {
        assert(matches.size() >= 2);
        std::string source = matches[1];
        return source;
    }
    else {
        return {};
    }
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

