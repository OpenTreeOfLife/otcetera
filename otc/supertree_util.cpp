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

std::string getSourceNodeName(const std::string& name) {
    std::regex e("(.*[ _])?(node\\d+)([ _].*)?");
    std::smatch matches;
    if (std::regex_match(name,matches,e))
    {
        assert(matches.size() >= 2);
        return matches[2];
    }
    else
        throw OTCError()<<"'"<<name<<"' does not contain a node name!";
}

bool culledAndCompleteIncompatWRTLeafSet(const OttIdSet & culled,
                                                const OttIdSet & complete,
                                                const OttIdSet & leafSet) {
    //TMP this could be more efficient. See areCompatibleDesIdSets
    const OttIdSet inter = set_intersection_as_set(culled, complete);
    if (inter.empty()) {
        return false;
    }
    if (inter == culled) {
        return false;
    }
    const OttIdSet compCulled = set_intersection_as_set(complete, leafSet);
    return (inter != compCulled);
}

}// namespace otc

