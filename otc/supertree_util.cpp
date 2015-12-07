#include "otc/supertree_util.h"
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
    const char* start = name.c_str();
    const char* end = name.c_str() + name.size();
    while(strchr(" \t_",*start) and start < end) {
        start++;
    }
    while(strchr(" \t_",*(end-1)) and start < end) {
        end--;
    }
    if (start >= end) {
        throw OTCError() << "Node name '" << name << "' contracted to nothing!";
    }
    const std::size_t offset = static_cast<std::size_t>(start - name.c_str());
    const std::size_t len = static_cast<std::size_t>(end - start);
    return name.substr(offset, len);
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

