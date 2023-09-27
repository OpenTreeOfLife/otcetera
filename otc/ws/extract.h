#ifndef EXTRACT_H
#define EXTRACT_H
#include <optional>
#include <vector>
#include "json.hpp"
#include <string>
#include "otc/ws/otc_web_error.h"
#include "otc/otc_base_includes.h"

template<typename T>
std::optional<T> convert_to(const nlohmann::json & j);

template <>
std::optional<bool> convert_to(const nlohmann::json & j) {
    if (j.is_boolean()) {
        return j.get<bool>();
    }
    return {};
}

template <>
std::optional<std::string> convert_to(const nlohmann::json & j) {
    if (j.is_string()) {
        return j.get<std::string>();
    }
    return {};
}

template <>
std::optional<int> convert_to(const nlohmann::json & j) {
    if (j.is_number()) {
        return j.get<int>();
    }
    return {};
}

template <>
std::optional<std::vector<std::string>> convert_to(const nlohmann::json & j) {
    if (not j.is_array()) {
        return {};
    }
    std::vector<std::string> v;
    for(auto& xj: j) {
        auto x = convert_to<std::string>(xj);
        if (not x) return {};
        v.push_back(*x);
    }
    return v;
}

#if defined(LONG_OTT_ID)
template <>
std::optional<otc::OttId> convert_to(const nlohmann::json & j) {
    return (j.is_number() ? j.get<OttId>() : {});
}
#endif

template <>
std::optional<otc::OttIdSet> convert_to(const nlohmann::json & j) {
    using namespace otc;
    if (not j.is_array()) {
        return {};
    }
    OttIdSet ids;
    for(auto& jid: j) {
        auto id = convert_to<OttId>(jid);
        if (not id) return {};
        ids.insert(*id);
    }
    return ids;
}



template <typename T>
constexpr const char* type_name_with_article();

template <> constexpr const char* type_name_with_article<bool>() {
    return "a boolean";
}
template <> constexpr const char* type_name_with_article<int>() {
    return "an integer";
}
template <> constexpr const char* type_name_with_article<std::string>() {
    return "a string";
}
#if defined(LONG_OTT_ID)
template <> constexpr const char* type_name_with_article<OttId>() {
    return "an OttId";
}
#endif
template <> constexpr const char* type_name_with_article<std::vector<std::string>>() {
    return "an array of strings";
}
template <> constexpr const char* type_name_with_article<otc::OttIdSet>() {
    return "an array of integers";
}

template<typename T>
std::optional<T> extract_argument(const nlohmann::json & j, const std::string& opt_name, bool required=false) {
    using namespace otc;
    auto opt = j.find(opt_name);
    if (opt == j.end()) {
        if (required) {
            throw OTCBadRequest("expecting ") << type_name_with_article<T>() << " argument called '" << opt_name << "'\n";
        }
        return {};
    }
    auto arg = convert_to<T>(*opt);
    if (not arg) {
        throw OTCBadRequest("expecting argument '") << opt_name << "' to be " << type_name_with_article<T>() <<"! Found '" << opt->dump() << "'\n";
    }
    return arg;
}

template<typename T>
T extract_argument_or_default(const nlohmann::json & j, const std::string& opt_name, const T& _default_) {
    auto arg = extract_argument<T>(j, opt_name);
    return (arg ? *arg : _default_);
}

template<typename T>
T extract_required_argument(const nlohmann::json & j, const std::string& opt_name) {
    auto arg = extract_argument<T>(j, opt_name, true);
    assert(arg);
    return *arg;
}


#endif
