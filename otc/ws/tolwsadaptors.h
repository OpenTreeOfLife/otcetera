#ifndef OTC_TOLWS_ADAPTORS_H
#define OTC_TOLWS_ADAPTORS_H
#include <restbed>
#include "otc/ws/tolws.h"
#include "otc/ws/trees_to_serve.h"
#include "otc/ws/parallelreadserialwrite.h"
#include "otc/otc_base_includes.h"
#include <optional>

inline std::optional<nlohmann::json> parse_body(const restbed::Bytes& body) {
    if (body.empty()) {
        return nlohmann::json();
    }
    try {
        return nlohmann::json::parse(body);
    }
    catch (...) {
        return {};
    }
}

inline std::optional<nlohmann::json> lookup(const nlohmann::json& j, const std::string& s) {
    auto x = j.find(s);
    if (x == j.end()) {
        return {};
    }
    return {*x};
}

inline std::optional<nlohmann::json> parse_body(const unsigned char* body) {
    if (not body) {
        return nlohmann::json();
    }
    try {
        return nlohmann::json::parse(body);
    }
    catch (...) {
        return {};
    }
}

inline nlohmann::json parse_body_or_throw(const restbed::Bytes& body) {
    auto oj = parse_body(body);
    if (not oj) {
        throw otc::OTCBadRequest("Could not parse body of call as JSON.\n");
    }
    return *oj;
}

inline nlohmann::json parse_body_or_throw(const unsigned char* body) {
    auto oj = parse_body(body);
    if (not oj) {
        throw otc::OTCBadRequest("Could not parse body of call as JSON.\n");
    }
    return *oj;
}

#endif
