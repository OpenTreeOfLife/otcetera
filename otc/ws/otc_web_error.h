#ifndef OTC_WEB_ERROR_H
#define OTC_WEB_ERROR_H

#include <string>
#include <stdexcept>
#include "otc/error.h"
#include "json.hpp"

namespace otc {


class OTCWebError : public std::exception
{
protected:
    int status_code_ = 500;
    nlohmann::json data;
public:
    int status_code() const {return status_code_;}

    const char * what() const noexcept {
        return data["message"].get<std::string>().c_str();
    }

    template <typename T> OTCWebError& operator<<(const T&);

    void prepend(const std::string& s);

    void append(const std::string& s);

    nlohmann::json& json() {return data;}

    OTCWebError() noexcept;
    OTCWebError(int c) noexcept;
    OTCWebError(const std::string & msg) noexcept;
    OTCWebError(int c, const std::string & msg) noexcept;
};

template <typename T>
OTCWebError& OTCWebError::operator<<(const T& t)
{
    std::ostringstream oss;
    oss << t;
    append(oss.str());
    return *this;
}

template <>
OTCWebError& OTCWebError::operator<<(const nlohmann::json& j);

template <>
OTCWebError& OTCWebError::operator<<(const std::string& j);

inline OTCWebError OTCBadRequest() {return OTCWebError(400);}
inline OTCWebError OTCBadRequest(const std::string& m) {
    return OTCWebError(400,m);
}


inline void OTCWebError::prepend(const std::string& s) {
    if (not data.count("message")) {
        data["message"] = s;
    } else {
        std::string tmp = data["message"];
        data["message"] = s + tmp;
    }
}

inline void OTCWebError::append(const std::string& s) {
    if (not data.count("message")) {
        data["message"] = s;
    } else {
        std::string tmp = data["message"];
        data["message"] = tmp + s;
    }
}

inline OTCWebError::OTCWebError() noexcept
    : OTCWebError(500,"") {
}

inline OTCWebError::OTCWebError(int c) noexcept
    :OTCWebError(c, "") {
}


inline OTCWebError::OTCWebError(const std::string & msg) noexcept
    :OTCWebError(500, msg) {
}

inline OTCWebError::OTCWebError(int c, const std::string & msg) noexcept
    :status_code_(c) {
    prepend(msg);
}

template<>
inline OTCWebError& OTCWebError::operator<<(const nlohmann::json& j) {
    for (auto & [key,value]: j.items()) {
        data[key] = value;
    }
    return *this;
}

template <>
inline OTCWebError& OTCWebError::operator<<(const std::string& j) {
    append(j);
    return *this;
}


} // namespace otc
#endif
