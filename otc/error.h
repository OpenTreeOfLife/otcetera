#ifndef OTCETERA_ERROR_H
#define OTCETERA_ERROR_H
#include <string>
#include <exception>
#include <sstream>

namespace otc {

class OTCError : public std::exception {
    protected:
    std::string message;
    public:
    const char * what() const noexcept {
        return message.c_str();
    }
    template <typename T> OTCError& operator<<(const T&);
    void prepend(const std::string& s) {
        message = s + message;
    }
    OTCError() noexcept {}
    OTCError(const std::string & msg) noexcept :message(msg) {}
};

template <typename T>
OTCError& OTCError::operator<<(const T& t) {
  std::ostringstream oss;
  oss << message << t;
  message = oss.str();
  return *this;
}

} //namespace otc
#endif
