#ifndef OTCETERA_ERROR_H
#define OTCETERA_ERROR_H
#include <string>
#include <exception>
#include "otc/otc_base_includes.h"
namespace otc {

struct OTCError : public std::exception {
	OTCError(const std::string & msg)
		:message(msg) {
	}
	const char * what () const noexcept {
		return message.c_str();
	}
	protected:
		std::string message;
};

} //namespace otc
#endif
