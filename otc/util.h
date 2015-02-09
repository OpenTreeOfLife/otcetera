#ifndef OTCETERA_UTIL_H
#define OTCETERA_UTIL_H

#include <iostream>
#include <string>
#include "otc/otcetera.h"
#include "otc/error.h"

namespace otc {

const std::string readStrContentOfUTF8File(const std::string &filepath);
const std::wstring readWStrContentOfUTF8File(const std::string &filepath);
bool openUTF8File(const std::string &filepath, std::wifstream & inp);

inline const std::wstring readWStrContentOfUTF8File(const std::string &filepath) {
	std::wifstream inp;
	if (!openUTF8File(filepath, inp)) {
		throw OTCError("Could not open \"" + filepath + "\"");
	}
	const std::wstring utf8content((std::istreambuf_iterator<wchar_t>(inp) ), (std::istreambuf_iterator<wchar_t>()));
	return utf8content;
}

} //namespace otc
#endif
