#ifndef OTCETERA_UTIL_H
#define OTCETERA_UTIL_H

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include "otc/otc_base_includes.h"
#include "otc/error.h"

namespace otc {

const std::string readStrContentOfUTF8File(const std::string &filepath);
const std::wstring readWStrContentOfUTF8File(const std::string &filepath);
bool openUTF8File(const std::string &filepath, std::ifstream & inp);
bool openUTF8WideFile(const std::string &filepath, std::wifstream & inp);


template<typename T>
bool isProperSubset(const T & small, const T & big);
template<typename T>
std::set<T> set_difference_as_set(const std::set<T> & small, const std::set<T> & big);
template<typename T, typename U>
bool contains(const T & container, const U & key);
template<typename T, typename U>
std::set<T> keys(const std::map<T, U> & container);
void writeOttSet(std::ostream & out, const char *indent, const std::set<long> &fir, const char * sep);
void writeOttSetDiff(std::ostream & out, const char *indent, const std::set<long> &fir, const char *firN, const std::set<long> & sec, const char *secN);

inline const std::wstring readWStrContentOfUTF8File(const std::string &filepath) {
	std::wifstream inp;
	if (!openUTF8WideFile(filepath, inp)) {
		throw OTCError("Could not open \"" + filepath + "\"");
	}
	const std::wstring utf8content((std::istreambuf_iterator<wchar_t>(inp) ), (std::istreambuf_iterator<wchar_t>()));
	return utf8content;
}

template<typename T, typename U>
bool contains(const T & container, const U & key) {
	return container.find(key) != container.end();
}

template<typename T, typename U>
std::set<T> keys(const std::map<T, U> & container) {
	std::set<T> k;
	for (auto x : container) {
		k.insert(x.first);
	}
	return k;
}


inline void writeOttSet(std::ostream & out,
						const char *indent,
						const std::set<long> &fir,
						const char * sep) {
	for (auto rIt = fir.begin(); rIt != fir.end(); ++rIt) {
		if (rIt != fir.begin()) {
			out << sep;
		}
		out << indent << "ott" << *rIt;
	}
}

inline void writeOttSetDiff(std::ostream & out,
							const char *indent,
							const std::set<long> &fir,
							const char *firN,
							const std::set<long> & sec,
							const char *secN) {
	for (auto rIt : fir) {
		if (sec.find(rIt) == sec.end()) {
			out << indent << "ott" << rIt << " is in " << firN << " but not " << secN << "\n";
		}
	}
	for (auto rIt : sec) {
		if (fir.find(rIt) == fir.end()) {
			out << indent << "ott" << rIt << " is in " << secN << " but not " << firN << "\n";
		}
	}
}

template<typename T>
inline bool isProperSubset(const T & small, const T & big) {
	if (big.size() <= small.size()) {
		return false;
	}
	for (auto rIt : small) {
		if (big.find(rIt) == big.end()) {
			return false;
		}
	}
	return true;
}

template<typename T>
std::set<T> set_difference_as_set(const std::set<T> & fir, const std::set<T> & sec) {
	std::set<T> d;
	set_difference(begin(fir), end(fir), begin(sec), end(sec), std::inserter(d, d.end()));
	return d;
}

} //namespace otc
#endif
