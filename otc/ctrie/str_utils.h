#ifndef OTC_CTRIE_STR_UTILS_H
#define OTC_CTRIE_STR_UTILS_H

#include <codecvt>
#include <locale>
#include "otc/error.h"
#include "otc/otc_base_includes.h"

namespace otc {

using stored_str_t = std::u32string;
using stored_char_t = char32_t;

// you MUST call set_global_conv_facet() before calling any of the functions that use glob_facet
int set_global_conv_facet();

std::string to_char_str(const stored_str_t & undecoded);
std::string to_char_str(const stored_char_t & undecoded);

// 
extern const std::ctype<char> * glob_facet;


//conversion from https://en.cppreference.com/w/cpp/locale/codecvt
// utility wrapper to adapt locale-bound facets for wstring/wbuffer convert
template <class Facet>
struct deletable_facet : Facet {
    template <class ...Args>
    deletable_facet(Args&& ...args) : Facet(std::forward<Args>(args)...) {}
    ~deletable_facet() {}
};

extern std::wstring_convert<deletable_facet<std::codecvt<char32_t, char, std::mbstate_t> >, char32_t> glob_conv32;
extern std::wstring_convert<std::codecvt_utf8_utf16<char32_t>, char32_t> glob_conv8;

inline std::u32string to_u32string(const std::string_view & undecoded) {
    return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
}

inline std::u32string to_u32string(const std::string & undecoded) {
    return glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
}


inline std::string to_char_str(const stored_char_t * undecoded) {
    return glob_conv8.to_bytes(undecoded);
}
 
inline std::string to_char_str(const stored_str_t & undecoded) {
    return glob_conv8.to_bytes(undecoded);
}

inline std::string to_char_str(const stored_char_t & undecoded) {
    return glob_conv8.to_bytes(undecoded);
}


inline std::u32string to_u32string_ci(const std::string_view & uncap_mod) {
    std::string undecoded{uncap_mod};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    std::u32string ret = glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    return ret;
}

inline std::string lower_case_version(const std::string & arg) {
    std::string undecoded{arg};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}

inline std::string upper_case_version(const std::string & arg) {
    std::string undecoded{arg};
    glob_facet->toupper(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}


inline bool starts_with(const stored_str_t & full, const stored_str_t & pref) {
    if (full.length() < pref.length()) {
        return false;
    }
    //std::cerr << "starts_with(" << to_char_str(full) << ", " << to_char_str(pref) << ")\n";
    return 0 == full.compare(0, pref.length(), pref);
}


} // namespace otc
#endif
