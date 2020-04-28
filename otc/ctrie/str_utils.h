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
    assert(glob_facet != nullptr);
    std::string undecoded{uncap_mod};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    std::u32string ret = glob_conv32.from_bytes(undecoded.data(), undecoded.data() + undecoded.length());
    return ret;
}

inline std::string lower_case_version(const std::string & arg) {
    assert(glob_facet != nullptr);
    std::string undecoded{arg};
    glob_facet->tolower(&undecoded[0], &undecoded[0] + undecoded.size());
    return undecoded;
}

inline std::string upper_case_version(const std::string & arg) {
    assert(glob_facet != nullptr);
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

// BDR: This does not handle all accented characters, but should get a fair number of them.
//      See https://en.wikipedia.org/wiki/Latin-1_Supplement_(Unicode_block)
//      See https://en.wikipedia.org/wiki/ISO/IEC_8859-1
//      See https://stackoverflow.com/questions/14094621/
inline unsigned char normalize_latin_char(unsigned char ch)
{
    static const char*
        //   "ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõö÷øùúûüýþÿ"
        tr = "AAAAAAECEEEEIIIIDNOOOOOx0UUUUYPsaaaaaaeceeeeiiiiOnooooo/0uuuuypy";

    if ( ch < 128)
        return ch;
    else if ( ch < 192)
    {
        if (ch == 171) return '"'; // << quote
        if (ch == 187) return '"'; // >> quote
        return ch;
    }

    return tr[ ch-192 ];
}

inline char32_t normalize_greek_or_coptic_uchar(char32_t ch)
{
    // See https://unicode.org/charts/PDF/U0370.pdf

    if (ch == 945)  return 'a'; // alpha
    if (ch == 946)  return 'b'; // beta
    if (ch == 947)  return 'g'; // gamma
    if (ch == 948)  return 'd'; // delta
    if (ch == 949)  return 'e'; // epsilon
    if (ch == 950)  return 'z'; // zeta
    if (ch == 951)  return 'e'; // eta
    if (ch == 952)  return 't'; // theta
    if (ch == 953)  return 'i'; // iota
    if (ch == 954)  return 'k'; // kappa
    if (ch == 955)  return 'l'; // lambda 
    if (ch == 956)  return 'm'; // mu

    return ch;
}

// https://www.codetable.net/decimal/8217
inline char32_t normalize_uchar(char32_t ch)
{
    if (ch < 256)
        return normalize_latin_char(ch);

    if (ch >= 880 and ch < 1024)
        return normalize_greek_or_coptic_uchar(ch);

    if (ch == 352)  return  'S'; // Š Latin capital S with Caron
    if (ch == 382)  return  'z'; // ž
    if (ch == 353)  return  's'; // ſ

    if (ch == 1086) return 'o';
    if (ch == 1089) return 'c';

    if (ch == 8201) return ' '; // thin space
    if (ch == 8220) return '"'; // left quote
    if (ch == 8221) return '"'; // right quote

    return ch;
}

// Currently we don't allow changing the number of chars, so we can't handle
// ligatures like fl or ae.
inline std::string normalize_query(const std::string_view & raw_query)
{
    auto uquery = to_u32string(raw_query);
    for (auto& c: uquery) {
        c = normalize_uchar(c);
        c = std::tolower(c);
    }
    return to_char_str(uquery);
}

} // namespace otc
#endif
