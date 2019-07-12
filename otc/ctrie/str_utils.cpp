
#include <locale>
#include "otc/ctrie/str_utils.h"

namespace otc {
const std::ctype<char> * glob_facet;
std::wstring_convert<deletable_facet<std::codecvt<char32_t, char, std::mbstate_t> >, char32_t> glob_conv32;
std::wstring_convert<std::codecvt_utf8_utf16<char32_t>, char32_t> glob_conv8;
std::locale global_locale;

int set_global_conv_facet() {
    try {
        global_locale = std::locale("en_US.utf8");
    } catch (const std::exception &) {
        try {
            global_locale = std::locale("en");
        } catch (const std::exception & x) {
             std::cerr << "locale \"en_US.utf8\" or \"en\" must be supported on the system to run this program.\n";
             std::cerr << x.what() << "\n";
             return 1;
        }
    }
    auto & f = std::use_facet<std::ctype<char> >(global_locale);
    glob_facet = &f;
    return 0;
}


const char PUNC_CHAR = '?';
const char UNKNOWN_CHAR = '@'; // some char that is not in any encoded indexed string
const char BRACE_CHAR = '|';
std::vector<char> ascii_char_to_norm;
std::map<char32_t, char> wide_to_norm;

std::string normalize_wide_query(const std::u32string & raw_query);


std::string normalize_wide_query(const std::u32string & wide_query) {
    std::string norm;
    norm.reserve(wide_query.size());
    std::cerr << "trying to normalize_query \"" << to_char_str(wide_query) << "\" ...";
    
    for (const auto & c: wide_query) {
        if (c < 127) {
            char tc = (char) c;
            std::cerr << "\n  norm ascii  to \"" << ascii_char_to_norm[tc] << "\" \n";
                
            norm.push_back(ascii_char_to_norm[tc]);
        } else {
            auto wcit = wide_to_norm.find(c);
            if (wcit == wide_to_norm.end()) {
                std::cerr << "\n  norm nonascii  to \"" << UNKNOWN_CHAR << "\" \n";
                norm.push_back(UNKNOWN_CHAR);
            } else {
                std::cerr << "\n  norm nonascii  to \"" << wcit->second << "\" \n";
                norm.push_back(wcit->second);
            }
        }
    }
    std::cerr << " converted  to \"" << norm << "\" \n";
    return norm;
}

// could use this for \"e  -> e as well
std::string normalize_query(const std::string & raw_query) {
    std::string norm;
    norm.reserve(raw_query.size());
    for (auto& c: raw_query) {
        if (c > 0 && c < 127) {
            norm.push_back(ascii_char_to_norm[c]);
        } else {
            auto wrq = to_u32string(raw_query);
            norm = normalize_wide_query(wrq);
            break;
        }
    }
    std::cerr << "normalized \"" << raw_query << "\" to \"" << norm << "\" \n";
    return norm;
}

std::string normalize_query(const std::string_view & raw_query) {
    std::string norm;
    norm.reserve(raw_query.size());
    for (auto& c: raw_query) {
        if (c > 0 && c < 127) {
            norm.push_back(ascii_char_to_norm[c]);
        } else {
            auto wrq = to_u32string(raw_query);
            norm = normalize_wide_query(wrq);
            break;
        }
    }
    std::cerr << "normalized \"" << raw_query << "\" to \"" << norm << "\" \n";
    return norm;
}

void init_char_maps(){
    set_global_conv_facet();
    std::vector<char> & c = ascii_char_to_norm;
    std::map<char32_t, char> & m = wide_to_norm;
    c.assign(128, PUNC_CHAR);
    const char SPACE_CHAR = ' ';
    const char SEP_CHAR = '-';
    c[' '] = SPACE_CHAR;
    c['-'] = SEP_CHAR;
    c[':'] = SEP_CHAR;
    c['_'] = SEP_CHAR;
    c['.'] = '.';
    c['!'] = PUNC_CHAR;
    c['#'] = PUNC_CHAR;
    c['%'] = PUNC_CHAR;
    c['&'] = PUNC_CHAR;
    c['*'] = PUNC_CHAR;
    c['+'] = PUNC_CHAR;
    c[','] = PUNC_CHAR;
    c['/'] = PUNC_CHAR;
    c[';'] = PUNC_CHAR;
    c['='] = PUNC_CHAR;
    c['?'] = PUNC_CHAR;
    c['@'] = PUNC_CHAR;
    c['^'] = PUNC_CHAR;
    c['{'] = BRACE_CHAR;
    c['}'] = BRACE_CHAR;
    c['\"'] = BRACE_CHAR;
    c['\''] = BRACE_CHAR;
    c['('] = BRACE_CHAR;
    c[')'] = BRACE_CHAR;
    c['<'] = BRACE_CHAR;
    c['>'] = BRACE_CHAR;
    c['['] = BRACE_CHAR;
    c[']'] = BRACE_CHAR;
    c['0'] = '0';
    c['1'] = '1';
    c['2'] = '2';
    c['3'] = '3';
    c['4'] = '4';
    c['5'] = '5';
    c['6'] = '6';
    c['7'] = '7';
    c['8'] = '8';
    c['9'] = '9';
    c['a'] = 'a';
    c['b'] = 'b';
    c['c'] = 'c';
    c['d'] = 'd';
    c['e'] = 'e';
    c['f'] = 'f';
    c['g'] = 'g';
    c['h'] = 'h';
    c['i'] = 'i';
    c['j'] = 'j';
    c['k'] = 'k';
    c['l'] = 'l';
    c['m'] = 'm';
    c['n'] = 'n';
    c['o'] = 'o';
    c['p'] = 'p';
    c['q'] = 'q';
    c['r'] = 'r';
    c['s'] = 's';
    c['t'] = 't';
    c['u'] = 'u';
    c['v'] = 'v';
    c['w'] = 'w';
    c['x'] = 'x';
    c['y'] = 'y';
    c['z'] = 'z';
    c['A'] = 'a';
    c['B'] = 'b';
    c['C'] = 'c';
    c['D'] = 'd';
    c['E'] = 'e';
    c['F'] = 'f';
    c['G'] = 'g';
    c['H'] = 'h';
    c['I'] = 'i';
    c['J'] = 'j';
    c['K'] = 'k';
    c['L'] = 'l';
    c['M'] = 'm';
    c['N'] = 'n';
    c['O'] = 'o';
    c['P'] = 'p';
    c['Q'] = 'q';
    c['R'] = 'r';
    c['S'] = 's';
    c['T'] = 't';
    c['U'] = 'u';
    c['V'] = 'v';
    c['W'] = 'w';
    c['X'] = 'x';
    c['Y'] = 'y';
    c['Z'] = 'z';

    m[L'\xa0'] = SPACE_CHAR;        //  160 ' '
    m[L'\xab'] = BRACE_CHAR; //  171 «
    m[L'\xc3'] = 'a';        //  195 Ã
    m[L'\xc4'] = 'a';  //  196 Ä
    m[L'\xc6'] = 'a';  //  198 Æ
    m[L'\xd6'] = 'o';  //  214 Ö
    m[L'\xd7'] = '*';  //  215 × maps to * which is not punc char so we have an ascii for hybrid
    m[L'\xdf'] = 's';  //  223 ß
    m[L'\xe0'] = 'a';  //  224 à
    m[L'\xe1'] = 'a';  //  225 á
    m[L'\xe3'] = 'a';  //  227 ã
    m[L'\xe4'] = 'a';  //  228 ä
    m[L'\xe5'] = 'a';  //  229 å
    m[L'\xe6'] = 'a';  //  230 æ
    m[L'\xe7'] = 'c';  //  231 ç
    m[L'\xe8'] = 'e';  //  232 è
    m[L'\xe9'] = 'e';  //  233 é
    m[L'\xea'] = 'e';  //  234 ê
    m[L'\xeb'] = 'e';  //  235 ë
    m[L'\xed'] = 'i';  //  237 í
    m[L'\xee'] = 'i';  //  238 î
    m[L'\xef'] = 'i';  //  239 ï
    m[L'\xf1'] = 'n';  //  241 ñ
    m[L'\xf3'] = 'o';  //  243 ó
    m[L'\xf4'] = 'o';  //  244 ô
    m[L'\xf6'] = 'o';  //  246 ö
    m[L'\xf8'] = 'o';  //  248 ø
    m[L'\xfa'] = 'u';  //  250 ú
    m[L'\xfc'] = 'u';  //  252 ü
    m[L'\xfd'] = 'y';  //  253 ý
    m[L'\xff'] = 'y';  //  255 ÿ
    m[L'\u0152'] = 'o';  //  338 Œ
    m[L'\u0153'] = 'o';  //  339 œ
    m[L'\u0160'] = 's';  //  352 Š
    m[L'\u0161'] = 's';  //  353 š
    m[L'\u017e'] = 'z';  //  382 ž
    m[L'\u017f'] = 's';  //  383 ſ
    m[L'\u03b1'] = 'a';  //  945 α
    m[L'\u03b2'] = 'b';  //  946 β
    m[L'\u03b3'] = 'g';  //  947 γ
    m[L'\u03b4'] = 'd';  //  948 δ
    m[L'\u03b5'] = 'e';  //  949 ε
    m[L'\u03b7'] = 'e';  //  951 η
    m[L'\u03b8'] = 't';  //  952 θ
    m[L'\u03bb'] = 'l';  //  955 λ
    m[L'\u03bc'] = 'm';  //  956 μ
    m[L'\u03ca'] = 'i';  //  970 ϊ
    m[L'\u2009'] = SPACE_CHAR;  //  8201 " "
    m[L'\u2019'] = BRACE_CHAR;  //  8217 ’
    m[L'\ufb02'] = 'f';  //  64258 ﬂ
    for (auto x : m) {
        std::cerr << "\"" << to_char_str(x.first) << "\" -> " << x.second << "\"\n"; 
    }
}

} // namespace otc
