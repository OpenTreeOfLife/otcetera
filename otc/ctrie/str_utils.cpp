
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

} // namespace otc
