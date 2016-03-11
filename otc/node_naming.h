#ifndef OTCETERA_NODE_NAMING_H
#define OTCETERA_NODE_NAMING_H

#include <string>
#include "otc/otc_base_includes.h"

namespace otc {

std::string makeName(const std::string& prefix, int number);
std::string makeMRCAName(int number1, int number2);

inline std::string makeName(const std::string& pre, int number) {
    return pre + std::to_string(number);
}

inline std::string makeMRCAName(int number1, int number2) {
    const static std::string mrca_prefix = "mrcaott";
    return mrca_prefix + std::to_string(number1) + "ott" + std::to_string(number2);
}

} //namespace otc
#endif
