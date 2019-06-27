#ifndef NOMENCLATURE_H
#define NOMENCLATURE_H

#include <string>

namespace Nomenclature
{
    struct Code
    {
        std::string name;
        std::string description;
        Code(const std::string&, const std::string&);
    };
    
    extern Code ICN;
    extern Code ICNP;
    extern Code ICZN;
    extern Code Undefined;
}

#endif
