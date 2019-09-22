#include "nomenclature.h"

using std::string;

namespace Nomenclature
{
    Code::Code(const string& n, const string& d):name(n),description(d) {}
    Code ICN ("ICN",  "plants, fungi, and some protists?");
    Code ICNP("ICNP", "bacteria");
    Code ICZN("ICZN", "animals");
    Code Undefined ("undefined", "governing code unclear, nonexistent, or multiple codes");
}
