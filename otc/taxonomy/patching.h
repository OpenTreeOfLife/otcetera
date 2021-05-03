#ifndef OTC_TAXONOMY_PATCHING_H
#define OTC_TAXONOMY_PATCHING_H

#include "otc/taxonomy/taxonomy.h"

namespace otc {

class PatchableTaxonomy: public RichTaxonomy {
    public:
    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    PatchableTaxonomy(const std::string& dir,
                      std::bitset<32> cf = std::bitset<32>(),
                      OttId kr = -1);
    PatchableTaxonomy(PatchableTaxonomy &&) = default;

};

PatchableTaxonomy load_patchable_taxonomy(const boost::program_options::variables_map& args);


} // namespace
#endif
