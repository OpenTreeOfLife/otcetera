#ifndef OTC_TAXONOMY_DIFF_MAKER_H
#define OTC_TAXONOMY_DIFF_MAKER_H

#include "otc/taxonomy/taxonomy.h"

namespace otc {

class TaxonomyDiffMaker: public RichTaxonomy {
    public:
    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    TaxonomyDiffMaker(const std::string& dir,
                      std::bitset<32> cf = std::bitset<32>(),
                      OttId kr = -1);
    TaxonomyDiffMaker(TaxonomyDiffMaker &&) = default;
    protected:
    std::map<std::string, std::vector<const RTRichTaxNode *> > synonym2node;
    
};

TaxonomyDiffMaker load_taxonomy_diff_maker(const boost::program_options::variables_map& args);


} // namespace
#endif
