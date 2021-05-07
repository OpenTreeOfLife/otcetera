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

    void write(const std::string& newdirname) const;
    void write_to_stream(std::ostream &) const;

    std::pair<bool, std::string> add_new_taxon(OttId oid,
                                               OttId parent_id,
                                               const std::string & name,
                                               const std::string & rank,
                                               const std::string & sourceinfo,
                                               const std::string & uniqname,
                                               const std::string & flags,
                                               OttId * homonym_of=nullptr) ;
    protected:
    void write_version_file_contents(std::ostream & out) const;
    void write_taxonomy_file_contents(std::ostream & tf) const;
    void write_synonyms_file_contents(std::ostream & sf) const;
    void write_forwards_file_contents(std::ostream & ff) const;

    //std::map<const RTRichTaxNode * , std::string> node_to_uniqname;
    std::map<std::string, std::vector<const RTRichTaxNode *> > synonym2node;

};

PatchableTaxonomy load_patchable_taxonomy(const boost::program_options::variables_map& args);


} // namespace
#endif
