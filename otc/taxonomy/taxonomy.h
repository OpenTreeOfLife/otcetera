#ifndef OTC_TAXONOMY_TAXONOMY_H
#define OTC_TAXONOMY_TAXONOMY_H
// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/utility/string_ref.hpp>
#include <bitset>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

#include "json.hpp"

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?


namespace otc {

enum TaxonomicRank {
    RANK_DOMAIN,
    RANK_SUPERKINGDOM,
    RANK_KINGDOM,
    RANK_SUBKINGDOM,
    RANK_INFRAKINGDOM,
    RANK_SUPERPHYLUM,
    RANK_PHYLUM,
    RANK_DIVISION,
    RANK_SUBPHYLUM,
    RANK_SUBDIVISION,
    RANK_INFRAPHYLUM,
    RANK_SUPERCLASS,
    RANK_CLASS,
    RANK_SUBCLASS,
    RANK_INFRACLASS,
    RANK_SUPERORDER,
    RANK_ORDER,
    RANK_SUBORDER,
    RANK_INFRAORDER,
    RANK_PARVORDER,
    RANK_SUPERFAMILY,
    RANK_FAMILY,
    RANK_SUBFAMILY,
    RANK_SUPERTRIBE,
    RANK_TRIBE,
    RANK_SUBTRIBE,
    RANK_GENUS,
    RANK_SUBGENUS,
    RANK_SECTION,
    RANK_SUBSECTION,
    RANK_SPECIES_GROUP,
    RANK_SPECIES_SUBGROUP,
    RANK_SPECIES,
    RANK_SUBSPECIES,
    RANK_INFRASPECIFICNAME,
    RANK_FORMA,
    RANK_SUBFORM,
    RANK_VARIETAS,
    RANK_VARIETY,
    RANK_SUBVARIETY,
    RANK_NO_RANK,
    RANK_NO_RANK_TERMINAL
};

extern const std::map<TaxonomicRank, std::string> rank_enum_to_name;
extern const std::string empty_string;
extern const std::set<std::string> indexed_source_prefixes;

inline const std::string & rank_to_string(const TaxonomicRank &r) {
    auto i = rank_enum_to_name.find(r);
    if (i == rank_enum_to_name.end()) {
        return empty_string;
    }
    return i->second;
}



inline std::vector<std::string> comma_separated_as_vec(const std::string & sourceinfo) {
    std::vector<std::string> mt;
    if (sourceinfo.empty()) {
        return mt;
    }
    std::string sis = std::string(sourceinfo);
    auto aslist = split_string(sis, ',');
    mt.assign(aslist.begin(), aslist.end());
    return mt;
}

inline nlohmann::json sources_vec_as_json(const std::vector<std::string> & vs) {
    nlohmann::json j = nlohmann::json::array();
    for (auto src_entry : vs) {
        j.push_back(src_entry);
    }
    return j;
}

struct TaxonomyRecord {
    std::string line;
    OttId id = 0;
    OttId parent_id = 0;
    int parent_index = 0;
    boost::string_ref name;
    boost::string_ref rank;
    boost::string_ref sourceinfo;
    boost::string_ref uniqname; // will point to name field, if empty in .tsv
    std::bitset<32> flags;
    int depth = 0;
    int out_degree = 0;
    TaxonomyRecord& operator=(TaxonomyRecord&& tr) = default;
    TaxonomyRecord& operator=(const TaxonomyRecord& tr) = delete;
    TaxonomyRecord(TaxonomyRecord&& tr) = default;
    TaxonomyRecord(TaxonomyRecord& tr) = delete;
    explicit TaxonomyRecord(const std::string& line);
    std::vector<std::string> sourceinfoAsVec() const {
        std::string si = std::string(sourceinfo);
        return comma_separated_as_vec(si);
    }
};

// returns a vector [0, 1, 2, ... sz-1]
inline std::vector<int> get_index_vec(std::size_t sz) {
    std::vector<int> iv(sz);
    for(int i = 0 ; i < (int)sz; i++) {
        iv[i] = i;
    }
    return iv;
}

class BaseTaxonomy {
    protected:
    std::unordered_map<OttId, OttId> forwards;
    OttId keep_root;
    std::bitset<32> cleaning_flags;
    std::string path;
    std::string version;
    std::string version_number;
    BaseTaxonomy(const std::string& dir, std::bitset<32> cf=std::bitset<32>(), OttId keep_root=-1);
    
    public:
    const std::string & get_version() const {
        return version;
    }

    const std::string & get_version_number() const {
        return version_number;
    }
};

class Taxonomy: public std::vector<TaxonomyRecord>, public BaseTaxonomy {
    protected:
    std::unordered_map<OttId, int> index;
    void read_forwards_file(std::string filepath);
    public:
    template <typename Tree_t> std::unique_ptr<Tree_t> get_tree(std::function<std::string(const TaxonomyRecord&)>) const;

    TaxonomyRecord& record_from_id(OttId id);
    
    const TaxonomyRecord& record_from_id(OttId id) const;
    
    TaxonomyRecord& record_from_unforwarded_id(OttId id) {
        return at(index.at(id));
    }
    
    const TaxonomyRecord& record_from_unforwarded_id(OttId id) const {
        return at(index.at(id));
    }

    /// Index of the root record
    constexpr static int root_index() {
        return 0;
    }


    /// Given an OTT ID, forward if necessary and return the current Ott ID.
    //    Returns -1 if not found and not known to be deprecated.
    //    Not implemented yet: return -2 if known to be deprecated.
    //  Note that relying on -1 and -2 requires having read in the deprecated.tsv
    //    file (which is optional)
    OttId map(OttId id) const;

    /// Write out a taxonomy to directory dirname
    void write(const std::string& dirname);

    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    Taxonomy(const std::string& dir, std::bitset<32> cf=std::bitset<32>(), OttId keep_root=-1);
    friend class RichTaxonomy;
};

class TaxonomicJuniorSynonym;

class RTRichTaxNodeData {
    public:
    //TaxonomicRank rank = TaxonomicRank::RANK_NO_RANK;
    //nlohmann::json sources;
    std::vector<const TaxonomicJuniorSynonym *> junior_synonyms;
    std::uint32_t trav_enter = UINT32_MAX;
    std::uint32_t trav_exit = UINT32_MAX;
    TaxonomicRank rank = TaxonomicRank::RANK_NO_RANK;
    std::bitset<32> flags;
    std::string source_info;
    boost::string_ref possibly_nonunique_name;
    boost::string_ref get_nonuniqname() const {
        return possibly_nonunique_name;
    }
    const std::string & get_rank() const {
        return rank_to_string(rank);
    }
    const std::bitset<32> get_flags() const {
        return flags;
    }
    std::vector<std::string> sourceinfoAsVec() const {
        return comma_separated_as_vec(source_info);
    }
    nlohmann::json get_sources_json() const {
        auto vs = this->sourceinfoAsVec();
        return sources_vec_as_json(vs);
    }
};

typedef RootedTreeNode<RTRichTaxNodeData> RTRichTaxNode;

class TaxonomicJuniorSynonym {
    public:
    TaxonomicJuniorSynonym(const std::string & namestring,
                           const RTRichTaxNode *senior_synonym,
                           const std::string & sources)
        :name(namestring),
        source_string(sources),
        primary(senior_synonym) {
    }
    const std::string name;
    const std::string source_string;
    const RTRichTaxNode * primary;
    const std::string & get_name() const {
        return name;
    }
    // deleting copy ctor because we are using string_refs, so it is important
    //    that the object not be copied.
    TaxonomicJuniorSynonym(const TaxonomicJuniorSynonym &) = delete;
};

#define MAP_FOREIGN_TO_POINTER
class RTRichTaxTreeData {
    public:
#if defined(MAP_FOREIGN_TO_POINTER)
    std::unordered_map<OttId, const RTRichTaxNode *> ncbi_id_map;
    std::unordered_map<OttId, const RTRichTaxNode *> gbif_id_map;
    std::unordered_map<OttId, const RTRichTaxNode *> worms_id_map;
    std::unordered_map<OttId, const RTRichTaxNode *> if_id_map;
    std::unordered_map<OttId, const RTRichTaxNode *> irmng_id_map;
#else
    std::unordered_map<OttId, OttId> ncbi_id_map;
    std::unordered_map<OttId, OttId> gbif_id_map;
    std::unordered_map<OttId, OttId> worms_id_map;
    std::unordered_map<OttId, OttId> if_id_map;
    std::unordered_map<OttId, OttId> irmng_id_map;
#endif
    std::unordered_map<std::bitset<32>, nlohmann::json> flags2json;
    std::map<boost::string_ref, const RTRichTaxNode *> name_to_node; // null if homonym, then check homonym2node
    std::map<boost::string_ref, const TaxonomyRecord *> name_to_record; // for filtered
    std::unordered_map<OttId, const RTRichTaxNode *> id_to_node;
    std::unordered_map<OttId, const TaxonomyRecord *> id_to_record;
    std::map<boost::string_ref, std::vector<const RTRichTaxNode *> > homonym_to_node;
    std::map<boost::string_ref, std::vector<const TaxonomyRecord *> > homonym_to_record;
    std::map<std::string, OttIdSet> non_unique_taxon_names;
    
};

typedef RootedTree<RTRichTaxNodeData, RTRichTaxTreeData> RichTaxTree;

class RichTaxonomy: public BaseTaxonomy {
    public:
    const RichTaxTree & getTaxTree() const {
        return *tree;
    }
    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    RichTaxonomy(const std::string& dir, std::bitset<32> cf = std::bitset<32>(), OttId kr = -1);
    RichTaxonomy(RichTaxonomy &&) = default;
    const RTRichTaxNode * included_taxon_from_id(OttId ott_id) const {
        //Returns node * or nullptr if not found.
        const auto & td = tree->get_data();
        auto i2n_it = td.id_to_node.find(ott_id);
        if (i2n_it == td.id_to_node.end()) {
            auto fit = forwards.find(ott_id);
            if (fit == forwards.end()) {
                return nullptr;
            }
            return included_taxon_from_id(fit->second);
        }
        return i2n_it->second;
    }
    void add_taxonomic_addition_string(const std::string &s);
    const OttIdSet & get_ids_to_suppress_from_tnrs() const {
        return ids_to_suppress_from_tnrs;
    }
    void set_ids_suppressed_from_summary_tree_alias(const OttIdSet * ott_id_set_ptr) {
        is_suppressed_from_synth = ott_id_set_ptr;
    }

    const OttIdSet * get_ids_suppressed_from_summary_tree_alias() const {
        return is_suppressed_from_synth;
    }
    private:
    std::vector<TaxonomyRecord> filtered_records;
    std::unique_ptr<RichTaxTree> tree;
    std::list<TaxonomicJuniorSynonym> synonyms;
    // flags: not_otu, environmental, environmental_inherited, viral, hidden, hidden_inherited, was_container
    //    are excluded from being returned in TNRS results.
    OttIdSet ids_to_suppress_from_tnrs;
    // If there is just one summary tree in memory, the web services code
    //    can pass a non-nullptr pointer in using the setter. This
    //    will allow the services to report "is_suppressed_from_synth" option.
    const OttIdSet * is_suppressed_from_synth = nullptr;
    void read_synonyms();
    void _fill_ids_to_suppress_set();
    RichTaxonomy(const RichTaxonomy &) = delete;
};



//class RTTaxNodeData {
//    public:
//    const TaxonomyRecord * taxonomy_line = nullptr;
//};

template <typename Node_t, typename TREE>
void populate_node_from_taxonomy_record(Node_t & nd,
                                    const TaxonomyRecord & line,
                                    std::function<std::string(const TaxonomyRecord&)> get_name,
                                    TREE & tree);


// default behavior is to set ID and Name from line
template <typename Node_t, typename TREE>
inline void populate_node_from_taxonomy_record(Node_t & nd,
                                           const TaxonomyRecord & line,
                                           std::function<std::string(const TaxonomyRecord&)> get_name,
                                           TREE & ) {
    nd.set_ott_id(line.id);
    nd.set_name(get_name(line));    
}

#if 0
No longer used?
// default behavior is to set ID and Name from line
template <typename TREE>
inline void populate_node_from_taxonomy_record(RootedTreeNode<RTTaxNodeData> & nd,
                                           const TaxonomyRecord & line,
                                           std::function<std::string(const TaxonomyRecord&)>,
                                           TREE &) {
    nd.set_ott_id(line.id);
    nd.get_data().taxonomy_line = &line;    
}
#endif


template <typename Tree_t>
std::unique_ptr<Tree_t> Taxonomy::get_tree(std::function<std::string(const TaxonomyRecord&)> get_name) const {
    const auto& taxonomy = *this;
    std::unique_ptr<Tree_t> tree(new Tree_t);
    vector<typename Tree_t::node_type*> node_ptr(size(), nullptr);
    for(auto i = 0U; i < taxonomy.size() ; i++) {
        const auto& line = taxonomy[i];
        // Make the tree
        typename Tree_t::node_type* nd = nullptr;
        if (i==0) {
            nd = tree->create_root();
        } else {
            auto parent_nd = node_ptr[line.parent_index];
            nd = tree->create_child(parent_nd);
        }
        populate_node_from_taxonomy_record(*nd, line, get_name, *tree);
        node_ptr[i] = nd;
    }
    return tree;
}
// formatting options.
std::string format_with_taxonomy(const std::string& orig, const std::string& format, const TaxonomyRecord& rec);
std::string format_without_taxonomy(const std::string& orig, const std::string& format);
char format_needs_taxonomy(const std::string& format);

OttId root_ott_id_from_file(const std::string& filename);
std::string get_taxonomy_dir(const boost::program_options::variables_map& args);
Taxonomy load_taxonomy(const boost::program_options::variables_map& args);
RichTaxonomy load_rich_taxonomy(const boost::program_options::variables_map& args);

} // namespace
#endif
