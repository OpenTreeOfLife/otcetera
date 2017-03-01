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

extern const std::map<TaxonomicRank, std::string> rankEnumToName;
extern const std::string emptyStringForMissingRank;
extern const std::set<std::string> indexed_source_prefixes;

inline const std::string & to_string(const TaxonomicRank &r) {
    auto i = rankEnumToName.find(r);
    if (i == rankEnumToName.end()) {
        return emptyStringForMissingRank;
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
    nlohmann::json j;
    for (auto src_entry : vs) {
        j.push_back(src_entry);
    }
    return j;

}

struct taxonomy_record {
    std::string line;
    long id = 0;
    long parent_id = 0;
    int parent_index = 0;
    boost::string_ref name;
    boost::string_ref rank;
    boost::string_ref sourceinfo;
    boost::string_ref uniqname; // will point to name field, if empty in .tsv
    std::bitset<32> flags;
    int depth = 0;
    int out_degree = 0;
    taxonomy_record& operator=(taxonomy_record&& tr) = default;
    taxonomy_record& operator=(const taxonomy_record& tr) = delete;
    taxonomy_record(taxonomy_record&& tr) = default;
    taxonomy_record(taxonomy_record& tr) = delete;
    explicit taxonomy_record(const std::string& line);
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


class Taxonomy: public std::vector<taxonomy_record> {
    protected:
    std::unordered_map<long, int> index;
    std::unordered_map<long, long> forwards;
    std::unordered_set<long> deprecated;
    long keep_root;
    std::bitset<32> cleaning_flags;
    std::string path;
    std::string version;
    std::string version_number;
    public:
    template <typename Tree_t> std::unique_ptr<Tree_t> getTree(std::function<std::string(const taxonomy_record&)>) const;

    const std::string & get_version() const {
        return version;
    }

    const std::string & get_version_number() const {
        return version_number;
    }

    taxonomy_record& record_from_id(long id);
    
    const taxonomy_record& record_from_id(long id) const;
    
    taxonomy_record& record_from_unforwarded_id(long id) {
        return at(index.at(id));
    }
    
    const taxonomy_record& record_from_unforwarded_id(long id) const {
        return at(index.at(id));
    }

    /// Index of the root record
    constexpr static int root_index() {
        return 0;
    }


    /// Given an OTT ID, forward if necessary and return the current Ott ID.
    //    Returns -1 if not found and not known to be deprecated.
    //    -2 if known to be deprecated.
    //  Note that relying on -1 and -2 requires having read in the deprecated.tsv
    //    file (which is optional)
    long map(long id) const;

    /// Write out a taxonomy to directory dirname
    void write(const std::string& dirname);

    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    Taxonomy(const std::string& dir, std::bitset<32> cf=std::bitset<32>(), long keep_root=-1);
    protected:
    void read_forwards_file(std::string filepath);
};

class TaxonomicJuniorSynonym;

class RTRichTaxNodeData {
    public:
    //TaxonomicRank rank = TaxonomicRank::RANK_NO_RANK;
    //nlohmann::json sources;
    const taxonomy_record * tax_record;
    std::vector<const TaxonomicJuniorSynonym *> junior_synonyms;
    boost::string_ref getName() const {
        return tax_record->name;
    }
    boost::string_ref getUniqname() const {
        return tax_record->uniqname;
    }
    boost::string_ref getRank() const {
        return tax_record->rank;
    }
    const std::bitset<32> getFlags() const {
        return tax_record->flags;
    }
    nlohmann::json getSourcesJSON() const {
        auto vs = tax_record->sourceinfoAsVec();
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
    const std::string & getName() const {
        return name;
    }
    // deleting copy ctor because we are using string_refs, so it is important
    //    that the object not be copied.
    TaxonomicJuniorSynonym(const TaxonomicJuniorSynonym &) = delete;
};


class RTRichTaxTreeData {
    public:
    std::unordered_map<unsigned long, const RTRichTaxNode *> ncbi_id_map;
    std::unordered_map<unsigned long, const RTRichTaxNode *> gbif_id_map;
    std::unordered_map<unsigned long, const RTRichTaxNode *> worms_id_map;
    std::unordered_map<unsigned long, const RTRichTaxNode *> if_id_map;
    std::unordered_map<unsigned long, const RTRichTaxNode *> irmng_id_map;
    std::unordered_map<std::bitset<32>, nlohmann::json> flags2json;
    std::map<boost::string_ref, const RTRichTaxNode *> name2node; // null if homonym, then check homonym2node

    std::unordered_map<OttId, const RTRichTaxNode *> id2node;
    std::map<boost::string_ref, std::vector<const RTRichTaxNode *> > homonym2node;
};

typedef RootedTree<RTRichTaxNodeData, RTRichTaxTreeData> RichTaxTree;

struct RichTaxonomy : public Taxonomy {
    const RichTaxTree & getTaxTree() const {
        return *tree;
    }
    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    RichTaxonomy(const std::string& dir, std::bitset<32> cf = std::bitset<32>(), long kr = -1);
    RichTaxonomy(RichTaxonomy &&) = default;
    const RTRichTaxNode * taxonFromId(OttId ott_id) const {
        //Returns node * or nullptr if not found.
        const auto & td = tree->getData();
        auto i2n_it = td.id2node.find(ott_id);
        if (i2n_it == td.id2node.end()) {
            auto fit = forwards.find(ott_id);
            if (fit == forwards.end()) {
                return nullptr;
            }
            return taxonFromId(fit->second);
        }
        return i2n_it->second;
    }
    void addTaxonomicAdditionString(const std::string &s);
    private:
    std::unique_ptr<RichTaxTree> tree;
    std::list<TaxonomicJuniorSynonym> synonyms;
    void readSynonyms();
    RichTaxonomy(const RichTaxonomy &) = delete;
};



class RTTaxNodeData {
    public:
    const taxonomy_record * taxonomy_line = nullptr;
};

template <typename Node_t, typename TREE>
void populateNodeFromTaxonomyRecord(Node_t & nd,
                                    const taxonomy_record & line,
                                    std::function<std::string(const taxonomy_record&)> getName,
                                    TREE & tree);


// default behavior is to set ID and Name from line
template <typename Node_t, typename TREE>
inline void populateNodeFromTaxonomyRecord(Node_t & nd,
                                           const taxonomy_record & line,
                                           std::function<std::string(const taxonomy_record&)> getName,
                                           TREE & ) {
    nd.setOttId(line.id);
    nd.setName(getName(line));    
}

// default behavior is to set ID and Name from line
template <typename TREE>
inline void populateNodeFromTaxonomyRecord(RootedTreeNode<RTTaxNodeData> & nd,
                                           const taxonomy_record & line,
                                           std::function<std::string(const taxonomy_record&)>,
                                           TREE &) {
    nd.setOttId(line.id);
    nd.getData().taxonomy_line = &line;    
}


template <typename Tree_t>
std::unique_ptr<Tree_t> Taxonomy::getTree(std::function<std::string(const taxonomy_record&)> getName) const {
    const auto& taxonomy = *this;
    std::unique_ptr<Tree_t> tree(new Tree_t);
    vector<typename Tree_t::node_type*> node_ptr(size(), nullptr);
    for(auto i = 0U; i < taxonomy.size() ; i++) {
        const auto& line = taxonomy[i];
        // Make the tree
        typename Tree_t::node_type* nd = nullptr;
        if (i==0) {
            nd = tree->createRoot();
        } else {
            auto parent_nd = node_ptr[line.parent_index];
            nd = tree->createChild(parent_nd);
        }
        populateNodeFromTaxonomyRecord(*nd, line, getName, *tree);
        node_ptr[i] = nd;
    }
    return tree;
}
// formatting options.
std::string format_with_taxonomy(const std::string& orig, const std::string& format, const taxonomy_record& rec);
std::string format_without_taxonomy(const std::string& orig, const std::string& format);
char format_needs_taxonomy(const std::string& format);

long root_ott_id_from_file(const std::string& filename);
std::string get_taxonomy_dir(const boost::program_options::variables_map& args);
Taxonomy load_taxonomy(const boost::program_options::variables_map& args);
RichTaxonomy load_rich_taxonomy(const boost::program_options::variables_map& args);

} // namespace
#endif
