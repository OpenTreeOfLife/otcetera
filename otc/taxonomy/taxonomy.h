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
#include <optional>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <string_view>
#include <bitset>
#include <variant>
#include "otc/taxonomy/flags.h"

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/flags.h"

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
    RANK_SUBTERCLASS,
    RANK_COHORT,
    RANK_SUBCOHORT,
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
    RANK_SERIES,
    RANK_SUBSERIES,
    RANK_SPECIES_GROUP,
    RANK_SPECIES_SUBGROUP,
    RANK_SPECIES,
    RANK_SUBSPECIES,
    RANK_INFRASPECIFICNAME,
    RANK_VARIETAS,
    RANK_VARIETY,
    RANK_SUBVARIETY,
    RANK_FORMA,
    RANK_SUBFORM,
    RANK_FORMA_SPECIALIS,
    RANK_MORPH,
    RANK_STRAIN,
    RANK_ISOLATE,
    RANK_BIOTYPE,
    RANK_PATHOGROUP,
    RANK_SEROGROUP,
    RANK_SEROTYPE,
    RANK_GENOTYPE,
    RANK_NO_RANK,
    RANK_NO_RANK_TERMINAL
};

bool rank_is_specific(TaxonomicRank rank);

extern const std::map<TaxonomicRank, std::string, std::less<>> rank_enum_to_name;
extern const std::map<std::string, TaxonomicRank, std::less<>> rank_name_to_enum;

extern const std::string empty_string;
extern const std::set<std::string> indexed_source_prefixes;

inline const std::string & rank_to_string(const TaxonomicRank &r) {
    auto i = rank_enum_to_name.find(r);
    if (i == rank_enum_to_name.end()) {
        return empty_string;
    }
    return i->second;
}


inline TaxonomicRank string_to_rank(const std::string_view& s, bool must_match_if_nonempty = false)
{
    // FIXME! 
    auto rank = rank_name_to_enum.find(s);
    if (rank == rank_name_to_enum.end()) {
        if (must_match_if_nonempty and !s.empty()) {
            throw OTCError() << "rank \"" << s << "\" not recognized.";
        }
        LOG(WARNING)<<"unknown rank '"<<s<<"'";
        return RANK_NO_RANK;
    }
    return rank->second;
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
    std::string_view name;
    std::string_view rank;
    std::string_view sourceinfo;
    std::string_view uniqname; // will point to name field, if empty in .tsv
    std::bitset<32> flags;
    int depth = 0; // 1 for root, 2 for root's children, etc
    int out_degree = 0;
    TaxonomyRecord& operator=(TaxonomyRecord&& tr) = default;
    TaxonomyRecord& operator=(const TaxonomyRecord& tr) = delete;
    TaxonomyRecord(TaxonomyRecord&& tr) = default;
    TaxonomyRecord(TaxonomyRecord& tr) = delete;
    bool is_extinct() const;
    explicit TaxonomyRecord(const std::string& line);
    explicit TaxonomyRecord(const std::string& line, bool is_short);
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

enum class reason_missing {unknown, not_an_id, never_minted_id, deprecated, pruned, forwarded, broken};

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

    // We should probably generalize this to record EITHER an OttId OR a reason why the OttId isn't found.
    virtual std::variant<OttId,reason_missing> get_unforwarded_id_or_reason(OttId id) const = 0;
    std::optional<OttId> get_unforwarded_id(OttId id) const;

    virtual ~BaseTaxonomy() = default;
};

class Taxonomy: public std::vector<TaxonomyRecord>, public BaseTaxonomy {
    protected:
    std::unordered_map<OttId, int> index;
    void read_forwards_file(std::string filepath);

    std::optional<int> maybe_index_from_id(OttId) const;
    int index_from_id(OttId) const;

public:
    static bool tolerate_synonyms_to_unknown_id;
    template <typename Tree_t> std::unique_ptr<Tree_t> get_tree(std::function<std::string(const TaxonomyRecord&)>,
                                                                bool process_source_maps = true) const;

    std::variant<OttId,reason_missing> get_unforwarded_id_or_reason(OttId id) const;

    TaxonomyRecord& record_from_id(OttId id);
    
    const TaxonomyRecord& record_from_id(OttId id) const;
    
    std::vector<std::string> path_from_id(OttId) const;

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
    void write(const std::string& dirname,
               bool copy_taxonomy_tsv_lines_raw,
               bool copy_synonyms_tsv_raw);
    void copy_relevant_synonyms(std::istream & inp, std::ostream & outp);

    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    Taxonomy(const std::string& dir, std::bitset<32> cf=std::bitset<32>(), OttId keep_root=-1);

    private:
    unsigned int read_input_taxonomy_stream(std::istream & taxonomy_stream);
    unsigned int read_ott_taxonomy_stream(std::istream & taxonomy_stream);
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
    std::string_view possibly_nonunique_name;
    uint32_t depth;
    std::string_view get_nonuniqname() const {
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
    std::string get_sources_as_fmt_str() const {
        return source_info;
    }
    nlohmann::json get_sources_json() const {
        auto vs = this->sourceinfoAsVec();
        return sources_vec_as_json(vs);
    }
    bool is_extinct() const;
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
    // deleting copy ctor because we are using string_views, so it is important
    //    that the object not be copied.
    TaxonomicJuniorSynonym(const TaxonomicJuniorSynonym &) = delete;
};

constexpr const char* sup_flag_comma = "not_otu,environmental,environmental_inherited,viral,hidden,hidden_inherited,was_container";

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
    std::map<std::string_view, const RTRichTaxNode *> name_to_node; // null if homonym, then check homonym2node
    std::map<std::string_view, const TaxonomyRecord *> name_to_record; // for filtered
    std::unordered_map<OttId, const RTRichTaxNode *> id_to_node;
    std::unordered_map<OttId, const TaxonomyRecord *> id_to_record;
    std::map<std::string_view, std::vector<const RTRichTaxNode *> > homonym_to_nodes;
    std::map<std::string_view, std::vector<const TaxonomyRecord *> > homonym_to_record;
    std::map<std::string, OttIdSet> non_unique_taxon_names;

    std::bitset<32> suppress_flags = flags_from_string(sup_flag_comma);
};

typedef RootedTree<RTRichTaxNodeData, RTRichTaxTreeData> RichTaxTree;

class ContextAwareCTrieBasedDB; // for fuzzy matching
class RichTaxonomy: public BaseTaxonomy {
    public:
    const RichTaxTree & get_tax_tree() const {
        return *tree;
    }
    /// Load the taxonomy from directory dir, and apply cleaning flags cf, and keep subtree below kr
    RichTaxonomy(const std::string& dir,
                 std::bitset<32> cf = std::bitset<32>(),
                 OttId kr = -1,
                 bool read_syn_type_as_src = false);
    RichTaxonomy(RichTaxonomy &&) = default;

    std::variant<OttId,reason_missing> get_unforwarded_id_or_reason(OttId id) const;

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
    
    const OttIdSet & get_ids_to_suppress_from_tnrs() const {
        return ids_to_suppress_from_tnrs;
    }
    
    void set_ids_suppressed_from_summary_tree_alias(const OttIdSet * ott_id_set_ptr) {
        is_suppressed_from_synth = ott_id_set_ptr;
    }

    bool node_is_suppressed_from_tnrs(const RTRichTaxNode * nd) const {
        auto& tree_data = tree->get_data();
        const auto& tax_record_flags = nd->get_data().get_flags();
        return (tree_data.suppress_flags & tax_record_flags).any();
    }

    const OttIdSet * get_ids_suppressed_from_summary_tree_alias() const {
        return is_suppressed_from_synth;
    }

    const std::list<TaxonomicJuniorSynonym> & get_synonyms_list() const {
        return synonyms;
    }

    const ContextAwareCTrieBasedDB * get_fuzzy_matcher() const {
        return fuzzy_match_db;
    }

    ContextAwareCTrieBasedDB * get_fuzzy_matcher() {
        return fuzzy_match_db;
    }

    void set_fuzzy_matcher(ContextAwareCTrieBasedDB * match_db) {
        fuzzy_match_db = match_db;
    }
    
    
    protected:
    RichTaxTree & get_mutable_tax_tree() const {
        return *tree;
    }
    void read_synonyms();
    void _fill_ids_to_suppress_set();
    
    bool read_synonym_type_as_src; // only relevant for source taxonomies with syntype info
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
    ContextAwareCTrieBasedDB * fuzzy_match_db = nullptr;
    RichTaxonomy(const RichTaxonomy &) = delete;
    private:
    void read_input_synonyms_stream(std::istream & synonyms_file);
    void read_ott_synonyms_stream(std::istream & synonyms_file);
};


//class RTTaxNodeData {
//    public:
//    const TaxonomyRecord * taxonomy_line = nullptr;
//};

template <typename Node_t, typename TREE>
void populate_node_from_taxonomy_record(Node_t & nd,
                                    const TaxonomyRecord & line,
                                    std::function<std::string(const TaxonomyRecord&)> get_name,
                                    TREE & tree,
                                    bool process_source_maps = true);


// default behavior is to set ID and Name from line
template <typename Node_t, typename TREE>
inline void populate_node_from_taxonomy_record(Node_t & nd,
                                           const TaxonomyRecord & line,
                                           std::function<std::string(const TaxonomyRecord&)> get_name,
                                           TREE & ,
                                           bool ) {
    nd.set_ott_id(line.id);
    nd.set_name(get_name(line));    
}

template<typename T>
inline void register_taxon_in_maps(std::map<std::string_view, const T *> & n2n,
                            std::map<std::string_view, std::vector<const T *> > & homonym_map,
                            std::string_view possibly_nonunique_name,
                            std::string_view uname,
                            const T * ti) {
    auto nit = n2n.lower_bound(possibly_nonunique_name);
    typedef std::pair<std::string_view, const T *> name_map_pair;
    if (nit == n2n.end() or nit->first != possibly_nonunique_name) {
        nit = n2n.insert(nit, name_map_pair(possibly_nonunique_name, ti));
    } else {
        if (nit->second != nullptr) {
            homonym_map[possibly_nonunique_name].push_back(nit->second);
            nit->second = nullptr;
        }
       homonym_map[possibly_nonunique_name].push_back(ti);
    }
    if (uname != possibly_nonunique_name) {
        auto r2 = n2n.insert(name_map_pair(uname, ti));
        assert(r2.second); // should be uniq.
    }
}

template<typename T>
void process_source_info_vec(const std::vector<std::string> & vs,
                             RTRichTaxTreeData & tree_data,
                             T & ,
                             const RTRichTaxNode * this_node);


template<typename T>
inline void process_source_info_vec(const std::vector<std::string> & vs,
                             RTRichTaxTreeData & tree_data,
                             T & ,
                             const RTRichTaxNode * this_node) {
    using std::string;
    for (auto src_entry : vs) {
        auto pref_id = split_string(src_entry, ':');
        if (pref_id.size() != 2) {
            throw OTCError() << "Expecting exactly 1 colon in a source ID string. Found: \"" << src_entry << "\".";
        }
        const string & prefix = *pref_id.begin();
        if (indexed_source_prefixes.count(prefix) == 0) {
            continue;
        }
        const string & id_str = *pref_id.rbegin();
        //data.sources.push_back(src_entry);
        std::size_t pos;
        //LOG(INFO) << src_entry;
        try {
            long raw_foreign_id  = std::stoul(id_str.c_str(), &pos);
            if (pos < id_str.length() || raw_foreign_id < 0) {
                throw OTCError() << "Could not convert ID to unsigned long \"" << src_entry << "\"";
            }
            OttId foreign_id = check_ott_id_size(raw_foreign_id);
#           if defined(MAP_FOREIGN_TO_POINTER)
                auto to_map_to = this_node;
#           else
                auto to_map_to = this_node->get_ott_id();
#           endif

            if (prefix == "ncbi") {
                tree_data.ncbi_id_map[foreign_id] = to_map_to;
            } else if (prefix == "gbif") {
                tree_data.gbif_id_map[foreign_id] = to_map_to;
            } else if (prefix == "worms") {
                tree_data.worms_id_map[foreign_id] = to_map_to;
            } else if (prefix == "if") {
                tree_data.if_id_map[foreign_id] = to_map_to;
            } else if (prefix == "irmng") {
                tree_data.irmng_id_map[foreign_id] = to_map_to;
            } else {
                assert(false);
            }
        } catch (OTCError & x) {
            throw;
        } catch (...) {
            LOG(WARNING) << "Could not convert ID to unsigned long \"" << src_entry << "\"";
        }
    }
}

template<typename T, typename U> 
inline void add_f_to_json_if_needed(T & flags2json,
                                    U & flags) {
    using std::string;
    using std::vector;
    using nlohmann::json;
    
    // If the flag combination is new, store the JSON representation
    if (flags2json.count(flags) == 0) {
        vector<string> vf = flags_to_string_vec(flags);
        flags2json[flags] = json();
        auto & fj = flags2json[flags];
        for (auto fs : vf) {
            fj.push_back(fs);
        }
    }
}

template <>
inline void populate_node_from_taxonomy_record(RTRichTaxNode & nd,
                                               const TaxonomyRecord & tr,
                                               std::function<std::string(const TaxonomyRecord&)> ,
                                               RichTaxTree & tree,
                                               bool process_source_maps) {
    using std::string;
    using std::vector;
    using std::string_view;
    RTRichTaxNode * this_node = &nd;
    nd.set_ott_id(tr.id);
    auto & data = nd.get_data();
    auto & tree_data = tree.get_data();
    nd.set_ott_id(tr.id);
    tree_data.id_to_node[tr.id] = this_node;
    this_node->set_name(string(tr.uniqname));
    const string & uname = this_node->get_name();
    if (tr.uniqname != tr.name) {
        string sn = string(tr.name);
        tree_data.non_unique_taxon_names[sn].insert(tr.id);
        auto nit = tree_data.non_unique_taxon_names.find(sn);
        assert(nit != tree_data.non_unique_taxon_names.end());
        data.possibly_nonunique_name = string_view(nit->first);
    } else {
        data.possibly_nonunique_name = string_view(nd.get_name());
    }
    data.flags = tr.flags;
    data.rank = string_to_rank(tr.rank);
    register_taxon_in_maps(tree_data.name_to_node,
                           tree_data.homonym_to_nodes,
                           data.possibly_nonunique_name,
                           uname,
                           this_node);
    auto flags = data.get_flags();
    //cout << "flags = " << flags << " name = " << this_node->get_name() << '\n';
    add_f_to_json_if_needed(tree_data.flags2json, flags);
    if (process_source_maps) {
        auto vs = tr.sourceinfoAsVec();
        data.source_info = string(tr.sourceinfo);
        process_source_info_vec(vs, tree_data, data, this_node);
    }
}

template <typename Tree_t>
std::unique_ptr<Tree_t> Taxonomy::get_tree(std::function<std::string(const TaxonomyRecord&)> get_name,
                                           bool process_source_maps) const {
    const auto& taxonomy = *this;
    std::unique_ptr<Tree_t> tree(new Tree_t);
    vector<typename Tree_t::node_type*> node_ptr(size(), nullptr);
    for(auto i = 0U; i < taxonomy.size() ; i++) {
        const auto& line = taxonomy[i];
        // std::cerr << "Taxonomy::get_tree line " << i << '\n';
        // Make the tree
        typename Tree_t::node_type* nd = nullptr;
        if (i==0) {
            nd = tree->create_root();
        } else {
            auto parent_nd = node_ptr[line.parent_index];
            nd = tree->create_child(parent_nd);
        }
        populate_node_from_taxonomy_record(*nd, line, get_name, *tree, process_source_maps);
        node_ptr[i] = nd;
    }
    return tree;
}

// formatting options.
std::string format_with_taxonomy(const std::string& orig, const std::string& format, const TaxonomyRecord& rec, const Taxonomy& taxonomy);
std::string format_without_taxonomy(const std::string& orig, const std::string& format);
char format_needs_taxonomy(const std::string& format);

OttId root_ott_id_from_file(const std::string& filename);
std::string get_taxonomy_dir(const boost::program_options::variables_map& args);
Taxonomy load_taxonomy(const boost::program_options::variables_map& args);
RichTaxonomy load_rich_taxonomy(const boost::program_options::variables_map& args);


const RTRichTaxNode* taxonomy_mrca(const std::vector<const RTRichTaxNode*>& nodes);
std::vector<const RTRichTaxNode*> exact_name_search(const RichTaxonomy& taxonomy,
                                                    const std::string& query,
                                                    bool include_suppressed);

std::vector<const RTRichTaxNode*> exact_name_search(const RichTaxonomy& taxonomy,
                                                    const RTRichTaxNode* context_root,
                                                    const std::string& query,
                                                    bool include_suppressed);

std::vector<const RTRichTaxNode*> exact_name_search(const RichTaxonomy& taxonomy,
                                                    const RTRichTaxNode* context_root,
                                                    const std::string& query,
                                                    std::function<bool(const RTRichTaxNode*)> ok = [](const RTRichTaxNode*){return true;});


template<typename N>
N * find_mrca_via_traversal_indices(N *f, N *s);

template<typename N>
inline N * find_mrca_via_traversal_indices(N *f, N *s) {
    const auto * fdata = &(f->get_data());
    const auto sec_ind = s->get_data().trav_enter;
    while (sec_ind < fdata->trav_enter || sec_ind > fdata->trav_exit) {
        f = f->get_parent();
        if (f == nullptr) {
            assert(false); 
            return nullptr;
        }
        fdata = &(f->get_data());
    }
    return f;
}


} // namespace
#endif
