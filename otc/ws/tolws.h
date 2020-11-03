
#ifndef OTC_TOLWS_H
#define OTC_TOLWS_H
#include <cstdint>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/filesystem/operations.hpp>
#include <stdexcept>
#include <memory>
#include "otc/newick_tokenizer.h"
#include "otc/newick.h"
#include "otc/tree.h"
#include "otc/error.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/ws/parallelreadserialwrite.h"
#include "otc/ws/otc_web_error.h"
#include "json.hpp"

#define REPORT_MEMORY_USAGE 1

#if defined(REPORT_MEMORY_USAGE)
#include "otc/memory_usage.h"
#endif

namespace otc {

typedef std::pair<const std::string *, const std::string *> src_node_id;
typedef std::vector<src_node_id> vec_src_node_id_mapper;

typedef RTRichTaxNode Taxon;

#define JOINT_MAPPING_VEC

#if defined(JOINT_MAPPING_VEC)
    enum SourceEdgeMappingType {
        CONFLICTS_WITH_MAPPING = 0,
        PARTIAL_PATH_OF_MAPPING = 1,
        RESOLVES_MAPPING = 2,
        SUPPORTED_BY_MAPPING = 3,
        TERMINAL_MAPPING = 4
    };
    typedef std::pair<SourceEdgeMappingType, std::uint32_t> semt_ind_t;
    typedef std::vector<semt_ind_t> vec_src_node_ids;
#else
    typedef std::vector<std::uint32_t> vec_src_node_ids;
#endif

class SumTreeNodeData {
    public:
    std::uint32_t trav_enter = UINT32_MAX;
    std::uint32_t trav_exit = UINT32_MAX;
#   if defined(JOINT_MAPPING_VEC)
        vec_src_node_ids source_edge_mappings;
#   else
        vec_src_node_ids supported_by;
        vec_src_node_ids conflicts_with;
        vec_src_node_ids resolves;
        vec_src_node_ids partial_path_of;
        vec_src_node_ids terminal;
#   endif
    bool was_uncontested = false;
    bool extinct_mark = false;  // extinctness means that the node has >= 1 descendant (including itself), and all descendants are extinct.
    bool is_extinct() const {return extinct_mark;}
    uint32_t num_tips = 0;
};

#if defined(REPORT_MEMORY_USAGE)
template<>
inline std::size_t calc_memory_used(const src_node_id &,
                                    MemoryBookkeeper &) {
    return  2 * sizeof(const std::string *); // we are aliasing stored strings, so we don't count len
}

template<>
inline std::size_t calc_memory_used(const SumTreeNodeData &d, MemoryBookkeeper &mb) {
    std::size_t total = 3 * sizeof(std::uint32_t); // traversals and num_tips
    mb["tree node data trav + num_tips"] += total;
    std::size_t x;
#if defined(JOINT_MAPPING_VEC)
    std::size_t v_el_size = sizeof(SourceEdgeMappingType) + sizeof(std::uint32_t); // 2 * sizeof(const std::string *);
    x = calc_memory_used_by_vector_eqsize(d.source_edge_mappings, v_el_size, mb);
    mb["tree node data source_edge_mappings"] += x; total += x;
#else
    std::size_t v_el_size = sizeof(std::uint32_t); // 2 * sizeof(const std::string *);
    x = calc_memory_used_by_vector_eqsize(d.supported_by, v_el_size, mb);
    mb["tree node data supported_by"] += x; total += x;
    x = calc_memory_used_by_vector_eqsize(d.conflicts_with, v_el_size, mb);
    mb["tree node data conflicts_with"] += x; total += x;
    x = calc_memory_used_by_vector_eqsize(d.resolves, v_el_size, mb);
    mb["tree node data resolves"] += x; total += x;
    x = calc_memory_used_by_vector_eqsize(d.partial_path_of, v_el_size, mb);
    mb["tree node data partial_path_of"] += x; total += x;
    x = calc_memory_used_by_vector_eqsize(d.terminal, v_el_size, mb);
    mb["tree node data terminal"] += x; total += x;
#endif
    total += sizeof(bool);
    return total;
}

#endif


typedef RootedTreeNode<SumTreeNodeData> SumTreeNode_t;
typedef std::vector<const SumTreeNode_t *> SumTreeNodeVec_t;
typedef std::pair<const SumTreeNode_t *, SumTreeNodeVec_t> BrokenMRCAAttachVec;

struct attachment_point_t
{
    std::string parent;
    std::vector<std::string> children_from_taxon;
};

struct contesting_tree_t
{
    std::string tree;
    std::vector<attachment_point_t> attachment_points;
};

class SumTreeData {
    public:
    // maps ottX or mrcaottXottY to node* if node lacks ottid.
    std::unordered_map<std::string, const SumTreeNode_t *> broken_name_to_node;
    std::unordered_map<OttId, const SumTreeNode_t *> id_to_node;
    
    // maps ottX to node* if taxon is broken.
    std::unordered_map<std::string, BrokenMRCAAttachVec> broken_taxa;

    // string -> [(Tree,[(parent,[children])])]
    std::unordered_map<std::string, std::vector<contesting_tree_t>> contesting_trees_for_taxon;
};

using SummaryTree_t = otc::RootedTree<SumTreeNodeData, SumTreeData>;


#if defined(REPORT_MEMORY_USAGE)

template<>
inline std::size_t calc_memory_used(const BrokenMRCAAttachVec &d, MemoryBookkeeper &mb) {
    std::size_t total = sizeof(const SumTreeNode_t *);
    total += calc_memory_used_by_vector_eqsize(d.second, sizeof(const SumTreeNode_t *), mb);
    return total;
}

template<>
inline std::size_t calc_memory_used(const SumTreeData &d, MemoryBookkeeper &mb) {
    const auto bn2nnum_unused_buckets = d.broken_name_to_node.bucket_count() - d.broken_name_to_node.size();
    std::size_t bn2nmem = bn2nnum_unused_buckets * (sizeof(std::string) + sizeof(SumTreeNode_t *));
    for (auto n : d.broken_name_to_node) {
        bn2nmem += calc_memory_used(n.first, mb) + sizeof(const SumTreeNode_t *);
    }
    const auto i2nnum_unused_buckets = d.id_to_node.bucket_count() - d.broken_name_to_node.size();
    std::size_t i2nmem = (d.id_to_node.size() + i2nnum_unused_buckets) * (sizeof(OttId) + sizeof(SumTreeNode_t *));
    const auto btnum_unused_buckets = d.broken_taxa.bucket_count() - d.broken_taxa.size();
    std::size_t btmem = btnum_unused_buckets * (sizeof(std::string) + sizeof(BrokenMRCAAttachVec ));
    for (auto n : d.broken_taxa) {
        btmem += calc_memory_used(n.first, mb) + calc_memory_used(n.second, mb);
    }
    mb["SumTreeData broken_name_to_node"] += bn2nmem;
    mb["SumTreeData id_to_node"] += i2nmem;
    mb["SumTreeData broken_taxa"] += btmem;
    return btmem + i2nmem + bn2nmem;
}
#endif

struct SourceTreeId {
    std::string tree_id;
    std::string study_id;
    std::string git_sha;
};

struct SummaryTreeAnnotation;
void from_json(const nlohmann::json &j, SummaryTreeAnnotation & sta);

struct SummaryTreeAnnotation {
    public:
        bool initialized = false;
        std::string date_completed;
        std::string filtered_flags;
        std::vector<std::string> filtered_flags_vec;
        std::string generated_by;
        OttId num_leaves_in_exemplified_taxonomy;
        OttId num_source_studies;
        OttId num_source_trees;
        OttId num_tips;
        OttId root_ott_id;
        std::string root_taxon_name;
        //source_id_map
        //sources
        std::string synth_id;
        std::string taxonomy_version;
        std::string tree_id;
        std::map<std::string, SourceTreeId> source_id_map;
        std::vector<std::string> sources;
        nlohmann::json full_source_id_map_json;
        OttIdSet suppressed_from_tree; // IDs not included in tree
        explicit SummaryTreeAnnotation() {
        }
        SummaryTreeAnnotation(const SummaryTreeAnnotation &) = delete;
        SummaryTreeAnnotation(const SummaryTreeAnnotation &&) = delete;
        SummaryTreeAnnotation &operator=(const SummaryTreeAnnotation &) = delete;
        SummaryTreeAnnotation &operator=(const nlohmann::json &j) {
            from_json(j, *this);
            return *this;
        }
};

inline std::string ott_id_to_idstr(OttId ott_id) {
    std::string ret;
    ret.reserve(12);
    ret = "ott";
    ret += std::to_string(ott_id);
    return ret;
}

enum NodeNameStyle {
    NNS_NAME_ONLY = 0,
    NNS_ID_ONLY = 1,
    NNS_NAME_AND_ID = 2
};

template<typename T>
void index_by_name_or_id(T & tree) {
    const std::string empty_string;
    auto & td = tree.get_data();
    auto & bm = td.broken_name_to_node;
    auto & im = td.id_to_node;
    for (auto nd : iter_pre(tree)) {
        if (nd->has_ott_id()) {
            nd->set_name(empty_string);
            im[nd->get_ott_id()] = nd;
        } else {
            bm[nd->get_name()] = nd;
        }
    }
}

/// The complete taxonomy of looking up is like:
///    ottX -> {too large, _ -> {never valid, _ -> {deprecated (previously valid), _ -> {pruned, _ -> {broken, OK!}}}}}
///     mrcaottXottY -> OK | BadOTT (ReasonOTTMissing) | BadMRCA (ReasonOTTMissing)
///       In this case if X or Y has a "broken" result, then we needn't fail.

struct node_lookup_t
{
    const SumTreeNode_t* node = nullptr;
    bool was_broken = false;

    bool broken() const {return was_broken;}
    bool present() const {return node and not was_broken;}
    // technically, this includes: {too large, never valid, deprecated (once valid), and PRUNED}
    bool invalid() const {return node == nullptr and not was_broken;}

    node_lookup_t() {};
    node_lookup_t(const SumTreeNode_t* n):node(n) {}
};

node_lookup_t find_node_by_id_str(const SummaryTree_t & tree, const RichTaxonomy&, const std::string & node_id);

class TreesToServe;

std::string available_trees_ws_method(const TreesToServe &tts);

std::string about_ws_method(const TreesToServe &tts,
                            const SummaryTree_t * tree_ptr,
                            const SummaryTreeAnnotation * sta,
                            bool include_sources);

std::string tax_about_ws_method(const RichTaxonomy &tts);

std::string node_info_ws_method(const TreesToServe & tts,
                                const SummaryTree_t * tree_ptr,
                                const SummaryTreeAnnotation * sta,
                                const std::string & node_id,
                                bool include_lineage);

std::string nodes_info_ws_method(const TreesToServe & tts,
                                 const SummaryTree_t * tree_ptr,
                                 const SummaryTreeAnnotation * sta,
                                 const std::vector<std::string> & node_id,
                                 bool include_lineage);

std::string mrca_ws_method(const TreesToServe & tts,
                           const SummaryTree_t * tree_ptr,
                           const SummaryTreeAnnotation * sta,
                           const std::vector<std::string> & node_id_vec,
                           const std::vector<std::string> & excluded_node_ids,
                           bool soft_exclude);

std::string induced_subtree_ws_method(const TreesToServe & tts,
                                      const SummaryTree_t * tree_ptr,
                                      const std::vector<std::string> & node_id_vec,
                                      NodeNameStyle label_format);

std::string newick_subtree_ws_method(const TreesToServe & tts,
                                     const SummaryTree_t * tree_ptr,
                                     const std::string & node_id,
                                     NodeNameStyle label_format,
                                     bool include_all_node_labels,
                                     int height_limit);

std::string arguson_subtree_ws_method(const TreesToServe & tts,
                                      const SummaryTree_t * tree_ptr,
                                      const SummaryTreeAnnotation * sta,
                                      const std::string & node_id,
                                      int height_limit);

std::string taxon_info_ws_method(const RichTaxonomy & taxonomy,
                                 const RTRichTaxNode * taxon_node,
                                 bool include_lineage,
                                 bool include_children,
                                 bool include_terminal_descendants);

std::string taxonomy_flags_ws_method(const RichTaxonomy & taxonomy);

std::string taxonomy_mrca_ws_method(const RichTaxonomy & taxonomy,
                                    const OttIdSet & ott_id_set);

std::string taxon_subtree_ws_method(const TreesToServe & tts,
                                    const RichTaxonomy & taxonomy,
                                    const RTRichTaxNode * taxon_node,
                                    NodeNameStyle label_format);

std::string tnrs_match_names_ws_method(const std::vector<std::string>& names,
                                       const std::optional<std::string>& context_name,
                                       bool do_approximate_matching,
                                       const std::optional<std::vector<std::string>>& ids,
                                       bool include_suppressed,
                                       const RichTaxonomy& taxonomy);
std::string tnrs_autocomplete_name_ws_method(const std::string& name,
                                             const std::string& context_name,
                                             bool include_suppressed,
                                             const RichTaxonomy& taxonomy);
std::string tnrs_contexts_ws_method();
std::string tnrs_infer_context_ws_method(const std::vector<std::string>& names,
                                         const RichTaxonomy& taxonomy);


std::string newick_conflict_ws_method(const SummaryTree_t & summary,
                                      const RichTaxonomy & taxonomy,
                                      const std::string& tree1s,
                                      const std::string& tree2s);

std::string phylesystem_conflict_ws_method(const SummaryTree_t & summary,
                                           const RichTaxonomy & taxonomy,
                                           const std::string& tree1s,
                                           const std::string& tree2s);

bool read_trees(const boost::filesystem::path & dirname, TreesToServe & tts);

void from_json(const nlohmann::json &j, SourceTreeId & sti);
void to_json(nlohmann::json &j, const SourceTreeId & sti);


inline std::string extract_string(const nlohmann::json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_string()) {
        return dc_el->get<std::string>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}

inline long extract_unsigned_long(const nlohmann::json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_number()) {
        return dc_el->get<long>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a non-negative integer.\n";
}


inline OttId extract_ott_id(const nlohmann::json &j, const char * field) {
    long ol = extract_unsigned_long(j, field);
    return check_ott_id_size(ol);
}

inline void to_json(nlohmann::json &j, const SourceTreeId & sti) {
    j = nlohmann::json{{"git_sha", sti.git_sha}, {"study_id", sti.study_id}, {"tree_id", sti.tree_id}};
}

inline void from_json(const nlohmann::json &j, SourceTreeId & sti) {
    sti.git_sha = extract_string(j, "git_sha");
    sti.study_id = extract_string(j, "study_id");
    sti.tree_id = extract_string(j, "tree_id");
}

nlohmann::json tax_about_json(const RichTaxonomy & taxonomy);
void tax_service_add_taxon_info(const RichTaxonomy & taxonomy, const RTRichTaxNode & nd_taxon, nlohmann::json & taxonrepr);
void tax_service_add_suppressed_taxon_info(const RichTaxonomy & taxonomy, const TaxonomyRecord & nd_taxon, nlohmann::json & taxonrepr);

std::string get_synth_node_label(const SumTreeNode_t* node);

inline const std::string & get_taxon_unique_name(const RTRichTaxNode & nd_taxon) {
    return nd_taxon.get_name();
}

inline void add_taxon_info(const RichTaxonomy & ,
                           const RTRichTaxNode & nd_taxon,
                           nlohmann::json & taxonrepr) {
    const auto & taxon_data = nd_taxon.get_data();
    taxonrepr["tax_sources"] = taxon_data.get_sources_json();
    taxonrepr["name"] = taxon_data.get_nonuniqname();
    taxonrepr["unique_name"] = get_taxon_unique_name(nd_taxon);
    taxonrepr["rank"] = taxon_data.get_rank();
    taxonrepr["ott_id"] = nd_taxon.get_ott_id();    
}

inline void add_taxon_record_info(const RichTaxonomy & ,
                           const TaxonomyRecord & record,
                           nlohmann::json & taxonrepr) {
    const std::string srcs{record.sourceinfo};
    const auto av = comma_separated_as_vec(srcs);
    taxonrepr["tax_sources"] = sources_vec_as_json(av);
    taxonrepr["name"] = std::string{record.name};
    taxonrepr["unique_name"] = std::string{record.uniqname};
    taxonrepr["rank"] = std::string{record.rank};
    taxonrepr["ott_id"] = record.id;    
}

} // namespace otc
#endif
