#ifndef OTC_TOLWS_H
#define OTC_TOLWS_H
#include <list>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/filesystem/operations.hpp>
#include "otc/newick_tokenizer.h"
#include "otc/newick.h"
#include "otc/tree.h"
#include "otc/error.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "json.hpp"

namespace otc {

typedef std::pair<const std::string *, const std::string *> src_node_id;
typedef std::vector<src_node_id> vec_src_node_ids;

template<typename T>
void index_nodes_by_name(T & tree) {
    auto & td = tree.get_data();
    auto & m = td.name_to_node;
    long ind = 0;
    for (auto nd : iter_pre(tree)) {
        m[nd->get_name()] = nd;
        nd->get_data().trav_enter = ind++;
    }
    for (auto pnd : iter_post(tree)) {
        auto fc = pnd->get_last_child();
        auto & d = pnd->get_data();
        if (fc == nullptr) {
            d.trav_exit = d.trav_enter;
            d.num_tips = 1;
        } else {
            d.trav_exit = fc->get_data().trav_exit;
            d.num_tips = 0;
            for (auto c : iter_child_const(*pnd)) {
                d.num_tips += c->get_data().num_tips;
            }
        }
    }
}

template<typename T>
void set_traversal_entry_exit(T & tree) {
    auto & td = tree.get_data();
    auto & m = td.name_to_node;
    long ind = 0;
    for (auto nd : iter_pre(tree)) {
        m[nd->get_name()] = nd;
        nd->get_data().trav_enter = ind++;
    }
    for (auto pnd : iter_post(tree)) {
        auto fc = pnd->get_last_child();
        auto & d = pnd->get_data();
        if (fc == nullptr) {
            d.trav_exit = d.trav_enter;
            d.num_tips = 1;
        } else {
            d.trav_exit = fc->get_data().trav_exit;
            d.num_tips = 0;
            for (auto c : iter_child_const(*pnd)) {
                d.num_tips += c->get_data().num_tips;
            }
        }
    }
}

class SumTreeNodeData {
    public:
    long trav_enter = -1;
    long trav_exit = -1;
    vec_src_node_ids supported_by;
    vec_src_node_ids conflicts_with;
    vec_src_node_ids resolves;
    vec_src_node_ids partial_path_of;
    vec_src_node_ids terminal;
    bool was_uncontested = false;
    long num_tips = 0;
};

typedef RootedTreeNode<SumTreeNodeData> SumTreeNode_t;
typedef std::vector<const SumTreeNode_t *> SumTreeNodeVec_t;
typedef std::pair<const SumTreeNode_t *, SumTreeNodeVec_t> BrokenMRCAAttachVec;

class SumTreeData {
    public:
    std::unordered_map<std::string, const SumTreeNode_t *> name_to_node;
    std::unordered_map<std::string, BrokenMRCAAttachVec> broken_taxa;
};
using TaxTree_t = otc::RootedTree<RTTaxNodeData, RTreeNoData>;
using SummaryTree_t = otc::RootedTree<SumTreeNodeData, SumTreeData>;

struct SourceTreeId {
    std::string tree_id;
    std::string study_id;
    std::string git_sha;
};

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
};

class TreesToServe {
        std::list< SummaryTreeAnnotation> annotation_list;
        std::list<std::unique_ptr<SummaryTree_t> > tree_list;
        std::map<std::string, const SummaryTree_t *> id_to_tree;
        std::map<std::string, const SummaryTreeAnnotation *> id_to_annotations;
        std::string default_synth_id;
        const RichTaxonomy * taxonomy_ptr = nullptr;
        std::map<std::string, const std::string *> stored_strings;
        std::list<std::string> stored_strings_list;
        const RichTaxTree * taxonomy_tree = nullptr;
    public:
        const std::string * getStoredString(const std::string & k) {
            const std::string * v = stored_strings[k];
            if (v == nullptr) {
                stored_strings_list.push_back(k);
                v = &(stored_strings_list.back());
                stored_strings[k] = v;
            }
            return v;
        }

        void setTaxonomy(const RichTaxonomy &taxonomy) {
            assert(taxonomy_ptr == nullptr);
            taxonomy_ptr = &taxonomy;
            taxonomy_tree = &(taxonomy.getTaxTree());
        }
        const RichTaxonomy & getTaxonomy() const {
            assert(taxonomy_ptr != nullptr);
            return *taxonomy_ptr;
        }
        void fillOttIdSet(const std::bitset<32> & flags, OttIdSet & ott_id_set) {
            ott_id_set.clear();
            for (const auto nd : iter_node_const(*taxonomy_tree)) {
                const auto & tax_record_flags = nd->get_data().get_flags();
                auto intersection = flags & tax_record_flags;
                if (!intersection.any()) {
                    ott_id_set.insert(nd->get_ott_id());
                }
            }
        }
        std::pair<SummaryTree_t &, SummaryTreeAnnotation &> getNewTreeAndAnnotations(const std::string & configfilename,
                                                                                     const std::string & filename) {
            
            OttIdSet ott_id_set;
            auto cleaning_flags = cleaning_flags_from_config_file(configfilename);
            fillOttIdSet(cleaning_flags, ott_id_set);

            assert(taxonomy_ptr != nullptr);
            ParsingRules parsingRules;
            parsingRules.ott_id_validator = &ott_id_set;
            parsingRules.include_internal_nodes_in_des_id_sets = true;
            parsingRules.set_ott_idForInternals = true;
            parsingRules.require_ott_ids = true;
            parsingRules.set_ott_ids = true;
            std::ifstream inp;
            if (!open_utf8_file(filename, inp)) {
                throw OTCError("Could not open \"" + filename + "\"");
            }
            LOG(INFO) << "reading \"" << filename << "\"...";
            ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
            FilePosStruct pos(filenamePtr);
            std::unique_ptr<SummaryTree_t> nt = read_next_newick<SummaryTree_t>(inp, pos, parsingRules);
            set_traversal_entry_exit(*nt);
            tree_list.push_back(move(nt));
            annotation_list.push_back(SummaryTreeAnnotation());
            return {*(tree_list.back()), annotation_list.back()};
        }
        void registerLastTreeAndAnnotations() {
            const SummaryTreeAnnotation & sta = annotation_list.back();
            const SummaryTree_t & tree = *(tree_list.back());
            default_synth_id = sta.synth_id; // @ TODO need a better system for deciding the default synth ID.
            id_to_tree[sta.synth_id] = &tree;
            id_to_annotations[sta.synth_id] = &sta;
        }
        void freeLastTreeAndAnnotations() {
            tree_list.back()->clear();
            tree_list.pop_back();
            annotation_list.pop_back();
        }

        const SummaryTreeAnnotation * getAnnotations(std::string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_annotations.find(key);
            return mit == id_to_annotations.end() ? nullptr : mit->second;
        }

        const SummaryTree_t * getSummaryTree(std::string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_tree.find(key);
            return mit == id_to_tree.end() ? nullptr : mit->second;
        }
        std::size_t getNumTrees() const {
            return id_to_tree.size();
        }
};

enum NodeNameStyle {
    NNS_NAME_ONLY = 0,
    NNS_ID_ONLY = 1,
    NNS_NAME_AND_ID = 2
};

void about_ws_method(const TreesToServe &tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     bool include_sources,
                     std::string & response_str,
                     int & status_code);
void tax_about_ws_method(const TreesToServe &tts,
                         std::string & response_str,
                         int & status_code);

void node_info_ws_method(const TreesToServe & tts,
                         const SummaryTree_t * tree_ptr,
                         const SummaryTreeAnnotation * sta,
                         const std::string & node_id,
                         bool include_lineage,
                         std::string & response_str,
                         int & status_code);

void mrca_ws_method(const TreesToServe & tts,
                    const SummaryTree_t * tree_ptr,
                    const SummaryTreeAnnotation * sta,
                    const std::vector<std::string> & node_id_vec,
                    std::string & response_str,
                    int & status_code);

void induced_subtree_ws_method(const TreesToServe & tts,
                               const SummaryTree_t * tree_ptr,
                               const SummaryTreeAnnotation * sta,
                               const std::vector<std::string> & node_id_vec,
                               NodeNameStyle label_format, 
                               std::string & response_str,
                               int & status_code);

void newick_subtree_ws_method(const TreesToServe & tts,
                              const SummaryTree_t * tree_ptr,
                              const SummaryTreeAnnotation * sta,
                              const std::string & node_id,
                              NodeNameStyle label_format, 
                              int height_limit,
                              std::string & response_str,
                              int & status_code);
void arguson_subtree_ws_method(const TreesToServe & tts,
                               const SummaryTree_t * tree_ptr,
                               const SummaryTreeAnnotation * sta,
                               const std::string & node_id,
                               int height_limit,
                               std::string & response_str,
                               int & status_code);
void taxon_info_ws_method(const TreesToServe & tts,
                          const RTRichTaxNode * taxon_node,
                          bool include_lineage,
                          bool include_children,
                          bool include_terminal_descendants,
                          std::string & response_str,
                          int & status_code);

bool read_trees(const boost::filesystem::path & dirname, TreesToServe & tts);

void from_json(const nlohmann::json &j, SummaryTreeAnnotation & sta);
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

inline OttId extract_unsigned_long(const nlohmann::json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_number()) {
        return dc_el->get<OttId>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a non-negative integer.\n";
}

inline void to_json(nlohmann::json &j, const SourceTreeId & sti) {
    j = nlohmann::json{{"git_sha", sti.git_sha}, {"study_id", sti.study_id}, {"tree_id", sti.tree_id}};
}

inline void from_json(const nlohmann::json &j, SourceTreeId & sti) {
    sti.git_sha = extract_string(j, "git_sha");
    sti.study_id = extract_string(j, "study_id");
    sti.tree_id = extract_string(j, "tree_id");
}


} // namespace otc
#endif