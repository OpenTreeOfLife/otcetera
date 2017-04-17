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

class OTCWebError : public std::exception {
protected:
    int status_code_ = 500;
    std::string message;
public:
    int status_code() const {return status_code_;}
    const char * what() const noexcept {
        return message.c_str();
    }
    template <typename T> OTCWebError& operator<<(const T&);
    void prepend(const std::string& s) {
        message = s + message;
    }
    OTCWebError() noexcept {}
    OTCWebError(int c) noexcept :status_code_(c) {}
    OTCWebError(const std::string & msg) noexcept :message(msg) {}
    OTCWebError(int c, const std::string & msg) noexcept :status_code_(c), message(msg) {}
};

template <typename T>
OTCWebError& OTCWebError::operator<<(const T& t) {
  std::ostringstream oss;
  oss << message << t;
  message = oss.str();
  return *this;
}

inline OTCWebError OTCBadRequest() {return OTCWebError(400);}
inline OTCWebError OTCBadRequest(const std::string& m) {return OTCWebError(400,m);}

class SumTreeNodeData {
    public:
    std::uint32_t trav_enter = UINT32_MAX;
    std::uint32_t trav_exit = UINT32_MAX;
    vec_src_node_ids supported_by;
    vec_src_node_ids conflicts_with;
    vec_src_node_ids resolves;
    vec_src_node_ids partial_path_of;
    vec_src_node_ids terminal;
    bool was_uncontested = false;
    uint32_t num_tips = 0;
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
        OttIdSet suppressed_from_tree; // IDs not included in tree
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
        const std::string * get_stored_string(const std::string & k) {
            const std::string * v = stored_strings[k];
            if (v == nullptr) {
                stored_strings_list.push_back(k);
                v = &(stored_strings_list.back());
                stored_strings[k] = v;
            }
            return v;
        }

        void set_taxonomy(const RichTaxonomy &taxonomy) {
            assert(taxonomy_ptr == nullptr);
            taxonomy_ptr = &taxonomy;
            taxonomy_tree = &(taxonomy.getTaxTree());
        }
        const RichTaxonomy & get_taxonomy() const {
            assert(taxonomy_ptr != nullptr);
            return *taxonomy_ptr;
        }
        void fill_ott_id_set(const std::bitset<32> & flags,
                              OttIdSet & ott_id_set,
                              OttIdSet & suppressed_from_tree) {
            ott_id_set.clear();
            for (const auto nd : iter_node_const(*taxonomy_tree)) {
                const auto & tax_record_flags = nd->get_data().get_flags();
                auto intersection = flags & tax_record_flags;
                const auto ott_id = nd->get_ott_id();
                if (!intersection.any()) {
                    ott_id_set.insert(ott_id);
                } else {
                    suppressed_from_tree.insert(ott_id);
                }
            }
        }
        typedef std::pair<SummaryTree_t &, SummaryTreeAnnotation &> SumTreeInitPair;
        SumTreeInitPair get_new_tree_and_annotations(const std::string & configfilename,
                                                     const std::string & filename) {
            
            OttIdSet ott_id_set, suppressed_id_set;
            auto cleaning_flags = cleaning_flags_from_config_file(configfilename);
            fill_ott_id_set(cleaning_flags, ott_id_set, suppressed_id_set);

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
            index_nodes_by_name(*nt);
            set_traversal_entry_exit_and_num_tips(*nt);
            tree_list.push_back(move(nt));
            annotation_list.push_back(SummaryTreeAnnotation());
            auto & sta = annotation_list.back();
            sta.suppressed_from_tree = suppressed_id_set;
            return {*(tree_list.back()), sta};
        }
        void register_last_tree_and_annotations() {
            const SummaryTreeAnnotation & sta = annotation_list.back();
            const SummaryTree_t & tree = *(tree_list.back());
            default_synth_id = sta.synth_id; // @ TODO need a better system for deciding the default synth ID.
            id_to_tree[sta.synth_id] = &tree;
            id_to_annotations[sta.synth_id] = &sta;
        }
        void free_last_tree_and_annotations() {
            tree_list.back()->clear();
            tree_list.pop_back();
            annotation_list.pop_back();
        }

        const SummaryTreeAnnotation * get_annotations(std::string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_annotations.find(key);
            return mit == id_to_annotations.end() ? nullptr : mit->second;
        }

        const SummaryTree_t * get_summary_tree(std::string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_tree.find(key);
            return mit == id_to_tree.end() ? nullptr : mit->second;
        }
        std::size_t get_num_trees() const {
            return id_to_tree.size();
        }
        void final_tree_added() {
            if (get_num_trees() == 1) {
                const auto & annot = annotation_list.back();
                const auto & sft = annot.suppressed_from_tree;
                assert(taxonomy_ptr != nullptr);
                // TODO: make taxonomy_ptr non-const or make the relevant field mutable. Probably the former.
                auto nct = const_cast<RichTaxonomy *>(taxonomy_ptr);
                nct->set_ids_suppressed_from_summary_tree_alias(&sft);
            }
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

std::string induced_subtree_ws_method(const TreesToServe & tts,
				      const SummaryTree_t * tree_ptr,
				      const SummaryTreeAnnotation * sta,
				      const std::vector<std::string> & node_id_vec,
				      NodeNameStyle label_format);

std::string newick_subtree_ws_method(const TreesToServe & tts,
				     const SummaryTree_t * tree_ptr,
				     const SummaryTreeAnnotation * sta,
				     const std::string & node_id,
				     NodeNameStyle label_format, 
				     int height_limit);

std::string arguson_subtree_ws_method(const TreesToServe & tts,
				      const SummaryTree_t * tree_ptr,
				      const SummaryTreeAnnotation * sta,
				      const std::string & node_id,
				      int height_limit);

std::string taxon_info_ws_method(const TreesToServe & tts,
				 const RTRichTaxNode * taxon_node,
				 bool include_lineage,
				 bool include_children,
				 bool include_terminal_descendants);

std::string taxonomy_mrca_ws_method(const TreesToServe & tts, const OttIdSet & ott_id_set);

void taxon_subtree_ws_method(const TreesToServe & tts,
                             const RTRichTaxNode * taxon_node,
                             NodeNameStyle label_format, 
                             std::string & response_str,
                             int & status_code);

std::string conflict_ws_method(const SummaryTree_t & summary,
			       const RichTaxonomy & taxonomy,
			       const std::string& tree1s,
			       const std::string& tree2s);

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
