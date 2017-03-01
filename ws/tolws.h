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
void setTravesalEntryExit(T & tree) {
    auto & td = tree.getData();
    auto & m = td.name2node;
    long ind = 0;
    for (auto nd : iter_pre(tree)) {
        m[nd->getName()] = nd;
        nd->getData().trav_enter = ind++;
    }
    for (auto pnd : iter_post(tree)) {
        auto fc = pnd->getLastChild();
        auto & d = pnd->getData();
        if (fc == nullptr) {
            d.trav_exit = d.trav_enter;
            d.num_tips = 1;
        } else {
            d.trav_exit = fc->getData().trav_exit;
            d.num_tips = 0;
            for (auto c : iter_child_const(*pnd)) {
                d.num_tips += c->getData().num_tips;
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
    std::unordered_map<std::string, const SumTreeNode_t *> name2node;
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
                const auto & tax_record_flags = nd->getData().getFlags();
                auto intersection = flags & tax_record_flags;
                if (!intersection.any()) {
                    ott_id_set.insert(nd->getOttId());
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
            parsingRules.ottIdValidator = &ott_id_set;
            parsingRules.includeInternalNodesInDesIdSets = true;
            parsingRules.setOttIdForInternals = true;
            parsingRules.requireOttIds = true;
            parsingRules.setOttIds = true;
            std::ifstream inp;
            if (!openUTF8File(filename, inp)) {
                throw OTCError("Could not open \"" + filename + "\"");
            }
            LOG(INFO) << "reading \"" << filename << "\"...";
            ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
            FilePosStruct pos(filenamePtr);
            std::unique_ptr<SummaryTree_t> nt = readNextNewick<SummaryTree_t>(inp, pos, parsingRules);
            setTravesalEntryExit(*nt);
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

inline void from_json(const nlohmann::json &j, SummaryTreeAnnotation & sta) {
    sta.date_completed = extract_string(j, "date_completed");
    sta.filtered_flags = extract_string(j, "filtered_flags");
    auto splitff = split_string(sta.filtered_flags, ',');
    sta.filtered_flags_vec.assign(splitff.begin(), splitff.end());
    // generated_by gets converted back to a string
    auto gb_el = j.find("generated_by");
    if (gb_el == j.end()) {
        throw OTCError() << "Missing generated_by field.\n";
    }
    sta.generated_by = gb_el->dump();
    sta.num_leaves_in_exemplified_taxonomy = extract_unsigned_long(j, "num_leaves_in_exemplified_taxonomy");
    sta.num_source_studies = extract_unsigned_long(j, "num_source_studies");
    sta.num_source_trees = extract_unsigned_long(j, "num_source_trees");
    sta.num_tips = extract_unsigned_long(j, "num_tips");
    sta.root_ott_id = extract_unsigned_long(j, "root_ott_id");
    sta.root_taxon_name = extract_string(j, "root_taxon_name");
    sta.synth_id = extract_string(j, "synth_id");
    sta.taxonomy_version = extract_string(j, "taxonomy_version");
    sta.tree_id = extract_string(j, "tree_id");
    auto sim_el = j.find("source_id_map");
    if (sim_el == j.end()) {
        throw OTCError() << "Missing source_id_map field.\n";
    }
    if (!sim_el->is_object()) {
        throw OTCError() << "Expected \"source_id_map\" field to be an object.\n";
    }
    try {
        for (nlohmann::json::const_iterator sim_it = sim_el->begin(); sim_it != sim_el->end(); ++sim_it) {
            sta.source_id_map[sim_it.key()] = sim_it.value(); 
        }
    } catch (OTCError & x) {
        throw OTCError() << "Error reading source_id_map field: " << x.what();
    }
    sta.full_source_id_map_json = *sim_el;
    auto s_el = j.find("sources");
    if (s_el == j.end()) {
        throw OTCError() << "Missing sources field.\n";
    }
    if (!s_el->is_array()) {
        throw OTCError() << "Expected \"sources\" field to be an array.\n";
    }
    sta.sources.resize(s_el->size());
    for (auto i = 0U; i < s_el->size(); ++i) {
        try {
            sta.sources[i] = s_el->at(i).get<std::string>();
        } catch (OTCError & x) {
            throw OTCError() << "Error expected each element of the sources array to be a string: " << x.what();
        }
    }
}


} // namespace otc
#endif