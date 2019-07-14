#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"
#include "ws/trees_to_serve.h"

using namespace std;
namespace fs = boost::filesystem;
using json = nlohmann::json;
typedef std::set<fs::path> fp_set;
typedef std::pair<bool, fp_set > bool_fp_set; 

bool_fp_set get_subdirs(const fs::path & dirname) {
    if (!fs::is_directory(dirname)) {
        LOG(ERROR) << "\"" << dirname << "\" is not a directory.";
        return bool_fp_set(false, fp_set());
    }
    fp_set fps;
    for (auto sub : fs::directory_iterator(dirname)) {
        if (fs::is_directory(sub)) {
            fps.insert(sub);
        }
    }
    return bool_fp_set(true, fps);
}

namespace otc {


void from_json(const nlohmann::json &j, SummaryTreeAnnotation & sta) {
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
    sta.root_ott_id = extract_ott_id(j, "root_ott_id");
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


bool read_tree_and_annotations(const fs::path & configpath,
                               const fs::path & treepath,
                               const fs::path & annotationspath,
                               const fs::path & brokentaxapath,
                               TreesToServe & tts);

// Globals. TODO: lock if we read twice
fp_set checked_dirs;
fp_set known_tree_dirs;

bool read_trees(const fs::path & dirname, TreesToServe & tts) {
    auto sdr = get_subdirs(dirname);
    if (!sdr.first) {
        return false;
    }
    const auto & subdir_set = sdr.second;
    for (auto p : subdir_set) {
        if (!contains(checked_dirs, p)) {
            checked_dirs.insert(p);
            fs::path configpath = p;
            configpath /= "config";
            fs::path treepath = p;
            treepath /= "labelled_supertree";
            treepath /= "labelled_supertree.tre";
            fs::path brokentaxapath = p;
            brokentaxapath /= "labelled_supertree";
            brokentaxapath /= "broken_taxa.json";
            fs::path annotationspath = p;
            annotationspath /= "annotated_supertree";
            annotationspath /= "annotations.json";
            bool was_tree_par = false;
            try {
                if (fs::is_regular_file(treepath)
                    && fs::is_regular_file(annotationspath)
                    && fs::is_regular_file(configpath)) {
                    if (read_tree_and_annotations(configpath, treepath, annotationspath, brokentaxapath, tts)) {
                        known_tree_dirs.insert(p);
                        was_tree_par = true;
                    }
                }
            } catch (const std::exception & x) {
                LOG(WARNING) << "Exception while reading summary tree directory:\n   ";
                LOG(WARNING) << x.what() << '\n';
            }
            if (!was_tree_par) {
                LOG(WARNING) << "Rejected \"" << p << "\" due to lack of " << treepath << " or lack of " << annotationspath << " or parsing error.\n";
            }
        }
    }
    tts.final_tree_added();
    return true;
}

#if defined(REPORT_MEMORY_USAGE)

template<>
inline std::size_t calc_memory_used(const RTRichTaxTreeData &d, MemoryBookkeeper &mb) {
    std::size_t sz_el_size = sizeof(OttId) + sizeof(const RTRichTaxNode *);
    std::size_t nm_sz = calc_memory_used_by_map_eqsize(d.ncbi_id_map, sz_el_size, mb);
    std::size_t gm_sz = calc_memory_used_by_map_eqsize(d.gbif_id_map, sz_el_size, mb);
    std::size_t wm_sz = calc_memory_used_by_map_eqsize(d.worms_id_map, sz_el_size, mb);
    std::size_t fm_sz = calc_memory_used_by_map_eqsize(d.if_id_map, sz_el_size, mb);
    std::size_t im_sz = calc_memory_used_by_map_eqsize(d.irmng_id_map, sz_el_size, mb);
    std::size_t f2j_sz = calc_memory_used_by_map_simple(d.flags2json, mb);
    std::size_t in_sz = calc_memory_used_by_map_simple(d.id_to_node, mb);
    std::size_t nn_sz = calc_memory_used_by_map_simple(d.name_to_node, mb);
    std::size_t nutn_sz = calc_memory_used_by_map_simple(d.non_unique_taxon_names, mb);
    std::size_t htn_sz = 0;
    for (auto el : d.homonym_to_node) {
        htn_sz += sizeof(std::string_view);
        htn_sz += calc_memory_used_by_vector_eqsize(el.second, sizeof(const RTRichTaxNode *), mb);
    }
    mb["taxonomy data ncbi map"] += nm_sz;
    mb["taxonomy data gbif map"] += gm_sz;
    mb["taxonomy data worms map"] += wm_sz;
    mb["taxonomy data indexfungorum map"] += fm_sz;
    mb["taxonomy data irmng map"] += im_sz;
    mb["taxonomy data flags2json"] += f2j_sz;
    mb["taxonomy data id_to_node"] += in_sz;
    mb["taxonomy data name_to_node"] += nn_sz;
    mb["taxonomy data non_unique_taxon_names"] += nutn_sz;
    mb["taxonomy data homonym_to_node"] += htn_sz;
    return nm_sz + gm_sz + wm_sz + fm_sz + im_sz + f2j_sz + in_sz + nn_sz + nutn_sz + htn_sz;
}

template<>
inline std::size_t calc_memory_used(const TaxonomicJuniorSynonym &d, MemoryBookkeeper &mb) {
    return calc_memory_used(d.name, mb) + calc_memory_used(d.source_string, mb) + sizeof(char *);
}

template<>
inline std::size_t calc_memory_used(const std::list<TaxonomicJuniorSynonym> &d, MemoryBookkeeper &mb) {
    std::size_t total = sizeof(size_t) * (2 + 2*d.size()) ; // start and capacity
    for (auto i = d.begin(); i != d.end(); ++i) {
        total += calc_memory_used(*i, mb);
    }
    return total;
}
template<>
inline std::size_t calc_memory_used(const RTRichTaxNodeData &rtn, MemoryBookkeeper &mb) {
    size_t total = 0;
    size_t x = calc_memory_used_by_vector_eqsize(rtn.junior_synonyms, sizeof(char *), mb);
    mb["taxonomy node data junior_synonyms"] += x; total += x;
    x = 2*sizeof(std::uint32_t);
    mb["taxonomy node data traversal"] += x; total += x;
    x = sizeof(TaxonomicRank);
    mb["taxonomy node data rank"] += x; total += x;
    x = sizeof(int32_t) + sizeof(std::bitset<32>);
    mb["taxonomy node data flags"] += x; total += x;
    x = calc_memory_used(rtn.source_info, mb);
    mb["taxonomy node data source_info"] += x; total += x;
    x = sizeof(std::string_view);
    mb["taxonomy node data nonunique name"] += x; total += x;
    return total;
}

template<>
inline std::size_t calc_memory_used(const RichTaxonomy &rt, MemoryBookkeeper &mb) {
    const auto & tax_tree = rt.get_tax_tree();
    std::size_t ttsz = calc_memory_used_by_tree(tax_tree, mb);
    const auto & syn_list = rt.get_synonyms_list();
    std::size_t slsz = calc_memory_used(syn_list, mb);
    const auto & suppress_trns_set = rt.get_ids_to_suppress_from_tnrs();
    std::size_t stssz = calc_memory_used(suppress_trns_set, mb);
    mb["taxonomy tree"] += ttsz;
    mb["taxonomy synonyms list"] += slsz;
    mb["taxonomy set of ids to suppress from tnrs"] += stssz;
    return ttsz + slsz + stssz;
}

#endif

bool read_tree_and_annotations(const fs::path & config_path,
                               const fs::path & tree_path,
                               const fs::path & annotations_path,
                               const fs::path & brokentaxa_path,
                               TreesToServe & tts) {
    std::string annot_str = annotations_path.native();
    std::ifstream annotations_stream(annot_str.c_str());
    json annotations_obj;
    try {
        annotations_stream >> annotations_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << annotations_path << "\" as JSON.\n";
        throw;
    }
    std::string bt_str = brokentaxa_path.native();
    std::ifstream brokentaxa_stream(bt_str.c_str());
    json brokentaxa_obj;
    try {
        brokentaxa_stream >> brokentaxa_obj;
    } catch (...) {
        LOG(WARNING) << "Could not read \"" << brokentaxa_path << "\" as JSON.\n";
        throw;
    }
    auto locked_taxonomy = tts.get_readable_taxonomy();
    const auto & taxonomy = locked_taxonomy.first;
#   if defined(REPORT_MEMORY_USAGE)
        MemoryBookkeeper tax_mem_b;
        std::size_t tree_mem = 0;
        auto tax_mem = calc_memory_used(taxonomy, tax_mem_b);
        write_memory_bookkeeping(LOG(INFO), tax_mem_b, "taxonomy", tax_mem);
#   endif
    auto tree_and_ann = tts.get_new_tree_and_annotations(config_path.native(), tree_path.native());
    try {
        SummaryTree_t & tree = tree_and_ann.first;
        SummaryTreeAnnotation & sta = tree_and_ann.second;
        sta = annotations_obj;
        json tref;
        tref["taxonomy"] = taxonomy.get_version();
        sta.full_source_id_map_json[taxonomy.get_version()] = tref;

        // read the node annotations and add them to the tree.
        auto node_obj = extract_obj(annotations_obj, "nodes");
        auto & sum_tree_data = tree.get_data();
        //const auto & n2n = sum_tree_data.name_to_node;
        for (json::const_iterator nit = node_obj.begin(); nit != node_obj.end(); ++nit) {
            string k = nit.key();
            bool was_broken = false;
            const SumTreeNode_t * stn = find_node_by_id_str(tree, k, was_broken);
            //auto stnit = n2n.find(k);
            if (stn == nullptr) {
                throw OTCError() << "Node " << k << " from annotations not found in tree.";
            }
            //const SumTreeNode_t * stn = stnit->second;
            SumTreeNode_t * mstn = const_cast<SumTreeNode_t *>(stn);
            SumTreeNodeData & mstnd = mstn->get_data();
            const auto & supportj = nit.value();
#           if defined(JOINT_MAPPING_VEC)
                vec_src_node_ids tmpv;
#           endif
            for (json::const_iterator sbit = supportj.begin(); sbit != supportj.end(); ++sbit) {
                const auto & sbk = sbit.key();
                const auto & sbv = sbit.value();
                if (sbk == "supported_by") {            
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::SUPPORTED_BY_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.supported_by = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "terminal") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::TERMINAL_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.terminal = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "conflicts_with") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::CONFLICTS_WITH_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.conflicts_with = extract_node_id_vec(tts, sbv)
#                   endif
                } else if (sbk == "partial_path_of") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::PARTIAL_PATH_OF_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.partial_path_of = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "resolves") {
#                   if defined(JOINT_MAPPING_VEC)
                        auto x = extract_node_id_vec(tts, sbv, SourceEdgeMappingType::RESOLVES_MAPPING);
                        tmpv.insert(end(tmpv), x.begin(), x.end());
#                   else
                        mstnd.resolves = extract_node_id_vec(tts, sbv);
#                   endif
                } else if (sbk == "was_uncontested") {
                    if (sbv.is_boolean()) {
                        mstnd.was_uncontested =  sbv.get<bool>();
                    } else {
                        throw OTCError() << "Expected was_uncontested to be a boolean.";
                    }
                } else {
                    if (sbk != "was_constrained") {
                        throw OTCError() << "Unrecognized annotations key " << sbit.key();
                    }
                }
            }
#           if defined(JOINT_MAPPING_VEC)
                mstnd.source_edge_mappings.clear();
                std::swap(mstnd.source_edge_mappings, tmpv);
                //LOG(INFO) << "mstnd.source_edge_mappings size = " << mstnd.source_edge_mappings.size() << " for k = " << k;
#           endif

        }
        auto & tree_broken_taxa = sum_tree_data.broken_taxa;
        // read the info from the broken taxa file
        if (brokentaxa_obj.count("non_monophyletic_taxa")
            && (!brokentaxa_obj["non_monophyletic_taxa"].is_null())) {
            auto & nmt_obj = extract_obj(brokentaxa_obj, "non_monophyletic_taxa");
            for (json::const_iterator btit = nmt_obj.begin(); btit != nmt_obj.end(); ++btit) {
                string broken_ott = btit.key();
                auto & dest_obj = btit.value();
                string mrca_id = extract_string(dest_obj, "mrca");
                auto & attach_obj = extract_obj(dest_obj, "attachment_points");
                list<string> attach_id_list;
                for (json::const_iterator ai_it = attach_obj.begin(); ai_it != attach_obj.end(); ++ai_it) {
                    attach_id_list.push_back(ai_it.key());
                }
                bool was_broken = false;
                const SumTreeNode_t * mrca_nd = find_node_by_id_str(tree, mrca_id, was_broken);
                vector<const SumTreeNode_t *> avec;
                avec.reserve(attach_id_list.size());
                for (auto attach_id : attach_id_list) {
                    auto anptr = find_node_by_id_str(tree, attach_id, was_broken);
                    assert(anptr != nullptr);
                    avec.push_back(anptr);
                }
                tree_broken_taxa[broken_ott] = BrokenMRCAAttachVec(mrca_nd, avec);
            }
        }
        sta.initialized = true;
        tts.register_last_tree_and_annotations();
#       if defined(REPORT_MEMORY_USAGE)
            MemoryBookkeeper tree_mem_b;
            tree_mem += calc_memory_used_by_tree(tree, tree_mem_b);
            write_memory_bookkeeping(LOG(INFO), tree_mem_b, "tree", tree_mem);
            LOG(INFO) << "tax + tree memory = " << tax_mem << " + " << tree_mem << " = " << tax_mem + tree_mem;
#       endif
    } catch (...) {
        tts.free_last_tree_and_annotations();
        throw;
    }
    return true;
}

}// namespace otc
