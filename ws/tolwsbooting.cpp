#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"

using namespace std;
namespace fs = boost::filesystem;
using json = nlohmann::json;
typedef std::set<fs::path> fp_set;
typedef std::pair<bool, fp_set > bool_fp_set; 

bool_fp_set get_subdirs(const fs::path & dirname) {
    if (!fs::is_directory(dirname)) {
        LOG(ERROR) << "\"" << dirname << "\" is not a directory.\n";
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
    const RichTaxonomy & taxonomy = tts.get_taxonomy();
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
        const auto & n2n = sum_tree_data.name_to_node;
        for (json::const_iterator nit = node_obj.begin(); nit != node_obj.end(); ++nit) {
            string k = nit.key();
            auto stnit = n2n.find(k);
            if (stnit == n2n.end()) {
                throw OTCError() << "Node " << k << " from annotations not found in tree.";
            }
            const SumTreeNode_t * stn = stnit->second;
            SumTreeNode_t * mstn = const_cast<SumTreeNode_t *>(stn);
            SumTreeNodeData & mstnd = mstn->get_data();
            const auto & supportj = nit.value();
            for (json::const_iterator sbit = supportj.begin(); sbit != supportj.end(); ++sbit) {
                const auto & sbk = sbit.key();
                const auto & sbv = sbit.value();
                if (sbk == "conflicts_with") {
                    mstnd.conflicts_with = extract_node_id_vec(tts, sbv);
                } else if (sbk == "supported_by") {
                    mstnd.supported_by = extract_node_id_vec(tts, sbv);
                } else if (sbk == "terminal") {
                    mstnd.terminal = extract_node_id_vec(tts, sbv);
                } else if (sbk == "partial_path_of") {
                    mstnd.partial_path_of = extract_node_id_vec(tts, sbv);
                } else if (sbk == "resolves") {
                    mstnd.resolves = extract_node_id_vec(tts, sbv);
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
                const SumTreeNode_t * mrca_nd = n2n.at(mrca_id);
                vector<const SumTreeNode_t *> avec;
                avec.reserve(attach_id_list.size());
                for (auto attach_id : attach_id_list) {
                    avec.push_back(n2n.at(attach_id));
                }
                tree_broken_taxa[broken_ott] = BrokenMRCAAttachVec(mrca_nd, avec);
            }
        }
        sta.initialized = true;
        tts.register_last_tree_and_annotations();
    } catch (...) {
        tts.free_last_tree_and_annotations();
        throw;
    }
    return true;
}

}// namespace otc