#include <memory>
#include <thread>
#include <cstdlib>
#include <restbed>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <regex>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem/operations.hpp>

#include "otc/util.h"
#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

#include "json.hpp"

using namespace std;
using namespace restbed;
using namespace otc;
namespace fs = boost::filesystem;
using json = nlohmann::json;

typedef std::set<fs::path> fp_set;
typedef std::pair<bool, fp_set > bool_fp_set; 
typedef std::pair<const std::string *, const std::string *> src_node_id;
typedef std::vector<src_node_id> vec_src_node_ids;

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
typedef vector<const SumTreeNode_t *> SumTreeNodeVec_t;
typedef pair<const SumTreeNode_t *, SumTreeNodeVec_t> BrokenMRCAAttachVec;
class SumTreeData {
    public:
    unordered_map<string, const SumTreeNode_t *> name2node;
    unordered_map<string, BrokenMRCAAttachVec> broken_taxa;
};

using TaxTree_t = RootedTree<RTTaxNodeData, RTreeNoData>;
using SummaryTree_t = RootedTree<SumTreeNodeData, SumTreeData>;
namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;


variables_map parse_cmd_line(int argc, char* argv[]) { 
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description output("Server options");
    output.add_options()
        ("tree-dir,D",value<string>(),"Filepath to directory that will hold synthetic tree output")
        ("port,P",value<long>(),"Port to bind to.")
        ("num-threads,t",value<long>(),"number of threads")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-tol-ws <taxonomy-dir> -D<dir> [OPTIONS]\n"
                                                    "Load taxonomy, check dir for synth tree outputs, and serve",
                                                    visible, invisible, p);

    return vm;
}

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
        json full_source_id_map_json;

};

void from_json(const json &j, SummaryTreeAnnotation & sta);
void from_json(const json &j, SourceTreeId & sti);
void to_json(json &j, const SourceTreeId & sti);

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
            d.trav_exit = fc->getData().trav_enter;
            d.num_tips = 0;
            for (auto c : iter_child_const(*pnd)) {
                d.num_tips += c->getData().num_tips;
            }
        }
    }
}

class TreesToServe {
        list< SummaryTreeAnnotation> annotation_list;
        list< unique_ptr<SummaryTree_t> > tree_list;
        map<string, const SummaryTree_t *> id_to_tree;
        map<string, const SummaryTreeAnnotation *> id_to_annotations;
        string default_synth_id;
        const Taxonomy * taxonomy_ptr = nullptr;
        OttIdSet ott_id_set;
        unique_ptr<TaxTree_t> taxonomy_tree;
        map<string, const string *> stored_strings;
        list<string> stored_strings_list;

    public:
        const string * getStoredString(const string & k) {
            const string * v = stored_strings[k];
            if (v == nullptr) {
                stored_strings_list.push_back(k);
                v = &(stored_strings_list.back());
                stored_strings[k] = v;
            }
            return v;
        }

        void setTaxonomy(const Taxonomy &taxonomy) {
            assert(taxonomy_ptr == nullptr);
            taxonomy_ptr = &taxonomy;
            ott_id_set.clear();
            auto nodeNamer = [](const auto&){return string();};
            taxonomy_tree = taxonomy.getTree<TaxTree_t>(nodeNamer);
        }
        const Taxonomy & getTaxonomy() const {
            assert(taxonomy_ptr != nullptr);
            return *taxonomy_ptr;
        }
        void fillOttIdSet(const std::bitset<32> & flags) {
            ott_id_set.clear();
            for (const auto nd : iter_node_const(*taxonomy_tree)) {
                const auto & tax_record = nd->getData().taxonomy_line;
                auto intersection = flags & tax_record->flags;
                if (!intersection.any()) {
                    ott_id_set.insert(tax_record->id);
                }
            }
        }
        pair<SummaryTree_t &, SummaryTreeAnnotation &> getNewTreeAndAnnotations(const string & configfilename,
                                                                                const string & filename) {
            
            auto cleaning_flags = cleaning_flags_from_config_file(configfilename);
            fillOttIdSet(cleaning_flags);

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

        const SummaryTreeAnnotation * getAnnotations(string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_annotations.find(key);
            return mit == id_to_annotations.end() ? nullptr : mit->second;
        }

        const SummaryTree_t * getSummaryTree(string synth_id) const {
            const auto & key = synth_id.empty() ? default_synth_id : synth_id;
            auto mit = id_to_tree.find(key);
            return mit == id_to_tree.end() ? nullptr : mit->second;
        }
        size_t getNumTrees() const {
            return id_to_tree.size();
        }
};

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
    return true;
}

vec_src_node_ids extract_node_id_vec(TreesToServe & tts,
                                     const json & sbv) {
    list<src_node_id> lsni;
    for (json::const_iterator jit = sbv.begin(); jit != sbv.end(); ++jit) {
        const string * kp = tts.getStoredString(jit.key());
        const auto & v = jit.value();
        for (json::const_iterator vit = v.begin(); vit != v.end(); ++vit) {
            const string * vp = tts.getStoredString(*vit);
            lsni.push_back(src_node_id(kp, vp));
        } 
    }
    return vec_src_node_ids(lsni.begin(), lsni.end());
}


const json & extract_obj(const json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_object()) {
        json & k = const_cast<json &>(j);
        return k[field];
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}

std::string extract_string(const json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_string()) {
        return dc_el->get<std::string>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}

OttId extract_unsigned_long(const json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_number()) {
        return dc_el->get<OttId>();
    }
    throw OTCError() << "Expected \"" << field << "\" field to be a non-negative integer.\n";
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
    const Taxonomy & taxonomy = tts.getTaxonomy();
    auto tree_and_ann = tts.getNewTreeAndAnnotations(config_path.native(), tree_path.native());
    try {
        SummaryTree_t & tree = tree_and_ann.first;
        SummaryTreeAnnotation & sta = tree_and_ann.second;
        sta = annotations_obj;
        json tref;
        tref["taxonomy"] = taxonomy.version;
        sta.full_source_id_map_json[taxonomy.version] = tref;

        // read the node annotations and add them to the tree.
        auto node_obj = extract_obj(annotations_obj, "nodes");
        auto & sum_tree_data = tree.getData();
        const auto & n2n = sum_tree_data.name2node;
        for (json::const_iterator nit = node_obj.begin(); nit != node_obj.end(); ++nit) {
            string k = nit.key();
            auto stnit = n2n.find(k);
            if (stnit == n2n.end()) {
                throw OTCError() << "Node " << k << " from annotations not found in tree.";
            }
            const SumTreeNode_t * stn = stnit->second;
            SumTreeNode_t * mstn = const_cast<SumTreeNode_t *>(stn);
            SumTreeNodeData & mstnd = mstn->getData();
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
        tts.registerLastTreeAndAnnotations();
    } catch (...) {
        tts.freeLastTreeAndAnnotations();
        throw;
    }
    return true;
}

void from_json(const json &j, SummaryTreeAnnotation & sta) {
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
        for (json::const_iterator sim_it = sim_el->begin(); sim_it != sim_el->end(); ++sim_it) {
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

void to_json(json &j, const SourceTreeId & sti) {
    j = json{{"git_sha", sti.git_sha}, {"study_id", sti.study_id}, {"tree_id", sti.tree_id}};
}

void from_json(const json &j, SourceTreeId & sti) {
    sti.git_sha = extract_string(j, "git_sha");
    sti.study_id = extract_string(j, "study_id");
    sti.tree_id = extract_string(j, "tree_id");
}


TreesToServe tts;
template<typename T>
bool extract_from_request(const json & j, string opt_name, T & setting, string & response, int & status_code);

template<>
bool extract_from_request(const json & j, string opt_name, bool & setting, string & response, int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_boolean()) {
            setting = opt->get<bool>();
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be a boolean.\n";
        status_code = 400;
    }
    return false;
}

template<>
bool extract_from_request(const json & j, string opt_name, string & setting, string & response, int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_string()) {
            setting = opt->get<string>();
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be a string.\n";
        status_code = 400;
    }
    return false;
}

template<>
bool extract_from_request(const json & j, string opt_name, vector<string> & setting, string & response, int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_array()) {
            for (json::const_iterator sIt = opt->begin(); sIt != opt->end(); ++sIt) {
                if (sIt->is_string()) {
                    setting.push_back(sIt->get<string>());
                } else {
                    response = "Expecting ";
                    response += opt_name;
                    response += " to be an array of strings.\n";
                    status_code = 400;
                    return false;
                }   
            }
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be an array of strings.\n";
        status_code = 400;
    }
    return false;
}

void add_basic_node_info(const Taxonomy & taxonomy, const SumTreeNode_t & nd, json & noderepr) {
    noderepr["node_id"] = nd.getName();
    noderepr["num_tips"] = nd.getData().num_tips;
    if (nd.hasOttId()) {
        auto nd_id = nd.getOttId();
        const auto & nd_taxon = taxonomy.record_from_id(nd_id);
        json taxon;
        taxon["tax_sources"] = nd_taxon.sourceinfoAsVec();
        auto un = string(nd_taxon.uniqname);
        auto n = string(nd_taxon.name);
        taxon["name"] = n;
        taxon["uniqname"] = un;
        auto r = string(nd_taxon.rank);
        taxon["rank"] = r;
        taxon["ott_id"] = nd_taxon.id;
        noderepr["taxon"] = taxon;
    }
}

inline void add_str_to_str_or_vec_string(json & o, const string * first, const string * second) {
    const char * studyc = first->c_str();
    if (o.count(studyc)) {
        if (o[studyc].is_array()) {
            o[studyc].push_back(*second);
        } else {
            string prev = o[studyc].get<string>();
            json v = {prev, *second};
            o[studyc] = v;
        }
    } else {
        o[studyc] = *second;
    }
}

void add_support_info_vec(const char * tag,
                          const vec_src_node_ids & v,
                          json & noderepr,
                          set<string> & usedSrcIds,
                          const string * extra_src=nullptr,
                          const string * extra_node_id=nullptr) {
    json o;
    for (const auto & study_node_pair : v ) {
        usedSrcIds.insert(*study_node_pair.first);
        add_str_to_str_or_vec_string(o, study_node_pair.first, study_node_pair.second);
    }
    if (extra_src && extra_node_id) {
        add_str_to_str_or_vec_string(o, extra_src, extra_node_id);
    }
    noderepr[tag] = o;
}

void add_node_support_info(const TreesToServe & tts,
                           const SumTreeNode_t & nd,
                           json & noderepr,
                           set<string> & usedSrcIds) {
    const Taxonomy & taxonomy = tts.getTaxonomy();
    const auto & d = nd.getData();
    const string * extra_src = nullptr;
    const string * extra_node_id = nullptr;
    if (nd.hasOttId()) {
        extra_src = &(taxonomy.version);
        extra_node_id = &(nd.getName());
        usedSrcIds.insert(taxonomy.version);
    }
    if (extra_src != nullptr || !d.supported_by.empty()) {
        add_support_info_vec("supported_by", d.supported_by, noderepr, usedSrcIds, extra_src, extra_node_id);
    }
    if (!d.conflicts_with.empty()) {
        add_support_info_vec("conflicts_with", d.conflicts_with, noderepr, usedSrcIds);
    }
    if (!d.resolves.empty()) {
        add_support_info_vec("resolves", d.resolves, noderepr, usedSrcIds);
    }
    if (!d.partial_path_of.empty()) {
        add_support_info_vec("partial_path_of", d.partial_path_of, noderepr, usedSrcIds);
    }
    if (!d.terminal.empty()) {
        add_support_info_vec("terminal", d.terminal, noderepr, usedSrcIds);
    }
    if (d.was_uncontested) {
        noderepr["was_uncontested"] = true;
        noderepr["was_constrained"] = true;
    }

}

const SumTreeNode_t * find_mrca(const SumTreeNode_t *f, const SumTreeNode_t *s) {
    const auto * fdata = &(f->getData());
    const auto sec_ind = s->getData().trav_enter;
    while (sec_ind < fdata->trav_enter || sec_ind > fdata->trav_exit) {
        f = f->getParent();
        if (f == nullptr) {
            assert(false); 
            return nullptr;
        }
        fdata = &(f->getData());
    }
    return f;
}

const std::regex mrca_id_pattern("^mrca(ott\\d+)(ott\\d+)$");
const std::regex ott_id_pattern("^ott(\\d+)$");

const SumTreeNode_t * find_node_by_id_str(const SummaryTree_t & tree,
                                          const string & node_id,
                                          bool & was_broken) {
    was_broken = false;
    const auto & tree_data = tree.getData();
    auto n2nit = tree_data.name2node.find(node_id);
    if (n2nit != tree_data.name2node.end()) {
        return n2nit->second;
    }
    std::smatch matches;
    if (std::regex_match(node_id, matches, mrca_id_pattern)) {
        assert(matches.size() >= 2);
        std::string first_id = matches[1];
        std::string second_id = matches[2];
        bool bogus = false;
        auto fir_nd = find_node_by_id_str(tree, first_id, bogus);
        if (fir_nd == nullptr) {
            return nullptr;
        }
        auto sec_nd = find_node_by_id_str(tree, second_id, bogus);
        if (sec_nd == nullptr) {
            return nullptr;
        }
        return find_mrca(fir_nd, sec_nd);
    }
    if (std::regex_match(node_id, matches, ott_id_pattern)) {
        auto bt_it = tree_data.broken_taxa.find(node_id);
        if (bt_it == tree_data.broken_taxa.end()) {
            return nullptr;
        }
        // if found we return the MRCA pointer in the first slot of the pair.
        was_broken = true;
        return bt_it->second.first;
    }
    return nullptr;
}


void about_ws_method(const TreesToServe &tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     bool include_sources,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["date_created"] = sta->date_completed;
    response["num_source_trees"] = sta->num_source_trees;
    response["taxonomy_version"] = sta->taxonomy_version;
    response["filtered_flags"] = sta->filtered_flags_vec;
    response["synth_id"] = sta->synth_id;
    if (include_sources) {
        response["source_id_map"] = sta->full_source_id_map_json;
        response["sources"] = sta->sources;
    }
    json root;
    auto root_node = tree_ptr->getRoot();
    add_basic_node_info(taxonomy, *root_node, root);
    response["root"] = root;
    response_str = response.dump(1);
}

void node_info_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     bool include_lineage,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_node_by_id_str(*tree_ptr, node_id, was_broken);
    bool is_broken = false;
    if (focal == nullptr) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return;
    }
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    add_basic_node_info(taxonomy, *focal, response);
    set<string> usedSrcIds;
    add_node_support_info(tts, *focal, response, usedSrcIds);
    if (include_lineage) {
        json lineage_arr;
        const SumTreeNode_t * anc = focal->getParent();
        while (anc) {
            json ancj;
            add_basic_node_info(taxonomy, *anc, ancj);
            add_node_support_info(tts, *anc, ancj, usedSrcIds);
            lineage_arr.push_back(ancj);
            anc = anc->getParent();
        }
        response["lineage"] = lineage_arr;
    }
    // now write source_id_map
    json sim;
    for (auto srcTag : usedSrcIds) {
        json jt;
        if (srcTag == taxonomy.version) {
            jt["taxonomy"] = taxonomy.version;
        } else {
            const auto & simentry = sta->source_id_map.at(srcTag);
            jt = simentry;
        }
        sim[srcTag] = jt;   
    }
    response["source_id_map"] = sim;
    response_str = response.dump(1);
}

void mrca_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const vector<string> & node_id_vec,
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_node_by_id_str(*tree_ptr, node_id, was_broken);
        if (n == nullptr) {
            response_str = "node_id \"";
            response_str += node_id;
            response_str += "\" was not recognized.\n";
            status_code = 400;
            return;
        }
        if (first) {
            first = false;
            focal = n;
        } else {
            focal = find_mrca(focal, n);
            if (focal == nullptr) {
                break;
            }
        }
    }
    bool is_broken = false;
    if (focal == nullptr) {
        response_str = "MRCA of taxa was not found.\n";
        status_code = 400;
        return;
    }
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    json mrcaj;
    add_basic_node_info(taxonomy, *focal, mrcaj);
    set<string> usedSrcIds;
    add_node_support_info(tts, *focal, mrcaj, usedSrcIds);
    response["mrca"] = mrcaj;
    // now write source_id_map
    json sim;
    for (auto srcTag : usedSrcIds) {
        json jt;
        if (srcTag == taxonomy.version) {
            jt["taxonomy"] = taxonomy.version;
        } else {
            const auto & simentry = sta->source_id_map.at(srcTag);
            jt = simentry;
        }
        sim[srcTag] = jt;   
    }
    response["source_id_map"] = sim;
    response_str = response.dump(1);
}



void about_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        stringstream id;
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        try {
            if (!body.empty()) {
                parsedargs = json::parse(body);
            }
        } catch (...) {
            rbody = "Could not parse body of call as JSON.\n";
            status_code = 400;
        }
        bool include_sources = false;
        string synth_id;
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_source_list", include_sources, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "synth_id", synth_id, rbody, status_code);
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.getAnnotations(synth_id);
            const SummaryTree_t * treeptr = tts.getSummaryTree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                about_ws_method(tts, treeptr, sta, include_sources, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}


void node_info_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        stringstream id;
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        try {
            if (!body.empty()) {
                parsedargs = json::parse(body);
            }
        } catch (...) {
            rbody = "Could not parse body of call as JSON.\n";
            status_code = 400;
        }
        bool include_lineage = false;
        string synth_id;
        string node_id;
        if (status_code == OK) {
            extract_from_request(parsedargs, "include_lineage", include_lineage, rbody, status_code);
        }
        if (status_code == OK) {
            extract_from_request(parsedargs, "synth_id", synth_id, rbody, status_code);
        }
        if (status_code == OK) {
            if (!extract_from_request(parsedargs, "node_id", node_id, rbody, status_code)) {
                rbody = "node_id is required.\n";
                status_code = 400;
            }
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.getAnnotations(synth_id);
            const SummaryTree_t * treeptr = tts.getSummaryTree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                node_info_ws_method(tts, treeptr, sta, node_id, include_lineage, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}

void mrca_method_handler( const shared_ptr< Session > session ) {
    const auto request = session->get_request( );
    size_t content_length = request->get_header( "Content-Length", 0 );
    session->fetch( content_length, [ request ]( const shared_ptr< Session > session, const Bytes & body ) {
        stringstream id;
        json parsedargs;
        std::string rbody;
        int status_code = OK;
        try {
            if (!body.empty()) {
                parsedargs = json::parse(body);
            }
        } catch (...) {
            rbody = "Could not parse body of call as JSON.\n";
            status_code = 400;
        }
        string synth_id;
        vector<string> node_id_vec;
        if (status_code == OK) {
            extract_from_request(parsedargs, "synth_id", synth_id, rbody, status_code);
        }
        if (status_code == OK) {
            if (!extract_from_request(parsedargs, "node_ids", node_id_vec, rbody, status_code)) {
                rbody = "node_ids is required.\n";
                status_code = 400;
            }
        }
        if (status_code == OK) {
            const SummaryTreeAnnotation * sta = tts.getAnnotations(synth_id);
            const SummaryTree_t * treeptr = tts.getSummaryTree(synth_id);
            if (sta == nullptr || treeptr == nullptr) {
               rbody = "Did not recognize the synth_id.\n";
               status_code = 400;
            } else {
                mrca_ws_method(tts, treeptr, sta, node_id_vec, rbody, status_code);
            }
        }
        session->close( OK, rbody, { { "Content-Length", ::to_string( rbody.length( ) ) } } );
    });
}


int main( const int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
        variables_map args = parse_cmd_line(argc,argv);
        int num_threads = 4;
        int port_number = 1984;
        if (args.count("num-threads")) {
            port_number = args["num-threads"].as<int>();
        }
        if (args.count("port")) {
            port_number = args["port"].as<int>();
        }
        if (!args.count("tree-dir")) {
            cerr << "Expecting a tree-dir argument for a path to a directory of synth outputs.\n";
            return 1;
        }
        const fs::path topdir{args["tree-dir"].as<string>()};
        // Must load taxonomy before trees
        auto taxonomy = load_taxonomy(args);
        tts.setTaxonomy(taxonomy);
        if (!read_trees(topdir, tts)) {
            return 2;
        }
        //
        if (tts.getNumTrees() == 0) {
            cerr << "No tree to serve. Exiting...\n";
            return 3;
        }
        ////// ROUTES
        auto r_about = make_shared< Resource >( );
        r_about->set_path( "/tree_of_life/about" );
        r_about->set_method_handler( "POST", about_method_handler );
        auto r_node_info = make_shared< Resource >( );
        r_node_info->set_path( "/tree_of_life/node_info" );
        r_node_info->set_method_handler( "POST", node_info_method_handler );
        auto r_mrca = make_shared< Resource >( );
        r_mrca->set_path( "/tree_of_life/mrca" );
        r_mrca->set_method_handler( "POST", mrca_method_handler );
        /////  SETTINGS
        auto settings = make_shared< Settings >( );
        settings->set_port( port_number );
        settings->set_worker_limit( num_threads );
        settings->set_default_header( "Connection", "close" );
        
        Service service;
        service.publish( r_about );
        service.publish( r_node_info );
        service.publish( r_mrca );
        LOG(INFO) << "starting service with " << num_threads << " on port " << port_number << "...\n";
        service.start( settings );
        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        LOG(ERROR) <<"otc-tol-ws: Error! " << e.what() << std::endl;
        exit(1);
    }
}