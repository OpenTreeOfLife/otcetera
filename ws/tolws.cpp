#include <regex>
#include "ws/tolws.h"
#include "ws/tolwsadaptors.h"
#if 0
#include <memory>
#include <thread>
#include <sstream>
#include <list>
#include <map>
#include <csignal>
#include <sstream>
#include <thread>
#include <ctime>
// PID reporting from https://github.com/Corvusoft/restbed/blob/master/example/signal_handling/source/example.cpp
#ifdef _WIN32
    #include <process.h>
#else
    #include <unistd.h>
#endif
#include "otc/util.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/config_file.h"
#endif

INITIALIZE_EASYLOGGINGPP


using namespace std;
using namespace boost::property_tree;


namespace otc {

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
    const RichTaxonomy & taxonomy = tts.getTaxonomy();
    auto tree_and_ann = tts.getNewTreeAndAnnotations(config_path.native(), tree_path.native());
    try {
        SummaryTree_t & tree = tree_and_ann.first;
        SummaryTreeAnnotation & sta = tree_and_ann.second;
        sta = annotations_obj;
        json tref;
        tref["taxonomy"] = taxonomy.get_version();
        sta.full_source_id_map_json[taxonomy.get_version()] = tref;

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



void add_taxon_info(const RichTaxonomy & , const RTRichTaxNode & nd_taxon, json & taxonrepr) {
    const auto & taxon_data = nd_taxon.getData();
    taxonrepr["tax_sources"] = taxon_data.getSourcesJSON();
    taxonrepr["name"] = string(taxon_data.getName());
    taxonrepr["uniqname"] = string(taxon_data.getUniqname());
    taxonrepr["rank"] = string(taxon_data.getRank());
    taxonrepr["ott_id"] = nd_taxon.getOttId();    
}

void add_basic_node_info(const RichTaxonomy & taxonomy, const SumTreeNode_t & nd, json & noderepr) {
    noderepr["node_id"] = nd.getName();
    noderepr["num_tips"] = nd.getData().num_tips;
    if (nd.hasOttId()) {
        auto nd_id = nd.getOttId();
        const auto * nd_taxon = taxonomy.taxonFromId(nd_id);
        if (nd_taxon == nullptr) {
            throw OTCError() << "OTT Id " << nd_id << " not found in taxonomy! Please report this bug";
        }
        json taxon;
        add_taxon_info(taxonomy, *nd_taxon, taxon);
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
    const auto & taxonomy = tts.getTaxonomy();
    const auto & d = nd.getData();
    const string * extra_src = nullptr;
    const string * extra_node_id = nullptr;
    if (nd.hasOttId()) {
        extra_src = &(taxonomy.get_version());
        extra_node_id = &(nd.getName());
        usedSrcIds.insert(taxonomy.get_version());
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

const SumTreeNode_t * find_mrca(const SumTreeNode_t *f,
                                const SumTreeNode_t *s) {
    const auto * fdata = &(f->getData());
    const auto sec_ind = s->getData().trav_enter;
    //cerr << "f->name = " << f->getName() <<" sec_ind = " << sec_ind << " enter, exit = " << fdata->trav_enter << " " << fdata->trav_exit << '\n';
    while (sec_ind < fdata->trav_enter || sec_ind > fdata->trav_exit) {
        f = f->getParent();
        if (f == nullptr) {
            assert(false); 
            return nullptr;
        }
        fdata = &(f->getData());
        //cerr << "f->name = " << f->getName() <<" sec_ind = " << sec_ind << " enter, exit = " << fdata->trav_enter << " " << fdata->trav_exit << '\n';
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

void tax_about_ws_method(const TreesToServe &tts,
                         string & response_str,
                         int & status_code) {
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    string weburl;
    weburl = "https://tree.opentreeoflife.org/about/taxonomy-version/ott";
    weburl += taxonomy.get_version_number();
    response["weburl"] = weburl;
    response["author"] = "open tree of life project";
    response["name"] = "ott";
    response["source"] = taxonomy.get_version();
    response["version"] = taxonomy.get_version_number();
    response_str = response.dump(1);
}


inline void add_lineage(json & j, const SumTreeNode_t * focal, const RichTaxonomy & taxonomy, set<string> & usedSrcIds) {
    json lineage_arr;
    const SumTreeNode_t * anc = focal->getParent();
    if (!anc) {
        vector<string> c;
        j["lineage"] = c;
        return;
    }
    while (anc) {
        json ancj;
        add_basic_node_info(taxonomy, *anc, ancj);
        add_node_support_info(tts, *anc, ancj, usedSrcIds);
        lineage_arr.push_back(ancj);
        anc = anc->getParent();
    }
    j["lineage"] = lineage_arr;
}

inline void add_source_id_map(json & j,
                              const set<string> & usedSrcIds,
                              const RichTaxonomy & taxonomy,
                              const SummaryTreeAnnotation * sta ) {
    json sim;
    for (auto srcTag : usedSrcIds) {
        json jt;
        if (srcTag == taxonomy.get_version()) {
            jt["taxonomy"] = taxonomy.get_version();
        } else {
            const auto & simentry = sta->source_id_map.at(srcTag);
            jt = simentry;
        }
        sim[srcTag] = jt;   
    }
    j["source_id_map"] = sim;
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
    if (focal == nullptr) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return;
    }
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    add_basic_node_info(taxonomy, *focal, response);
    add_node_support_info(tts, *focal, response, usedSrcIds);
    if (include_lineage) {
        add_lineage(response, focal, taxonomy, usedSrcIds);
    }
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    if (was_broken) {
        response["response_for_mrca_of_broken_taxon"] = true; //@TODO: discuss and document
    }
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
    if (!focal->hasOttId()) {
        auto anc = focal->getParent();
        assert(anc != nullptr);
        while (!anc->hasOttId()) {
            anc = anc->getParent();
            if (anc == nullptr) {
                response_str = "No ancestors were taxa. That is odd.\n";
                status_code = 500;
                return;
            }
        }
        json nt;
        const RTRichTaxNode * anc_taxon = taxonomy.taxonFromId(anc->getOttId());
        if (anc_taxon == nullptr) {
            throw OTCError() << "Ancd OTT Id " << anc->getOttId() << " not found in taxonomy! Please report this bug";
        }
        add_taxon_info(taxonomy, *anc_taxon, nt);
        response["nearest_taxon"] = nt;
    }
    response["mrca"] = mrcaj;
    add_source_id_map(response, usedSrcIds, taxonomy, sta);
    response_str = response.dump(1);
}


const SumTreeNode_t * get_node_for_subtree(const SummaryTree_t * tree_ptr,
                                           const string & node_id, 
                                           int height_limit,
                                           int tip_limit,
                                           string & response_str,
                                           int & status_code) {
    assert(tree_ptr != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = find_node_by_id_str(*tree_ptr, node_id, was_broken);
    if (focal == nullptr) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return nullptr;
    }
    if (was_broken) {
        response_str = "node_id was not found.\n";
        status_code = 400;
        return nullptr;
    }
    if (focal->getData().num_tips > tip_limit
        && height_limit < 0) {
        response_str = "The requested subtree is too large to be returned via the API. Download the entire tree.\n";
        status_code = 400;
        return nullptr;  
    }
    return focal;
}

class NodeNamerSupportedByStasher {
    public:
        mutable set<const string *> study_id_set;
        NodeNameStyle nns;
        const RichTaxonomy & taxonomy;
        NodeNamerSupportedByStasher(NodeNameStyle in_nns, const RichTaxonomy &tax)
            :nns(in_nns),
            taxonomy(tax) {
        }

        string operator()(const SumTreeNode_t *nd) const {
            const string & id_str = nd->getName(); // in this tree the "name" is really the "ott###" string
            const SumTreeNodeData & d = nd->getData();
            for (auto & p : d.supported_by) {
                study_id_set.insert(p.first);
            }
            if (nns != NodeNameStyle::NNS_ID_ONLY && nd->hasOttId()) {
                const auto * tr = taxonomy.taxonFromId(nd->getOttId());
                if (tr == nullptr) {
                    throw OTCError() << "OTT Id " << nd->getOttId() << " in namer not found in taxonomy! Please report this bug";
                }
                string taxon_name = string(tr->getData().getUniqname());
                if (nns == NodeNameStyle::NNS_NAME_AND_ID) {
                    string ret;
                    ret.reserve(taxon_name.length() + 1 + id_str.length());
                    ret = taxon_name;
                    ret += ' ';
                    ret += id_str;
                    return ret;    
                } else {
                    return taxon_name;
                }
            }
            return id_str;
        }
};

template<typename C, typename T, typename Y>
inline void writeVisitedNewick(std::ostream & out,
                               const C & visited,
                               T nd,
                               Y & nodeNamer) {
    assert(nd != nullptr);
    if (!(nd->isTip())) {
        bool first = true;
        for (auto c : iter_child_const(*nd)) {
            if (visited.count(c)) {
                if (first) {
                    out << '(';
                    first = false;
                } else {
                    out << ',';
                }
                writeVisitedNewick<C, T, Y>(out, visited, c, nodeNamer);
            }
        }
        if (!first) {
            out << ')';
        }
    }
    writeEscapedForNewick(out, nodeNamer(nd));
}


void induced_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const vector<string> & node_id_vec,
                     NodeNameStyle label_format, 
                     string & response_str,
                     int & status_code) {
    assert(tree_ptr != nullptr);
    assert(sta != nullptr);
    bool was_broken = false;
    const SumTreeNode_t * focal = nullptr;
    bool first = true;
    set<const SumTreeNode_t *> tip_nodes;
    for (auto node_id : node_id_vec) {
        const SumTreeNode_t * n = find_node_by_id_str(*tree_ptr, node_id, was_broken);
        tip_nodes.insert(n);
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
    set<const SumTreeNode_t *> visited;
    visited.insert(focal);
    for (auto tni : tip_nodes) {
        auto cnd = tni;
        while (visited.count(cnd) == 0) {
            visited.insert(cnd);
            cnd = cnd->getParent(); 
        } 
    }
    const auto & taxonomy = tts.getTaxonomy();
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    writeVisitedNewick(out, visited, focal, nnsbs);
    json response;
    response["newick"] = out.str();
    json ss_arr;
    for (auto study_it_ptr : nnsbs.study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    response["supporting_studies"] = ss_arr;
    response_str = response.dump(1);
}

void newick_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     NodeNameStyle label_format, 
                     int height_limit,
                     string & response_str,
                     int & status_code) {
    const int NEWICK_TIP_LIMIT = 25000;
    const SumTreeNode_t * focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT, response_str, status_code);
    if (focal == nullptr) {
        return;
    }
    assert(sta != nullptr);
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    NodeNamerSupportedByStasher nnsbs(label_format, taxonomy);
    ostringstream out;
    writeNewickGeneric<const SumTreeNode_t *, NodeNamerSupportedByStasher>(out, focal, nnsbs, height_limit);
    response["newick"] = out.str();
    json ss_arr;
    for (auto study_it_ptr : nnsbs.study_id_set) {
        ss_arr.push_back(*study_it_ptr);
    }
    response["supporting_studies"] = ss_arr;
    response_str = response.dump(1);
}


template<typename T>
inline void write_arguson(json & j,
                          const TreesToServe & tts,
                          const SummaryTreeAnnotation * sta,
                          const RichTaxonomy & taxonomy,
                          T nd,
                          long height_limit,
                          set<string> & usedSrcIds) {
    assert(nd != nullptr);
    if (!(nd->isTip()) && height_limit != 0) {
        json c_array;
        const long nhl = height_limit - 1;
        for (auto c : iter_child_const(*nd)) {
            json cj;
            write_arguson<T>(cj, tts, sta, taxonomy, c, nhl, usedSrcIds);
            c_array.push_back(cj);
        }
        j["children"] = c_array;
    }
    add_basic_node_info(taxonomy, *nd, j);
    add_node_support_info(tts, *nd, j, usedSrcIds);
}

void arguson_subtree_ws_method(const TreesToServe & tts,
                     const SummaryTree_t * tree_ptr,
                     const SummaryTreeAnnotation * sta,
                     const string & node_id,
                     int height_limit,
                     string & response_str,
                     int & status_code) {
    const int NEWICK_TIP_LIMIT = 25000;
    auto focal = get_node_for_subtree(tree_ptr, node_id, height_limit, NEWICK_TIP_LIMIT, response_str, status_code);
    if (focal == nullptr) {
        return;
    }
    assert(sta != nullptr);
    const auto & taxonomy = tts.getTaxonomy();
    status_code = OK;
    json response;
    response["synth_id"] = sta->synth_id;
    set<string> usedSrcIds;
    json a;
    write_arguson(a, tts, sta, taxonomy, focal, height_limit, usedSrcIds);
    add_lineage(a, focal, taxonomy, usedSrcIds);
    response["arguson"] = a;
    response_str = response.dump(1);
}


void taxon_info_ws_method(const TreesToServe & tts,
                          const RTRichTaxNode * taxon_node,
                          bool include_lineage,
                          bool include_children,
                          bool include_terminal_descendants,
                          string & response_str,
                          int & status_code) {
    const auto & taxonomy = tts.getTaxonomy();
    assert(taxon_node != nullptr);
    const auto & taxonomy_tree = taxonomy.getTaxTree();
    const auto & taxonomy_tree_data = taxonomy_tree.getData();
    const auto & node_data = taxon_node->getData();
    json response;
    response["source"] = taxonomy.get_version();
    response["ott_id"] = taxon_node->getOttId();
    response["name"] = string(node_data.getName());
    response["uniqname"] = string(node_data.getUniqname());
    response["tax_sources"] = node_data.getSourcesJSON();
    response["flags"] = flags_to_string_vec(node_data.getFlags());
    json syn_list;
    for (auto tjs : node_data.junior_synonyms) {
        syn_list.push_back(tjs->getName());
    }
    response["synonyms"] = syn_list;
    //add_lineage(a, focal, taxonomy, usedSrcIds);
    //response["arguson"] = a;
    response_str = response.dump(1);
}

} //namespace otc

int main( const int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
        po::variables_map args = parse_cmd_line(argc,argv);
        return run_server(args);
    } catch (std::exception& e) {
        LOG(ERROR) <<"otc-tol-ws: Error! " << e.what() << std::endl;
        return 1;
    }
}
