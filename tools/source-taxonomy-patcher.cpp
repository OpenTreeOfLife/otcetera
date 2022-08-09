#include "otc/taxonomy/taxonomy-diff.h"
#include <boost/program_options.hpp>
#include "otc/taxonomy/diff_maker.h"
#include <boost/filesystem/operations.hpp>
INITIALIZE_EASYLOGGINGPP

using namespace otc;

namespace fs = boost::filesystem;
using std::string;
using std::list;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;
using std::map;
using std::string_view;
using std::ofstream;
using json = nlohmann::json;
using vec_strv_t = std::vector<std::string_view>;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

unsigned focal_id = 993561;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("oldtaxonomy", value<string>(),"Filename for the old taxonomy")
        ("edits", value<string>(),"Filename for source edits from otc-taxonomy-diff-maker")
        ("outdir", value<string>(),"Filepath for directory to write the output to.")
        ;

    options_description output("Output options");
    output.add_options()
        ("write-to-stdout","Primarily for debugging. Writes contents of taxonomy output to stdout. Only used if write-taxonomy is not used.")
        ;

    options_description visible;
    visible.add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("oldtaxonomy", 1);
    p.add("edits", 1);
    p.add("outdir", 1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-source-taxonomy-patcher <taxonomy-dir> <edit.json> <outdir>\n"
                                                    "Read a taxonomy and edit JSON files for source taxonomies",
                                                    visible, invisible, p);

    return vm;
}

const json * find_object_ptr(const json & j,
                             const std::string & prop_name,
                             bool required) {
    auto aIt = j.find(prop_name);
    if (aIt == j.end()) {
        if (required) {
            throw OTCError() << "Expecting a \"" << prop_name << "\" property in JSON object";
        }
        return nullptr;
    }
    return &(*aIt);
}

std::pair<bool, std::string> get_string_property(const json & j,
                                                 const std::string & prop_name,
                                                 bool required=false) {
    auto jp = find_object_ptr(j, prop_name, required);
    if (jp != nullptr) {
        try {
            std::string str_form = jp->get<std::string>();
            return std::pair<bool, std::string>(true, str_form);
        } catch (...) {
            throw OTCError() << "Expecting \"" << prop_name << "\" property to be a string";
        }
    }
    return std::pair<bool, std::string>(false, "");
}

std::pair<bool, const json *> get_object_property(const json & j,
                                                  const std::string & prop_name,
                                                  bool required) {
    auto jp = find_object_ptr(j, prop_name, required);
    if (jp != nullptr) {
        if (!jp->is_object()) {
            throw OTCError() << "Expecting \"" << prop_name << "\" property to be an object";
        }
        return std::pair<bool, const json *>(true, jp);
    }
    return std::pair<bool, const json *>(false, nullptr);
}

std::pair<bool, const json *> get_array_property(const json & j,
                                                  const std::string & prop_name,
                                                  bool required) {
    auto jp = find_object_ptr(j, prop_name, required);
    if (jp != nullptr) {
        if (!jp->is_array()) {
            throw OTCError() << "Expecting \"" << prop_name << "\" property to be an array";
        }
        return std::pair<bool, const json *>(true, jp);
    }
    return std::pair<bool, const json *>(false, nullptr);
}

const unsigned parse_as_unsigned(const json & j) {
    if (!j.is_number_unsigned()) {
        throw OTCError() << "Expecting a non negative integer";
    }
    return j.get<unsigned>();
}

std::pair<bool, unsigned> get_unsigned_property(const json & j,
                                                const std::string & prop_name,
                                                bool required) {
    auto jp = find_object_ptr(j, prop_name, required);
    if (jp != nullptr) {
        if (!jp->is_number_unsigned()) {
            throw OTCError() << "Expecting \"" << prop_name << "\" property to be a non negative integer";
        }
        unsigned uint_form = jp->get<unsigned>();
        return std::pair<bool, unsigned>(false, uint_form);
    }
    return std::pair<bool, unsigned>(false, UINT_MAX);
}


inline AlphaEdit parse_alpha(const json & edit_obj) {
    AlphaEdit ed;
    auto op_str = get_string_property(edit_obj, "operation", true).second;
    auto s2aIt = str2aeo.find(op_str);
    if (s2aIt == str2aeo.end()) {
        throw OTCError() << "Unrecognized operation " << op_str << " in alpha edit.";
    }
    AlphaEditOp op_enum = s2aIt->second;
    ed.operation = op_enum;
    if (op_enum == AlphaEditOp::NO_CHANGE) {
        // no op
    } else if (op_enum == AlphaEditOp::CHANGED_ID) {
        ed.first_id = get_unsigned_property(edit_obj, "from", true).second;
        ed.second_id = get_unsigned_property(edit_obj, "to", true).second;
    } else {
        ed.first_id = get_unsigned_property(edit_obj, "taxon_id", true).second;
        if (op_enum == AlphaEditOp::CHANGED_NAME) {
            ed.first_str = get_string_property(edit_obj, "from", true).second;
            ed.second_str = get_string_property(edit_obj, "to", true).second;
        } else if (op_enum == AlphaEditOp::DELETED_SYN || op_enum == AlphaEditOp::ADDED_SYN) {
            ed.first_str = get_string_property(edit_obj, "synonym", true).second;
            std::string x = get_string_property(edit_obj, "type", false).second;
            if (!x.empty()) {
                ed.second_str = x;
            }
        } else if (op_enum == AlphaEditOp::DELETE_TAXON) {
            ed.first_str = get_string_property(edit_obj, "name", true).second;
        } else if (op_enum == AlphaEditOp::CHANGED_RANK) {
            string x = get_string_property(edit_obj, "from", true).second;
            ed.first_rank = string_to_rank(x, true);
            x = get_string_property(edit_obj, "to", true).second;
            ed.second_rank = string_to_rank(x, true);
        } else if (op_enum == AlphaEditOp::CHANGED_FLAGS) {
            string x = get_string_property(edit_obj, "from", true).second;
            ed.first_flags = flags_from_string(x);
            x = get_string_property(edit_obj, "to", true).second;
            ed.second_flags = flags_from_string(x);
        } else {
            assert(op_enum == AlphaEditOp::ADD_TAXON);
            ed.first_str = get_string_property(edit_obj, "name", true).second;
            string x = get_string_property(edit_obj, "rank", false).second;
            ed.first_rank = string_to_rank(x, true);
            x = get_string_property(edit_obj, "flags", false).second;
            if (!x.empty()) {
                ed.first_flags = flags_from_string(x);
            }
        }
    }
    return ed;
}

inline AlphaGroupEdit parse_alpha_group(const json & edit_obj) {
    AlphaGroupEdit ed;
        auto op_str = get_string_property(edit_obj, "operation", true).second;
    auto s2aIt = str2ageo.find(op_str);
    if (s2aIt == str2ageo.end()) {
        throw OTCError() << "Unrecognized operation " << op_str << " in group edit.";
    }
    AlphaGroupEditOp op_enum = s2aIt->second;
    ed.operation = op_enum;
    if (op_enum == AlphaGroupEditOp::NO_GR_CHANGE) {
        // no op
    } else if (op_enum == AlphaGroupEditOp::GR_CHANGED_ID) {
        ed.first_id = get_unsigned_property(edit_obj, "from", true).second;
        ed.second_id = get_unsigned_property(edit_obj, "to", true).second;
    } else {
        ed.first_id = get_unsigned_property(edit_obj, "taxon_id", true).second;
        if (op_enum == AlphaGroupEditOp::GR_CHANGED_NAME) {
            ed.first_str = get_string_property(edit_obj, "from", true).second;
            ed.second_str = get_string_property(edit_obj, "to", true).second;
        }
        if (op_enum == AlphaGroupEditOp::ADD_TAXA
           // || op_enum == AlphaGroupEditOp::ADD_DEL_TAXA
            || op_enum == AlphaGroupEditOp::NEW_GROUPING) {
            auto ada = get_array_property(edit_obj, "added", true).second;
            for (auto aed : *ada) {
                unsigned atax_id = parse_as_unsigned(aed);
                if (atax_id == focal_id) {
                    LOG(DEBUG) << "found " << focal_id << " in added for taxon " << ed.first_id ;   
                }
                ed.newChildIds.insert(atax_id);
            }
        }
        // if (op_enum == AlphaGroupEditOp::DEL_TAXA || op_enum == AlphaGroupEditOp::ADD_DEL_TAXA) {
        //     auto dda = get_array_property(edit_obj, "deleted", true).second;
        //     for (auto aed : *dda) {
        //         ed.delIds.insert(parse_as_unsigned(aed));
        //     }
        // }
        if (op_enum == AlphaGroupEditOp::NEW_GROUPING) {
            ed.first_str = get_string_property(edit_obj, "name", true).second;
        }
    }
    return ed;    
}

inline AlphaGroupEdit parse_higher_taxon(const json & edit_obj) {
    return parse_alpha_group(edit_obj); 
}

void parse_source_edits_json(std::istream & inp,
                             list<AlphaEdit> & alpha_taxa,
                             list<AlphaGroupEdit> & alpha_groups,
                             list<AlphaGroupEdit> & higher) {
    json edits_obj = json::parse(inp);
    if (!edits_obj.is_object()) {
        throw OTCError() << "Expecting top level entity in edit JSON to be an object.";
    }
    auto alpha_array = get_array_property(edits_obj, "alpha", true).second;
    unsigned obj_num = 0;
    for (auto aed : *alpha_array) {
        if (aed.is_object()) {
            alpha_taxa.push_back(parse_alpha(aed));
        } else {
            throw OTCError() << "Expecting \"alpha\" to be an array object, but element " << obj_num << " of array was not an object.";
        }
        ++obj_num;
    }
    auto alpha_groups_array = get_array_property(edits_obj, "alpha_groups", true).second;
    obj_num = 0;
    for (auto aged : *alpha_groups_array) {
        if (aged.is_object()) {
            alpha_groups.push_back(parse_alpha_group(aged));
        } else {
            throw OTCError() << "Expecting \"alpha_groups\" to be an array object, but element " << obj_num << " of array was not an object.";
        }
        ++obj_num;
    }
    auto higher_array = get_array_property(edits_obj, "higher_taxa", true).second;
    obj_num = 0;
    for (auto hed : *higher_array) {
        if (hed.is_object()) {
            higher.push_back(parse_higher_taxon(hed));
        } else {
            throw OTCError() << "Expecting \"higher_taxa\" to be an array object, but element " << obj_num << " of array was not an object.";
        }
        ++obj_num;
    }
    LOG(DEBUG) << "alpha_taxa.size() = " << alpha_taxa.size();
    LOG(DEBUG) << "alpha_groups.size() = " << alpha_groups.size();
    LOG(DEBUG) << "higher.size() = " << higher.size();
}

std::list<TaxonomicJuniorSynonym> new_synonyms;
std::set<const RTRichTaxNode *> deleted_nodes;
std::set<const RTRichTaxNode *> detached_nodes;
std::set<const RTRichTaxNode *> attached_nodes;

inline void handle_alpha(const AlphaEdit & aed, RichTaxTree & tree,
                         RTRichTaxTreeData & tree_data) {
    if (aed.operation == AlphaEditOp::NO_CHANGE) {
        return;
    }
    OttId tax_id = aed.first_id;
    RTRichTaxNode * nd = nullptr;
    RTRichTaxNodeData * nd_data = nullptr;
    if (aed.operation != AlphaEditOp::ADD_TAXON) {
            try {
            nd = const_cast<RTRichTaxNode *>(tree_data.id_to_node.at(tax_id));
        } catch (...) {
            throw OTCError() << "Could not find ID " << tax_id << " in taxonomy id_to_node";
        }
        nd_data = &(nd->get_data());
    }
    if (aed.operation == AlphaEditOp::CHANGED_ID) {
        nd->set_ott_id(aed.second_id);
        tree_data.id_to_node[aed.second_id] = nd;
    } else if (aed.operation == AlphaEditOp::CHANGED_NAME) {
        nd->set_name(aed.second_str);
    } else if (aed.operation == AlphaEditOp::DELETED_SYN) {
        while (true) {
            auto jsIt = nd_data->junior_synonyms.begin();
            for (; jsIt != nd_data->junior_synonyms.end(); ++jsIt) {
                if ((*jsIt)->name == aed.first_str) {
                    break;
                }
            }
            if (jsIt != nd_data->junior_synonyms.end()) {
                nd_data->junior_synonyms.erase(jsIt);
            } else {
                break;
            }
        }
    } else if (aed.operation == AlphaEditOp::ADDED_SYN) {
        string mt;
        new_synonyms.emplace_back(aed.first_str, nd, mt);
        TaxonomicJuniorSynonym & js = *new_synonyms.rbegin();
        nd_data->junior_synonyms.push_back(&js);
    } else if (aed.operation == AlphaEditOp::DELETE_TAXON) {
        deleted_nodes.insert(nd);
        collapse_split_dont_del_node(nd);
    } else if (aed.operation == AlphaEditOp::CHANGED_RANK) {
        nd_data->rank = aed.second_rank;
    } else if (aed.operation == AlphaEditOp::CHANGED_FLAGS) {
        nd_data->flags = aed.second_flags;
    } else if (aed.operation == AlphaEditOp::ADD_TAXON) {
        RTRichTaxNode * newnd = tree.create_node(nullptr);
        newnd->set_ott_id(aed.first_id);
        newnd->set_name(aed.first_str);
        auto & newnd_data = newnd->get_data();
        newnd_data.rank = aed.first_rank;
        newnd_data.flags = aed.first_flags;
        tree_data.id_to_node[aed.first_id] = newnd;
    }
}

void handle_alpha_group(const AlphaGroupEdit & aed, RichTaxTree & tree, RTRichTaxTreeData & tree_data) {
    if (aed.operation == AlphaGroupEditOp::NO_GR_CHANGE) {
        return;
    }
    OttId tax_id = aed.first_id;
    RTRichTaxNode * nd = nullptr;
    RTRichTaxNodeData * nd_data = nullptr;
    auto & id2nd = tree_data.id_to_node;
    if (aed.operation == AlphaGroupEditOp::NEW_GROUPING) {
        nd = tree.create_node(nullptr);
        nd->set_ott_id(aed.first_id);
        nd->set_name(aed.first_str);
        nd_data = &(nd->get_data());
        //newnd_data.rank = aed.first_rank;
        //newnd_data.flags = aed.first_flags;
        id2nd[aed.first_id] = nd;
    } else {
            try {
            nd = const_cast<RTRichTaxNode *>(id2nd.at(tax_id));
        } catch (...) {
            throw OTCError() << "Could not find ID " << tax_id << " in taxonomy id_to_node";
        }
        nd_data = &(nd->get_data());
    }
    if (aed.operation == AlphaGroupEditOp::GR_CHANGED_ID) {
        nd->set_ott_id(aed.second_id);
        tree_data.id_to_node[aed.second_id] = nd;
    } else if (aed.operation == AlphaGroupEditOp::GR_CHANGED_NAME) {
        nd->set_name(aed.second_str);
    } else if (aed.operation == AlphaGroupEditOp::ADD_TAXA
        // || aed.operation == AlphaGroupEditOp::ADD_DEL_TAXA
        || aed.operation == AlphaGroupEditOp::NEW_GROUPING
        ) {
        bool added_focal = false;
        for (auto add_id : aed.newChildIds) {
            RTRichTaxNode * nc = const_cast<RTRichTaxNode *>(id2nd.at(add_id));
            if (add_id == focal_id) {
                added_focal = true;
                LOG(DEBUG) << "Adding " << focal_id << " to " << tax_id;
            }
            if (nc->get_parent() != nullptr) {
                nc->detach_this_node();
                detached_nodes.insert(nc);
            }
            nd->add_child(nc);
            attached_nodes.insert(nc);
        }
    }
    //if (aed.operation == AlphaGroupEditOp::ADD_DEL_TAXA
    //    || aed.operation == AlphaGroupEditOp::DEL_TAXA) {
        // no-op // for (auto add_id : aed.delIds) {
        //     RTRichTaxNode * nc = const_cast<RTRichTaxNode *>(id2nd.at(add_id));
        //     if (nc->get_parent() != nullptr) {
        //         nc->detach_this_node();
        //         detached_nodes.insert(nc);
        //     }
        // }
    //} else 
    if (aed.operation == AlphaGroupEditOp::DELETED_GROUPING) {
        if (!contains(deleted_nodes, nd)) {
            deleted_nodes.insert(nd);
            collapse_split_dont_del_node(nd);
        }
    }
}

void handle_higher(const AlphaGroupEdit & aed, RichTaxTree & tree, RTRichTaxTreeData & tree_data) {
    handle_alpha_group(aed, tree, tree_data);
}

void write_patched(const RichTaxTree & tree, const std::string outdir) {
    fs::path new_dir = outdir;
    if (! fs::exists(new_dir)) {
        fs::create_directories(new_dir);
    }
    const char * sep = "\t|\t";
    ofstream tf;
    tf.open((new_dir/"taxonomy.tsv").string());
    if (!tf.good()) {
        throw OTCError() << "could not open " << (new_dir/"taxonomy.tsv").string() ;
    }
    tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\t|\t\n";
    
    const auto sfname = "synonyms.tsv";
    ofstream soutpf;
    soutpf.open((new_dir/sfname).string());
    if (!soutpf.good()) {
        throw OTCError() << "could not open " << (new_dir/sfname).string() ;
    }
    soutpf << "uid\t|\tname\t|\ttype\t|\t\n ";

    for(auto nd: iter_pre(tree)) {
        tf << nd->get_ott_id() << sep;
        auto par = nd->get_parent();
        if (par != nullptr) {
            tf << par->get_ott_id();
        }
        tf << sep << nd->get_name() << sep;
        const auto & nd_data = nd->get_data();
        tf << rank_enum_to_name.at(nd_data.rank) << sep;
        tf << flags_to_string(nd_data.flags) << sep << '\n';
        auto jsIt = nd_data.junior_synonyms.begin();
        for (; jsIt != nd_data.junior_synonyms.end(); ++jsIt) {
            soutpf << nd->get_ott_id() << sep 
                   << (*jsIt)->name << sep 
                   << (*jsIt)->source_string << sep << '\n';
        }
    }

}

void patch_source_taxonomy(RichTaxonomy & otaxonomy,
                           const list<AlphaEdit> & alpha_taxa,
                           const list<AlphaGroupEdit> & alpha_groups,
                           const list<AlphaGroupEdit> & higher,
                           const std::string outdir) {
    RichTaxTree & tree = const_cast<RichTaxTree &>(otaxonomy.get_tax_tree());
    RTRichTaxTreeData & tree_data = tree.get_data();
    for (auto aed : alpha_taxa) {
        handle_alpha(aed, tree, tree_data);
    }
    for (auto aged : alpha_groups) {
        handle_alpha_group(aged, tree, tree_data);
    }
    for (auto aged : higher) {
        handle_higher(aged, tree, tree_data);
    }
    LOG(DEBUG) << "deleted_nodes.size() = " << deleted_nodes.size();
    LOG(DEBUG) << "detached_nodes.size() = " << detached_nodes.size();
    LOG(DEBUG) << "attached_nodes.size() = " << attached_nodes.size();
    write_patched(tree, outdir);
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        if (!args.count("oldtaxonomy")) {
            cerr << "oldtaxonomy expected as first unnamed argument\n";
            return 1;
        }
        if (!args.count("edits")) {
            cerr << "edits expected as second unnamed argument\n";
            return 1;
        }
        if (!args.count("outdir")) {
            cerr << "outdir expected as third unnamed argument\n";
            return 1;
        }
        string otd = args["oldtaxonomy"].as<string>();
        string edit_fp = args["edits"].as<string>();
        string outdir = args["outdir"].as<string>();

        std::ifstream edit_stream(edit_fp);
        if (!edit_stream.good()) {
            cerr << "otc-source-taxonomy-patcher: Could not open \"" << edit_fp << "\"" << std::endl;
            return 1;
        }
        list<AlphaEdit> alpha_taxa;
        list<AlphaGroupEdit> alpha_groups;
        list<AlphaGroupEdit> higher;

        try {
            parse_source_edits_json(edit_stream, alpha_taxa, alpha_groups, higher);
        } catch (...) {
            cerr <<  "otc-source-taxonomy-patcher: Could not parse \"" << edit_fp << "\"" << std::endl;
            throw;
        }
        OttId keep_root = -1;
        bitset<32> cleaning_flags = 0;
        LOG(INFO) << "loading old taxonomy\n";
        Taxonomy::tolerate_synonyms_to_unknown_id = true;
        RichTaxonomy otaxonomy = {otd, cleaning_flags, keep_root, true};
        patch_source_taxonomy(otaxonomy, alpha_taxa, alpha_groups, higher, outdir);
    } catch (std::exception& e) {
        cerr << "otc-source-taxonomy-patcher: Error! " << e.what() << std::endl;
        return 1;
    }
}

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

