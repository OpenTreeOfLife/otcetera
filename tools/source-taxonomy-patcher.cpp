#include "otc/taxonomy/taxonomy-diff.h"
#include <boost/program_options.hpp>
#include "otc/taxonomy/diff_maker.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

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
using json = nlohmann::json;
using vec_strv_t = std::vector<std::string_view>;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

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


// std::list<TaxonomyAmendmentPtr> parse_taxon_amendments_json(std::istream & inp) {
//     json edits_obj = json::parse(inp);
//     std::list<TaxonomyAmendmentPtr> edit_list;
//     if (edits_obj.is_object()) {
//         edit_list.push_back(parse_taxon_amendment_obj(edits_obj));
//     } else if (edits_obj.is_array()) {
//         unsigned obj_num = 0;
//         for (auto eo : edits_obj) {
//             if (eo.is_object()) {
//                 edit_list.push_back(parse_taxon_amendment_obj(eo));
//             } else {
//                 throw OTCError() << "Expecting amendment object, but element " << obj_num << " of array was not an object.";
//             }
//             ++obj_num;
//         }
//     } else {
//         throw OTCError() << "Expecting amendment array or object.";
//     }
//     return edit_list;
// }

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


AlphaEdit parse_alpha(const json & edit_obj) {
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

AlphaGroupEdit parse_alpha_group(const json & edit_obj) {
    AlphaGroupEdit ed;
        auto op_str = get_string_property(edit_obj, "operation", true).second;
    auto s2aIt = str2ageo.find(op_str);
    if (s2aIt == str2ageo.end()) {
        throw OTCError() << "Unrecognized operation " << op_str << " in group edit.";
    }
    AlphaGroupEditOp op_enum = s2aIt->second;
    ed.operation = op_enum;
    // if (op_enum == AlphaEditOp::NO_CHANGE) {
    //     // no op
    // } else if (op_enum == AlphaEditOp::CHANGED_ID) {
    //     ed.first_id = get_unsigned_property(edit_obj, "from", true).second;
    //     ed.second_id = get_unsigned_property(edit_obj, "to", true).second;
    // } else {
    //     ed.first_id = get_unsigned_property(edit_obj, "taxon_id", true).second;
    //     if (op_enum == AlphaEditOp::CHANGED_NAME) {
    //         ed.first_str = get_string_property(edit_obj, "from", true).second;
    //         ed.second_str = get_string_property(edit_obj, "to", true).second;
    //     } else if (op_enum == AlphaEditOp::DELETED_SYN || op_enum == AlphaEditOp::ADDED_SYN) {
    //         ed.first_str = get_string_property(edit_obj, "synonym", true).second;
    //         std::string x = get_string_property(edit_obj, "type", false).second;
    //         if (!x.empty()) {
    //             ed.second_str = x;
    //         }
    //     } else if (op_enum == AlphaEditOp::DELETE_TAXON) {
    //         ed.first_str = get_string_property(edit_obj, "name", true).second;
    //     } else if (op_enum == AlphaEditOp::CHANGED_RANK) {
    //         string x = get_string_property(edit_obj, "from", true).second;
    //         ed.first_rank = string_to_rank(x, true);
    //         x = get_string_property(edit_obj, "to", true).second;
    //         ed.second_rank = string_to_rank(x, true);
    //     } else if (op_enum == AlphaEditOp::CHANGED_FLAGS) {
    //         string x = get_string_property(edit_obj, "from", true).second;
    //         ed.first_flags = flags_from_string(x);
    //         x = get_string_property(edit_obj, "to", true).second;
    //         ed.second_flags = flags_from_string(x);
    //     } else {
    //         assert(op_enum == AlphaEditOp::ADD_TAXON);
    //         ed.first_str = get_string_property(edit_obj, "name", true).second;
    //         string x = get_string_property(edit_obj, "rank", false).second;
    //         ed.first_rank = string_to_rank(x, true);
    //         x = get_string_property(edit_obj, "flags", false).second;
    //         if (!x.empty()) {
    //             ed.first_flags = flags_from_string(x);
    //         }
    //     }
    // }
    return ed;    
}

AlphaGroupEdit parse_higher_taxon(const json & edit_obj) {
    AlphaGroupEdit ed;
    return ed;    
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
            higher.push_back(parse_alpha_group(aged));
        } else {
            throw OTCError() << "Expecting \"alpha_groups\" to be an array object, but element " << obj_num << " of array was not an object.";
        }
        ++obj_num;
    }
    auto higher_array = get_array_property(edits_obj, "higher_taxa", true).second;
    obj_num = 0;
    for (auto hed : *alpha_groups_array) {
        if (hed.is_object()) {
            alpha_groups.push_back(parse_higher_taxon(hed));
        } else {
            throw OTCError() << "Expecting \"higher_taxa\" to be an array object, but element " << obj_num << " of array was not an object.";
        }
        ++obj_num;
    }
}

void patch_source_taxonomy(RichTaxonomy & otaxonomy,
                           list<AlphaEdit> & alpha_taxa,
                           list<AlphaGroupEdit> & alpha_groups,
                           list<AlphaGroupEdit> & higher,
                           std::string outdir) {

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
        RichTaxonomy otaxonomy = {otd, cleaning_flags, keep_root};
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

