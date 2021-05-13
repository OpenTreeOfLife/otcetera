#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <bitset>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/patching.h"
#include "otc/config_file.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;
using json = nlohmann::json;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("taxonomy", value<string>(),"Filename for the taxonomy")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("edits", value<string>(), "filepath of JSON file with terse taxonomy edits (this is only relevant when the --write-taxonomy option in effect)")
        ;

    options_description output("Output options");
    output.add_options()
        ("amend-status-to-stdout", "Primarily for debugging. Writes amendment status info to stdout rather than stderr.")
        ("write-to-stdout","Primarily for debugging. Writes contents of taxonomy output to stdout. Only used if write-taxonomy is not used.")
        ("write-taxonomy",value<string>(),"Write out the result as a taxonomy to directory 'arg'")
        ;

    options_description visible;
    visible.add(taxonomy).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-patcher <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy and edit JSON files",
                                                    visible, invisible, p);

    return vm;
}


////////////////////////////////////////////////////////////////////////////////
// possible content for a JSON parsing helpers
////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////
// end possible content for a JSON parsing helpers
////////////////////////////////////////////////////////////////////////////////
// possible content for an amendments header file
////////////////////////////////////////////////////////////////////////////////

class TaxonAmendment {
    public:
        virtual std::pair<bool, std::string> patch(PatchableTaxonomy &) = 0;
        virtual ~TaxonAmendment(){
        }
};

class TaxonAdditionAmendment: public TaxonAmendment {
    public:
    TaxonAdditionAmendment(const json & taxon_obj)
        :rank(TaxonomicRank::RANK_NO_RANK) {
        this->taxon_id = get_unsigned_property(taxon_obj, "ott_id", true).second;
        this->parent_id = get_unsigned_property(taxon_obj, "parent", true).second;
        this->name = get_string_property(taxon_obj, "name", true).second;
        auto si = get_string_property(taxon_obj, "sourceinfo", false);
        if (si.first) {
            this->source_info = si.second;
        }
        auto ri = get_string_property(taxon_obj, "rank", false);
        if (ri.first) {
            this->rank = string_to_rank(ri.second, true);
        }
    }
    
    virtual ~TaxonAdditionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(PatchableTaxonomy &);
    
    private:
        OttId taxon_id;
        OttId parent_id;
        std::string source_info;
        TaxonomicRank rank;
        std::string name;
};


inline std::pair<bool, std::string> TaxonAdditionAmendment::patch(PatchableTaxonomy &t) {
    std::string empty;
    auto rank_str = rank_enum_to_name.at(rank);
    return t.add_new_taxon(taxon_id, parent_id, name, rank_str, source_info, empty, empty);
}

class ForwardAdditionAmendment : public TaxonAmendment {
    public:
    ForwardAdditionAmendment(const json & taxon_obj) {
        this->former_id = get_unsigned_property(taxon_obj, "former", true).second;
        this->redirect_to_id = get_unsigned_property(taxon_obj, "redirect_to", true).second;
    }
    
    virtual ~ForwardAdditionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(PatchableTaxonomy &t) {
        return t.add_forward(former_id, redirect_to_id);
    }
    
    private:
        OttId former_id;
        OttId redirect_to_id;
};

typedef std::shared_ptr<TaxonAmendment> TaxonAmendmentPtr;

TaxonAmendmentPtr parse_taxon_amendment_obj(const json & edit_obj) {
    if (! edit_obj.is_object()) {
        throw OTCError() << "Expecting a taxon amendment object";
    }
    auto action = get_string_property(edit_obj, "action", true).second;
    if (action == "add") {
        auto taxon_j = get_object_property(edit_obj, "taxon", false);
        if (taxon_j.first) {
            return std::make_shared<TaxonAdditionAmendment>(*(taxon_j.second));
        }
        auto forward_j = get_object_property(edit_obj, "forward", false);
        if (forward_j.first) {
           return std::make_shared<ForwardAdditionAmendment>(*(forward_j.second));
        }
        
    } else {
        throw OTCError() << "Taxon amendment with action \"" << action << "\" not implemented.";
    }
    TaxonAmendmentPtr ta;
    return ta;
}


std::list<TaxonAmendmentPtr> parse_taxon_amendments_json(std::istream & inp) {
    json edits_obj = json::parse(inp);
    std::list<TaxonAmendmentPtr> edit_list;
    if (edits_obj.is_object()) {
        edit_list.push_back(parse_taxon_amendment_obj(edits_obj));
    } else if (edits_obj.is_array()) {
        unsigned obj_num = 0;
        for (auto eo : edits_obj) {
            if (eo.is_object()) {
                edit_list.push_back(parse_taxon_amendment_obj(eo));
            } else {
                throw OTCError() << "Expecting amendment object, but element " << obj_num << " of array was not an object.";
            }
            ++obj_num;
        }
    } else {
        throw OTCError() << "Expecting amendment array or object.";
    }
    return edit_list;
}

////////////////////////////////////////////////////////////////////////////////
// end possible content for an amendments header file
////////////////////////////////////////////////////////////////////////////////




void edit_taxonomy(PatchableTaxonomy & taxonomy,
                   const std::list<TaxonAmendmentPtr> & edit_list,
                   bool amend_status_to_stdout) {
    std::size_t num_attempts = 0;
    std::size_t num_applied = 0;
    std::ostream & outp = (amend_status_to_stdout ? std::cout : std::cerr);
    outp << edit_list.size() << " amendments to be processed." << std::endl;
    for (auto tap : edit_list) {
        ++num_attempts;
        auto x = tap->patch(taxonomy);
        if (x.first) {
            ++num_applied;
        } else {
            outp << "amendement #" << num_attempts << " not applied: " << x.second << std::endl;
        }
    }
    outp << num_applied << "/" << num_attempts << " amendments applied." << std::endl;
}


int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        const bool do_json_edits = bool(args.count("edits"));
        const bool amend_status_to_stdout = bool(args.count("amend-status-to-stdout"));
        std::ostream & out = (amend_status_to_stdout ? std::cout : std::cerr);
        std::list<TaxonAmendmentPtr> edits;
        if (do_json_edits) {
            string edit_fp = args["edits"].as<string>();
            std::ifstream edit_stream(edit_fp);
            if (!edit_stream.good()) {
                out << "otc-taxonomy-parser: Could not open \"" << edit_fp << "\"" << std::endl;
                return 1;
            }
            try {
                edits = parse_taxon_amendments_json(edit_stream);
            } catch (...) {
                out <<  "otc-taxonomy-parser: Could not parse \"" << edit_fp << "\"" << std::endl;
                throw;
            }
        }
        out << "loading taxonomy" << std::endl;
        auto taxonomy = load_patchable_taxonomy(args);
        out << "loaded" << std::endl;
        if (do_json_edits) {
            string edit_fp = args["edits"].as<string>();
            edit_taxonomy(taxonomy, edits, amend_status_to_stdout);
        }
        if (args.count("write-taxonomy")) {
            taxonomy.write(args["write-taxonomy"].as<string>());
        } else if (args.count("write-to-stdout")) {
            taxonomy.write_to_stream(std::cout);
        }
    } catch (std::exception& e) {
        cerr << "otc-taxonomy-parser: Error! " << e.what() << std::endl;
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

