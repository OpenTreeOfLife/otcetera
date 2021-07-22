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
#include "otc/taxonomy/diff_maker.h"
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
        ("oldtaxonomy", value<string>(),"Filename for the old taxonomy")
        ("newtaxonomy", value<string>(),"Filename for the new taxonomy")
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
    p.add("newtaxonomy", 1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-diff-maker <taxonomy-dir> [OPTIONS]\n"
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

class TaxonomyAmendment {
    public:
        virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &) = 0;
        virtual ~TaxonomyAmendment(){
        }
};

class BaseForwardAmendment : public TaxonomyAmendment {
    public:
    BaseForwardAmendment(const json & taxon_obj) {
        this->former_id = get_unsigned_property(taxon_obj, "former", true).second;
        this->redirect_to_id = get_unsigned_property(taxon_obj, "redirect_to", true).second;
    }
    
    virtual ~BaseForwardAmendment(){
    }

    protected:
        OttId former_id;
        OttId redirect_to_id;
};


class BaseSynonymAmendment : public TaxonomyAmendment {
    public:
    BaseSynonymAmendment(const json & taxon_obj) {
        this->ott_id = get_unsigned_property(taxon_obj, "ott_id", true).second;
        this->name = get_string_property(taxon_obj, "name", true).second;
        this->source_info = get_string_property(taxon_obj, "sourceinfo", false).second;
    }
    
    virtual ~BaseSynonymAmendment(){
    }
    
    protected:
        OttId ott_id;
        std::string name;
        std::string source_info;
};

class BaseTaxonAmendment: public TaxonomyAmendment {
    public:
    BaseTaxonAmendment(const json & taxon_obj, bool all_req)
        :parent_id(UINT_MAX),
        rank(TaxonomicRank::RANK_NO_RANK),
        flags_set(false) {
        this->taxon_id = get_unsigned_property(taxon_obj, "ott_id", true).second;
        this->parent_id = get_unsigned_property(taxon_obj, "parent", all_req).second;
        this->name = get_string_property(taxon_obj, "name", all_req).second;
        auto si = get_string_property(taxon_obj, "sourceinfo", false);
        if (si.first) {
            this->source_info = si.second;
        }
        auto ri = get_string_property(taxon_obj, "rank", false);
        if (ri.first) {
            this->rank = string_to_rank(ri.second, true);
        }
        auto fi = get_string_property(taxon_obj, "flags", false);
        if (fi.first) {
            this->flags = flags_from_string(fi.second);
            flags_set = true;
        }
    }
    
    virtual ~BaseTaxonAmendment(){
    }

    
    protected:
    OttId taxon_id;
    OttId parent_id;
    std::string source_info;
    TaxonomicRank rank;
    std::string name;
    tax_flags flags;
    bool flags_set;
};


class TaxonAdditionAmendment: public BaseTaxonAmendment {
    public:
    TaxonAdditionAmendment(const json & taxon_obj)
        :BaseTaxonAmendment(taxon_obj,true) {
    }
    
    virtual ~TaxonAdditionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        std::string empty;
        std::string fs = flags_to_string(flags);
        auto rank_str = rank_enum_to_name.at(rank);
        return t.add_new_taxon(taxon_id, parent_id, name, rank_str, source_info, empty, fs);
    }
};

class TaxonEditAmendment: public BaseTaxonAmendment {
    public:
    TaxonEditAmendment(const json & taxon_obj)
        :BaseTaxonAmendment(taxon_obj, false) {
    }
    
    virtual ~TaxonEditAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        std::string empty;
        std::string fs;
        auto rank_str = rank_enum_to_name.at(rank);
        if (flags_set) {
            fs = flags_to_string(flags);
        } else {
            fs = empty;
        }
        return t.edit_taxon(taxon_id, parent_id, name, rank_str, source_info, empty, fs, flags_set);
    }
};

class TaxonDeletionAmendment: public BaseTaxonAmendment {
    public:
    TaxonDeletionAmendment(const json & taxon_obj)
        :BaseTaxonAmendment(taxon_obj, false) {
    }
    
    virtual ~TaxonDeletionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        return t.delete_taxon(taxon_id);
    }
};

class ForwardAdditionAmendment : public BaseForwardAmendment {
    public:
    ForwardAdditionAmendment(const json & taxon_obj)
        :BaseForwardAmendment(taxon_obj) {
    }
    
    virtual ~ForwardAdditionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        return t.add_forward(former_id, redirect_to_id);
    }
};

class ForwardDeletionAmendment : public BaseForwardAmendment {
    public:
    ForwardDeletionAmendment(const json & taxon_obj)
        :BaseForwardAmendment(taxon_obj) {
    }
    
    virtual ~ForwardDeletionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        return t.delete_forward(former_id, redirect_to_id);
    }
};


class SynonymAdditionAmendment : public BaseSynonymAmendment {
    public:
    SynonymAdditionAmendment(const json & taxon_obj)
        :BaseSynonymAmendment(taxon_obj) {
    }
    
    virtual ~SynonymAdditionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        return t.add_synonym(name, ott_id, source_info);
    }
};

class SynonymDeletionAmendment : public BaseSynonymAmendment {
    public:
    SynonymDeletionAmendment(const json & taxon_obj)
        :BaseSynonymAmendment(taxon_obj) {
    }
    
    virtual ~SynonymDeletionAmendment(){
    }

    virtual std::pair<bool, std::string> patch(TaxonomyDiffMaker &t) {
        return t.delete_synonym(name, ott_id);
    }
};

typedef std::shared_ptr<TaxonomyAmendment> TaxonomyAmendmentPtr;

TaxonomyAmendmentPtr parse_taxon_amendment_obj(const json & edit_obj) {
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
        auto syn_j = get_object_property(edit_obj, "synonym", false);
        if (syn_j.first) {
           return std::make_shared<SynonymAdditionAmendment>(*(syn_j.second));
        }
        throw OTCError() << "Expecting add action to contain taxon, forward, or synonym object.";
    } else if (action == "delete") {
        auto taxon_j = get_object_property(edit_obj, "taxon", false);
        if (taxon_j.first) {
            return std::make_shared<TaxonDeletionAmendment>(*(taxon_j.second));
        }
        auto forward_j = get_object_property(edit_obj, "forward", false);
        if (forward_j.first) {
           return std::make_shared<ForwardDeletionAmendment>(*(forward_j.second));
        }
        auto syn_j = get_object_property(edit_obj, "synonym", false);
        if (syn_j.first) {
           return std::make_shared<SynonymDeletionAmendment>(*(syn_j.second));
        }
        throw OTCError() << "Expecting delete action to contain taxon, forward, or synonym object.";
    } else if (action == "edit") {
        auto taxon_j = get_object_property(edit_obj, "taxon", false);
        if (taxon_j.first) {
            return std::make_shared<TaxonEditAmendment>(*(taxon_j.second));
        }
        throw OTCError() << "Expecting edit action to contain taxon object.";
    } else {
        throw OTCError() << "Taxon amendment with action \"" << action << "\" not implemented.";
    }
    TaxonomyAmendmentPtr ta;
    return ta;
}


std::list<TaxonomyAmendmentPtr> parse_taxon_amendments_json(std::istream & inp) {
    json edits_obj = json::parse(inp);
    std::list<TaxonomyAmendmentPtr> edit_list;
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




bool edit_taxonomy(TaxonomyDiffMaker & taxonomy,
                   const std::list<TaxonomyAmendmentPtr> & edit_list,
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
    return num_applied == num_attempts;
}


int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc, argv);
        std::ostream & out = std::cout;
        using list_amend_t = std::list<TaxonomyAmendmentPtr>;
        if (!args.count("oldtaxonomy")) {
            cerr << "oldtaxonomy expected as first unnamed argument\n";
            return 1;
        }
        if (!args.count("newtaxonomy")) {
            cerr << "newtaxonomy expected as second unnamed argument\n";
            return 1;
        }
        string otd = args["oldtaxonomy"].as<string>();
        string ntd = args["newtaxonomy"].as<string>();
        OttId keep_root = -1;
        bitset<32> cleaning_flags = 0;
        out << "loading old taxonomy" << std::endl;
        TaxonomyDiffMaker otaxonomy = {otd, cleaning_flags, keep_root};
        TaxonomyDiffMaker ntaxonomy = {ntd, cleaning_flags, keep_root};
        
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

