// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <bitset>
#include <fstream>
#include <regex>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
namespace fs = boost::filesystem;

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
#include "otc/config_file.h"
#include "otc/util.h"
#include "otc/otc_base_includes.h"
#include "otc/ctrie/context_ctrie_db.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::ofstream;
using std::map;
using std::set;
using std::optional;
using nlohmann::json;
using std::ifstream;
using std::unordered_set;

using std::string_view;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

namespace otc
{

bool rank_is_specific(TaxonomicRank rank)
{
    // taxomachine includes "species", "subspecies", "variety", "varietas", "forma", "form"
    if (rank == TaxonomicRank::RANK_SPECIES) return true;
    if (rank == TaxonomicRank::RANK_SUBSPECIES) return true;
    if (rank == TaxonomicRank::RANK_VARIETY) return true;
    if (rank == TaxonomicRank::RANK_VARIETAS) return true;
    if (rank == TaxonomicRank::RANK_FORMA) return true;
    // There is no RANK_FORM?

    return false;
}

const map<string, TaxonomicRank, std::less<>> rank_name_to_enum = 
    {   {"domain", RANK_DOMAIN},
        {"superkingdom", RANK_SUPERKINGDOM},
        {"kingdom", RANK_KINGDOM},
        {"subkingdom", RANK_SUBKINGDOM},
        {"infrakingdom", RANK_INFRAKINGDOM},
        {"superphylum", RANK_SUPERPHYLUM},
        {"phylum", RANK_PHYLUM},
        {"division", RANK_DIVISION},
        {"subphylum", RANK_SUBPHYLUM},
        {"subdivision", RANK_SUBDIVISION},
        {"infraphylum", RANK_INFRAPHYLUM},
        {"superclass", RANK_SUPERCLASS},
        {"class", RANK_CLASS},
        {"subclass", RANK_SUBCLASS},
        {"infraclass", RANK_INFRACLASS},
        {"subterclass", RANK_SUBTERCLASS},
        {"cohort", RANK_COHORT},
        {"subcohort", RANK_SUBCOHORT},
        {"superorder", RANK_SUPERORDER},
        {"order", RANK_ORDER},
        {"suborder", RANK_SUBORDER},
        {"infraorder", RANK_INFRAORDER},
        {"parvorder", RANK_PARVORDER},
        {"superfamily", RANK_SUPERFAMILY},
        {"family", RANK_FAMILY},
        {"subfamily", RANK_SUBFAMILY},
        {"supertribe", RANK_SUPERTRIBE},
        {"tribe", RANK_TRIBE},
        {"subtribe", RANK_SUBTRIBE},
        {"genus", RANK_GENUS},
        {"subgenus", RANK_SUBGENUS},
        {"section", RANK_SECTION},
        {"subsection", RANK_SUBSECTION},
        {"species group", RANK_SPECIES_GROUP},
        {"species subgroup", RANK_SPECIES_SUBGROUP},
        {"species", RANK_SPECIES},
        {"subspecies", RANK_SUBSPECIES},
        {"infraspecificname", RANK_INFRASPECIFICNAME},
        {"forma", RANK_FORMA},
        {"subform", RANK_SUBFORM},
        {"varietas", RANK_VARIETAS},
        {"variety", RANK_VARIETY},
        {"subvariety", RANK_SUBVARIETY},
        {"no rank", RANK_NO_RANK},
        {"no rank - terminal", RANK_NO_RANK_TERMINAL},
        {"natio", RANK_INFRASPECIFICNAME} // not really a rank, should go in subsequent version of OTT
    };

const map<TaxonomicRank, string, std::less<>> rank_enum_to_name = 
    {   {RANK_DOMAIN, "domain"},
        {RANK_SUPERKINGDOM, "superkingdom"},
        {RANK_KINGDOM, "kingdom"},
        {RANK_SUBKINGDOM, "subkingdom"},
        {RANK_INFRAKINGDOM, "infrakingdom"},
        {RANK_SUPERPHYLUM, "superphylum"},
        {RANK_PHYLUM, "phylum"},
        {RANK_DIVISION, "division"},
        {RANK_SUBPHYLUM, "subphylum"},
        {RANK_SUBDIVISION, "subdivision"},
        {RANK_INFRAPHYLUM, "infraphylum"},
        {RANK_SUPERCLASS, "superclass"},
        {RANK_CLASS, "class"},
        {RANK_SUBCLASS, "subclass"},
        {RANK_INFRACLASS, "infraclass"},
        {RANK_SUBTERCLASS, "subterclass"},
        {RANK_COHORT, "cohort"},
        {RANK_SUBCOHORT, "subcohort"},
        {RANK_SUPERORDER, "superorder"},
        {RANK_ORDER, "order"},
        {RANK_SUBORDER, "suborder"},
        {RANK_INFRAORDER, "infraorder"},
        {RANK_PARVORDER, "parvorder"},
        {RANK_SUPERFAMILY, "superfamily"},
        {RANK_FAMILY, "family"},
        {RANK_SUBFAMILY, "subfamily"},
        {RANK_SUPERTRIBE, "supertribe"},
        {RANK_TRIBE, "tribe"},
        {RANK_SUBTRIBE, "subtribe"},
        {RANK_GENUS, "genus"},
        {RANK_SUBGENUS, "subgenus"},
        {RANK_SECTION, "section"},
        {RANK_SUBSECTION, "subsection"},
        {RANK_SPECIES_GROUP, "species group"},
        {RANK_SPECIES_SUBGROUP, "species subgroup"},
        {RANK_SPECIES, "species"},
        {RANK_SUBSPECIES, "subspecies"},
        {RANK_INFRASPECIFICNAME, "infraspecificname"},
        {RANK_FORMA, "forma"},
        {RANK_SUBFORM, "subform"},
        {RANK_VARIETAS, "varietas"},
        {RANK_VARIETY, "variety"},
        {RANK_SUBVARIETY, "subvariety"},
        {RANK_NO_RANK, "no rank"},
        {RANK_NO_RANK_TERMINAL, "no rank - terminal"},
        {RANK_INFRASPECIFICNAME, "natio"} // not really a rank, should go in subsequent version of OTT
    };

const std::string empty_string;
const set<string> indexed_source_prefixes = {"ncbi", "gbif", "worms", "if", "irmng"};
std::set<std::string> rank_strings;

bool TaxonomyRecord::is_extinct() const
{
    return ::is_extinct(flags);
}

TaxonomyRecord::TaxonomyRecord(const string& line_)
    :line(line_) {
    // parse the line
    // also see boost::make_split_iterator
    const char* start[8];
    const char* end[8];
    start[0] = line.c_str();
    for(int i=0; i<7; i++) {
        end[i] = std::strstr(start[i],"\t|\t");
        start[i+1] = end[i] + 3;
    }
    char *temp;
    id = std::strtoul(start[0], &temp, 10);
    parent_id = std::strtoul(start[1], &temp, 10);
    name = string_view(start[2], end[2] - start[2]);
    rank = string_view(start[3], end[3] - start[3]);
    sourceinfo = string_view(start[4], end[4] - start[4]);
    uniqname = string_view(start[5], end[5] - start[5]);
    flags = flags_from_string(start[6], end[6]);
    if (not uniqname.size()) {
        uniqname = name;
    }
    rank_strings.insert(string(rank));
}

// is_input_form will be true when the headers lack sourceinfo and uniqname
TaxonomyRecord::TaxonomyRecord(const string& line_, bool /* is_input_form */)
    :line(line_) {
    // parse the line
    // also see boost::make_split_iterator
    const char* start[6];
    const char* end[6];
    start[0] = line.c_str();
    for(int i = 0; i < 5; i++) {
        end[i] = std::strstr(start[i],"\t|\t");
        start[i + 1] = end[i] + 3;
        // std::cerr << "start,end[" << i << "] = " << start[i] << ", " << end[i] << std::endl;
    }
    char *temp;
    id = std::strtoul(start[0], &temp, 10);
    parent_id = std::strtoul(start[1], &temp, 10);
    name = string_view(start[2], end[2] - start[2]);
    rank = string_view(start[3], end[3] - start[3]);
    sourceinfo = string_view();
    uniqname = string_view();
    flags = flags_from_string(start[4], end[4]);
    if (not uniqname.size()) {
        uniqname = name;
    }
    rank_strings.insert(string(rank));
}

optional<int> Taxonomy::maybe_index_from_id(OttId id) const
{
    auto loc = index.find(id);
    if (loc == index.end()) {
        auto loc2 = forwards.find(id);

        // not in taxonomy or forwarding list";
        if (loc2 == forwards.end()) return {};

        OttId newid = loc2->second;
        loc = index.find(newid);
        // If id is in the forwarding table, then newid should be in the taxonomy
        assert(loc != index.end());
    }
    return loc->second;
}

int Taxonomy::index_from_id(OttId id) const
{
    if (auto index = maybe_index_from_id(id))
        return *index;
    else
        throw OTCError() << "ID " << id << " not in taxonomy or forwarding list";
}

const TaxonomyRecord& Taxonomy::record_from_id(OttId id) const {
    return (*this)[index_from_id(id)];
}

TaxonomyRecord& Taxonomy::record_from_id(OttId id) {
    return (*this)[index_from_id(id)];
}

OttId Taxonomy::map(OttId old_id) const {
    if (index.count(old_id)) {
        return old_id;
    }
    auto loc = forwards.find(old_id);
    if (loc != forwards.end()) {
        return loc->second;
    }
    //if (deprecated.count(old_id) != 0) {
    //    return -2;
    //}
    return -1;
}


void Taxonomy::write(const std::string& newdirname, bool copy_taxonomy_tsv_lines_raw) {
    fs::path old_dir = path;
    fs::path new_dir = newdirname;
    if (! fs::exists(new_dir)) {
        fs::create_directories(new_dir);
    }
    
    // Copy the other files.
    for(const auto& name: {"about.json", "conflicts.tsv", "deprecated.tsv",
                "log.tsv", "otu_differences.tsv", "synonyms.tsv", "weaklog.csv"}) {
        if (fs::exists(old_dir/name)) {
            fs::copy_file(old_dir/name,new_dir/name);
        }
    }
    // Write the new version file.
    {
        ofstream version_file((new_dir/"version.txt").string());
        version_file << version;
        if (keep_root != -1 or cleaning_flags.any()) {
            version_file << "modified: ";
            if (keep_root != -1) {
                version_file << "root=" << keep_root <<"  ";
            }
            if (cleaning_flags.any()) {
                version_file << "cleaning_flags=" << flags_to_string(cleaning_flags);
            }
            version_file << "\n";
        }
        version_file.close();
    }
    // Write the new taxonomy file.
    {
        ofstream tf ((new_dir/"taxonomy.tsv").string());
        tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t" << std::endl;
        
        if (copy_taxonomy_tsv_lines_raw) {
            for(const auto& r: *this) {
                tf << r.line <<"\n";
            }
        } else {
           string sep = "\t|\t";
           for (auto& rec: *this) {
               tf << rec.id << sep;
               if (rec.parent_id > 0) {
                   tf << rec.parent_id;
               }
               tf << sep;
               tf << rec.name << sep;
               tf << rec.rank << sep;
               tf << rec.sourceinfo << sep;
               if (rec.uniqname != rec.name) {
                   tf << rec.uniqname;
               }
               tf << sep;
               tf << flags_to_string(rec.flags) << sep;
               tf << '\n';
           }

           
        }
        tf.close();
    }
    // Write the new forwards file.
    {
        ofstream ff((new_dir/"forwards.tsv").string());
        ff << "id\treplacement\n";
        for(const auto& p: forwards) {
            ff << p.first << '\t' << p.second << '\n';
        }
        ff.close();
    }
}

vector<string> Taxonomy::path_from_id(OttId I) const
{
    vector<string> path;
    int i = index_from_id(I);
    while(true)
    {
        auto& rec = (*this)[i];
        path.push_back( (string)rec.name );
        if (i == root_index()) break;
        i = rec.parent_index;
    }
    return path;
}


const std::regex ott_version_pattern("^([0-9.]+)draft.*");

BaseTaxonomy::BaseTaxonomy(const string& dir,
                   bitset<32> cf,
                   OttId kr)
    :keep_root(kr),
     cleaning_flags(cf),
     path(dir),
     version(strip_trailing_whitespace(read_str_content_of_utf8_file(dir + "/version.txt"))) {
    std::smatch matches;
    if (std::regex_match(version, matches, ott_version_pattern)) {
        assert(matches.size() == 2);
        version_number = matches[1];
    } else {
        throw OTCError() << "Could not parse version number out of ott version string " << version;
    }
}

std::optional<OttId> BaseTaxonomy::get_unforwarded_id(OttId id) const
{
    auto id_or_reason = get_unforwarded_id_or_reason(id);
    if (std::holds_alternative<OttId>(id_or_reason))
        return std::get<OttId>(id_or_reason);
    else
        return {};
}

std::variant<OttId,reason_missing> Taxonomy::get_unforwarded_id_or_reason(OttId id) const
{
    if (index.count(id))
        return id;
    if (auto iter = forwards.find(id); iter != forwards.end())
        return iter->second;
    return reason_missing::unknown;
}

Taxonomy::Taxonomy(const string& dir,
                   bitset<32> cf,
                   OttId kr)
    :BaseTaxonomy(dir, cf, kr) {
    string filename = path + "/taxonomy.tsv";
    // 1. Open the file.
    ifstream taxonomy_stream(filename);
    if (not taxonomy_stream) {
        throw OTCError() << "Could not open file '" << filename << "'.";
    }
    // 2. Read and check the first line
    string line;
    std::getline(taxonomy_stream,line);
    unsigned int count;
    if (line == "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\t|\t") {
        count = read_input_taxonomy_stream(taxonomy_stream);
    } else if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t") {
        throw OTCError() << "First line of file '" << filename << "' is not a taxonomy header.";
    } else {
        count = read_ott_taxonomy_stream(taxonomy_stream);
    }
    LOG(TRACE) << "records read = " << count;
    LOG(TRACE) << "records kept = " << size();
    taxonomy_stream.close();
    /*
    if (read_deprecated) {
        read_deprecated_file(path + "/deprecated.tsv");
    }
    */
    read_forwards_file(path + "/forwards.tsv");
}

unsigned int Taxonomy::read_input_taxonomy_stream(std::istream & taxonomy_stream) {
    // 3. Read records up to the record containing the root.
    unsigned int count = 0;
    string line;

    if (keep_root != -1) {
        while(std::getline(taxonomy_stream, line)) {
            count++;
            // Add line to vector
            emplace_back(TaxonomyRecord{line, true});
            if (back().id == keep_root) {
                break;
            }
            pop_back();
        }
        if (empty()) {
            throw OTCError() << "Root id '" << keep_root << "' not found.";
        }
    } else {
        std::getline(taxonomy_stream, line);
        count++;
        // Add line to vector
        emplace_back(TaxonomyRecord{line, true});
    }
    back().depth = 1;
    if ((back().flags & cleaning_flags).any()) {
        throw OTCError() << "Root taxon (ID = " << back().id << ") removed according to cleaning flags!";
    }
    index[back().id] = size() - 1;
    // 4. Read the remaining records
    while(std::getline(taxonomy_stream, line)) {
        count++;
        if (count % 100 == 0) {
            std::cerr << "line " << count << ": " << line << '\n';
        }
        // Add line to vector
        emplace_back(TaxonomyRecord{line, true});
        // Eliminate records that match the cleaning flags
        if ((back().flags & cleaning_flags).any()) {
            pop_back();
            continue;
        }
        // Eliminate records whose parents have been eliminated, or are not found.
        auto loc = index.find(back().parent_id);
        if (loc == index.end()) {
            pop_back();
            continue;
        }
        back().parent_index = loc->second;
        back().depth = (*this)[back().parent_index].depth + 1;
        (*this)[back().parent_index].out_degree++;
        index[back().id] = size() - 1;
    }
    return count;
}

unsigned int Taxonomy::read_ott_taxonomy_stream(std::istream & taxonomy_stream) {
    // 3. Read records up to the record containing the root.
    unsigned int count = 0;
    string line;
    if (keep_root != -1) {
        while(std::getline(taxonomy_stream, line)) {
            count++;
            // Add line to vector
            emplace_back(line);
            if (back().id == keep_root) {
                break;
            }
            pop_back();
        }
        if (empty()) {
            throw OTCError() << "Root id '" << keep_root << "' not found.";
        }
    } else {
        std::getline(taxonomy_stream, line);
        count++;
        // Add line to vector
        emplace_back(line);
    }
    back().depth = 1;
    if ((back().flags & cleaning_flags).any()) {
        throw OTCError() << "Root taxon (ID = " << back().id << ") removed according to cleaning flags!";
    }
    index[back().id] = size() - 1;
    // 4. Read the remaining records
    while(std::getline(taxonomy_stream, line)) {
        count++;
        // Add line to vector
        emplace_back(line);
        // Eliminate records that match the cleaning flags
        if ((back().flags & cleaning_flags).any()) {
            pop_back();
            continue;
        }
        // Eliminate records whose parents have been eliminated, or are not found.
        auto loc = index.find(back().parent_id);
        if (loc == index.end()) {
            pop_back();
            continue;
        }
        back().parent_index = loc->second;
        back().depth = (*this)[back().parent_index].depth + 1;
        (*this)[back().parent_index].out_degree++;
        index[back().id] = size() - 1;
    }
    return count;
}

std::variant<OttId,reason_missing> RichTaxonomy::get_unforwarded_id_or_reason(OttId id) const
{
    const auto & td = tree->get_data();
    if (td.id_to_node.count(id)) {
        return id;
    }
    if (auto iter = forwards.find(id); iter != forwards.end()) {
        return iter->second;
    }
    return reason_missing::unknown;
}


RichTaxonomy::RichTaxonomy(const std::string& dir, std::bitset<32> cf, OttId kr)
    :BaseTaxonomy(dir, cf, kr) {
    { //braced to reduce scope of light_taxonomy to reduced memory
        Taxonomy light_taxonomy(dir, cf, kr); 
        auto nodeNamer = [](const auto&){return string();};
        cerr << "light_taxonomy.get_tree<RichTaxTree>(nodeNamer)..." << std::endl;
        tree = light_taxonomy.get_tree<RichTaxTree>(nodeNamer);
        cerr << "... tree returned" << std::endl;
        auto & tree_data = tree->get_data();
        std::swap(forwards, light_taxonomy.forwards);
        //std::swap(deprecated, light_taxonomy.deprecated);
        std::swap(keep_root, light_taxonomy.keep_root);
        std::swap(cleaning_flags, light_taxonomy.cleaning_flags);
        std::swap(path, light_taxonomy.path);
        std::swap(version, light_taxonomy.version);
        std::swap(version_number, light_taxonomy.version_number);
        for (auto tr_it = light_taxonomy.begin(); tr_it != light_taxonomy.end(); ++tr_it) {
            const auto & tr = *tr_it;
            auto ott_id = tr.id;
            if (tree_data.id_to_node.count(ott_id) == 0) {
                filtered_records.push_back(TaxonomyRecord(tr.line));
            }
        }
        for (auto tr_it = filtered_records.begin(); tr_it != filtered_records.end(); ++tr_it) {
            const auto & tr = *tr_it;
            register_taxon_in_maps(tree_data.name_to_record,
                                   tree_data.homonym_to_record,
                                   tr.name,
                                   tr.uniqname,
                                   &tr);
            tree_data.id_to_record[tr.id] = &tr;
        }
    }
    compute_depth(*tree);
    set_traversal_entry_exit(*tree);
    _fill_ids_to_suppress_set();
    this->read_synonyms();
    const auto & td = tree->get_data();
    // LOG(INFO) << "# of taxa stored in taxonomy, but filtered from taxonomy tree = " << filtered_records.size();
    // LOG(INFO) << "last # in ncbi_id_map = " << (td.ncbi_id_map.empty() ? 0 : max_numeric_key(td.ncbi_id_map));
    // LOG(INFO) << "last # in gbif_id_map = " <<  (td.gbif_id_map.empty() ? 0 : max_numeric_key(td.gbif_id_map));
    // LOG(INFO) << "last # in worms_id_map = " <<  (td.worms_id_map.empty() ? 0 : max_numeric_key(td.worms_id_map));
    // LOG(INFO) << "last # in if_id_map = " <<  (td.if_id_map.empty() ? 0 : max_numeric_key(td.if_id_map));
    // LOG(INFO) << "last # in irmng_id_map = " <<  (td.irmng_id_map.empty() ? 0 : max_numeric_key(td.irmng_id_map));
}


void Taxonomy::read_forwards_file(string filepath)
{
    // 1. Read forwards file and create id -> forwarded_id map
    ifstream forwards_stream(filepath);
    int i = 1;
    if (forwards_stream) {
        string line;
        std::getline(forwards_stream, line);
        while(std::getline(forwards_stream, line)) {
            char* temp;
            long old_id = std::strtoul(line.c_str(), &temp, 10);
            if (*temp != '\t') {
                throw OTCError() << "Expecting a tab after first Id in forwards file, " << filepath << "\nGot:\n" << line ;
            }
            const char* temp2 = temp+1;
            long new_id = std::strtoul(temp2, &temp, 10);
            cerr << i++ << ": " << old_id << " -> " << new_id << '\n';
            forwards[check_ott_id_size(old_id)] = check_ott_id_size(new_id);
        }
    }


    // 2. walk through the full forwards table, and try to find any that are multiple
    //    step paths of forwards.
    unordered_set<OttId> failed_forwards;
    for (auto& [old_id,new_id] : forwards)
    {
        unordered_set<OttId> visited;
        visited.insert(old_id);
        visited.insert(new_id);

        // Iterate the new_id
        assert(new_id > 0);
        while(forwards.count(new_id))
        {
            new_id = forwards.at(new_id);
            assert(new_id > 0);
            if (visited.count(new_id))
                throw OTCError()<<"forwarding loop from id "<<old_id<<"!";
            visited.insert(new_id);
        }

        if (not index.count(new_id))
        {
            LOG(DEBUG) << "OTT id "<<old_id<<" forwarded to non-existent (possibly pruned) id "<<new_id;
            failed_forwards.insert(old_id);
        }

        assert(new_id != -1);
        assert(new_id >= 0);
    }

    // 3. Remove forwards that forward to ids that are non-existent for any reason
    for(auto id: failed_forwards)
        forwards.erase(id);
}

std::string get_taxonomy_dir(const variables_map& args) {
    if (args.count("taxonomy")) {
        return args["taxonomy"].as<string>();
    }
    OTCError E;
    E << "Taxonomy dir not specified on command line.";
    vector<string> config_files;
    if (args.count("config")) {
        config_files.push_back(args["config"].as<string>());
    }
    auto dot_file = dot_opentree();
    if (not dot_file and not std::getenv("HOME")){
        E << "\n  Not looking in ~/.opentree: $HOME is not set";
    } else if (not dot_file) {
        E << "\n  Not looking in ~/.opentree: cannot open file";
    } else {
        config_files.push_back(*dot_file);
    }
    auto dir = load_config(config_files,"opentree","ott");
    if (not dir) {
        if (config_files.empty()) {
            E << "\n  No config files specified";
        } else {
            for(const auto& f: config_files) {
                E << "\n  '" << f << "': No variable ott in section [opentree]";
            }
        }
        throw E;
    }
    return *dir;
}

OttId root_ott_id_from_file(const string& filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    try {
        return pt.get<OttId>("synthesis.root_ott_id");
    } catch (...) {
        return -1;
    }
}

Taxonomy load_taxonomy(const variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    OttId keep_root = -1;
    if (args.count("root")) {
        keep_root = args["root"].as<OttId>();
    } else if (args.count("xroot")) {
        keep_root = args["xroot"].as<OttId>();
    } else if (args.count("config")) {
        keep_root = root_ott_id_from_file(args["config"].as<string>());
    }
    bitset<32> cleaning_flags = 0;
    if (args.count("config")) {
        cleaning_flags = cleaning_flags_from_config_file(args["config"].as<string>());
    }
    if (args.count("clean")) {
        cleaning_flags = flags_from_string(args["clean"].as<string>());
    }
    return {taxonomy_dir, cleaning_flags, keep_root};
}

RichTaxonomy load_rich_taxonomy(const variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    OttId keep_root = -1;
    if (args.count("root")) {
        keep_root = args["root"].as<OttId>();
    } else if (args.count("config")) {
        keep_root = root_ott_id_from_file(args["config"].as<string>());
    }
    bitset<32> cleaning_flags = 0;
    if (args.count("config")) {
        cleaning_flags = cleaning_flags_from_config_file(args["config"].as<string>());
    }
    if (args.count("clean")) {
        cleaning_flags = flags_from_string(args["clean"].as<string>());
    }
    return {taxonomy_dir, cleaning_flags, keep_root};
}

bool RTRichTaxNodeData::is_extinct() const {
    return ::is_extinct(flags);
}


void RichTaxonomy::read_synonyms() {
    string filename = path + "/synonyms.tsv";
    ifstream synonyms_file(filename);
        if (not synonyms_file) {
        throw OTCError() << "Could not open file '" << filename << "'.";
    }
    // 2. Read and check the first line
    string line;
    std::getline(synonyms_file, line);
    if (line == "uid\t|\tname\t|\ttype\t|\t") {
        read_input_synonyms_stream(synonyms_file);
    } else  if (line != "name\t|\tuid\t|\ttype\t|\tuniqname\t|\tsourceinfo\t|\t") {
        throw OTCError() << "First line of file '" << filename << "' is not a synonym header.";
    } else {
        read_ott_synonyms_stream(synonyms_file);
    }
}

void RichTaxonomy::read_input_synonyms_stream(std::istream & synonyms_file) {
    RTRichTaxTreeData & tree_data = this->tree->get_data();
    string line;
    while(std::getline(synonyms_file, line)) {
        const char* start[3];
        const char* end[3];
        start[0] = line.c_str();
        for(int i=0; i < 2; i++) {
            end[i] = std::strstr(start[i],"\t|\t");
            start[i + 1] = end[i] + 3;
        }
        end[2] = start[0] + line.length() - 3 ; // -3 for the \t|\t
        char *temp;
        string name = string(start[1], end[1] - start[1]);
        unsigned long raw_id = std::strtoul(start[0], &temp, 10);
        OttId ott_id = check_ott_id_size(raw_id);
        const RTRichTaxNode * primary = tree_data.id_to_node.at(ott_id);
        string sourceinfo;
        
        this->synonyms.emplace_back(name, primary, sourceinfo);
        TaxonomicJuniorSynonym & tjs = *(this->synonyms.rbegin());
        
        auto vs = comma_separated_as_vec(sourceinfo);
        process_source_info_vec(vs, tree_data, tjs, primary);
        RTRichTaxNode * mp = const_cast<RTRichTaxNode *>(primary);
        mp->get_data().junior_synonyms.push_back(&tjs);
    }
}

void RichTaxonomy::read_ott_synonyms_stream(std::istream & synonyms_file) {
    RTRichTaxTreeData & tree_data = this->tree->get_data();
    string line;
    while(std::getline(synonyms_file, line)) {
        const char* start[5];
        const char* end[5];
        start[0] = line.c_str();
        for(int i=0; i < 4; i++) {
            end[i] = std::strstr(start[i],"\t|\t");
            start[i + 1] = end[i] + 3;
        }
        end[4] = start[0] + line.length() - 3 ; // -3 for the \t|\t
        char *temp;
        string name = string(start[0], end[0] - start[0]);
        unsigned long raw_id = std::strtoul(start[1], &temp, 10);
        OttId ott_id = check_ott_id_size(raw_id);
        const RTRichTaxNode * primary = tree_data.id_to_node.at(ott_id);
        string sourceinfo = string(start[4], end[4] - start[4]);
        
        this->synonyms.emplace_back(name, primary, sourceinfo);
        TaxonomicJuniorSynonym & tjs = *(this->synonyms.rbegin());
        
        auto vs = comma_separated_as_vec(sourceinfo);
        process_source_info_vec(vs, tree_data, tjs, primary);
        RTRichTaxNode * mp = const_cast<RTRichTaxNode *>(primary);
        mp->get_data().junior_synonyms.push_back(&tjs);
    }
}

void RichTaxonomy::_fill_ids_to_suppress_set() {
    for (const auto nd : iter_node_const(*tree)) {
        if (node_is_suppressed_from_tnrs(nd)) {
            ids_to_suppress_from_tnrs.insert(nd->get_ott_id());
        }
    }
}

string format_with_taxonomy(const string& orig, const string& format, const TaxonomyRecord& rec, const Taxonomy& taxonomy) {
    string result;
    int pos = 0;
    do {
        auto loc = format.find('%', pos);
        if (loc == string::npos) {
            result += format.substr(pos);
            break;
        }
        result += format.substr(pos, loc - pos);
        loc++;
        const auto nc = format[loc];
        if (nc == 0) {
            std::abort();
        } 
        if (nc == 'I') {
            result += std::to_string(rec.id);
        } else if (nc == 'N') {
            result += rec.name;
        } else if (nc == 'U') {
            result += rec.uniqname;
        } else if (nc == 'R') {
            result += rec.rank;
        } else if (nc == 'P') {
            vector<string> path = taxonomy.path_from_id(rec.id);
            for(std::size_t i=0;i < path.size(); i++) {
                result += path[i];
                if (i != path.size()-1) {
                    result += " < ";
                }
            }
        } else if (nc == 'F') {
            result += flags_to_string(rec.flags);
        } else if (format[loc] == 'S') {
            result += rec.sourceinfo;
        } else if (format[loc] == 'L') {
            result += orig;
        } else if (format[loc] == '%') {
            result += '%';
        } else {
            throw OTCError() << "Invalid format specification '%" << nc << "' in taxonomy-based format string '" << format <<"'";
        }
        pos = loc + 1;
    } while (pos < static_cast<int>(format.size()));
    return result;
}

char format_needs_taxonomy(const string& format) {
    int pos = 0;
    do {
        auto loc = format.find('%', pos);
        if (loc == string::npos) {
            break;
        }
        loc++;
        const auto nc = format[loc];
        if (nc == 0) {
            std::abort();
        }
        if (nc == 'I' 
            || nc == 'N'
            || nc == 'U'
            || nc == 'R'
            || nc == 'P'
            || nc == 'F'
            || nc == 'S') {
            return nc;
        } else if (nc == 'L' || nc == '%') {
            ; // pass
        }
        else {
            throw OTCError() << "Invalid format specification '%" << nc << "' in format string '" << format << "'";
        }
        pos = loc + 1;
    } while (pos < static_cast<int>(format.size()));
    return false;
}


string format_without_taxonomy(const string& orig, const string& format) {
    string result;
    int pos = 0;
    do {
        auto loc = format.find('%', pos);
        if (loc == string::npos) {
            result += format.substr(pos);
            break;
        }
        result += format.substr(pos, loc - pos);
        loc++;
        const auto nc = format[loc];
        if (nc == 0) {
            std::abort();
        }
        if (nc == 'L') {
            result += orig;
        } else if (format[loc] == '%') {
            result += '%';
        } else {
            throw OTCError() << "Invalid format specification '%" << nc << "' in non-taxonomy-based format string '" << format <<"'";
        }
        pos = loc + 1;
    } while (pos < static_cast<int>(format.size()));
    return result;
}


string extract_long_as_string(const json & j, string opt_name) {
    long r;
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_number()) {
            r = opt->get<long>();
            string rs = std::to_string(r);
            return rs;
        }
    }
    throw OTCError() << "Expecting a \"" << opt_name << "\" property that refers to an integer.";
}

string extract_string(const json & j, string opt_name) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_number()) {
            return opt->get<string>();
        }
    }
    throw OTCError() << "Expecting a \"" << opt_name << "\" property that refers to an string.";
}


// BDR: factored this code out of taxonomy_mrca_ws_method below for use in tnrs
const RTRichTaxNode* taxonomy_mrca(const std::vector<const RTRichTaxNode*>& nodes)
{
    if (nodes.empty()) {
        return nullptr;
    }
    auto focal = nodes[0];
    for(auto& node: nodes) {
        focal = find_mrca_via_traversal_indices(focal, node);
        if (not focal) {
            throw OTCError() << "MRCA of taxa was not found. Please report this bug!\n";
        }
    }
    return focal;
}

// FIXME: move this out of here

std::vector<const RTRichTaxNode*> exact_name_search(const RichTaxonomy& taxonomy,
                                                    const RTRichTaxNode* context_root,
                                                    const std::string& query,
                                                    bool include_suppressed)
{
    if (include_suppressed) {
        return exact_name_search(taxonomy, context_root, query);
    }
    std::function<bool(const RTRichTaxNode*)> ok = [&](const RTRichTaxNode* taxon) {
        return not taxonomy.node_is_suppressed_from_tnrs(taxon);
    };
    return exact_name_search(taxonomy, context_root, query, ok);
}


std::vector<const RTRichTaxNode*> exact_name_search(const RichTaxonomy& taxonomy,
                                                    const std::string& query,
                                                    bool include_suppressed)
{
    const RTRichTaxNode* context_root = taxonomy.get_tax_tree().get_root();
    return exact_name_search(taxonomy, context_root, query, include_suppressed);
}


vector<const RTRichTaxNode *> exact_name_search_slow(const RichTaxonomy& /*taxonomy*/,
                                                     const RTRichTaxNode* context_root,
                                                     const std::string&  query_ref,
                                                     std::function<bool(const RTRichTaxNode*)> ok)
{
    std::string query{query_ref};
    for (auto& c: query) {
        c = std::tolower(c);
    }
    vector<const RTRichTaxNode*> hits;
    for(auto taxon: iter_post_n_const(*context_root)) {
        if (not ok(taxon)) {
            continue;
        }
        if (lcase_string_equals(query, taxon->get_data().get_nonuniqname())) {
            hits.push_back(taxon);
        }
    }
    return hits;
}

vector<const RTRichTaxNode *> exact_name_search(const RichTaxonomy& taxonomy,
                                                const RTRichTaxNode* context_root,
                                                const std::string&  query_ref,
                                                std::function<bool(const RTRichTaxNode*)> ok)
{
    auto ctp = taxonomy.get_fuzzy_matcher();
    assert(ctp);

    auto results = ctp->to_taxa(ctp->exact_query(query_ref), context_root, taxonomy, true);
    vector<const RTRichTaxNode*> hits;
    for(auto& result: results)
    {
        if (not result.is_synonym())
        {
            auto t = result.get_taxon();
            if (ok(t))
                hits.push_back(t);
        }
    }

#ifdef DEBUG_NAME_SEARCH
    {
        auto hits2 = exact_name_search_slow(taxonomy, context_root, query_ref, ok);
        std::sort(hits.begin(), hits.end());
        std::sort(hits2.begin(), hits2.end());
        LOG(INFO)<<"exact_name_search: query = '"<<query_ref<<"'  context_id = "<<context_root->get_ott_id();
        if (hits != hits2)
        {
            LOG(INFO)<<"ctrie match:";
            for(int i=0;i<hits.size();i++)
                LOG(INFO)<<"   "<<hits[i]->get_data().get_nonuniqname();
            LOG(INFO)<<"lcase match:";
            for(int i=0;i<hits2.size();i++)
                LOG(INFO)<<"   "<<hits2[i]->get_data().get_nonuniqname();
        }
        else
            LOG(INFO)<<"exact name search: "<<hits.size()<<" names agree";
    }
#endif

    return hits;
}


} //namespace otc
