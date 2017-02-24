// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
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
using nlohmann::json;

using boost::string_ref;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

namespace otc
{

std::set<std::string> rank_strings;

taxonomy_record::taxonomy_record(const string& line_)
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
//    boost::container::small_vector<string_ref,10> words;
//    for (string_ref&& r : iter_split(v, line, token_finder(is_any_of(","))) |
//             transformed([](R const& r){return boost::string_ref(&*r.begin(), r.size());})
//        )
//        words.push_back(r);

    char *temp;
    id = std::strtoul(start[0], &temp, 10);
    parent_id = std::strtoul(start[1], &temp, 10);
    name = string_ref(start[2], end[2] - start[2]);
    rank = string_ref(start[3], end[3] - start[3]);
    sourceinfo = string_ref(start[4], end[4] - start[4]);
    uniqname = string_ref(start[5], end[5] - start[5]);
    flags = flags_from_string(start[6], end[6]);
    if (not uniqname.size()) {
        uniqname = name;
    }
    rank_strings.insert(string(rank));
}

const taxonomy_record& Taxonomy::record_from_id(long id) const {
    auto loc = index.find(id);
    if (loc == index.end()) {
        auto loc2 = forwards.find(id);
        if (loc2 == forwards.end()) {
            throw OTCError() << "ID " << id << " not in taxonomy or forwarding list";
        }
        long newid = loc2->second;
        loc = index.find(newid);
        // If id is in the forwarding table, then newid should be in the taxonomy
        assert(loc != index.end());
    }
    int index = loc->second;
    return (*this)[index];
}

taxonomy_record& Taxonomy::record_from_id(long id) {
    auto loc = index.find(id);
    if (loc == index.end()) {
        auto loc2 = forwards.find(id);
        if (loc2 == forwards.end()) {
            throw OTCError() << "ID " << id <<" not in taxonomy or forwarding list";
        }
        long newid = loc2->second;
        loc = index.find(newid);
        // If id is in the forwarding table, then newid should be in the taxonomy
        assert(loc != index.end());
    }
    int index = loc->second;
    return (*this)[index];
}

long Taxonomy::map(long old_id) const {
    if (index.count(old_id)) {
        return old_id;
    }
    auto loc = forwards.find(old_id);
    if (loc != forwards.end()) {
        return loc->second;
    }
    return -1;
}

void Taxonomy::write(const std::string& newdirname) {
    fs::path old_dir = path;
    fs::path new_dir = newdirname;
    if (fs::exists(new_dir)) {
        throw OTCError() << "File '" << newdirname << "' already exists!";
    }
    fs::create_directories(new_dir);
    // Copy the other files.
    for(const auto& name: {"about.json", "conflicts.tsv", "deprecated.tsv",
                "log.tsv", "otu_differences.tsv", "synonyms.tsv", "weaklog.csv"}) {
        fs::copy_file(old_dir/name,new_dir/name);
    }
    // Write the new version file.
    {
        ofstream version_file((new_dir/"version.txt").string());
        version_file << version;
        if (keep_root != -1 or cleaning_flags.any()) {
            version_file<<"modified: ";
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
        tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t"<<std::endl;
        string sep = "\t|\t";
        for(const auto& r: *this) {
            tf << r.line <<"\n";
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


const std::regex ott_version_pattern("^([.0-9]+)draft.*");

Taxonomy::Taxonomy(const string& dir, bitset<32> cf, long kr)
    :keep_root(kr),
     cleaning_flags(cf),
     path(dir),
     version(strip_trailing_whitespace(readStrContentOfUTF8File(dir + "/version.txt"))) {
    std::smatch matches;
    if (std::regex_match(version, matches, ott_version_pattern)) {
        assert(matches.size() == 2);
        version_number = matches[1];
    } else {
        throw OTCError() << "Could not parse version number out of ott version string " << version;
    }

    string filename = path + "/taxonomy.tsv";
    // 1. Open the file.
    std::ifstream taxonomy_stream(filename);
    if (not taxonomy_stream) {
        throw OTCError() << "Could not open file '" << filename << "'.";
    }
    // 2. Read and check the first line
    string line;
    std::getline(taxonomy_stream,line);
    if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t") {
        throw OTCError() << "First line of file '" << filename << "' is not a taxonomy header.";
    }
    // 3. Read records up to the record containing the root.
    int count = 0;
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
        std::getline(taxonomy_stream,line);
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
    while(std::getline(taxonomy_stream,line)) {
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
    LOG(TRACE)<<"records read = "<<count;
    LOG(TRACE)<<"records kept = "<<size();
    taxonomy_stream.close();
    std::ifstream forwards_stream(path + "/forwards.tsv");
    if (forwards_stream) {
        std::getline(forwards_stream, line);
        while(std::getline(forwards_stream, line)) {
            char* temp;
            long old_id = std::strtoul(line.c_str(), &temp, 10);
            assert(*temp == '\t');
            const char* temp2 = temp+1;
            long new_id = std::strtoul(temp2, &temp, 10);
            if (index.count(new_id)) {
                forwards[old_id] = new_id;
            }
        }
    }

    //for (auto rs : rank_strings) {
    //    std::cerr << "rankstring: \"" << rs << "\"\n"; 
    //}
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
                E << "\n  '"<<f<<"': No variable ott in section [opentree]";
            }
        }
        throw E;
    }
    return *dir;
}

long root_ott_id_from_file(const string& filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    try {
        return pt.get<long>("synthesis.root_ott_id");
    } catch (...) {
        return -1;
    }
}

Taxonomy load_taxonomy(const variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    long keep_root = -1;
    if (args.count("root")) {
        keep_root = args["root"].as<long>();
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
    long keep_root = -1;
    if (args.count("root")) {
        keep_root = args["root"].as<long>();
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

const map<string, TaxonomicRank> rankNameToEnum = 
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

const map<TaxonomicRank, string> rankEnumToName = 
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
const std::string emptyStringForMissingRank;

const set<string> indexed_source_prefixes = {"ncbi", "gbif", "worms", "if", "irmng"};

template<typename T>
void process_source_info_vec(const std::vector<std::string> & vs,
                             RTRichTaxTreeData & tree_data,
                             T & data,
                             const RTRichTaxNode * this_node) {
    for (auto src_entry : vs) {
        auto pref_id = split_string(src_entry, ':');
        if (pref_id.size() != 2) {
            throw OTCError() << "Expecting exactly 1 colon in a source ID string. Found: \"" << src_entry << "\".";
        }
        const string & prefix = *pref_id.begin();
        if (indexed_source_prefixes.count(prefix) == 0) {
            continue;
        }
        const string & id_str = *pref_id.rbegin();
        //data.sources.push_back(src_entry);
        std::size_t pos;
        //LOG(INFO) << src_entry;
        try {
            unsigned long foreign_id  = std::stoul(id_str.c_str(), &pos);
            if (pos < id_str.length()) {
                throw OTCError() << "Could not convert ID to unsigned long \"" << src_entry << "\"";
            }
            if (prefix == "ncbi") {
                tree_data.ncbi_id_map[foreign_id] = this_node;
            } else if (prefix == "gbif") {
                tree_data.gbif_id_map[foreign_id] = this_node;
            } else if (prefix == "worms") {
                tree_data.worms_id_map[foreign_id] = this_node;
            } else if (prefix == "if") {
                tree_data.if_id_map[foreign_id] = this_node;
            } else if (prefix == "irmng") {
                tree_data.irmng_id_map[foreign_id] = this_node;
            } else {
                assert(false);
            }
        } catch (OTCError & x) {
            throw;
        } catch (...) {
            LOG(WARNING) << "Could not convert ID to unsigned long \"" << src_entry << "\"";
        }
    }
}


// default behavior is to set ID and Name from line
template <>
inline void populateNodeFromTaxonomyRecord(RTRichTaxNode & nd,
                                           const taxonomy_record & line,
                                           std::function<std::string(const taxonomy_record&)>,
                                           RichTaxTree & tree) {
    RTRichTaxNode * this_node = &nd;
    nd.setOttId(line.id);
    auto & data = nd.getData();
    auto & tree_data = tree.getData();
    nd.setOttId(line.id);
    tree_data.id2node[nd.getOttId()] = this_node;
    data.tax_record = &line;
    auto name = data.getName();
    auto nit = tree_data.name2node.lower_bound(name);
    typedef std::pair<boost::string_ref, const RTRichTaxNode *> name_map_pair;
    if (nit->first != name) {
        nit = tree_data.name2node.insert(nit, name_map_pair(name, this_node));
    } else {
        if (nit->second != nullptr) {
            tree_data.homonym2node[name].push_back(nit->second);
            nit->second = nullptr;
        }
        tree_data.homonym2node[name].push_back(this_node);
    }
    auto uname = data.getUniqname();
    if (uname != name) {
        auto r2 = tree_data.name2node.insert(name_map_pair(uname, this_node));
        assert(r2.second); // should be uniq.
    }
    auto flags = data.getFlags();
    // If the flag combination is new, store the JSON representation
    if (tree_data.flags2json.count(flags) == 0) {
        vector<string> vf = flags_to_string_vec(flags);
        tree_data.flags2json[flags] = json();
        auto & fj = tree_data.flags2json[flags];
        for (auto fs : vf) {
            fj.push_back(fs);
        }
    }
    /*
    const string rank = string(line.rank);
    if (rank == "natio") {
        LOG(WARNING) << "Converting rank natio to RANK_INFRASPECIFICNAME";
    }
    auto reit = rankNameToEnum.find(rank);
    if (reit == rankNameToEnum.end()) {
        throw OTCError() << "Rank string not recognized: \"" << rank << "\"";
    }
    data.rank = reit->second;
    */
    auto vs = line.sourceinfoAsVec();
    process_source_info_vec(vs, tree_data, data, this_node);
}



void RichTaxonomy::readSynonyms() {
    RTRichTaxTreeData & tree_data = this->tree->getData();
    string filename = path + "/synonyms.tsv";
    std::ifstream synonyms_file(filename);
        if (not synonyms_file) {
        throw OTCError() << "Could not open file '" << filename << "'.";
    }
    // 2. Read and check the first line
    string line;
    std::getline(synonyms_file, line);
    if (line != "name\t|\tuid\t|\ttype\t|\tuniqname\t|\tsourceinfo\t|\t") {
        throw OTCError() << "First line of file '" << filename << "' is not a synonym header.";
    }
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
        unsigned long id = std::strtoul(start[1], &temp, 10);
        const RTRichTaxNode * primary = tree_data.id2node.at(id);
        string sourceinfo = string(start[4], end[4] - start[4]);
        
        this->synonyms.emplace_back(name, primary, sourceinfo);
        TaxonomicJuniorSynonym & tjs = *(this->synonyms.rbegin());
        auto nit = tree_data.name2node.lower_bound(name);
        boost::string_ref name_ref = tjs.name;
        typedef std::pair<boost::string_ref, const RTRichTaxNode *> name_map_pair;
        if (nit->first != name_ref) {
            nit = tree_data.name2node.insert(nit, name_map_pair(name_ref, primary));
        } else {
            if (nit->second != nullptr) {
                tree_data.homonym2node[name_ref].push_back(nit->second);
                nit->second = nullptr;
            }
            tree_data.homonym2node[name_ref].push_back(primary);
        }
        
        auto vs = comma_separated_as_vec(sourceinfo);
        process_source_info_vec(vs, tree_data, tjs, primary);
        RTRichTaxNode * mp = const_cast<RTRichTaxNode *>(primary);
        mp->getData().junior_synonyms.push_back(&tjs);
    }
}

RichTaxonomy::RichTaxonomy(const std::string& dir, std::bitset<32> cf, long kr)
    :Taxonomy(dir, cf, kr) {
    auto nodeNamer = [](const auto&){return string();};
    this->tree = getTree<RichTaxTree>(nodeNamer);
    this->readSynonyms();
    // Could call:
    // index.clear(); 
    // to save about 8M RAM, but this disables some Taxonomy functionality! DANGEROUS move
}

string format_with_taxonomy(const string& orig, const string& format, const taxonomy_record& rec) {
    string result;
    int pos = 0;
    do {
        auto loc = format.find('%',pos);
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
            result += rec.name.to_string();
        } else if (nc == 'U') {
            result += rec.uniqname.to_string();
        } else if (nc == 'R') {
            result += rec.rank.to_string();
        } else if (nc == 'F') {
            result += flags_to_string(rec.flags);
        } else if (format[loc] == 'S') {
            result += rec.sourceinfo.to_string();
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
        auto loc = format.find('%',pos);
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

void RichTaxonomy::addTaxonomicAdditionString(const std::string &s) {
    json j;
    try {
        j = json::parse(s);
    } catch (std::exception & x) {
        throw OTCError() << "Error parsing JSON for taxonomic addition: " << x.what();
    }
    auto tax_arr_it = j.find("taxa");
    if (tax_arr_it == j.end() || !tax_arr_it->is_array()) {
        throw OTCError() << "Expecting a taxonomic addition to have \"taxa\" property referring to an array.";
    }
    auto tax_arr = j["taxa"];
    string amend_id = extract_string(j, "id");
    for (json::const_iterator tax_add_it = tax_arr.begin(); tax_add_it != tax_arr.end(); ++tax_add_it) {
        vector<string> elements;
        const json & taxon_addition = *tax_add_it;
        elements.reserve(8);
        string ott_id = extract_long_as_string(taxon_addition, "ott_id");
        string source = amend_id;
        source += ":";
        source += ott_id;
        elements.push_back(ott_id);
        elements.push_back(extract_long_as_string(taxon_addition, "parent"));
        elements.push_back(extract_string(taxon_addition, "name"));
        elements.push_back(extract_string(taxon_addition, "rank"));
        elements.push_back(source);
        elements.push_back(string());
        elements.push_back(string());
        elements.push_back(string());
        string fake_line = boost::algorithm::join(elements, "\t|\t");
        //process_taxonomy_line(fake_line);
    }
}

} //namespace otc
