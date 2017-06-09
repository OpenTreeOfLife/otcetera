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
using std::ifstream;
using std::unordered_set;

using boost::string_ref;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

namespace otc
{


const map<string, TaxonomicRank> rank_name_to_enum = 
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

const map<TaxonomicRank, string> rank_enum_to_name = 
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
const std::string empty_string;
const set<string> indexed_source_prefixes = {"ncbi", "gbif", "worms", "if", "irmng"};
std::set<std::string> rank_strings;

template<typename T>
void register_taxon_in_maps(std::map<boost::string_ref, const T *> & n2n,
                            std::map<boost::string_ref, std::vector<const T *> > & homonym_map,
                            boost::string_ref possibly_nonunique_name,
                            boost::string_ref uname,
                            const T * ti) {
    auto nit = n2n.lower_bound(possibly_nonunique_name);
    typedef std::pair<boost::string_ref, const T *> name_map_pair;
    if (nit->first != possibly_nonunique_name) {
        nit = n2n.insert(nit, name_map_pair(possibly_nonunique_name, ti));
    } else {
        if (nit->second != nullptr) {
            homonym_map[possibly_nonunique_name].push_back(nit->second);
            nit->second = nullptr;
        }
       homonym_map[possibly_nonunique_name].push_back(ti);
    }
    if (uname != possibly_nonunique_name) {
        auto r2 = n2n.insert(name_map_pair(uname, ti));
        assert(r2.second); // should be uniq.
    }
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

const TaxonomyRecord& Taxonomy::record_from_id(OttId id) const {
    auto loc = index.find(id);
    if (loc == index.end()) {
        auto loc2 = forwards.find(id);
        if (loc2 == forwards.end()) {
            throw OTCError() << "ID " << id << " not in taxonomy or forwarding list";
        }
        OttId newid = loc2->second;
        loc = index.find(newid);
        // If id is in the forwarding table, then newid should be in the taxonomy
        assert(loc != index.end());
    }
    return (*this)[loc->second];
}

TaxonomyRecord& Taxonomy::record_from_id(OttId id) {
    auto loc = index.find(id);
    if (loc == index.end()) {
        auto loc2 = forwards.find(id);
        if (loc2 == forwards.end()) {
            throw OTCError() << "ID " << id <<" not in taxonomy or forwarding list";
        }
        OttId newid = loc2->second;
        loc = index.find(newid);
        // If id is in the forwarding table, then newid should be in the taxonomy
        assert(loc != index.end());
    }
    return (*this)[loc->second];
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

void Taxonomy::write(const std::string& newdirname) {
    fs::path old_dir = path;
    fs::path new_dir = newdirname;
    if (fs::exists(new_dir)) {
        throw OTCError() << "File '" << newdirname << "' already exists!";
    }
    fs::create_directories(new_dir);
    // Copy the other files.
    for(const auto& name: {"about.json", "conflicts.tsv", "deprecated.tsv",
                "log.tsv", "otu_differences.tsv", "synonyms.tsv", "weaklog.csv"})
    {
	if (fs::exists(old_dir/name))
	    fs::copy_file(old_dir/name,new_dir/name);
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

RichTaxonomy::RichTaxonomy(const std::string& dir, std::bitset<32> cf, OttId kr)
    :BaseTaxonomy(dir, cf, kr) {
    { //braced to reduce scope of light_taxonomy to reduced memory
        Taxonomy light_taxonomy(dir, cf, kr); 
        auto nodeNamer = [](const auto&){return string();};
        tree = light_taxonomy.get_tree<RichTaxTree>(nodeNamer);
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
    LOG(INFO) << "# of taxa stored in taxonomy, but filtered from taxonomy tree = " << filtered_records.size();
    LOG(INFO) << "last # in ncbi_id_map = " << (td.ncbi_id_map.empty() ? 0 : max_numeric_key(td.ncbi_id_map));
    LOG(INFO) << "last # in gbif_id_map = " <<  (td.gbif_id_map.empty() ? 0 : max_numeric_key(td.gbif_id_map));
    LOG(INFO) << "last # in worms_id_map = " <<  (td.worms_id_map.empty() ? 0 : max_numeric_key(td.worms_id_map));
    LOG(INFO) << "last # in if_id_map = " <<  (td.if_id_map.empty() ? 0 : max_numeric_key(td.if_id_map));
    LOG(INFO) << "last # in irmng_id_map = " <<  (td.irmng_id_map.empty() ? 0 : max_numeric_key(td.irmng_id_map));
}


void Taxonomy::read_forwards_file(string filepath) {
    ifstream forwards_stream(filepath);
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
            forwards[check_ott_id_size(old_id)] = check_ott_id_size(new_id);
        }
    }
    // walk through the full forwards table, and try to find any that are multiple
    //    step paths of forwards.
    unordered_set<OttId> need_iterating;
    for (auto old_new : forwards) {
        auto new_id = old_new.second;
        if (new_id >= 0) {
            OttId nnid = this->map(new_id);
            if (nnid != new_id) {
                need_iterating.insert(old_new.first);
            }
        }
    }
    while (!need_iterating.empty()) {
        unordered_set<OttId> scratch;
        for (auto old_id: need_iterating) {
            auto fm_it = forwards.find(old_id);
            assert(fm_it != forwards.end());
            auto curr_new_id = fm_it->second;
            OttId nnid = this->map(fm_it->second);
            assert(nnid != fm_it->second);
            fm_it->second = nnid;
            if (nnid > 0 && nnid != this->map(nnid)) {
                scratch.insert(old_id);
            }
        }
        need_iterating = scratch;
    }
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

template<typename T>
void process_source_info_vec(const std::vector<std::string> & vs,
                             RTRichTaxTreeData & tree_data,
                             T & ,
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
            long raw_foreign_id  = std::stoul(id_str.c_str(), &pos);
            if (pos < id_str.length() || raw_foreign_id < 0) {
                throw OTCError() << "Could not convert ID to unsigned long \"" << src_entry << "\"";
            }
            OttId foreign_id = check_ott_id_size(raw_foreign_id);
#           if defined(MAP_FOREIGN_TO_POINTER)
                auto to_map_to = this_node;
#           else
                auto to_map_to = this_node->get_ott_id();
#           endif

            if (prefix == "ncbi") {
                tree_data.ncbi_id_map[foreign_id] = to_map_to;
            } else if (prefix == "gbif") {
                tree_data.gbif_id_map[foreign_id] = to_map_to;
            } else if (prefix == "worms") {
                tree_data.worms_id_map[foreign_id] = to_map_to;
            } else if (prefix == "if") {
                tree_data.if_id_map[foreign_id] = to_map_to;
            } else if (prefix == "irmng") {
                tree_data.irmng_id_map[foreign_id] = to_map_to;
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
inline void populate_node_from_taxonomy_record(RTRichTaxNode & nd,
                                           const TaxonomyRecord & tr,
                                           std::function<std::string(const TaxonomyRecord&)>,
                                           RichTaxTree & tree) {
    RTRichTaxNode * this_node = &nd;
    nd.set_ott_id(tr.id);
    auto & data = nd.get_data();
    auto & tree_data = tree.get_data();
    nd.set_ott_id(tr.id);
    tree_data.id_to_node[tr.id] = this_node;
    this_node->set_name(string(tr.uniqname));
    const string & uname = this_node->get_name();
    if (tr.uniqname != tr.name) {
        string sn = string(tr.name);
        tree_data.non_unique_taxon_names[sn].insert(tr.id);
        auto nit = tree_data.non_unique_taxon_names.find(sn);
        assert(nit != tree_data.non_unique_taxon_names.end());
        data.possibly_nonunique_name = string_ref(nit->first);
    } else {
        data.possibly_nonunique_name = string_ref(nd.get_name());
    }
    data.flags = tr.flags;
    data.rank = rank_name_to_enum.at(string(tr.rank));
    register_taxon_in_maps(tree_data.name_to_node,
                           tree_data.homonym_to_node,
                           data.possibly_nonunique_name,
                           uname,
                           this_node);
    auto flags = data.get_flags();
    // If the flag combination is new, store the JSON representation
    if (tree_data.flags2json.count(flags) == 0) {
        vector<string> vf = flags_to_string_vec(flags);
        tree_data.flags2json[flags] = json();
        auto & fj = tree_data.flags2json[flags];
        for (auto fs : vf) {
            fj.push_back(fs);
        }
    }
    auto vs = tr.sourceinfoAsVec();
    data.source_info = string(tr.sourceinfo);
    process_source_info_vec(vs, tree_data, data, this_node);
}



void RichTaxonomy::read_synonyms() {
    RTRichTaxTreeData & tree_data = this->tree->get_data();
    string filename = path + "/synonyms.tsv";
    ifstream synonyms_file(filename);
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
        unsigned long raw_id = std::strtoul(start[1], &temp, 10);
        OttId ott_id = check_ott_id_size(raw_id);
        const RTRichTaxNode * primary = tree_data.id_to_node.at(ott_id);
        string sourceinfo = string(start[4], end[4] - start[4]);
        
        this->synonyms.emplace_back(name, primary, sourceinfo);
        TaxonomicJuniorSynonym & tjs = *(this->synonyms.rbegin());
        auto nit = tree_data.name_to_node.lower_bound(name);
        boost::string_ref name_ref = tjs.name;
        typedef std::pair<boost::string_ref, const RTRichTaxNode *> name_map_pair;
        if (nit->first != name_ref) {
            nit = tree_data.name_to_node.insert(nit, name_map_pair(name_ref, primary));
        } else {
            if (nit->second != nullptr) {
                tree_data.homonym_to_node[name_ref].push_back(nit->second);
                nit->second = nullptr;
            }
            tree_data.homonym_to_node[name_ref].push_back(primary);
        }
        
        auto vs = comma_separated_as_vec(sourceinfo);
        process_source_info_vec(vs, tree_data, tjs, primary);
        RTRichTaxNode * mp = const_cast<RTRichTaxNode *>(primary);
        mp->get_data().junior_synonyms.push_back(&tjs);
    }
}

void RichTaxonomy::_fill_ids_to_suppress_set() {
    const string sup_flag_comma = "not_otu,environmental,environmental_inherited,viral,hidden,hidden_inherited,was_container";
    const auto suppress_flags = flags_from_string(sup_flag_comma);
    for (const auto nd : iter_node_const(*tree)) {
        const auto & tax_record_flags = nd->get_data().get_flags();
        auto intersection = suppress_flags & tax_record_flags;
        if (intersection.any()) {
            const auto ott_id = nd->get_ott_id();
            ids_to_suppress_from_tnrs.insert(ott_id);
        }
    }
}

string format_with_taxonomy(const string& orig, const string& format, const TaxonomyRecord& rec) {
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

void RichTaxonomy::add_taxonomic_addition_string(const std::string &s) {
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
