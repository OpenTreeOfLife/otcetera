#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#include "otc/taxonomy/patching.h"

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
// namespace po = boost::program_options;
// using po::variables_map;

namespace otc
{
PatchableTaxonomy load_patchable_taxonomy(const boost::program_options::variables_map& args) {
    string taxonomy_dir = get_taxonomy_dir(args);
    OttId keep_root = -1;
    bitset<32> cleaning_flags = 0;
    return {taxonomy_dir, cleaning_flags, keep_root};
}

PatchableTaxonomy::PatchableTaxonomy(const std::string& dir,
                                     std::bitset<32> cf,
                                     OttId kr)
    :RichTaxonomy(dir, cf, kr) {
    std::cerr << filtered_records.size() << " filtered_records" << std::endl;
}

std::pair<bool, std::string> PatchableTaxonomy::add_new_taxon(OttId oid,
                                                         OttId parent_id,
                                                         const std::string & name,
                                                         const std::string & rank,
                                                         const std::string & sourceinfo,
                                                         const std::string & uniqname,
                                                         const std::string & flags) {
    vector<string> elements;
    elements.reserve(8);
    elements.push_back(std::to_string(oid));
    elements.push_back(std::to_string(parent_id));
    elements.push_back(name);
    elements.push_back(rank);
    elements.push_back(sourceinfo);
    elements.push_back(uniqname);
    elements.push_back(flags);
    elements.push_back(string());
    string fake_line = boost::algorithm::join(elements, "\t|\t");
    return std::pair<bool, std::string>{false, "not implemented"}; //add_taxon_record(fake_line);
}

void PatchableTaxonomy::write_version_file_contents(std::ostream & out) const {
    out << version << std::endl;
}

void PatchableTaxonomy::write_taxonomy_file_contents(std::ostream & tf) const {
    tf << "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t" << std::endl;

    // string sep = "\t|\t";
    // for (auto& rec: *this) {
    //    tf << rec.id << sep;
    //    if (rec.parent_id > 0) {
    //        tf << rec.parent_id;
    //    }
    //    tf << sep;
    //    tf << rec.name << sep;
    //    tf << rec.rank << sep;
    //    tf << rec.sourceinfo << sep;
    //    if (rec.uniqname != rec.name) {
    //        tf << rec.uniqname;
    //    }
    //    tf << sep;
    //    tf << flags_to_string(rec.flags) << sep;
    //    tf << '\n';
    // }
}

void PatchableTaxonomy::write_synonyms_file_contents(std::ostream & sf) const {
    sf << "name\t|\tuid\t|\ttype\t|\tuniqname\t|\tsourceinfo\t|\t" << std::endl;
    // string sep = "\t|\t";
    // for (auto& rec: *this) {
    //    tf << rec.id << sep;
    //    if (rec.parent_id > 0) {
    //        tf << rec.parent_id;
    //    }
    //    tf << sep;
    //    tf << rec.name << sep;
    //    tf << rec.rank << sep;
    //    tf << rec.sourceinfo << sep;
    //    if (rec.uniqname != rec.name) {
    //        tf << rec.uniqname;
    //    }
    //    tf << sep;
    //    tf << flags_to_string(rec.flags) << sep;
    //    tf << '\n';
    // }
}

void PatchableTaxonomy::write_forwards_file_contents(std::ostream & ff) const {
    ff << "id\treplacement\n";
    for(const auto& p: forwards) {
        ff << p.first << '\t' << p.second << '\n';
    }
}


void PatchableTaxonomy::write(const std::string& newdirname) const{
    fs::path old_dir = path;
    fs::path new_dir = newdirname;
    if (! fs::exists(new_dir)) {
        fs::create_directories(new_dir);
    }
    {
        ofstream version_file((new_dir/"version.txt").string());
        write_version_file_contents(version_file);
        version_file.close();
    }
    {
        ofstream tf ((new_dir/"taxonomy.tsv").string());
        write_taxonomy_file_contents(tf);
        tf.close();
    }
    {
        ofstream sf((new_dir/"synonyms.tsv").string());
        write_synonyms_file_contents(sf);
        sf.close();
    }
    {
        ofstream ff((new_dir/"forwards.tsv").string());
        write_forwards_file_contents(ff);
        ff.close();
    }
}

void PatchableTaxonomy::write_to_stream(std::ostream & out) const {
    out << "version.txt\n";
    write_version_file_contents(out);
    out << "taxonomy.tsv\n";
    write_taxonomy_file_contents(out);
    out << "synonyms.tsv\n";
    write_synonyms_file_contents(out);
    out << "forwards.tsv";
    write_forwards_file_contents(out);
}

} //namespace otc
