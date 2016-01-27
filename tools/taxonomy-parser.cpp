#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <bitset>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

// 1. Write a parser to read the lines faster
// 2. Avoid memory allocation -- by mmapping the taxonomy file?
// 3. Convert the flags into a bitmask
// 4. Should the Rank be a converted to an integer?
// 5. Can we assign OTT IDs to internal nodes of a tree while accounting for Incertae Sedis taxa?
// * What are the triplet-inference rules for the Incertae Sedis problem?

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) 
{ 
  using namespace po;

  // named options
  options_description invisible("Invisible options");
  invisible.add_options()
    ("taxonomy", value<string>(),"Filename for the taxonomy")
    ;

  options_description visible("All options");
  visible.add_options()
      ("help,h", "Produce help message")
      ("config,c",value<string>(),"Config file containing flags to filter")
      ("clean",value<string>(),"Comma-separated string of flags to filter")
      ("write-tree,T","Write out the result as a tree")
      ("root,r", value<int>(), "OTT id of root node of subtree to keep")
//    ("quiet,q","QUIET mode (all logging disabled)")
//    ("trace,t","TRACE level debugging (very noisy)")
//    ("verbose,v","verbose")
    ;

  options_description all("All options");
  all.add(invisible).add(visible);

  // positional options
  positional_options_description p;
  p.add("taxonomy", -1);

  variables_map args;     
  store(command_line_parser(argc, argv).options(all).positional(p).run(), args);
  notify(args);    

  if (args.count("help")) {
    cout<<"Usage: taxonomy-parser <taxonomy-dir> [OPTIONS]\n";
    cout<<"Select columns from a Tracer-format data file.\n\n";
    cout<<visible<<"\n";
    exit(0);
  }

  return args;
}

// https://github.com/OpenTreeOfLife/taxomachine/blob/master/src/main/java/org/opentree/taxonomy/OTTFlag.java

// http://www.boost.org/doc/libs/1_50_0/libs/spirit/doc/html/spirit/qi/reference/string/symbols.html

auto get_symbols()
{
    symbols<char, int> sym;

    sym.add
        ("not_otu", 0)
        ("environmental", 1)
        ("environmental_inherited", 2)
        ("viral", 3)
        ("hidden", 4)
        ("hidden_inherited", 5)
//        unclassified_direct
        ("was_container", 6)
        ("barren", 8)
        ("extinct", 9)
//        extinct_direct)
        ("extinct_inherited", 11)
        ("major_rank_conflict", 12)
//        major_rank_conflict_direct
        ("major_rank_conflict_inherited", 14)
        ("unclassified", 15)
        ("unclassified_inherited", 16)
        ("edited", 17)
        ("hybrid", 18)
        ("incertae_sedis", 19)
        ("incertae_sedis_inherited", 20)
//	      incertae_sedis_direct
        ("infraspecific", 22)
        ("sibling_lower", 23)
        ("sibling_higher", 24)
        ("tattered", 25)
        ("tattered_inherited", 26)
        ("forced_visible", 27)
        ("unplaced", 28)
        ("unplaced_inherited", 29)
        ("inconsistent", 30)
        ("merged", 31)
        ;
    return sym;
}

auto flag_symbols = get_symbols();

enum treemachine_prune_flags
{
	not_otu = 0,
    environmental = 1,
    environmental_inherited = 2,
	viral = 3,
    hidden = 4,
    hidden_inherited = 5,
//    unclassified_direct
    was_container = 6,
    barren = 8,
    extinct = 9,
//    extinct_direct,
    extinct_inherited = 11,
    major_rank_conflict = 12,
// major_rank_conflict_direct
	major_rank_conflict_inherited = 14,
    unclassified = 15,
	unclassified_inherited = 16,
    edited = 17,
    hybrid = 18,
    incertae_sedis = 19,
	incertae_sedis_inherited = 20,
//	incertae_sedis_direct
	infraspecific = 22,
	sibling_lower = 23,
    sibling_higher = 24,
	tattered = 25,
	tattered_inherited = 26,
	forced_visible = 27,
	unplaced = 28,
    unplaced_inherited = 29,
	inconsistent = 30,
    merged = 31
};

auto flag_from_string(const char* start, const char* end)
{
    int n = end - start;
    assert(n >= 0);
    bitset<32> flags;
    if (n > 0)
    {
        int flag = 0;
        boost::spirit::qi::parse(start, end, flag_symbols, flag);
        if (start != end)
        {
            std::cout<<"fail!";
            std::abort();
        }
        flags.set(flag);
    }
    return flags;
}


auto flags_from_string(const char* start, const char* end)
{
    assert(start <= end);

    bitset<32> flags;
    while (start < end)
    {
        assert(start <= end);
        const char* sep = std::strchr(start, ',');
        if (not sep) sep = end;
        flags |= flag_from_string(start, sep);
        start = sep + 1;
    }
    return flags;
}

auto flags_from_string(const string& s)
{
    const char* start = s.c_str();
    const char* end = start + s.length();
    return flags_from_string(start, end);
}

struct taxonomy_record
{
    long id = 0;
    long parent_id = 0;
    int parent_index = 0;
    string name;
    string rank;
    string sourceinfo;
    string uniqname;
    bitset<32> flags;
    bitset<32> marks;
    Tree_t::node_type* node_ptr = nullptr;
    taxonomy_record(taxonomy_record&& tr) = default;
    taxonomy_record(long i1, long i2, string&& s, bitset<32> f): id(i1), parent_id(i2), name(std::move(s)), flags(f) {}
    taxonomy_record(long i1, long i2, const string& s, bitset<32> f): id(i1), parent_id(i2), name(s), flags(f) {}
    explicit taxonomy_record(const string& line);
};

taxonomy_record::taxonomy_record(const string& line)
{
    // parse the line
    const char* start = line.c_str();
    char* end = nullptr;
    id = std::strtoul(line.c_str(),&end,10);
    start = end + 3;
    parent_id = std::strtoul(start,&end,10);
//            taxonomy.push_back(line);
    start = end + 3;
    const char* end3 = std::strstr(start,"\t|\t");
    name = line.substr(start - line.c_str(),end3 - start);
//            cerr<<id<<"\t"<<parent_id<<"\t'"<<name<<"'\n";
    const char* end4 = std::strstr(end3+3,"\t|\t");
    const char* end5 = std::strstr(end4+3,"\t|\t");
    const char* end6 = std::strstr(end5+3,"\t|\t");
    const char* end7 = std::strstr(end6+3,"\t|\t");
    flags = flags_from_string(end6+3,end7);
}

void propagate_marks_to_descendants(vector<taxonomy_record>& taxonomy)
{
    for(auto& record: taxonomy)
        if (record.parent_id)
            record.marks |= taxonomy[record.parent_index].marks;
}

void mark_taxonomy_with_cleaning_flags(vector<taxonomy_record>& taxonomy, bitset<32> cleaning_flags, int bit)
{
    int matched = 0;
    for(auto& record: taxonomy)
        if ((record.flags & cleaning_flags).any())
        {
            matched++;
            record.marks.set(bit);
        }
    std::cerr<<"#lines directly matching cleaning flags = "<<matched<<std::endl;
}

std::unique_ptr<Tree_t> tree_from_taxonomy(vector<taxonomy_record>& taxonomy)
{
    std::unique_ptr<Tree_t> tree(new Tree_t);
    int nodes = 0;
    for(int i=0;i<taxonomy.size();i++)
    {
        const auto& line = taxonomy[i];

        if (line.marks.test(1)) continue;
        if (not line.marks.test(2)) continue;

        // Make the tree
        Tree_t::node_type* nd = nullptr;
        if (not line.parent_id)
            nd = tree->createRoot();
        else if (not taxonomy[line.parent_index].marks.test(2))
            nd = tree->createRoot();
        else
        {
            auto parent_nd = taxonomy[line.parent_index].node_ptr;
            nd = tree->createChild(parent_nd);
        }
        nd->setOttId(line.id);
        nd->setName(line.name);
        taxonomy[i].node_ptr = nd;
        nodes++;
    }
    cerr<<"#tree nodes = "<<nodes<<std::endl;
    return tree;
}

auto cleaning_flags_from_config_file(const string& filename)
{
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    string cleaning_flags_string = pt.get<std::string>("taxonomy.cleaning_flags");
    return flags_from_string(cleaning_flags_string);
}


int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
//        if (argc != 2)
//            throw OTCError()<<"Expecting exactly 1 argument, but got "<<argc-1<<".";
        variables_map args = parse_cmd_line(argc,argv);

        unsigned keep_root = -1;
        if (args.count("root"))
            keep_root = args["root"].as<int>();
        
        bitset<32> cleaning_flags = 0;
        if (args.count("config"))
            cleaning_flags |= cleaning_flags_from_config_file(args["config"].as<string>());
        if (args.count("clean"))
            cleaning_flags |= flags_from_string(args["clean"].as<string>());
        
        string taxonomy_dir = args["taxonomy"].as<string>();
        string filename = taxonomy_dir + "/taxonomy.tsv";

        std::ifstream taxonomy_stream(filename);
        if (not taxonomy_stream)
            throw OTCError()<<"Could not open file '"<<filename<<"'.";
    
        extern std::string foobar;
        string line;
        int count = 0;
        vector<taxonomy_record> taxonomy;
        std::unordered_map<long,int> index;
        std::getline(taxonomy_stream,line);
        if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t")
            throw OTCError()<<"First line of file '"<<filename<<"' is not a taxonomy header.";

        int root = -1;
        while(std::getline(taxonomy_stream,line))
        {
            // Add line to vector
            int my_index = taxonomy.size();
            taxonomy.emplace_back(line);
            auto& record = taxonomy.back();
            index[record.id] = my_index;
            if (record.parent_id)
                record.parent_index = index.at(record.parent_id);
            else
                root = my_index;
            count++;
        }
        cerr<<"#lines = "<<count<<std::endl;

        // Mark records with bit #1 if they match the cleaning flags
        mark_taxonomy_with_cleaning_flags(taxonomy, cleaning_flags, 1);

        // Mark the root
        if (keep_root == -1)
            taxonomy[root].marks.set(2);
        else if (not index.count(keep_root))
            throw OTCError()<<"Can't find root id '"<<keep_root<<"'";
        else
            taxonomy[index.at(keep_root)].marks.set(2);

        // Propagate marks
        propagate_marks_to_descendants(taxonomy);

        if (args.count("write-tree"))
            writeTreeAsNewick(cout, *tree_from_taxonomy(taxonomy));
        std::cout<<std::endl;
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
