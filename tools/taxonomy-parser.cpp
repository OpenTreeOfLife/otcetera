#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

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
      ("write-tree,t","Write out the result as a tree")
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

struct taxonomy_record
{
    long id = 0;
    long parent_id = 0;
    int parent_index = 0;
    string name;
    string rank;
    string sourceinfo;
    string uniqname;
    unsigned flags = 0;
    unsigned marks = 0;
    Tree_t::node_type* node_ptr = nullptr;
    taxonomy_record(long i1, long i2, string&& s, unsigned f): id(i1), parent_id(i2), name(std::move(s)), flags(f) {}
    taxonomy_record(long i1, long i2, const string& s, unsigned f): id(i1), parent_id(i2), name(s), flags(f) {}
};

// http://www.boost.org/doc/libs/1_50_0/libs/spirit/doc/html/spirit/qi/reference/string/symbols.html

unsigned flag_from_string(const char* start, const char* end)
{
    int n = end - start;
    assert(n >= 0);
    if (n == 0) return 0;
    int flag = 0;
//    for(auto i = start;i<end;i++)
//        std::cerr<<*i;
//    std::cerr<<" = ";
    boost::spirit::qi::parse(start, end, flag_symbols, flag);
    if (start != end)
    {
        std::cout<<"fail!";
        return 0;
    }
//    std::cerr<<flag<<std::endl;
    return (1<<flag);
}


unsigned flags_from_string(const char* start, const char* end)
{
    assert(start <= end);

    unsigned flags = 0;
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

int flags_from_string(const string& s)
{
    const char* start = s.c_str();
    const char* end = start + s.length();
    return flags_from_string(start, end);
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
        
        unsigned cleaning_flags = 0;
        if (args.count("config"))
        {
            boost::property_tree::ptree pt;
            boost::property_tree::ini_parser::read_ini(args["config"].as<string>(), pt);
            string cleaning_flags_string = pt.get<std::string>("taxonomy.cleaning_flags");
            std::cerr<<cleaning_flags_string<<std::endl;
            cleaning_flags = flags_from_string(cleaning_flags_string);
            std::cerr<<cleaning_flags<<std::endl;
        }
        
        string taxonomy_dir = args["taxonomy"].as<string>();
        string filename = taxonomy_dir + "/taxonomy.tsv";

        std::ifstream taxonomy(filename);
        if (not taxonomy)
            throw OTCError()<<"Could not open file '"<<filename<<"'.";
    
        extern std::string foobar;
        string line;
        int count = 0;
        vector<taxonomy_record> lines;
        std::unordered_map<long,Tree_t::node_type*> id_to_node;
        std::unordered_map<long,int> index;
        std::getline(taxonomy,line);
        if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t")
            throw OTCError()<<"First line of file '"<<filename<<"' is not a taxonomy header.";

        std::unique_ptr<Tree_t> tree(new Tree_t());
        int matched = 0;
        int nodes = 0;
        while(std::getline(taxonomy,line))
        {
            // parse the line
            int b1 = line.find("\t|\t",0);
            const char* start = line.c_str();
            char* end = nullptr;
            long id = std::strtoul(line.c_str(),&end,10);
            start = end + 3;
            long parent_id = std::strtoul(start,&end,10);
//            lines.push_back(line);
            start = end + 3;
            const char* end3 = std::strstr(start,"\t|\t");
            string name = line.substr(start - line.c_str(),end3 - start);
//            cerr<<id<<"\t"<<parent_id<<"\t'"<<name<<"'\n";
            const char* end4 = std::strstr(end3+3,"\t|\t");
            const char* end5 = std::strstr(end4+3,"\t|\t");
            const char* end6 = std::strstr(end5+3,"\t|\t");
            const char* end7 = std::strstr(end6+3,"\t|\t");
            unsigned flags = flags_from_string(end6+3,end7);

            // Add line to vector
            int my_index = lines.size();
            index[id] = my_index;
            lines.push_back({id,parent_id,name,flags});
            if (parent_id)
            {
                int parent_index = index.at(parent_id);
                lines.back().parent_index = parent_index;
                lines.back().marks |= lines[parent_index].marks;
            }
            count++;

            // Compare with cleaning flags
            if ((lines.back().flags & cleaning_flags) != 0)
            {
                matched++;
                lines.back().marks |= 1;
            }
            if (lines.back().id == keep_root or (keep_root == -1 and not parent_id))
                lines.back().marks |= 2;

        }
        
        cerr<<"#lines = "<<count<<std::endl;
        cerr<<"#matched lines = "<<matched<<std::endl;

        if (args.count("write-tree"))
        {
            for(int i=0;i<lines.size();i++)
            {
                const auto& line = lines[i];

                if ((line.marks & 1) != 0) continue;
                if ((line.marks & 2) == 0) continue;

                nodes++;
                // Make the tree
                Tree_t::node_type* nd = nullptr;
                if (not line.parent_id)
                    nd = tree->createRoot();
                else if ((lines[line.parent_index].marks & 2) == 0)
                    nd = tree->createRoot();
                else
                {
                    auto parent_nd = lines[line.parent_index].node_ptr;
                    nd = tree->createChild(parent_nd);
                }
                nd->setOttId(line.id);
                nd->setName(line.name);
                lines[i].node_ptr = nd;
            }
            writeTreeAsNewick(cout, *tree);
        }
        cerr<<"#tree nodes = "<<nodes<<std::endl;
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
