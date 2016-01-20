#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/tree_operations.h"

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

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

struct taxonomy_record
{
    int id = 0;
    int parent_id = 0;
    string name;
    taxonomy_record(int i1, int i2, string&& s): id(i1), parent_id(i2), name(std::move(s)) {}
    taxonomy_record(int i1, int i2, const string& s): id(i1), parent_id(i2), name(s) {}
};

int main(int argc, char* argv[])
{
    try
    {
//        if (argc != 2)
//            throw OTCError()<<"Expecting exactly 1 argument, but got "<<argc-1<<".";
        variables_map args = parse_cmd_line(argc,argv);

        if (args.count("config"))
        {
            
            boost::property_tree::ptree pt;
            boost::property_tree::ini_parser::read_ini(args["config"].as<string>(), pt);
            string cleaning_flags = pt.get<std::string>("taxonomy.cleaning_flags");
            std::cout<<cleaning_flags<<std::endl;
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
        std::getline(taxonomy,line);
        if (line != "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t")
            throw OTCError()<<"First line of file '"<<filename<<"' is not a taxonomy header.";

        std::unique_ptr<Tree_t> tree(new Tree_t());
        while(std::getline(taxonomy,line))
        {
            int b1 = line.find("\t|\t",0);
            const char* start = line.c_str();
            char* end = nullptr;
            int id = std::strtoul(line.c_str(),&end,10);
            start = end + 3;
            int parent_id = std::strtoul(start,&end,10);
//            lines.push_back(line);
            start = end + 3;
            const char* end3 = std::strstr(start,"\t|\t");
            string name = line.substr(start - line.c_str(),end3 - start);
//            cerr<<id<<"\t"<<parent_id<<"\t'"<<name<<"'\n";

            Tree_t::node_type* nd = nullptr;
            if (not parent_id)
            {
                nd = tree->createRoot();
            }
            else
            {
                auto parent_nd = id_to_node.at(parent_id);
                nd = tree->createChild(parent_nd);
            }
            nd->setOttId(id);
            nd->setName(name);
            id_to_node[id] = nd;
            
            lines.push_back({id,parent_id,name});
            count++;
        }
        cerr<<"#lines = "<<count<<std::endl;
        if (args.count("write-tree"))
            writeTreeAsNewick(cout, *tree);
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
