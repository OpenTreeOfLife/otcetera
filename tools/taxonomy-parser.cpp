// TODO: mmap via BOOST https://techoverflow.net/blog/2013/03/31/mmap-with-boost-iostreams-a-minimalist-example/
// TODO: write out a reduced taxonomy

#include <iostream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/utility/string_ref.hpp>
#include <bitset>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"

INITIALIZE_EASYLOGGINGPP

using namespace otc;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;

using boost::spirit::qi::symbols;
using boost::string_ref;
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
      ("write-taxonomy",value<string>(),"Write out the result as a taxonomy to directory 'arg'")
      ("root,r", value<long>(), "OTT id of root node of subtree to keep")
      ("name,N", value<long>(), "Return name of the given ID")
      ("uniqname,U", value<long>(), "Return unique name for the given ID")
      ("report-lost-taxa",value<string>(), "A tree to report missing taxa for")
      ("version,V","Taxonomy version")
      ("quiet,q","QUIET mode (all logging disabled)")
      ("trace,t","TRACE level debugging (very noisy)")
      ("verbose,v","verbose")
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

  el::Configurations defaultConf;
  if (args.count("quiet"))
      defaultConf.set(el::Level::Global,  el::ConfigurationType::Enabled, "false");
  else
  {
      defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");
      defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");

      if (args.count("trace"))
          defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "true");

      if (args.count("verbose"))
          defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
  }
  el::Loggers::reconfigureLogger("default", defaultConf);

  return args;
}

long n_nodes(const Tree_t& T) {
#pragma clang diagnostic ignored  "-Wunused-variable"
#pragma GCC diagnostic ignored  "-Wunused-variable"
    long count = 0;
    for(auto nd: iter_post_const(T)){
        count++;
    }
    return count;
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        if (not args.count("taxonomy"))
            throw OTCError()<<"Please specify the taxonomy directory!";
        
        long keep_root = -1;
        if (args.count("root"))
            keep_root = args["root"].as<long>();
        
        bitset<32> cleaning_flags = 0;
        if (args.count("config"))
            cleaning_flags |= cleaning_flags_from_config_file(args["config"].as<string>());
        if (args.count("clean"))
            cleaning_flags |= flags_from_string(args["clean"].as<string>());

        string taxonomy_dir = args["taxonomy"].as<string>();

        Taxonomy taxonomy(taxonomy_dir, cleaning_flags, keep_root);

        if (args.count("write-tree"))
        {
            auto nodeNamer = [](const auto& record){return string(record.name)+"_ott"+std::to_string(record.id);};
            writeTreeAsNewick(cout, *taxonomy.getTree<Tree_t>(nodeNamer));
            std::cout<<std::endl;
        }
        if (args.count("write-taxonomy"))
            taxonomy.write(args["write-taxonomy"].as<string>());
        if (args.count("name"))
        {
            long id = args["name"].as<long>();
            std::cout<<taxonomy[taxonomy.index.at(id)].name<<std::endl;
        }
        if (args.count("uniqname"))
        {
            long id = args["uniqname"].as<long>();
            std::cout<<taxonomy[taxonomy.index.at(id)].uniqname<<std::endl;
        }
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
