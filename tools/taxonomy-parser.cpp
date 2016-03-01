#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
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
using std::unique_ptr;

using boost::spirit::qi::symbols;
using boost::string_ref;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

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

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<long>(), "OTT id of root node of subtree to keep")
        ;

    options_description output("Output options");
    output.add_options()
        ("write-tree,T","Write out the result as a tree")
        ("write-taxonomy",value<string>(),"Write out the result as a taxonomy to directory 'arg'")
        ("name,N", value<long>(), "Return name of the given ID")
        ("uniqname,U", value<long>(), "Return unique name for the given ID")
        ("report-lost-taxa",value<string>(), "A tree to report missing taxa for")
        ("version,V","Taxonomy version")
        ;

    options_description visible;
    visible.add(taxonomy).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: taxonomy-parser <taxonomy-dir> [OPTIONS]\n"
                                                    "Read a taxonomy, clean it, and then make a tree or some other operation.",
                                                    visible, invisible, p);
    return vm;
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

long root_ott_id_from_file(const string& filename)
{
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(filename, pt);
    try {
        return pt.get<long>("synthesis.root_ott_id");
    }
    catch (...)
    {
        return -1;
    }
}

void report_lost_taxa(const Taxonomy& taxonomy, const string& filename)
{
    vector<unique_ptr<Tree_t>> trees;
    std::function<bool(unique_ptr<Tree_t>)> a = [&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    ParsingRules rules;
    rules.requireOttIds = false;
    otc::processTrees(filename,rules,a);//[&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;});
    const auto& T =  trees[0];

    std::unordered_map<long, const Tree_t::node_type*> ottid_to_node;
    for(auto nd: iter_pre_const(*T))
        if (nd->hasOttId())
            ottid_to_node[nd->getOttId()] = nd;

    vector<const taxonomy_record*> records;
    for(const auto& rec: taxonomy)
        records.push_back(&rec);
    
    std::sort(records.begin(), records.end(), [](const auto& a, const auto& b) {return a->depth < b->depth;});
    for(const auto& rec: records)
        if (not ottid_to_node.count(rec->id))
            std::cout<<"depth="<<rec->depth<<"   id="<<rec->id<<"   uniqname='"<<rec->uniqname<<"'\n";
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        if (not args.count("taxonomy"))
            throw OTCError()<<"Please specify the taxonomy directory!";
        
        string taxonomy_dir = args["taxonomy"].as<string>();

        long keep_root = -1;
        if (args.count("root"))
            keep_root = args["root"].as<long>();
        else if (args.count("config"))
            keep_root = root_ott_id_from_file(args["config"].as<string>());
        
        bitset<32> cleaning_flags = 0;
        if (args.count("config"))
            cleaning_flags |= cleaning_flags_from_config_file(args["config"].as<string>());
        if (args.count("clean"))
            cleaning_flags |= flags_from_string(args["clean"].as<string>());

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
        if (args.count("report-lost-taxa"))
        {
            string treefile = args["report-lost-taxa"].as<string>();
            report_lost_taxa(taxonomy,treefile);
        }
        if (args.count("version"))
        {
            std::cout<<taxonomy.version<<std::endl;
        }
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
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

