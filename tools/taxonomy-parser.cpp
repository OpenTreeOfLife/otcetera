#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <bitset>
#include <regex>

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
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

using boost::spirit::qi::symbols;
using boost::string_ref;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;
using namespace boost::property_tree;

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
	("any-flags",value<string>(),"Show nodes with one of these flags")
	("all-flags",value<string>(),"Show nodes with all of these flags")
        ("show-root,R","Show the ottid of the root node")
        ("find,S",value<string>(),"Show taxa whose names match regex <arg>")
        ("degree,D",value<long>(),"Show out the degree of node <arg>")
        ("children,C",value<long>(),"Show the children of node <arg>")
        ("parent,P",value<long>(),"Show the parent taxon of node <arg>")
        ("high-degree-nodes",value<int>(),"Show the top <arg> high-degree nodes")
        ("write-tree,T","Write out the result as a tree")
        ("write-taxonomy",value<string>(),"Write out the result as a taxonomy to directory 'arg'")
        ("name,N", value<long>(), "Return name of the given ID")
        ("uniqname,U", value<long>(), "Return unique name for the given ID")
        ("report-lost-taxa",value<string>(), "A tree to report missing taxa for")
        ("version,V","Taxonomy version")
        ;

    options_description formatting("Formatting options");
    formatting.add_options()
	("format",value<string>()->default_value("%I"),"Form of line to write for each taxonomy record");

    options_description visible;
    visible.add(taxonomy).add(output).add(formatting).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("taxonomy", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-taxonomy-parser <taxonomy-dir> [OPTIONS]\n"
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

void show_rec(const taxonomy_record& rec)
{
    std::cout<<rec.id<<"   '"<<rec.uniqname<<"'   '"<<rec.rank<<"'   depth = "<<rec.depth<<"   out-degree = "<<rec.out_degree<<"    flags = "<<flags_to_string(rec.flags)<<"\n";
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        auto taxonomy = load_taxonomy(args);

        if (args.count("show-root"))
        {
            show_rec(taxonomy[0]);
            exit(0);
        }
	else if (args.count("any-flags") or args.count("all-flags"))
	{
	    // Currently this option is mostly used to print out the incertae sedis ott numbers.

	    std::bitset<32> all_flags;
	    if (args.count("all-flags"))
		all_flags = flags_from_string(args["all-flags"].as<string>());

	    std::bitset<32> any_flags = all_flags;
	    if (args.count("any-flags"))
		any_flags = flags_from_string(args["any-flags"].as<string>());

	    bool id_only = args.count("id-only");

	    string format=args["format"].as<string>();

            for(const auto& rec: taxonomy)
		if ((rec.flags&any_flags).any() and (all_flags == (rec.flags&all_flags)))
		    std::cout<<format_with_taxonomy("No original label",format,rec)<<"\n";
	}
        if (args.count("find"))
        {
            string s = args["find"].as<string>();
            std::regex e(s);
            for(const auto& rec: taxonomy)
            {
                std::cmatch m;
                if (std::regex_match(rec.name.data(), rec.name.data()+rec.name.size(), m, e))
                    show_rec(rec);
            }
            exit(0);
        }
        else if (args.count("degree"))
        {
            long id = args["degree"].as<long>();
            std::cout<<"degree = "<<taxonomy[taxonomy.index.at(id)].out_degree<<std::endl;
            exit(0);
        }
        else if (args.count("children"))
        {
            long id = args["children"].as<long>();

            for(const auto& rec: taxonomy)
            {
                long parent_id = taxonomy[rec.parent_index].id;
                if (parent_id == id)
                    show_rec(rec);
            }
            exit(0);
        }
        else if (args.count("parent"))
        {
            long id = args["parent"].as<long>();
            auto parent_index = taxonomy[taxonomy.index.at(id)].parent_index;
            show_rec(taxonomy[parent_index]);
            exit(0);
        }
        else if (args.count("high-degree-nodes"))
        {
            int n = args["high-degree-nodes"].as<int>();
            vector<int> index(taxonomy.size());
            for(int i=0;i<index.size();i++)
                index[i] = i;
            std::sort(index.begin(), index.end(), [&taxonomy](int i, int j){return taxonomy[i].out_degree > taxonomy[j].out_degree;});
            for(int i=0;i<n;i++)
                show_rec(taxonomy[index[i]]);
            exit(0);
        }
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

