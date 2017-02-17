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

    options_description selection("Selection options");
    selection.add_options()
	("any-flags",value<string>(),"Show records with one of these flags")
	("all-flags",value<string>(),"Show records with all of these flags")
	("in-tree",value<string>(),"Show records from OTT ids in tree")
	("in-file",value<string>(),"Show records of integer ids in file");

    options_description output("Output options");
    output.add_options()
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
	("format",value<string>()->default_value("ott%I\t'%U'\tflags=%F"),"Form of line to write for each taxonomy record");

    options_description visible;
    visible.add(taxonomy).add(selection).add(output).add(formatting).add(otc::standard_options());

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

vector<long> get_ids_from_file(const string& filename)
{
    vector<long> ids;
    std::ifstream file(filename);
    while (file)
    {
	long i;
	file >> i;
	ids.push_back(i);
    }
    return ids;
}

vector<long> get_ids_from_tree(const string& filename)
{
    vector<long> ids;
    auto tree = get_tree<Tree_t>(filename);
    for(auto nd: iter_post_const(*tree))
    {
	auto id = nd->getOttId();
	if (id != -1)
	    ids.push_back(id);
    }
    return ids;
}

vector<long> get_ids_matching_regex(const Taxonomy& taxonomy, const string& rgx)
{
    std::regex e(rgx);
    vector<long> ids;
    for(const auto& rec: taxonomy)
    {
	std::cmatch m;
	if (std::regex_match(rec.name.data(), rec.name.data()+rec.name.size(), m, e))
	    ids.push_back(rec.id);
    }
    return ids;
}

bool has_flags(tax_flags flags, tax_flags any_flags, tax_flags all_flags)
{
    if ((flags&all_flags) != all_flags) return false;
    if (any_flags.any() and (flags&any_flags).none()) return false;
    return true;
}

void show_taxonomy_ids(const Taxonomy& taxonomy, const string& format, const vector<long>& ids,
		       std::function<bool(tax_flags)> flags_match)
{
    for(auto id: ids)
    {
	try
	{
	    auto& rec = taxonomy.record_from_id(id);
	    if (flags_match(rec.flags))
		std::cout<<format_with_taxonomy("No original label",format,rec)<<"\n";
	}
	catch (...) {
	    std::cerr<<"id="<<id<<": not in taxonomy!\n";
	}
    }
}

std::function<bool(tax_flags)> get_flags_match(variables_map& args)
{
    tax_flags all_flags;
    if (args.count("all-flags"))
	all_flags = flags_from_string(args["all-flags"].as<string>());

    tax_flags any_flags;
    if (args.count("any-flags"))
	any_flags = flags_from_string(args["any-flags"].as<string>());

    
    if (any_flags.any() and all_flags.any())
	return [all_flags,any_flags](tax_flags flags) {return (flags&any_flags).any() and
		                                              (flags&all_flags)==all_flags; };
    else if (any_flags.any())
	return [any_flags](tax_flags flags) { return (flags&any_flags).any(); };
    else if (all_flags.any())
	return [all_flags](tax_flags flags) { return (flags&all_flags)==all_flags; };
    else
	return [](tax_flags){return true;};
}


int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        auto args = parse_cmd_line(argc,argv);

	auto format = args["format"].as<string>();

	auto flags_match = get_flags_match(args);

	auto taxonomy = load_taxonomy(args);

        if (args.count("show-root"))
        {
            show_rec(taxonomy[0]);
            exit(0);
        }
	else if (args.count("in-tree"))
	{
	    auto ids = get_ids_from_tree(args["in-tree"].as<string>());

	    show_taxonomy_ids(taxonomy, format, ids, flags_match);
	    exit(0);
	}
	else if (args.count("in-file"))
	{
	    auto ids = get_ids_from_file(args["in-file"].as<string>());

	    show_taxonomy_ids(taxonomy, format, ids, flags_match);
	    exit(0);
	}
        else if (args.count("find"))
        {
	    vector<long> ids = get_ids_matching_regex(taxonomy, args["find"].as<string>());

	    show_taxonomy_ids(taxonomy, format, ids, flags_match);
            exit(0);
        }
	else if (args.count("any-flags") or args.count("all-flags"))
	{
	    string format=args["format"].as<string>();

            for(const auto& rec: taxonomy)
		if (flags_match(rec.flags))
		    std::cout<<format_with_taxonomy("No original label",format,rec)<<"\n";
	    exit(0);
	}

        if (args.count("degree"))
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

