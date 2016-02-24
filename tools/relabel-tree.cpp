#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/tokenizer.hpp>
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
        ("tree", value<string>(),"Filename for the taxonomy")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<long>(), "OTT id of root node of subtree to keep")
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;

    options_description output("Relabelling options");
    output.add_options()
        ("format-tax",value<string>()->default_value("%L"),"Form of labels to write for taxonomy nodes.")
        ("format-unknown",value<string>()->default_value("%L"),"Form of labels to write for non-taxonomy nodes.")
//        ("label-regex", value<string>(), "Return name of the given ID")
        ;

    options_description tree("Tree options");
    tree.add_options()
        ("del-monotypic","Remove monotypic nodes.")
        ;

    options_description visible;
    visible.add(taxonomy).add(output).add(tree).add(otc::standard_options());

    // positional options
    positional_options_description p;
    p.add("tree", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-relabel-tree <newick-file> [OPTIONS]\n"
                                                    "Rewrite node labels for a tree.\n\n"
                                                    "Format strings have the following interpretation:\n"
                                                    "  %I=id %N=name %U=uniqname %R=rank %S=sourceinfo %%=% %L=original label",
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

unique_ptr<Tree_t> get_tree(const string& filename)
{
    vector<unique_ptr<Tree_t>> trees;
    std::function<bool(unique_ptr<Tree_t>)> a = [&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    ParsingRules rules;
    rules.requireOttIds = false;
    otc::processTrees(filename,rules,a);//[&](unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;});
    return std::move(trees[0]);
}

string format_with_taxonomy(const string& orig, const string& format, const taxonomy_record& rec)
{
    string result;
    
    int pos = 0;
    do {
        auto loc = format.find('%',pos);
        if (loc == string::npos)
        {
            result += format.substr(pos);
            break;
        }
        result += format.substr(pos,loc-pos);
        loc++;
        if (format[loc] == 0)
            std::abort();
        else if (format[loc] == 'I')
            result += std::to_string(rec.id);
        else if (format[loc] == 'N')
            result += rec.name.to_string();
        else if (format[loc] == 'U')
            result += rec.uniqname.to_string();
        else if (format[loc] == 'R')
            result += rec.rank.to_string();
        else if (format[loc] == 'S')
            result += rec.sourceinfo.to_string();
        else if (format[loc] == 'L')
            result += orig;
        else if (format[loc] == '%')
            result += '%';
        else
            throw OTCError()<<"Invalid format specification %"<<format[loc]<<" in taxonomy-based format string '"<<format<<"'";
        pos = loc + 1;
    }
    while (pos < format.size());

    return result;
}

string format_without_taxonomy(const string& orig, const string& format)
{
    string result;
    
    int pos = 0;
    do {
        auto loc = format.find('%',pos);
        if (loc == string::npos)
        {
            result += format.substr(pos);
            break;
        }
        result += format.substr(pos,loc-pos);
        loc++;
        if (format[loc] == 0)
            std::abort();
        else if (format[loc] == 'L')
            result += orig;
        else if (format[loc] == '%')
            result += '%';
        else
            throw OTCError()<<"Invalid format specification %"<<format[loc]<<" in non-taxonomy-based format string '"<<format<<"'";
        pos = loc + 1;
    }
    while (pos < format.size());

    return result;
}

void suppressMonotypicFast(Tree_t& tree)
{
    std::vector<Tree_t::node_type*> remove;
    for(auto nd:iter_pre(tree))
        if (nd->isOutDegreeOneNode())
            remove.push_back(nd);

    for(auto nd: remove)
    {
        auto child = nd->getFirstChild();
        child->detachThisNode();
        nd->addSibOnRight(child);
        nd->detachThisNode();
        delete nd;
    }
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        if (not args.count("tree"))
            throw OTCError()<<"Please specify the newick tree to be relabelled!";
        
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

        bool keep_non_ott = not args.count("unlabel-non-taxa");

        auto tree = get_tree(args["tree"].as<string>());

        Taxonomy taxonomy(taxonomy_dir, cleaning_flags, keep_root);

        string format_tax = args["format-tax"].as<string>();
        string format_unknown = args["format-unknown"].as<string>();

        if (args.count("del-monotypic"))
            suppressMonotypicFast(*tree);
        
        for(auto nd: iter_pre(*tree))
        {
            if (nd->hasOttId())
            {
                int id = nd->getOttId();
                const auto& record = taxonomy.record_from_id(id);
                string name = format_with_taxonomy(nd->getName(), format_tax, record);
                nd->setName(std::move(name));
            }
            else
            {
                string name = format_without_taxonomy(nd->getName(), format_unknown);
                nd->setName(std::move(name));
            }

        }

        writeTreeAsNewick(cout, *tree);
    }
    catch (std::exception& e)
    {
        cerr<<"otc-taxonomy-parser: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
