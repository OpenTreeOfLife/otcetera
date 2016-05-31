#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <regex>
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
        ("replace,R",value<string>(),"Perform a regex replacement on all labels")
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

char format_needs_taxonomy(const string& format)
{
    int pos = 0;
    do {
        auto loc = format.find('%',pos);
        if (loc == string::npos)
            break;
        loc++;
        if (format[loc] == 0)
            std::abort();
        else if (format[loc] == 'I')
            return 'I';
        else if (format[loc] == 'N')
            return 'N';
        else if (format[loc] == 'U')
            return 'U';
        else if (format[loc] == 'R')
            return 'R';
        else if (format[loc] == 'S')
            return 'S';
        else if (format[loc] == 'L')
            ;
        else if (format[loc] == '%')
            ;
        else
            throw OTCError()<<"Invalid format specification '%"<<format[loc]<<"' in format string '"<<format<<"'";
        pos = loc + 1;
    }
    while (pos < static_cast<int>(format.size()));

    return false;
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
            throw OTCError()<<"Invalid format specification '%"<<format[loc]<<"' in taxonomy-based format string '"<<format<<"'";
        pos = loc + 1;
    }
    while (pos < static_cast<int>(format.size()));

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
            throw OTCError()<<"Invalid format specification '%"<<format[loc]<<"' in non-taxonomy-based format string '"<<format<<"'";
        pos = loc + 1;
    }
    while (pos < static_cast<int>(format.size()));

    return result;
}

int main(int argc, char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        variables_map args = parse_cmd_line(argc,argv);

        if (not args.count("tree"))
            throw OTCError()<<"Please specify the newick tree to be relabelled!";
        
        string format_tax = args["format-tax"].as<string>();
        string format_unknown = args["format-unknown"].as<string>();

        boost::optional<Taxonomy> taxonomy = boost::none;
        
        if (format_needs_taxonomy(format_tax))
            taxonomy = load_taxonomy(args);

        if (char c = format_needs_taxonomy(format_unknown))
            throw OTCError()<<"Cannot use taxonomy-based specifier '%"<<c<<"' in non-taxonomy format string '"<<format_unknown<<"'";
        
        auto tree = get_tree<Tree_t>(args["tree"].as<string>());

        if (args.count("del-monotypic"))
            suppressMonotypicFast(*tree);

        for(auto nd: iter_pre(*tree))
        {
            if (nd->hasOttId())
            {
                int id = nd->getOttId();
                if (taxonomy)
                {
                    const auto& record = (*taxonomy).record_from_id(id);
                    nd->setName( format_with_taxonomy(nd->getName(), format_tax, record) );
                }
                else
                    nd->setName( format_without_taxonomy(nd->getName(), format_tax ) );
            }
            else
            {
                string name = format_without_taxonomy(nd->getName(), format_unknown);
                nd->setName(std::move(name));
            }
        }

        if (args.count("replace"))
        {
            string match_replace = args["replace"].as<string>();
            if (match_replace.empty())
                throw OTCError()<<"Empty pattern for argument 'replace'!";
            char sep = match_replace[0];

            if (std::count(match_replace.begin(), match_replace.end(), sep) !=3)
                throw OTCError()<<"Delimiter '"<<sep<<"' does not occur exactly three times in '"<<match_replace<<"'!";
            if (match_replace.back() != sep)
                throw OTCError()<<"Pattern '"<<match_replace<<"' does not end with delimiter '"<<sep<<"'";

            int loc2 = match_replace.find(sep,1);
            int loc3 = match_replace.find(sep,loc2+1);
            
            std::regex match (match_replace.substr(1,loc2-1));
            string replace = match_replace.substr(loc2+1,loc3-loc2-1);

            for(auto nd: iter_pre(*tree))
            {
                string name = nd->getName();
                name = std::regex_replace(name,match,replace);
                nd->setName(name);
            }
        }

        writeTreeAsNewick(cout, *tree);
        std::cout<<std::endl;
    }
    catch (std::exception& e)
    {
        cerr<<"otc-relabel-tree: Error! "<<e.what()<<std::endl;
        exit(1);
    }
}
