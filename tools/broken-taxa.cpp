#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <iomanip>
#include <regex>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <bitset>
#include <set>
#include <map>
#include "json.hpp"

#include "otc/error.h"
#include "otc/tree.h"
#include "otc/otcli.h"
#include "otc/tree_operations.h"
#include "otc/taxonomy/taxonomy.h"
#include "otc/taxonomy/flags.h"
INITIALIZE_EASYLOGGINGPP

using namespace otc;
namespace fs = boost::filesystem;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::unique_ptr;
using std::set;
using std::map;
using json = nlohmann::json;

using boost::spirit::qi::symbols;
using namespace boost::spirit;

using Tree_t = RootedTree<RTNodeNoData, RTreeNoData>;

namespace po = boost::program_options;
using po::variables_map;

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;

    // named options
    options_description invisible("Invisible options");
    invisible.add_options()
        ("broken_taxa_json", value<string>(),"broken_taxa_json")
        ("tree", value<string>(),"Filename for synthesis tree")
        ;

    options_description taxonomy("Taxonomy options");
    taxonomy.add_options()
        ("config,c",value<string>(),"Config file containing flags to filter")
        ("clean",value<string>(),"Comma-separated string of flags to filter")
        ("root,r", value<OttId>(), "OTT id of root node of subtree to keep")
        ("taxonomy", value<string>(),"Directory name for the taxonomy")
        ;

    options_description output("Output options");
    output.add_options()
        ("json,j", value<string>(), "filepath to an output JSON log")
        ;

    options_description visible;
    visible.add(taxonomy).add(output).add(otc::standard_options());

    // positional options
    positional_options_description p;
//    p.add("tree", -1);
    p.add("broken_taxa_json", -1);

    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-broken-taxa <broken-taxa-json> [OPTIONS]\n"
                                                    "Manipulate JSON with info on broken taxa.\n",
                                                    visible, invisible, p);
    return vm;
}

OttId get_ott_id_from_string(const string& s) {
    if (s.size() <= 3 or s[0] != 'o' or s[1] != 't' or s[2] != 't') {
        throw OTCError() << "String '" << s << "' does not begin with 'ott'";
    }
    return check_ott_id_size(std::stol(s.substr(3)));
}

void add_name_and_rank(json& broken_taxa, const Taxonomy& taxonomy)
{
    for(auto& [ott_id_string, tax]: broken_taxa["non_monophyletic_taxa"].items())
    {
        auto ott_id = get_ott_id_from_string( ott_id_string );
        auto& record = taxonomy.record_from_id(ott_id);
        tax["name"] = record.name;
        tax["rank"] = record.rank;
        tax["depth"] = record.depth;
    }
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        variables_map args = parse_cmd_line(argc,argv);

        json broken_taxa;
        if (args.count("broken_taxa_json")) {
            string json_filename = args["broken_taxa_json"].as<string>();
            std::ifstream json_stream( json_filename );
            if (!json_stream.good()) {
                throw OTCError() << "Could not open JSON broken taxa file at \"" <<fs::absolute(json_filename) <<"\"";
            }
            json_stream >> broken_taxa;
        } else {
            throw OTCError() << "Broken taxa JSON file not specified!";
        }
        auto taxonomy = load_taxonomy(args);
        add_name_and_rank(broken_taxa, taxonomy);
        std::cout << std::setw(1) << broken_taxa << std::endl;
    }
    catch (std::exception& e) {
        cerr << "otc-broken-taxa: Error! " << e.what() << std::endl;
        exit(1);
    }
}
