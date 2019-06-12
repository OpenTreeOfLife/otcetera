#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <boost/program_options.hpp>
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

variables_map parse_cmd_line(int argc,char* argv[]) {
    using namespace po;
    options_description invisible("Invisible options");
    invisible.add_options()("taxonomy", value<string>(), "Filename for the taxonomy");
    options_description visible;
    visible.add(otc::standard_options());
    positional_options_description p;
    p.add("taxonomy", -1);
    variables_map vm = otc::parse_cmd_line_standard(argc, argv,
                                                    "Usage: otc-tnrs-cli <taxonomy-dir> [OPTIONS]\n"
                                                    "Build data structures for name matching, and allow interactive testing.",
                                                    visible, invisible, p);

    return vm;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    try {
        auto args = parse_cmd_line(argc,argv);
        auto taxonomy = load_taxonomy(args);
    } catch (std::exception& e) {
        cerr << "otc-tnrs-cli: Error! " << e.what() << std::endl;
        return 1;
    }
}
