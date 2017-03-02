#include "otc/otcli.h"
#include <cstring>
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;
template<typename T>
bool cleanNamesAndWrite(OTCLI & otCLI, std::unique_ptr<T> tree);

bool hasProblematicCharacter(const std::string & name);
std::string generateCleanedName(const std::string & name);


static const char * problematicChars = ",():;[]{}'\t\n";
static const char replacementChar = ' ';
static bool logToStderr = true;

bool hasProblematicCharacter(const std::string & name) {
    for (const auto & c : name) {
        if (strchr(problematicChars, c) != nullptr) {
            return true;
        }
    }
    return false;
}

std::string generateCleanedName(const std::string & name) {
    std::string x = name;
    std::size_t pos = 0U;
    for (const auto & c : name) {
        if (strchr(problematicChars, c) != nullptr) {
            x[pos] = replacementChar;
        }
        ++pos;
    }
    if (logToStderr) {
        std::cerr << "Replacing \"" << name << "\" with \"" << x << "\"\n";
    }
    return x;
}


template<typename T>
bool cleanNamesAndWrite(OTCLI & otCLI, std::unique_ptr<T> tree) {
    for (auto nd : iter_node(*tree)) {
        const auto name = nd->get_name();
        if (hasProblematicCharacter(name)) {
            nd->setName(generateCleanedName(name));
        }
    }
    writeTreeAsNewick(std::cout, *tree);
    return true;
}

int main(int argc, char *argv[]) {
    OTCLI otCLI("otc-munge-names",
                 "takes a filepath to a newick file. Writes the tree to standard output after replacing each " \
                 "\ninstance of a problematic character in a name with space." \
                 "\nThe problematic characters are considered to be:\n" \
                 "   , ( ) : ; [ ] { } ' <TAB> <NEWLINE>",
                 {"some.tre"});
    std::function<bool (OTCLI &, std::unique_ptr<Tree_t>)> wnl = cleanNamesAndWrite<Tree_t>;
    return tree_processing_main<Tree_t>(otCLI, argc, argv, wnl, nullptr, 1);
}

