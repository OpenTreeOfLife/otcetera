#if !defined OTCLI_H
#define OTCLI_H
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "otc/otc_base_includes.h"
#include "otc/newick.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <memory>
namespace otc {
bool get_bool(const std::string& arg, const std::string& context);

class OTCLI {
    public:
        OTCLI(const char *title,
              const char *descrip,
              const char *usage,
              bool silentExecution=false);
        int exitCode;
        bool verbose;
        bool currReadingDotTxtFile;
        std::string currentFilename;
        std::string prefixForFiles;
        std::string currTmpFilepath;
        void * blob;

        void addFlag(char flag, const std::string & help, bool (*cb)(OTCLI &, const std::string &), bool argNeeded) {
            if (clientDefFlagHelp.find(flag) != clientDefFlagHelp.end()) {
                throw OTCError("Clashing flag assignment");
            }
            clientDefFlagHelp[flag] = help;
            clientDefFlagCallbacks[flag] = cb;
            if (argNeeded) {
                clientDefArgNeeded.insert(flag);
            }
        }
        bool handleFlag(const std::string & flagWithoutDash);
        bool parseArgs(int argc, char *argv[], std::vector<std::string> & args);
        void printHelp(std::ostream & outStream);
        
        /*int readFilepath(const std::string &fp,
                          ProcessedTreeValidationFunction func=0L,
                          void * blob=0L);
        // reads a NexSON v1.2 at filepath and returns a NxsSimpleTree with the associated treeID
        NxsSimpleTree * readTreeFromNexSONv_1_2(const std::string &filepath, const std::string & tree_id); */
        bool isDotTxtFile(const std::string &fp);
        auto getTitle() const {
            return this->titleStr;
        }
        ParsingRules & getParsingRules() {
            return parsingRules;
        }
        const ParsingRules & getParsingRules() const {
            return parsingRules;
        }
        void turnOnVerboseMode();
        void turnOffVerboseMode();
    private:
        std::string titleStr;
        std::string descriptionStr;
        std::string usageStr;
        std::map<char, std::string> clientDefFlagHelp;
        std::map<char, bool (*)(OTCLI &, const std::string &)> clientDefFlagCallbacks;
        std::set<char> clientDefArgNeeded;
        ParsingRules parsingRules;
        std::list<std::string> extraArgs;
    public:
        std::ostream & out;
        std::ostream & err;
};

template<typename T>
inline bool processTrees(const std::string& filename,
                         const ParsingRules& parsingRules,
                         std::function<bool (std::unique_ptr<T>)> treePtr)
{
    std::ifstream inp;
    if (!openUTF8File(filename, inp)) {
        throw OTCError("Could not open \"" + filename + "\"");
    }
    LOG(INFO) << "reading \"" << filename << "\"...";
    
    ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
    FilePosStruct pos(filenamePtr);
    unsigned treeNumInThisFile = 1;
    for (;;) {
        std::unique_ptr<T> nt = readNextNewick<T>(inp, pos, parsingRules);
        if (nt == nullptr) {
            break;
        }
        if (parsingRules.pruneUnrecognizedInputTips) {
            pruneTipsWithoutIds(*nt);
            if (nt->getRoot() == nullptr) {
                continue;
            }
        }
        std::string treeName = std::string("tree ") + std::to_string(treeNumInThisFile++);
        treeName.append(" from ");
        treeName.append(filepathToFilename(filename));
        nt->setName(treeName);
        const auto cbr = treePtr(std::move(nt));
        if (not cbr) {
            return false;
        }
    }
    return true;
}

template<typename Tree_t>
std::vector<std::unique_ptr<Tree_t>> get_trees(const std::vector<std::string>& filenames, const ParsingRules& rules)
{
    std::vector<std::unique_ptr<Tree_t>> trees;
    std::function<bool(std::unique_ptr<Tree_t>)> proc = [&](std::unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    for(auto& filename: filenames)
	otc::processTrees(filename,rules,proc);
    return trees;
}

template<typename Tree_t>
std::vector<std::unique_ptr<Tree_t>> get_trees(const std::string& filename, const ParsingRules& rules)
{
    const std::vector<std::string> filenameVec = {filename};
    return get_trees<Tree_t>(filenameVec, rules);
}

template<typename Tree_t>
std::unique_ptr<Tree_t> get_tree(const std::string& filename, const ParsingRules& rules)
{
    auto trees = get_trees<Tree_t>(filename, rules);
    return std::move(trees[0]);
}

template<typename Tree_t>
std::unique_ptr<Tree_t> get_tree(const std::string& filename)
{
    ParsingRules rules;
    rules.requireOttIds = false;
    return get_tree<Tree_t>(filename, rules);
}

template<typename T>
int treeProcessingMain(OTCLI & otCLI,
                          int argc,
                          char * argv[],
                          std::function<bool (OTCLI &, std::unique_ptr<T>)> treePtr,
                          int (*summarizePtr)(OTCLI &),
                          unsigned minNumTrees);

template<typename T>
inline int treeProcessingMain(OTCLI & otCLI,
                                 int argc,
                                 char * argv[],
                                 std::function<bool (OTCLI &, std::unique_ptr<T>)> treePtr,
                                 int (*summarizePtr)(OTCLI &),
                                 unsigned minNumTrees) {
    std::vector<std::string> filenameVec;
    if (!otCLI.parseArgs(argc, argv, filenameVec)) {
        otCLI.exitCode = 1;
        return otCLI.exitCode;
    }
    if (filenameVec.size() < minNumTrees) {
        otCLI.printHelp(otCLI.err);
        otCLI.err << otCLI.getTitle() << ": Expecting at least " << minNumTrees << " tree filepath(s).\n";
        otCLI.exitCode = 1;
        return otCLI.exitCode;
    }
    try {
        if (treePtr) {
            std::function<bool(std::unique_ptr<T>)> treePtr1 = 
                [&otCLI,&treePtr](std::unique_ptr<T> t)
                {return treePtr(otCLI,std::move(t));};
            
            for (const auto & filename : filenameVec) {
                otCLI.currentFilename = filepathToFilename(filename);
                const auto cbr = processTrees(filename, otCLI.getParsingRules(), treePtr1);
                if (not cbr) {
                    otCLI.exitCode = 2;
                    return otCLI.exitCode;
                }
            }
        }
        if (summarizePtr) {
            return summarizePtr(otCLI);
        }
    } catch (std::exception & x) {
        std::cerr << "ERROR. Exiting due to an exception:\n" << x.what() << std::endl;
        otCLI.exitCode = 3;
        return otCLI.exitCode;
    }
    return otCLI.exitCode;
}

template<typename Tree>
inline std::set<long> getAllOTTIds(const Tree& taxonomy) {
    std::set<long> o;
    for (auto nd : iter_node_const(taxonomy)) {
        if (nd->hasOttId()) {
            o.insert(nd->getOttId());
        }
    }
    return o;
}

template<>
inline std::set<long> getAllOTTIds(const TreeMappedWithSplits &taxonomy) {
    return taxonomy.getRoot()->getData().desIds;
}
 
template<typename T>
class TaxonomyDependentTreeProcessor {
    public:
        using node_type = typename T::node_type;
        using node_data_type = typename T::node_type::data_type;
        using tree_data_type = typename T::data_type;
        using tree_type = T;

        std::unique_ptr<T> taxonomy;
        std::set<long> ottIds;

        virtual bool processTaxonomyTree(OTCLI & otCLI) {
            ottIds = getAllOTTIds(*taxonomy);
            if (not taxonomy->getRoot()->hasOttId())
                throw OTCError()<<"Taxonomy root does not have an OTT ID!";
            otCLI.getParsingRules().ottIdValidator = &ottIds;
            otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
            return true;
        }
        virtual bool processSourceTree(OTCLI & , std::unique_ptr<T> tree) {
            assert(tree != nullptr);
            assert(taxonomy != nullptr);
            return true;
        }
        virtual bool summarize(OTCLI &) {
            return true;
        }
        TaxonomyDependentTreeProcessor()
            :taxonomy(nullptr) {
        }
        virtual ~TaxonomyDependentTreeProcessor(){}
    
};

template<typename T>
inline bool taxDependentProcessNextTree(OTCLI & otCLI, std::unique_ptr<T> tree) {
    TaxonomyDependentTreeProcessor<T> * tdtp = static_cast<TaxonomyDependentTreeProcessor<T> *>(otCLI.blob);
    assert(tdtp != nullptr);
    assert(tree != nullptr);
    if (tdtp->taxonomy == nullptr) {
        tdtp->taxonomy = std::move(tree);
        return tdtp->processTaxonomyTree(otCLI);
    }
    return tdtp->processSourceTree(otCLI, std::move(tree));
}

template<typename T>
int taxDependentTreeProcessingMain(OTCLI & otCLI,
                                   int argc,
                                   char *argv[],
                                   TaxonomyDependentTreeProcessor<T> & proc,
                                   unsigned int numTrees,
                                   bool includeInternalNodesInDesIdSets) {
    assert(otCLI.blob == nullptr);
    otCLI.blob = static_cast<void *>(&proc);
    otCLI.getParsingRules().includeInternalNodesInDesIdSets = includeInternalNodesInDesIdSets;
    std::function<bool (OTCLI &, std::unique_ptr<T>)> pcb = taxDependentProcessNextTree<T>;
    auto rc = treeProcessingMain<T>(otCLI, argc, argv, pcb, nullptr, numTrees);
    if (rc == 0) {
        return (proc.summarize(otCLI) ? 0 : 1);
    }
    return rc;
}

inline bool get_bool(const std::string& arg, const std::string& context) {
    if (arg == "true" or arg == "yes" or arg == "True" or arg == "Yes") {
        return true;
    } else if (arg == "false" or arg == "no" or arg == "False" or arg == "No") {
        return false;
    }
    throw OTCError() << context << "'" << arg << "' is not a recognized boolean value.";
}

boost::program_options::options_description standard_options();

boost::program_options::variables_map cmd_line_set_logging(const boost::program_options::variables_map& vm);

std::vector<std::string> cmd_line_response_file_contents(const boost::program_options::variables_map& vm);

boost::program_options::variables_map parse_cmd_line_response_file(int argc, char* argv[],
                                                                   const std::string& message,
                                                                   boost::program_options::options_description visible_,
                                                                   boost::program_options::options_description invisible,
                                                                   boost::program_options::positional_options_description p);

boost::program_options::variables_map parse_cmd_line_standard(int argc, char* argv[],
                                                              const std::string& message,
                                                              boost::program_options::options_description visible_,
                                                              boost::program_options::options_description invisible,
                                                              boost::program_options::positional_options_description p);

} // namespace otc
#endif
