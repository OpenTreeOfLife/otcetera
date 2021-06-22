#if !defined OTCLI_H
#define OTCLI_H
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include "otc/otc_base_includes.h"
#include "otc/newick.h"
#include "otc/util.h"
#include "otc/tree_iter.h"
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <memory>
#include "assert.hh"
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

        void add_flag(char flag, const std::string & help, bool (*cb)(OTCLI &, const std::string &), bool argNeeded) {
            if (clientDefFlagHelp.find(flag) != clientDefFlagHelp.end()) {
                throw OTCError("Clashing flag assignment");
            }
            clientDefFlagHelp[flag] = help;
            clientDefFlagCallbacks[flag] = cb;
            if (argNeeded) {
                clientDefArgNeeded.insert(flag);
            }
        }
        bool handle_flag(const std::string & flagWithoutDash);
        bool parse_args(int argc, char *argv[], std::vector<std::string> & args);
        void print_help(std::ostream & outStream);
        
        /*int readFilepath(const std::string &fp,
                          ProcessedTreeValidationFunction func=0L,
                          void * blob=0L);
        // reads a NexSON v1.2 at filepath and returns a NxsSimpleTree with the associated treeID
        NxsSimpleTree * readTreeFromNexSONv_1_2(const std::string &filepath, const std::string & tree_id); */
        bool is_dot_txt_file(const std::string &fp);
        auto get_title() const {
            return this->titleStr;
        }
        ParsingRules & get_parsing_rules() {
            return parsingRules;
        }
        const ParsingRules & get_parsing_rules() const {
            return parsingRules;
        }
        void turn_on_verbose_mode();
        void turn_off_verbose_mode();
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
inline bool process_trees(const std::string& filename,
                         const ParsingRules& parsingRules,
                         std::function<bool (std::unique_ptr<T>)> treePtr) {
    std::ifstream inp;
    if (!open_utf8_file(filename, inp)) {
        throw OTCError("Could not open \"" + filename + "\"");
    }
    LOG(INFO) << "reading \"" << filename << "\"...";
    
    ConstStrPtr filenamePtr = ConstStrPtr(new std::string(filename));
    FilePosStruct pos(filenamePtr);
    unsigned treeNumInThisFile = 1;
    for (;;) {
        std::unique_ptr<T> nt = read_next_newick<T>(inp, pos, parsingRules);
        if (nt == nullptr) {
            break;
        }
        if (parsingRules.prune_unrecognized_input_tips) {
            prune_tips_without_ids(*nt);
            if (nt->get_root() == nullptr) {
                continue;
            }
        }
        std::string treeName = std::string("tree ") + std::to_string(treeNumInThisFile++);
        treeName.append(" from ");
        treeName.append(filepath_to_filename(filename));
        nt->set_name(treeName);
        const auto cbr = treePtr(std::move(nt));
        if (not cbr) {
            return false;
        }
    }
    return true;
}

template<typename Tree_t>
std::vector<std::unique_ptr<Tree_t>> get_trees(const std::vector<std::string>& filenames, const ParsingRules& rules) {
    std::vector<std::unique_ptr<Tree_t>> trees;
    std::function<bool(std::unique_ptr<Tree_t>)> proc = [&](std::unique_ptr<Tree_t> t) {trees.push_back(std::move(t));return true;};
    for(auto& filename: filenames) {
        otc::process_trees(filename,rules,proc);
    }
    return trees;
}

template<typename Tree_t>
std::vector<std::unique_ptr<Tree_t>> get_trees(const std::string& filename, const ParsingRules& rules) {
    const std::vector<std::string> filenameVec = {filename};
    return get_trees<Tree_t>(filenameVec, rules);
}

template<typename Tree_t>
std::unique_ptr<Tree_t> get_tree(const std::string& filename, const ParsingRules& rules) {
    auto trees = get_trees<Tree_t>(filename, rules);
    return std::move(trees[0]);
}

template<typename Tree_t>
std::unique_ptr<Tree_t> get_tree(const std::string& filename) {
    ParsingRules rules;
    rules.require_ott_ids = false;
    return get_tree<Tree_t>(filename, rules);
}

template<typename T>
int tree_processing_main(OTCLI & otCLI,
                          int argc,
                          char * argv[],
                          std::function<bool (OTCLI &, std::unique_ptr<T>)> treePtr,
                          int (*summarizePtr)(OTCLI &),
                          std::function<bool(OTCLI &)>,
                          unsigned minNumTrees);

template<typename T>
inline int tree_processing_main(OTCLI & otCLI,
                                 int argc,
                                 char * argv[],
                                 std::function<bool (OTCLI &, std::unique_ptr<T>)> treePtr,
                                 int (*summarizePtr)(OTCLI &),
                                 std::function<bool (OTCLI &)> preTreeHook,
                                 unsigned minNumTrees) {
    std::vector<std::string> filenameVec;
    if (!otCLI.parse_args(argc, argv, filenameVec)) {
        otCLI.exitCode = 1;
        return otCLI.exitCode;
    }
    if (preTreeHook) {
        if (!preTreeHook(otCLI)) {
            return -1;
        }
    }
    
    if (filenameVec.size() < minNumTrees) {
        otCLI.print_help(otCLI.err);
        otCLI.err << otCLI.get_title() << ": Expecting at least " << minNumTrees << " tree filepath(s).\n";
        otCLI.exitCode = 1;
        return otCLI.exitCode;
    }
    try {
        if (treePtr) {
            std::function<bool(std::unique_ptr<T>)> treePtr1 = 
                [&otCLI,&treePtr](std::unique_ptr<T> t) {
                    return treePtr(otCLI,std::move(t));
                };
            
            for (const auto & filename : filenameVec) {
                otCLI.currentFilename = filepath_to_filename(filename);
                const auto cbr = process_trees(filename, otCLI.get_parsing_rules(), treePtr1);
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
inline OttIdSet get_all_ott_ids(const Tree& taxonomy) {
    OttIdSet o;
    for (auto nd : iter_node_const(taxonomy)) {
        if (nd->has_ott_id()) {
            o.insert(nd->get_ott_id());
        }
    }
    return o;
}

template<>
inline OttIdSet get_all_ott_ids(const TreeMappedWithSplits &taxonomy) {
    return taxonomy.get_root()->get_data().des_ids;
}
 
template<typename T>
class TaxonomyDependentTreeProcessor {
    public:
        using node_type = typename T::node_type;
        using node_data_type = typename T::node_type::data_type;
        using tree_data_type = typename T::data_type;
        using tree_type = T;

        std::unique_ptr<T> taxonomy;
        OttIdSet ottIds;
        virtual bool pretree_read_hook(OTCLI & ) {
            return true;
        }

        virtual bool process_taxonomy_tree(OTCLI & otCLI) {
            ottIds = get_all_ott_ids(*taxonomy);
            if (not taxonomy->get_root()->has_ott_id()) {
                throw OTCError() << "Taxonomy root does not have an OTT ID!";
            }
            otCLI.get_parsing_rules().ott_id_validator = &ottIds;
            otCLI.get_parsing_rules().include_internal_nodes_in_des_id_sets = false;
            return true;
        }
        virtual bool process_source_tree(OTCLI & , std::unique_ptr<T> tree) {
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
inline bool tax_dependent_process_next_tree(OTCLI & otCLI, std::unique_ptr<T> tree) {
    TaxonomyDependentTreeProcessor<T> * tdtp = static_cast<TaxonomyDependentTreeProcessor<T> *>(otCLI.blob);
    assert(tdtp != nullptr);
    assert(tree != nullptr);
    if (tdtp->taxonomy == nullptr) {
        tdtp->taxonomy = std::move(tree);
        return tdtp->process_taxonomy_tree(otCLI);
    }
    return tdtp->process_source_tree(otCLI, std::move(tree));
}


template<typename T>
int tax_dependent_tree_processing_main(OTCLI & otCLI,
                                   int argc,
                                   char *argv[],
                                   TaxonomyDependentTreeProcessor<T> & proc,
                                   unsigned int num_trees,
                                   bool include_internal_nodes_in_des_id_sets) {
    assert(otCLI.blob == nullptr);
    otCLI.blob = static_cast<void *>(&proc);
    otCLI.get_parsing_rules().include_internal_nodes_in_des_id_sets = include_internal_nodes_in_des_id_sets;
    std::function<bool (OTCLI &, std::unique_ptr<T>)> pcb = tax_dependent_process_next_tree<T>;
    using namespace std::placeholders;
    //auto prh = std::bind(&TaxonomyDependentTreeProcessor<T>::pretree_read_hook, proc, _2);
    std::function<bool(OTCLI &) > prh = [&proc] (OTCLI & o) {return proc.pretree_read_hook(o);};
    auto rc = tree_processing_main<T>(otCLI, argc, argv, pcb, nullptr, prh, num_trees);
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
