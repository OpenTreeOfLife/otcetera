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
namespace otc {

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
			for (const auto & filename : filenameVec) {
				std::ifstream inp;
				if (!openUTF8File(filename, inp)) {
					throw OTCError("Could not open \"" + filename + "\"");
				}
				LOG(INFO) << "reading \"" << filename << "\"...";
				otCLI.currentFilename = filepathToFilename(filename);
				for (;;) {
					std::unique_ptr<T> nt = readNextNewick<T>(inp, filename, otCLI.getParsingRules());
					if (nt == nullptr) {
						break;
					}
					const auto cbr = treePtr(otCLI, std::move(nt));
					if (!cbr) {
						otCLI.exitCode = 2;
						return otCLI.exitCode;
					}
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
			ottIds = taxonomy->getRoot()->getData().desIds;
			otCLI.getParsingRules().ottIdValidator = &ottIds;
			otCLI.getParsingRules().includeInternalNodesInDesIdSets = false;
			return true;
		}
		virtual bool processSourceTree(OTCLI & , std::unique_ptr<T> tree) {
			assert(tree != nullptr);
			assert(taxonomy != nullptr);
			return true;
		}
		virtual bool summarize(const OTCLI &) {
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

} // namespace otc
#endif
