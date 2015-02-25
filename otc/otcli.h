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
						  bool (*treePtr)(OTCLI &, std::unique_ptr<T>),
						  int (*summarizePtr)(OTCLI &));

template<typename T>
inline int treeProcessingMain(OTCLI & otCLI,
								 int argc,
								 char * argv[],
								 bool (*treePtr)(OTCLI &, std::unique_ptr<T>),
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
			for (auto filename : filenameVec) {
				std::ifstream inp;
				if (!openUTF8File(filename, inp)) {
					throw OTCError("Could not open \"" + filename + "\"");
				}
				LOG(INFO) << "reading \"" << filename << "\"...";
				for (;;) {
					std::unique_ptr<T> nt = readNextNewick<T>(inp, filename, otCLI.getParsingRules());
					if (nt == nullptr) {
						break;
					}
					auto cbr = treePtr(otCLI, std::move(nt));
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

} // namespace otc
#endif
