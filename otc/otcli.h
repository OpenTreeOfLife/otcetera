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

		bool handleFlag(const std::string & flagWithoutDash);
		bool parseArgs(int argc, char *argv[], std::vector<std::string> & args);
		void printHelp(std::ostream & outStream);
		
		/*int readFilepath(const std::string &fp,
						  ProcessedTreeValidationFunction func=0L,
						  void * blob=0L);
		// reads a NexSON v1.2 at filepath and returns a NxsSimpleTree with the associated treeID
		NxsSimpleTree * readTreeFromNexSONv_1_2(const std::string &filepath, const std::string & tree_id); */
		bool isDotTxtFile(const std::string &fp);
	private:
		std::string titleStr;
		std::string descriptionStr;
		std::string usageStr;
	public:
		std::ostream & out;
		std::ostream & err;
};

template<typename T>
int treeProcessingMain(OTCLI & otCLI,
						  int argc,
						  char * argv[],
						  bool (*treePtr)(OTCLI &, std::unique_ptr<RootedTree<T> >),
						  int (*summarizePtr)(OTCLI &));

template<typename T>
inline int treeProcessingMain(OTCLI & otCLI,
								 int argc,
								 char * argv[],
								 bool (*treePtr)(OTCLI &, std::unique_ptr<RootedTree<T> >),
								 int (*summarizePtr)(OTCLI &)) {
	std::vector<std::string> filenameVec;
	if (!otCLI.parseArgs(argc, argv, filenameVec)) {
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
				for (;;) {
					std::unique_ptr<RootedTree<T> > nt = readNextNewick<T>(inp, filename);
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
