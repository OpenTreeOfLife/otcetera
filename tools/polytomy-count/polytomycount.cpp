#include "otcli.h"
using namespace otc;
bool countPolytomies(OTCLI & otCLI, std::unique_ptr<RootedTree<bool> > tree);

bool countPolytomies(OTCLI & otCLI, std::unique_ptr<RootedTree<bool> > tree) {
	otCLI.out << countPolytomies(*tree) << std::endl;
	return true;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("polytomycount",
				 "takes a filepath to a newick file and reports the number of polytomies in each tree (one line per tree)",
				 {"some.tre"});
	return treeProcessingMain<bool>(otCLI, argc, argv, countPolytomies, nullptr);
}

