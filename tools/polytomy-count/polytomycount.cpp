#include "pharse.h"
OTCLI gOTCli("polytomycount",
			 "takes a filepath to a newick file and reports the number of polytomies in each tree (one line per tree)",
			 {"some.tre"});

bool countPolytomies(std::unique_ptr<RootedTree<bool> tree) {
	gOTCli.out << tree.countPolytomies() << std::endl;
}

int main(int argc, char *argv[]) {
	return treeProcessingMain<bool>(gOTCli, argc, argv, countPolytomies, nullptr);
}

