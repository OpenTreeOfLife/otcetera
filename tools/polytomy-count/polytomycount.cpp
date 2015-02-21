#include "otc/otcli.h"
#include "otc/tree_operations.h"
using namespace otc;
template<typename T, typename U>
bool writeNumPolytomies(OTCLI & otCLI, std::unique_ptr<RootedTree<T, U> > tree);

template<typename T, typename U>
inline bool writeNumPolytomies(OTCLI & otCLI, std::unique_ptr<RootedTree<T, U> > tree) {
	otCLI.out << countPolytomies(*tree) << std::endl;
	return true;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("polytomycount",
				 "takes a filepath to a newick file and reports the number of polytomies in each tree (one line per tree)",
				 {"some.tre"});
	return treeProcessingMain<RTNodeNoData, RTreeNoData>(otCLI,
														 argc,
														 argv,
														 writeNumPolytomies,
														 nullptr);
}

