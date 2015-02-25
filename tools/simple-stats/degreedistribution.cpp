#include "otc/otcli.h"
#include "otc/tree_operations.h"
using namespace otc;
typedef RootedTree<RTNodeNoData, RTreeNoData> Tree_t;

template<typename T>
bool writeDegreeDistribution(OTCLI & , std::unique_ptr<T> tree);

template<typename T>
inline bool writeDegreeDistribution(OTCLI & , std::unique_ptr<T> tree) {
	std::map<unsigned long, unsigned long> degreeDistribution;
	for (auto nd : ConstPreorderIter<T>(*tree)) {
		auto od = nd->getOutDegree();
		degreeDistribution[od] += 1;
	}
	std::cout << "Out-degree\tCount\n";
	for (auto p: degreeDistribution) {
		std::cout << p.first << "\t" << p.second << "\n";
	}
	return true;
}

int main(int argc, char *argv[]) {
	OTCLI otCLI("otcdegreedistribution",
				 "takes a filepath to a newick file and reports the number of nodes of each out-degree",
				 {"some.tre"});
	return treeProcessingMain<Tree_t>(otCLI, argc, argv, writeDegreeDistribution, nullptr);
}

