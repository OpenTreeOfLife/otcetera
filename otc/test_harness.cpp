#include "otc/test_harness.h"
#include "otc/otcli.h"
namespace otc {

TestHarness::TestHarness(int argc, char *argv[])
	:initFailed(true) {
	OTCLI otCLI(argv[0], "a test harness", "path/to/a/data/dir");
	std::vector<std::string> dataDirArgVec;
	if (otCLI.parseArgs(argc, argv, dataDirArgVec)) {
		if (dataDirArgVec.size() != 1) {
			std::cerr << "path to data file dir is required as the only argument.\n";
		} else {
			dataDir = dataDirArgVec[0];
			initFailed = false;
		}
	}
}


int TestHarness::runTests(const TestsVec & tests) {
	if (this->initFailed) {
		LOG(ERROR) << tests.size() << " unavaible due to incorrect initializaion of the TestHarness.\n";
		return (int)tests.size();
	}
	std::vector<std::string> failures;
	for (auto tpIt : tests) {
		auto fn = tpIt.second;
		char resp;
		try {
			resp = fn(*this);
		} catch (...) {
			resp = 'E';
		}
		std::cerr << resp;
		if (resp != '.') {
			failures.push_back(tpIt.first);
		}
	}
	std::cerr << '\n';
	for (auto f : failures) {
		std::cerr << f << '\n';
	}
	if (!failures.empty()) {
		std::cerr << "FAILED " << failures.size() << " out of " << tests.size() << "tests.\n";
	}
	return (int) failures.size();
}

} // namespace otc
