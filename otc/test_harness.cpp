#include "otc/test_harness.h"
#include "otc/otcli.h"
namespace otc {

TestHarness::TestHarness(int argc, char *argv[])
    :initFailed(true) {
    OTCLI otCLI(argv[0], "a test harness", "path/to/a/data/dir", true);
    //otCLI.turn_on_verbose_mode();
    std::vector<std::string> dataDirArgVec;
    if (otCLI.parse_args(argc, argv, dataDirArgVec)) {
        if (dataDirArgVec.size() != 1) {
            std::cerr << "path to data file dir is required as the only argument.\n";
        } else {
            dataDir = dataDirArgVec[0];
            initFailed = false;
        }
    }
}

int TestHarness::run_tests(const TestsVec & tests) {
    if (this->initFailed) {
        std::cerr << "ERROR: " << tests.size() << " unavaible due to incorrect initializaion of the TestHarness.\n";
        return static_cast<int>(tests.size());
    }
    std::vector<std::string> failures;
    for (const auto & tpIt : tests) {
        const auto & fn = tpIt.second;
        char resp;
#if defined(SUPPRESS_TEST_EXCEPTIONS)
        try {
#endif
            resp = fn(*this);
#if defined(SUPPRESS_TEST_EXCEPTIONS)
        } catch (std::exception & x) {
            std::cerr << "exception in test: " << x.what() << '\n';
            resp = 'E';
        } catch (...) {
            std::cerr << "unrecognized exception in test.\n";
            resp = 'E';
        }
#endif
        std::cerr << resp;
        if (resp != '.') {
            failures.push_back(tpIt.first);
        }
    }
    std::cerr << '\n';
    for (const auto & f : failures) {
        std::cerr << f << '\n';
    }
    if (!failures.empty()) {
        std::cerr << "FAILED " << failures.size() << "/" << tests.size() << " tests.\n";
    }
    return static_cast<int>(failures.size());
}

} // namespace otc
