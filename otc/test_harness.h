#if !defined (TEST_HARNESS_H)
#define TEST_HARNESS_H
#include <functional>
#include "otc/otc_base_includes.h"
namespace otc {

template<typename T>
bool test_vec_element_equality(const std::vector<T> & expected, const std::vector<T> & obtained) noexcept;
template<typename T>
void test_complete_diff_message(const T & expected, const T &obtained) noexcept;

template<typename T>
void test_complete_diff_message(const T & expected, const T &obtained) noexcept {
    std::cerr << expected << " != " << obtained << std::endl;
}
template<>
void test_complete_diff_message<std::string>(const std::string & expected, const std::string &obtained) noexcept {
    std::cerr << '\"' << expected << "\" != \"" << obtained << '\"' << std::endl;
}

template<typename T>
bool test_vec_element_equality(const std::vector<T> & expected, const std::vector<T> & obtained) noexcept {
    bool differed = false;
    try {
        const auto minlen = std::min(expected.size(), obtained.size());
        for (auto i = 0U; i < minlen; ++i) {
            const auto & e = expected[i];
            const auto & o = obtained[i];
            if (e != o) {
                std::cerr << "  Element with index=" << i << " differ [obtained != expected] ";
                test_complete_diff_message(o, e);
                differed = true;
            }
        }
        const auto lendiff = static_cast<long>(expected.size()) - static_cast<long>(obtained.size());
        if (lendiff > 0) {
            std::cerr << lendiff << " too few elements obtained.\n";
        } else if (lendiff < 0) {
            std::cerr << -lendiff << " too many elements obtained.\n";
        }
        differed = (differed || (lendiff != 0));
    } catch (...) {
        return false;
    }
    return !differed;
}

class TestHarness;
typedef std::function<char(const TestHarness &)> TestCallBack;
typedef std::pair<const std::string, TestCallBack> TestFn;
typedef std::vector<TestFn> TestsVec;

class TestHarness {
    private:
        std::string dataDir; //dir that holds the tests
        bool initFailed;
    public:
        TestHarness(int argc, char * argv []);
        int run_tests(const TestsVec & tests);
        std::string get_filepath(const std::string & fn) const {
            auto fp = dataDir;
            fp += "/";
            fp += fn;
            return fp;
        }
        bool open_test_file(const std::string & fn, std::ifstream &inp) const {
            inp.open(get_filepath(fn));
            return inp.good();
        }
};

} //namespace otc
#endif
