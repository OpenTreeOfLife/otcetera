#include "git_version.h"
#include "config.h"
#include "config_file.h"
#include <boost/version.hpp>
#include <iostream>
#include <string>
int main() {
    std::cout << "otcetera version reporter\nSee https://github.com/OpenTreeOfLife/otcetera" << std::endl;
    std::cout << "  version = " << PACKAGE_VERSION << std::endl;
#ifdef GIT_MESSAGE
    std::cout <<"  git_commit = " << GIT_MESSAGE << std::endl;
#endif

    std::cout << "  git_sha = " << GIT_SHAID << std::endl;
#ifdef GIT_COMMIT_DATE
    std::cout << "  git_commit_date = " << GIT_COMMIT_DATE << std::endl;
#endif
    std::cout << "Compiled with:" << std::endl;
    std::cout << "  build_date = " << __DATE__ << " " << __TIME__ << std::endl;
    std::cout << "  Architecture = " << _ARCH_ << std::endl;
#ifdef __GNUC__
    std::cout << "  Compiler = GCC " << __VERSION__ << std::endl;
#endif
    std::cout << "  BOOST_LIB_VERSION = " << BOOST_LIB_VERSION << std::endl;
#ifdef CONFIG_FLAGS
    std::cout << "  CXXFLAGS = " << CONFIG_FLAGS << std::endl;
#endif
    auto configpath = otc::dot_opentree();
    if (configpath) {
        std::cout << "Default config filepath = \"" << *configpath << '\"' << std::endl;
    } else {
        std::cout << "No default config filepath found." << std::endl;
    }
    return 0;
}
