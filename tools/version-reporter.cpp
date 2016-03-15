#include "gitversion.h"
#include "config.h"
#include <boost/version.hpp>
#include <iostream>
int main() {
    std::cout << "otcetera version reporter\nSee https://github.com/OpenTreeOfLife/otcetera" << std::endl;
    std::cout << "version = " << VERSION << std::endl;
    std::cout << "git_sha = " << gitversion << std::endl;
    std::cout << "Compiled with:" << std::endl;
    std::cout << "BOOST_LIB_VERSION = " << BOOST_LIB_VERSION << std::endl;
    return 0;
}