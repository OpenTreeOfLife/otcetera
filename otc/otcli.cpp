#include "otc/otcli.h"
#include <cstring>
#include <string>

///////////////////////////////////////////////////////////////
// pragmas are MTH mods to silence clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored  "-Wglobal-constructors"
#pragma clang diagnostic ignored  "-Wexit-time-destructors"
INITIALIZE_EASYLOGGINGPP

static el::Configurations defaultConf;

#pragma clang diagnostic pop
///////////////////////////////////////////////////////////////

namespace otc {
bool debuggingOutputEnabled = false;
long ottIDBeingDebugged = -1;
OTCLI::OTCLI(const char *title,
          const char *descrip,
          const char *usage,
          bool quietExecution)
        :exitCode(0),
        verbose(false),
        currReadingDotTxtFile(false),
        blob(nullptr),
        titleStr(title),
        descriptionStr(descrip),
        usageStr(usage),
        out(std::cout),
        err(std::cerr) {
    defaultConf.setToDefault();
    if (quietExecution) {
        defaultConf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
    } else {
        defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");
        defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
    }
    el::Loggers::reconfigureLogger("default", defaultConf);
}

void OTCLI::printHelp(std::ostream & outStream) {
    outStream << this->titleStr << ": " << this->descriptionStr << ".\n";
    outStream << "The most common usage is simply:\n";
    outStream << "    " << this->titleStr << " " << this->usageStr << "\n";
    outStream << "Standard command-line flags:\n";
    outStream << "    -h on the command line shows this help message\n";
    outStream << "    -fFILE treat each line of FILE as an arg\n";
    outStream << "    -q QUIET mode (all logging disabled)\n";
    outStream << "    -t TRACE level debugging (very noisy)\n";
    outStream << "    -v verbose\n";
    if (!clientDefFlagHelp.empty()) {
        outStream << "Exe-specific command-line flags:\n";
    }
    for (const auto & c : clientDefFlagHelp) {
        if (contains(clientDefArgNeeded, c.first)) {
            outStream << "    -" << c.first << "ARG " << c.second << '\n';
        } else {
            outStream << "    -" << c.first << ' ' << c.second << '\n';
        }
    }
}

void OTCLI::turnOnVerboseMode() {
    this->verbose = true;
    debuggingOutputEnabled = true;
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

void OTCLI::turnOffVerboseMode() {
    this->verbose = false;
    debuggingOutputEnabled = false;
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

bool OTCLI::handleFlag(const std::string & flagWithoutDash) {
    bool recursionNeeded = false;
    auto f = flagWithoutDash[0];
    if (clientDefFlagCallbacks.find(f) != clientDefFlagCallbacks.end()) {
        const auto cb = clientDefFlagCallbacks[f];
        try {
            const auto n = flagWithoutDash.substr(1);
            if (contains(clientDefArgNeeded, f) && n.empty()) {
                this->err << "Expecting an argument value after the  -" << f << " flag.\n";
                return false;
            }
            const auto rc = cb(*this, n);
            return rc;
        } catch (std::exception & x) {
            this->err << x.what() << '\n';
            return false;
        }
    }
    if (f == 'h') {
        this->printHelp(this->out);
        this->exitCode = 1;
        return false;
    } else if (f == 'v') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        turnOnVerboseMode();
    } else if (f == 't') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        this->verbose = true;
        defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "true");
        el::Loggers::reconfigureLogger("default", defaultConf);
    } else if (f == 'q') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        this->verbose = false;
        debuggingOutputEnabled = false;
        defaultConf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
        el::Loggers::reconfigureLogger("default", defaultConf);
    } else if (f == 'f') {
        if (flagWithoutDash.length() == 1) {
            this->err << "Expecting an argument value after the  -f flag.\n";
            return false;
        }
        extraArgs = readLinesOfFile(flagWithoutDash.substr(1));
    } else {
        std::cerr  << "Unrecognized flag: -" + flagWithoutDash << '\n';
        return false;
    }
    if (recursionNeeded) {
        return this->handleFlag(flagWithoutDash.substr(1));
    }
    return true;
}
bool OTCLI::parseArgs(int argc, char *argv[], std::vector<std::string> & args) {
    this->exitCode = 0;
    std::list<std::string> allArgs;
    for (int i = 1; i < argc; ++i) {
        const std::string filepath = argv[i];
        allArgs.push_back(filepath);
    }
    for (auto aIt = allArgs.begin(); aIt != allArgs.end(); ++aIt) {
        const auto & filepath = *aIt;
        const auto slen = filepath.length();
        if (slen > 1U && filepath[0] == '-') {
            extraArgs.clear();
            const std::string flagWithoutDash = filepath.substr(1);
            if (!this->handleFlag(flagWithoutDash)) {
                return false;
            }
            if (!extraArgs.empty()) {
                auto b = aIt;
                ++b;
                allArgs.insert(b, begin(extraArgs), end(extraArgs));
            }
        } else {
            const std::string filepathstr(filepath);
            args.push_back(filepathstr);
        }
    }
    return true;
}

bool OTCLI::isDotTxtFile(const std::string &fp) {
    const size_t fnl = fp.length();
    return (fp.substr(fnl - 4) == std::string(".txt"));
}

} // namespace otc
