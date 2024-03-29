#include "otc/otcli.h"
#include <cstring>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include <g3log/logworker.hpp>
#include <g3log/loglevels.hpp>

///////////////////////////////////////////////////////////////
// pragmas are MTH mods to silence clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored  "-Wglobal-constructors"
#pragma clang diagnostic ignored  "-Wexit-time-destructors"

#pragma clang diagnostic pop
///////////////////////////////////////////////////////////////

namespace po = boost::program_options;
using po::variables_map;
using std::string;
using std::vector;
using std::optional;
namespace fs = std::filesystem;

namespace otc {
bool debugging_output_enabled = false;

std::unique_ptr<g3::LogWorker>  default_worker;

class ConsoleSink
{
public:
    void logToStderr(g3::LogMessageMover logEntry) {
        std::cerr << logEntry.get().timestamp()<<" ["<<logEntry.get().level()<<"]  "<<logEntry.get().message() << std::endl;
    }
};

class FileSink
{
    std::ofstream out;
public:
    void logToFile(g3::LogMessageMover logEntry) {
        out << logEntry.get().timestamp()<<" ["<<logEntry.get().level()<<"]  "<<logEntry.get().message() << std::endl;
    }

    FileSink(const fs::path& filepath)
        :out(filepath, std::ios::app)
        {
            out<<"======================== BEGIN ========================\n";
        }
};

void initialize_logging(const string& name, std::optional<fs::path> logfile_dir)
{
    fs::path logfile = fs::path(name+".log.txt").filename();

    // This is for backward compatibility.  After all users switch to the command-line flag we can remove it.
    if (auto logfile_dir_env = std::getenv("OTCETERA_LOGFILE"); logfile_dir_env and not logfile_dir)
        logfile_dir = string(logfile_dir_env);

    if (logfile_dir)
        logfile = *logfile_dir / logfile;

    auto worker = g3::LogWorker::createLogWorker();
    // Log to console
    auto handle1 = worker->addSink(std::make_unique<ConsoleSink>(), &ConsoleSink::logToStderr);
    // Log to file
    auto handle2 = worker->addSink(std::make_unique<FileSink>(logfile), &FileSink::logToFile);

    // Initialize Logging
    g3::only_change_at_initialization::addLogLevel(ERROR);
    g3::only_change_at_initialization::addLogLevel(TRACE);
    g3::initializeLogging(worker.get());
    default_worker = std::move(worker);
}

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
        err(std::cerr)
{
    if (quietExecution)
        g3::log_levels::disableAll();
    else
    {
        g3::log_levels::enable(TRACE);
        g3::log_levels::enable(DEBUG);
    }
}

void OTCLI::print_help(std::ostream & outStream) {
    outStream << this->titleStr << ": " << this->descriptionStr << ".\n";
    outStream << "The most common usage is simply:\n";
    outStream << "    " << this->titleStr << " " << this->usageStr << "\n";
    outStream << "Standard command-line flags:\n";
    outStream << "    -h on the command line shows this help message\n";
    outStream << "    -fFILE treat each line of FILE as an arg\n";
    outStream << "    -l logfile directory.\n";
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

void OTCLI::turn_on_verbose_mode() {
    this->verbose = true;
    debugging_output_enabled = true;
    g3::log_levels::enable(DEBUG);
}

void OTCLI::turn_off_verbose_mode() {
    this->verbose = false;
    debugging_output_enabled = false;
    g3::log_levels::disable(DEBUG);
}

bool OTCLI::handle_flag(const std::string & flagWithoutDash) {
    bool recursionNeeded = false;
    auto f = flagWithoutDash[0];
    if (clientDefFlagCallbacks.find(f) != clientDefFlagCallbacks.end()) {
        const auto cb = clientDefFlagCallbacks[f];
        try {
            const auto n = flagWithoutDash.substr(1);
            if (contains(clientDefArgNeeded, f)) {
                if (n.empty()) {
                    this->err << "Expecting an argument value after the  -" << f << " flag.\n";
                    return false;
                }
            } else if (not n.empty()) {
                this->err << "Unexpected argument value after the  -" << f << " flag.\n";
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
        this->print_help(this->out);
        this->exitCode = 1;
        return false;
    } else if (f == 'v') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        turn_on_verbose_mode();
    } else if (f == 'l') {
        if (flagWithoutDash.length() == 1) {
            this->err << "Expecting an argument value after the  -l flag.\n";
            return false;
        }
        logfile_dir_arg = flagWithoutDash.substr(1);
    } else if (f == 't') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        this->verbose = true;
        g3::log_levels::enable(TRACE);
    } else if (f == 'q') {
        if (flagWithoutDash.length() > 1) {
            recursionNeeded = true;
        }
        this->verbose = false;
        debugging_output_enabled = false;
        g3::log_levels::disableAll();
    } else if (f == 'f') {
        if (flagWithoutDash.length() == 1) {
            this->err << "Expecting an argument value after the  -f flag.\n";
            return false;
        }
        extraArgs = read_lines_of_file(flagWithoutDash.substr(1));
    } else {
        std::cerr  << "Unrecognized flag: -" + flagWithoutDash << '\n';
        return false;
    }
    if (recursionNeeded) {
        return this->handle_flag(flagWithoutDash.substr(1));
    }
    return true;
}

bool OTCLI::parse_args(int argc, char *argv[], std::vector<std::string> & args) {
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
            if (!this->handle_flag(flagWithoutDash)) {
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
    initialize_logging(argv[0], logfile_dir_arg);
    return true;
}

bool OTCLI::is_dot_txt_file(const std::string &fp) {
    const size_t fnl = fp.length();
    return (fp.substr(fnl - 4) == std::string(".txt"));
}

po::options_description standard_options() {
    using namespace po;
    options_description standard("Standard command-line flags");
    standard.add_options()
      ("help,h", "Produce help message")
      ("response-file,f", value<string>(), "Treat contents of file <arg> as a command line.")
      ("quiet,q","QUIET mode (all logging disabled)")
      ("trace,t","TRACE level debugging (very noisy)")
      ("logdir,l", value<string>(), "Directory to put log file in")
      ("verbose,v","verbose")
    ;
    return standard;
}

variables_map cmd_line_set_logging(const string& name, const po::variables_map& vm)
{
    optional<fs::path> logdir;
    if (vm.count("logdir"))
        logdir = vm["logdir"].as<string>();

    initialize_logging(name, logdir);

    if (vm.count("quiet")) {
        g3::log_levels::disableAll();
    } else {
        g3::log_levels::disable(TRACE);
        g3::log_levels::disable(DEBUG);
        if (vm.count("trace")) {
            g3::log_levels::enable(TRACE);
            g3::log_levels::enable(DEBUG);
        }
        if (vm.count("verbose")) {
            g3::log_levels::enable(DEBUG);
        }
    }
    return vm;
}

vector<string> read_response_file_content(const char * fn) {
    vector<string> args;
    std::ifstream ifs(fn);
    if (not ifs) {
        throw OTCError() << "Could not open the response file \"" << fn << "\" !";
    }
    // Read the whole file into a string
    std::stringstream ss;
    ss << ifs.rdbuf();
    // Split the file content
    boost::char_separator<char> sep(" \n\r");
    std::string ResponsefileContents( ss.str() );
    boost::tokenizer<boost::char_separator<char> > tok(ResponsefileContents, sep);
    copy(tok.begin(), tok.end(), back_inserter(args));
    return args;
}

vector<string> cmd_line_response_file_contents(const po::variables_map& vm) {
    vector<string> args;
    if (vm.count("response-file")) {
        return read_response_file_content(vm["response-file"].as<string>().c_str());
    }
    return args;
}

vector<string> expand_for_response_file(int argc, char* argv[]) {
    vector<string> args;
    for (int i = 1; i < argc; ++i) {
        std::string na = argv[i];
        bool insert_file_flag = false;
        std::string rfp;
        if (na.size() > 1) {
            if (na[0] == '-') {
                if (na[1] == 'f') {
                    insert_file_flag = true;
                    if (na.size() > 2) {
                        rfp = na.substr(2);
                    } else {
                        ++i;
                        if (i >= argc) {
                            throw OTCError() << "command line cannot end with -f\n";
                        }
                        rfp = string(argv[i]);
                    }
                } else {
                    auto pos = na.find("--response-file");
                    if (pos == 0) {
                        insert_file_flag = true;
                        if (pos != std::string::npos) {
                            std::size_t fchar = 15;
                            if (na.size() > 15) {
                                if (na[fchar] == '=') {
                                    ++fchar;
                                }
                            }
                            if (na.size() > fchar) {
                                rfp = na.substr(fchar);
                            } else {
                                ++i;
                                if (i >= argc) {
                                    throw OTCError() << "command line cannot end with --response-file\n";
                                }
                                rfp = string(argv[i]);                   
                            }
                        }
                    }
                }
            }
        }
        if (insert_file_flag) {
            if (rfp.size() < 1) {
                throw OTCError() << "cannot have an empty string as a --response-file\n";
            }
            auto a = read_response_file_content(rfp.c_str());
            args.insert(args.end(), a.begin(), a.end());
        } else {
            args.push_back(na);
        }
    }
    return args;
}


variables_map parse_cmd_line_response_file(int argc, char* argv[],
                                           po::options_description visible,
                                           po::options_description invisible,
                                           po::positional_options_description p) {
    using namespace po;
    auto expandedargs = expand_for_response_file(argc, argv);
    variables_map vm;
    options_description all;
    all.add(invisible).add(visible);
    store(command_line_parser(expandedargs).options(all).positional(p).run(), vm);
    notify(vm);
    return vm;
}

variables_map parse_cmd_line_standard(int argc, char* argv[],
                                      const string& message,
                                      po::options_description visible,
                                      po::options_description invisible,
                                      po::positional_options_description p) {
    using namespace po;
    variables_map vm = parse_cmd_line_response_file(argc, argv, visible, invisible, p);
    if (vm.count("help")) {
        std::cout << message << "\n";
        std::cout << visible << "\n";
        if (vm.count("verbose")) {
            std::cout << invisible << "\n";
        }
        exit(0);
    }
    cmd_line_set_logging(argv[0], vm);
    return vm; 
}


} // namespace otc
