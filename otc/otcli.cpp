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


OTCLI::OTCLI(const char *title,
		  const char *descrip,
		  const char *usage)
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
	defaultConf.set(el::Level::Trace, 
					el::ConfigurationType::Enabled, "false");
	defaultConf.set(el::Level::Debug, 
					el::ConfigurationType::Enabled, "false");
	el::Loggers::reconfigureLogger("default", defaultConf);
}


void OTCLI::printHelp(std::ostream & outStream) {
	outStream << this->titleStr << ": " << this->descriptionStr << ".\n";
	outStream << "\nThe most common usage is simply:\n";
	outStream << this->titleStr << " " << this->usageStr << "\n";
	outStream << "\nCommand-line flags:\n\n";
	outStream << "    -h on the command line shows this help message\n\n";
	outStream << "    -v verbose outStreamput\n\n";
}

bool OTCLI::handleFlag(const std::string & flagWithoutDash) {
	bool recursionNeeded = false;
	if (flagWithoutDash[0] == 'h') {
		this->printHelp(this->out);
		this->exitCode = 1;
		return false;
	} else if (flagWithoutDash[0] == 'v') {
		if (flagWithoutDash.length() > 1) {
			recursionNeeded = true;
		}
		this->verbose = true;
		defaultConf.set(el::Level::Debug, 
					el::ConfigurationType::Enabled, "true");
		el::Loggers::reconfigureLogger("default", defaultConf);
	}
	// this is where, we'd insert a procedure for generic flag handling...
	// IF WE HAD ONE!
	if (recursionNeeded) {
		return this->handleFlag(flagWithoutDash.substr(1));
	}
	return true;
}
bool OTCLI::parseArgs(int argc, char *argv[], std::vector<std::string> & args) {
	this->exitCode = 0;
	for (int i = 1; i < argc; ++i) {
		auto filepath = argv[i];
		auto slen = strlen(filepath);
		if (slen > 1U && filepath[0] == '-') {
			const std::string flagWithoutDash(filepath + 1);
			if (!this->handleFlag(flagWithoutDash)) {
				return false;
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
