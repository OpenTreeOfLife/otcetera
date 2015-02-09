	std::vector<std::string> args;
	if (!gOTCli.parseArgs(argc, argv, args)) {
		return 1;
	}
	if (args.size() != 2) {
		cerr << "Expecting a tree file and taxonomy tree.\n";
		return 1;
	}
	try{
		gOTCli.readFilepath(args[0], newTreeHook);
		gOTCli.readFilepath(args[1], newTreeHook);
	} catch (NxsException &x) {
		std::cerr << x.what() << "\n";
		return 1;
	}
	if (gOTCli.exitCode == 0) {
		summarize(std::cout);
	}
	return gOTCli.exitCode;
