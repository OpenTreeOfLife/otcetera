This is a directory of the expected output for test cases of the otcetera "tools."

The tests are invoked by the "make check" target in the tools directory. The
`test_otc_tools.py` looks for subdirectories that have a `*.json` file, and 
reads them to determine an invocation to use (with filepaths that are relative
to the `data` dir) and where to look for the expected output.

By convention the subdirectories of this dir are named to correspond to the 
subdirs of the tools directory which they test.

## Syntax for describing the invocation.
Each JSON just holds a list of objects (Look at the `*.json` files in the subdirs of this directory) with the invocations and a pointer to the location
of the expected output. There are 2 styles supported:

single file invocations like this:

    "invocation" : ["otcpolytomycount", "<INFILE>"],
    "expected": "polytomycount-expected"

where any subdir of `polytomycount-expected` will be used to provide the filename to be used in place of `<INFILE>`. Each of these subdirs should
have the expected output. So adding a new test is just a matter of adding
a new subdir. This invocation only works for scripts that just take 1 filepath
as an argument.

Or with a list of filenames:

    "invocation" : ["otcprunetaxonomy", "<INFILELIST>"],
    "infile_list": ["three-genus-taxonomy.tre", "three-genus-subsample.tre"],
    "expected": "prunesingletons-expected"

Where the "expected" dir holds the output for this single invocation.

## results checked.
A test need not check every aspect of the run. The check is determined by the existence of a file in the expected output location. The files are:

  * `output` text file. Each line will be checked against standard output.
        Test fails if they are not identical
  * `out.tre`. Standard output will be checked for any difference in rooted
        tree sense from this tree. The `tools/rooted-tree-diff.py` script will 
        be used
  * `exit` should hold a non-zero # that is the expected exit code (the lack 
        of the file means that 0 is expected)
