#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
#include <unordered_map>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
	CharString outputFileName;
	CharString referenceFileName;
	CharString insertFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("fasta_insert");
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the reference fasta file that your other file will be inserted into. This might be a reference genome.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "reference-file");
	addOption(parser, seqan::ArgParseOption("i", "insert-file", "Path to the fasta file you wish to insert into your reference file.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "insert-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

	addOption(parser, seqan::ArgParseOption("e", "exclude", "Remove reads in fasta file that have IDs in the text file."));

	setShortDescription(parser, "Methylation tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "August 2017");
	addUsageLine(parser, "-r reference.fasta -i insert.fasta -o result.fasta [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Allows you to insert a fasta file into another fasta file at a given position.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.referenceFileName, parser, "reference-file");
	getOptionValue(options.insertFileName, parser, "insert-file");
	getOptionValue(options.outputFileName, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;
}

/*
Current progress: It compiles!
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	return 0;
}
