#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
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
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("fast_extract");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

	setShortDescription(parser, "unique_read_extractor");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i input.bam -o output.bam [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Exclude any mapped read that occurs more than one.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.outputFileName, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;
}

/*
Aim: Ronseal, do what it says on the tin

Current progress: It compiles!

*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//open the fasta/fastq file	
	BamFileIn bamFileIn;
	if (!open(bamFileIn, toCString(options.inputFileName)))
	{
		cerr << "ERROR: Could not open file." << endl;
		return 1;
	}

	map<string,int> counter;

	//output stream
	BamFileOut bamFileOut(context(bamFileIn), toCString(options.outputFileName));
	
	try
	{
		//Read the header
		BamHeader header;
		readHeader(header, bamFileIn);
		BamAlignmentRecord record;

        	while (!atEnd(bamFileIn))
        	{
			readRecord(record, bamFileIn);
			counter[toCString(record.qName)]++;
		}

		//at the begining and reread the header to get back to the first read position
		setPosition(bamFileIn,0);
		BamHeader outputheader;
		readHeader(outputheader, bamFileIn);
		writeHeader(bamFileOut, outputheader);

		//now cycl through
		while (!atEnd(bamFileIn))
                {
			readRecord(record, bamFileIn);
			if(counter[toCString(record.qName)] < 2)
			{
				writeRecord(bamFileOut, record);	
			}
		}
	}
	catch (Exception const & e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		return 1;
	}

	

	return 0;
}
