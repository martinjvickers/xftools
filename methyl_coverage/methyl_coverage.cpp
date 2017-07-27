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
#include "boost/lexical_cast.hpp"
using boost::bad_lexical_cast;
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
	CharString cgInputFile;
	CharString chgInputFile;
	CharString chhInputFile;
	CharString referenceFile;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("single-c_combine");
	addOption(parser, seqan::ArgParseOption("cg", "cg-input-file", "Path to the CG input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
	setRequired(parser, "cg-input-file");
	addOption(parser, seqan::ArgParseOption("chg", "chg-input-file", "Path to the CHG input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "chg-input-file");
	addOption(parser, seqan::ArgParseOption("chh", "chh-input-file", "Path to the CHH input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "chh-input-file");
	addOption(parser, seqan::ArgParseOption("r", "ref", "Path to the Genome Reference input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "ref");
	setShortDescription(parser, "XFTOOLS");
	setVersion(parser, "0.0.1");
	setDate(parser, "Apr 2017");
	addUsageLine(parser, "--cg sample.cg.w1.gff --chg sample.chg.w1.gff --chh sample.chh.w1.gff -r ref.fasta  [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Calculate coverage");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.cgInputFile, parser, "cg-input-file");
	getOptionValue(options.chgInputFile, parser, "chg-input-file");
	getOptionValue(options.chhInputFile, parser, "chh-input-file");
	getOptionValue(options.referenceFile, parser, "ref");

	return seqan::ArgumentParser::PARSE_OK;

}

//count
int count_cs(ModifyStringOptions options, SeqFileIn &referenceFileIn, map<string, int> &map)
{
	CharString refid;
	IupacString refseq;
	SeqFileIn refFileIn;

	while(!atEnd(referenceFileIn))
	{
		int counter = 0;
		readRecord(refid, refseq, referenceFileIn);

		for(int i = 0; i < length(refseq); i++)
			if(refseq[i] == 'C' || refseq[i] == 'G')
				counter++;

		map.insert(pair<string,int> (toCString(refid), counter));
	}

	return 0;
}

void histo(GffFileIn &gffCGFileIn, GffFileIn &gffCHGFileIn, GffFileIn &gffCHHFileIn, map<int,int> &counter)
{
	// value (as in c+t) followed by counter
//	map<int, int> counter;

	while(!atEnd(gffCGFileIn))
        {
		GffRecord record;
		readRecord(record, gffCGFileIn);
		int c;
		int t;
		
		// go through tags
		for(int i = 0; i < length(record.tagNames); i++)
		{
			if(record.tagNames[i] == "c")
				c = stoi(toCString(record.tagValues[i]));
			if(record.tagNames[i] == "t")
				t = stoi(toCString(record.tagValues[i]));
		}

		counter[c+t]++;
	}
}

int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//open first file
	GffFileIn gffCGFileIn, gffCHGFileIn, gffCHHFileIn;
	if ( !open(gffCGFileIn, toCString(options.cgInputFile)) || !open(gffCHGFileIn, toCString(options.chgInputFile)) || !open(gffCHHFileIn, toCString(options.chhInputFile)) )
        {
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

	SeqFileIn referenceFileIn;
	if ( !open(referenceFileIn, toCString(options.referenceFile)) )
	{
		std::cerr << "ERROR: Could not open the file.\n";
                return 1;
	}

	map<string, int> contig_c_counts;

	// create a map that is <chr,total_c's> from our reference
	count_cs(options, referenceFileIn, contig_c_counts);

//	for(auto i : contig_c_counts)
//		cout << i.first << " " << i.second << endl;

	//create histo
	map<int, int> counter;
	histo(gffCGFileIn, gffCHGFileIn, gffCHHFileIn, counter);

	for(auto i : counter)
		cout << i.first << " " << i.second << endl;

	close(gffCGFileIn), close(gffCHGFileIn), close(gffCHHFileIn), close(referenceFileIn);

	return 0;
}
