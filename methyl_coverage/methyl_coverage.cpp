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
	bool byContig = false;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("methyl_coverage");
	addOption(parser, seqan::ArgParseOption("cg", "cg-input-file", "Path to the CG input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
	setRequired(parser, "cg-input-file");
	addOption(parser, seqan::ArgParseOption("chg", "chg-input-file", "Path to the CHG input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "chg-input-file");
	addOption(parser, seqan::ArgParseOption("chh", "chh-input-file", "Path to the CHH input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "chh-input-file");
	addOption(parser, seqan::ArgParseOption("r", "ref", "Path to the Genome Reference input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
        setRequired(parser, "ref");
	setShortDescription(parser, "XFTOOLS");
	setVersion(parser, "0.0.4");
	setDate(parser, "March 2019");
	addUsageLine(parser, "--cg sample.cg.w1.gff --chg sample.chg.w1.gff --chh sample.chh.w1.gff -r ref.fasta  [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Calculates methylation coverage histogram");
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

void histo(GffFileIn &gffCGFileIn, GffFileIn &gffCHGFileIn, GffFileIn &gffCHHFileIn, map<string, map<int,int>> &counter, map<int,int> &all)
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

		counter[toCString(record.ref)][c+t]++;
		all[c+t]++;
	}

        while(!atEnd(gffCHGFileIn))
        {
                GffRecord record;
                readRecord(record, gffCHGFileIn);
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

                counter[toCString(record.ref)][c+t]++;
                all[c+t]++;
        }

        while(!atEnd(gffCHHFileIn))
        {
                GffRecord record;
                readRecord(record, gffCHHFileIn);
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

                counter[toCString(record.ref)][c+t]++;
                all[c+t]++;
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

	long long int total = 0;
	for(auto i : contig_c_counts)
		total = total + i.second;

	//create histo
	map<string, map<int, int>> counter;
	map<int, int> all;
	histo(gffCGFileIn, gffCHGFileIn, gffCHHFileIn, counter, all);

	if(options.byContig == true)
	{
		for(auto ref : counter)
		{
			cout << "Histogram for " << ref.first << endl;
			for(auto i : ref.second)
			{
				cout << i.first << " " << i.second << endl;
			}
		}
	}

	cout << "Combined histogram " << endl;
	long long int runningtotal = 0;
	for (auto iter = all.rbegin(); iter != all.rend(); ++iter) {
		runningtotal = runningtotal + iter->second;
		cout << iter->first << " " << iter->second << " " << runningtotal << " " << (double)iter->second / (double)total << " " << (double)runningtotal / (double)total << endl;
	}

        //now we have the running total, we can calculate the median coverage.
        int stop_pos = runningtotal / 2;
        int current_pos = 0;
        int last_coverage = 0;
        for(auto i : all) {
           if(current_pos > stop_pos)
           {
              cout << "Median Cytosine Coverage = " << last_coverage << endl;
              break;
           }
           current_pos = current_pos + i.second;
           last_coverage = i.first;
        }

	close(gffCGFileIn), close(gffCHGFileIn), close(gffCHHFileIn), close(referenceFileIn);

	return 0;
}
