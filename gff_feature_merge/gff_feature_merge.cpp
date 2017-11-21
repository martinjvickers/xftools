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
#include <cassert>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
using namespace seqan;
using namespace std;


#include <unordered_map>
#include <set>

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString outputFileName;
	int size;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("gff_feature_merge");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	addOption(parser, seqan::ArgParseOption("s", "size", "Merge Size", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "size", "4");

	setShortDescription(parser, "GFF Feature Merge");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i input.gff-o output.gff [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Given a GFF of features, this program will merge features into single features if they lay within a defined (user supplied) range.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	getOptionValue(options.size, parser, "size");

	return seqan::ArgumentParser::PARSE_OK;
}

/*
Aim:
Current progress: It compiles!
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//get our value and add value to map
	GffFileIn gffIn;
	GffRecord record;
        if (!open(gffIn, toCString(options.inputFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }

	//output GFF to cout
	GffFileOut gffOut(std::cout, Gff());

	GffRecord last_record;
	GffRecord merging;
	int count = 0;
	int nvalue = 0;

	while (!atEnd(gffIn)) // loop through GFF file
        {
		readRecord(record, gffIn); // read record

		if(count==0) // if it's the first record
		{
			merging = record; // mark as our merging read.
			for(int i = 0; i < length(record.tagNames); i++)
                                {
                                        if(record.tagNames[i] == "n")
                                        {
                                                nvalue = nvalue + stoi(toCString(record.tagValues[i]));
                                        }
                        }

		} else {
			if( (record.ref == merging.ref) && ((record.beginPos - (merging.endPos+1)) < options.size) )
	               	{
				merging.endPos = record.endPos;

				//get the n= tag and add the current one to it.
				for(int i = 0; i < length(record.tagNames); i++)
				{
					if(record.tagNames[i] == "n")
					{
						nvalue = nvalue + stoi(toCString(record.tagValues[i]));
					}
				}
			} else {
				//this is where we have finished the merging
				//cout << merging << endl;
				clear(merging.tagNames);
				clear(merging.tagValues);
				merging.score = GffRecord::INVALID_SCORE();
				appendValue(merging.tagNames, "windowSize");
				appendValue(merging.tagValues, to_string((merging.endPos+1)-merging.beginPos));
				appendValue(merging.tagNames, "n");
				appendValue(merging.tagValues, to_string(nvalue));
				merging.score = (float)nvalue / (float) ((merging.endPos+1)-merging.beginPos);
				writeRecord(gffOut, merging);
				clear(merging);
				merging = record;
				nvalue = 0;
				
				for(int i = 0; i < length(record.tagNames); i++)
				{
					//now get the current value
					if(record.tagNames[i] == "n")
                	                {
                                		nvalue = nvalue + stoi(toCString(record.tagValues[i]));
                                	}
				}
			}
		}
		//move along here, nothing to see
		count++;
	}

	//there is always one last result left over
	clear(merging.tagNames);
	clear(merging.tagValues);
	merging.score = GffRecord::INVALID_SCORE();
	appendValue(merging.tagNames, "windowSize");
	appendValue(merging.tagValues, to_string((merging.endPos+1)-merging.beginPos));
	appendValue(merging.tagNames, "n");
	appendValue(merging.tagValues, to_string(nvalue));
	merging.score = (float)nvalue / (float) ((merging.endPos+1)-merging.beginPos);
	writeRecord(gffOut, merging);

	return 0;
}
