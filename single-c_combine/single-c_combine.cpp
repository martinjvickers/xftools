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
	vector<CharString> inputFileNameList;
	CharString outputFileNameList;
	CharString label;
};

struct WindowValues
{
	std::map<string, string> tagmap;
	char strand;
};

typedef std::map<pair<int,int>, WindowValues> resultsmap;
typedef std::map<string, resultsmap> completemap;

resultsmap newresults;
GffFileOut gffFileOut;

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("single-c_combine");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	addOption(parser, seqan::ArgParseOption("l", "label", "Column 3 GFF output label. Useful if using SignalMap", seqan::ArgParseArgument::STRING, "TEXT"));
	setRequired(parser, "label");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "October 2017");
	addUsageLine(parser, "-i input_w1.gff [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Merge multiple single-c gff files.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	for(int i = 0; i < getOptionValueCount(parser, "input-file"); i++)
	{
		CharString tmpVal;
		getOptionValue(tmpVal, parser, "input-file", i);
		options.inputFileNameList.push_back(tmpVal);
	}

	getOptionValue(options.label, parser, "label");
	getOptionValue(options.outputFileNameList, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;

}

bool is_digits(const std::string &str)
{
	return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}

std::map<string, string> addtags(StringSet<CharString> tagNames, StringSet<CharString> tagValues, std::map<string, string> _tagmap)
{
        for(int i = 0; i < length(tagNames); i++)
        {
                //for each tag, find it in the _tagmap
                std::map<string, string>::iterator it;
                it = _tagmap.find(toCString(tagNames[i]));

                //if it exists, do something
                if(it != _tagmap.end())
                {
			int toadd, current;
                        current = 0;
                        string currentstr = _tagmap[toCString(tagNames[i])];
                        try {
				if(!is_digits(toCString(tagValues[i])))
				{
					cerr << "ERROR: " << tagValues[i] << " isn't just a number. This will give you inconsistant results." << endl;
				}
                                toadd = boost::lexical_cast<int>(toCString(tagValues[i]));
                                current = boost::lexical_cast<int>(_tagmap[toCString(tagNames[i])]);
                        }
                        catch (bad_lexical_cast &) {
                                //don't actually do anything then
                        }

                        _tagmap[toCString(tagNames[i])] = std::to_string(toadd + current);

                }
                else //just add it
                {
                        _tagmap[toCString(tagNames[i])] = toCString(tagValues[i]);
                }
        }

        return _tagmap;
}

//given a filename and a chromosome, add those results to the newresults map
int add_new_data(CharString filename, CharString chromosome)
{
	GffFileIn gffFileIn;
	if (!open(gffFileIn, toCString(filename)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

	//read in the GFF file
	GffRecord record;
	while (!atEnd(gffFileIn))
	{
		try
		{
			readRecord(record, gffFileIn);
			pair <int, int> feature_position;
			feature_position = make_pair(record.beginPos, record.endPos);

			if(chromosome == record.ref)
                        {
                                WindowValues tag = newresults[feature_position];
                                tag.tagmap = addtags(record.tagNames, record.tagValues, tag.tagmap);
				tag.strand = record.strand;
                                newresults[feature_position] = tag;
			}
		}
		catch(Exception const & e)
		{
			cerr << "ERROR: " << e.what() << endl;
		}
	}

	close(gffFileIn);

	return 0;
}

int write_data(CharString ref, ModifyStringOptions options)
{
        CharString source = ".";
        CharString type = ".";
        float score = 0;
        char strand = '.';
        char phase = '.';

	for(auto& location : newresults)
        {
		GffRecord outputrecord;
		outputrecord.ref = ref;
		outputrecord.source = source;
		outputrecord.type = options.label;
		outputrecord.beginPos = location.first.first;
		outputrecord.endPos = location.first.second;
		outputrecord.score = GffRecord::INVALID_SCORE();
		outputrecord.strand = location.second.strand;

		int c,t;

		for(auto& tag : location.second.tagmap)
		{
			appendValue(outputrecord.tagNames, tag.first);
			appendValue(outputrecord.tagValues, tag.second);
			if(tag.first == "c")
				c = atoi(toCString(tag.second));
			if(tag.first == "t")
				t = atoi(toCString(tag.second));
		}
		score = (float)c / (float)(c + t);
		outputrecord.score = score;

		writeRecord(gffFileOut, outputrecord);
	}

	return 0;
}

int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;


        //open file for writing out
        if (!open(gffFileOut, toCString(options.outputFileNameList)))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

	//open first file
	GffFileIn gffFirstFileIn;
	if (!open(gffFirstFileIn, toCString(options.inputFileNameList[0])))
        {
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

	int count = 0;
	CharString lastref;
	GffRecord firstrecord;

	try
	{
		while (!atEnd(gffFirstFileIn))
		{

			//read record
			readRecord(firstrecord, gffFirstFileIn);
			pair <int, int> feature_position;
			feature_position = make_pair(firstrecord.beginPos, firstrecord.endPos);

			//populate initial value
			if(count == 0)
			{
				lastref = firstrecord.ref;
			}

			//populate we're working on the same chromosome as the last
			if(lastref == firstrecord.ref)
			{
				WindowValues tag = newresults[feature_position];
				tag.tagmap = addtags(firstrecord.tagNames, firstrecord.tagValues, tag.tagmap);
				tag.strand = firstrecord.strand;
				newresults[feature_position] = tag;
//				cout << firstrecord.ref << " " << firstrecord.beginPos << " " << firstrecord.endPos << " " << tag.strand << endl;
			} else {

				//now do the same for every other file
				for(int i = 1; i < options.inputFileNameList.size(); i++)
				{
					add_new_data(options.inputFileNameList[i], lastref);
				}
				write_data(lastref, options);
				lastref = firstrecord.ref; // now update
				clear(newresults);

				//now add the value we are sat at
				WindowValues tag = newresults[feature_position];
				tag.strand = firstrecord.strand;
				tag.tagmap = addtags(firstrecord.tagNames, firstrecord.tagValues, tag.tagmap);
				newresults[feature_position] = tag;
			}
			count++;
		}
	}
	catch(Exception const & e)
	{
		cerr << "ERROR: " << e.what() << endl;
	}

	for(int i = 1; i < options.inputFileNameList.size(); i++)
	{
		add_new_data(options.inputFileNameList[i], lastref);
	}
	write_data(lastref, options); //write out final chromosome

        close(gffFileOut);
	close(gffFirstFileIn);

	return 0;
}
