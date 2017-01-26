#include "common.h"
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
	setVersion(parser, "0.0.1");
	setDate(parser, "Jan 2017");
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
	vector<GffRecord> merged;
	int count = 0;
	bool meh = true;
	//put all of our annotations into a vector and store the pairs in an interval map
	while (!atEnd(gffIn))
        {
		readRecord(record, gffIn);
		if(count==0)
		{
			merging = record;
			meh = false;
		} else if(meh == true) {
			merging = record;
			meh = false;
		} else {
			if( (record.ref == merging.ref) && ((record.beginPos - merging.endPos) < options.size) )
	               	{
		//		cout << "Extending " << merging.ref << " " << merging.beginPos << " " << merging.endPos << endl;
				merging.endPos = record.endPos;
				meh = false;
			} else {
				//this is where we have finished the merging
				//cout << merging << endl;
				clear(merging.tagNames);
				clear(merging.tagValues);
				merging.score = GffRecord::INVALID_SCORE();
				appendValue(merging.tagNames, "windowSize");
				appendValue(merging.tagValues, to_string(merging.endPos-merging.beginPos));
				//merging.score = 0;
				writeRecord(gffOut, merging);
				clear(merging);
				merging = record;
		//		meh = true;
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
	appendValue(merging.tagValues, to_string(merging.endPos-merging.beginPos));
	writeRecord(gffOut, merging);

	return 0;
}
