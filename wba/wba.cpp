#include "common.h"

//http://marknelson.us/2011/09/03/hash-functions-for-c-unordered-containers/
struct ModifyStringOptions
{
        CharString inputFileName;
	CharString inputAnnotationFileName;
	CharString outputFileName;
	bool exclude;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("wba");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("a", "input-annotation-file", "Path to the input filter file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-annotation-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

	setShortDescription(parser, "Window by Annotation");
	setVersion(parser, "0.0.1");
	setDate(parser, "Jan 2017");
	addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Extract or exclude reads based on an input list text file of IDs.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.inputAnnotationFileName, parser, "input-annotation-file");
	getOptionValue(options.outputFileName, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;
}

struct Feature
{
	int endPos;
	char strand;
	CharString ref;
	int score;
};

/*
Aim: Given an annotation file and an input w1 file, this will add all of the w1 records that are within that annotation file.
Current progress: It compiles!
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	multimap<int,Feature> meh;
	multimap<int,Feature>::iterator it, itlow, itup;

	GffFileIn gffFileIn;
	if (!open(gffFileIn, toCString(options.inputAnnotationFileName)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

        GffRecord record;
        while (!atEnd(gffFileIn))
        {
                try
                {
                        readRecord(record, gffFileIn);

			//lets make an annotation object
			Feature feat;
			feat.endPos = record.endPos;
			feat.strand = record.strand;
			feat.ref = record.ref;
			feat.score = 0;

			//insert
			meh.insert(make_pair(record.beginPos+1, feat));
                }
                catch (Exception const & e)
                {
                        std::cerr << "ERROR: " << e.what() << std::endl;
                        return 1;
                }

        }
/*
	cout << "Input : " << meh.size() << endl;

	for(auto& a : meh)
		cout << a.first << " " << a.second.endPos << " " << a.second.strand << " " << a.second.ref << endl;
*/

	GffFileIn gffFileInput;
        if (!open(gffFileInput, toCString(options.inputFileName)))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

	GffRecord inputrecord;
        while (!atEnd(gffFileInput))
        {
		try
		{
			readRecord(inputrecord, gffFileInput);
			//get the first one
			//Returns an iterator pointing to the first element in the container whose key is considered to go after k.
			itup = meh.upper_bound(inputrecord.endPos);
		//	itup--;

			//Returns an iterator pointing to the first element in the container whose key is not considered to go before k (i.e., either it is equivalent or goes after).
			itlow = meh.lower_bound(inputrecord.beginPos+1);

			//now as i iterate, check to see if the iterator is larger than the beginPos
			//of our input
			//IN WORDS
			//I have a w1, and i have found the first key that is above the current test
			//then iterate backwards
			//until the iterator is 

			cout << "checking " << inputrecord.beginPos+1 << " is less than either of these " << (*itup).first << endl;
			for(it = itup; ((inputrecord.beginPos+1 > (*it).first)); --it)
			{
				cout << "Checking : " << inputrecord.beginPos+1 << " against " << (*it).first << " " << (*it).second.endPos << endl;
			}

/*
			for(it = itup; ((*it).first < inputrecord.beginPos+1) && ((*it).second.endPos < inputrecord.beginPos+1) ; it--)
			{
				cout << "Checking : " << inputrecord.beginPos+1 << " against " << (*it).first << " " << (*it).second.endPos << endl;
				if((inputrecord.beginPos+1 >= (*it).first) && (inputrecord.endPos <= (*it).second.endPos))
				{
					cout << inputrecord.beginPos+1 << " " << inputrecord.endPos << " is within " << (*it).first << " " << (*it).second.endPos << endl;
					(*it).second.score++;
				} else {
					cout << inputrecord.beginPos+1 << " " << inputrecord.endPos << " is NOT within " << (*it).first << " " << (*it).second.endPos << endl;
				}
			}
*/
		}
		catch (Exception const & e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl;
                        return 1;
		}
	}


        for(auto& a : meh)
                cout << a.second.ref << " " << a.first << " " << a.second.endPos << " " << a.second.strand << " " << a.second.score <<  endl;

/*
        multimap<int,int>::iterator it, itlow, itup;
        itlow = meh.lower_bound(10000);
        itup = meh.upper_bound(20000);

        for (it=itlow; it!=itup; ++it)
        {
                cout << (*it).first << " " << (*it).second <<  endl;
        }
*/
	return 0;
}
