#include "common.h"
#include <unordered_map>
#include "key.cpp"

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString inputAnnotationFileName;
	CharString outputFileName;
	bool exclude;
};

struct WindowValues
{
        int c;
        int t;
	int n;
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

/*
Aim: Given an annotation file and an input w1 file, this will add all of the w1 records that are within that annotation file.

Current progress: It compiles!

*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	map<MyPair,int> annotation;

	GffFileIn gffAnnotationIn;
	GffRecord annotationrecord;
        if (!open(gffAnnotationIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }
	while (!atEnd(gffAnnotationIn))
        {
		try
		{
			readRecord(annotationrecord, gffAnnotationIn);
			MyPair temp(annotationrecord.ref, annotationrecord.beginPos, annotationrecord.endPos);
			annotation[temp] = annotation[temp]+1;
		} catch (int e) {
			cout << "Huh " << e <<  endl;
			break;
		}
	}

	///NOW START READING OUR INPUT FILE
	GffFileIn gffInputIn;
	GffRecord inputrecord;
	if (!open(gffInputIn, toCString(options.inputFileName)))
	{
        	std::cerr << "ERROR: Could not open example.gff" << std::endl;
	        return 1;
	}

	double sum = 0.0;

	while (!atEnd(gffInputIn))
        {

		readRecord(inputrecord, gffInputIn);
		for(auto const& p: annotation)
        	{
			//cout << p.first.getBegin()<< " " << p.first.getEnd() << " " << p.second << endl;
			if(p.first.within(inputrecord.ref, inputrecord.beginPos))
			{
				cout << inputrecord.ref << " " << inputrecord.beginPos << " is within  " << p.first.getID() << " " << p.first.getBegin()<< " " << p.first.getEnd() << " " << p.second << endl;
				break;
			}
        	} 
		sum = sum + inputrecord.score;

	}

	cout << "Sum " << sum << endl;
	return 0;
}
