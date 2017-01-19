#include "common.h"
#include <unordered_map>
#include "key.cpp"
#include <functional>

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

typedef pair<string,pair<int,int>> Name;

namespace std {
    template <>
        class hash<Name>{
        public :
            size_t operator()(const Name &name ) const
            {
		//can i somehow return a hash 
                return hash<string>()(name.first) ^ hash<int>()(name.second.first) ^ hash<int>()(name.second.second);
            }
    };
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

	unordered_map<Name,int> ids;

	GffFileIn gffAnnotationIn;
	GffRecord annotationrecord;
        if (!open(gffAnnotationIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }
	while (!atEnd(gffAnnotationIn))
        {
		readRecord(annotationrecord, gffAnnotationIn);
		ids[Name(toCString(annotationrecord.ref), make_pair(annotationrecord.beginPos, annotationrecord.endPos))] = 0;
	}

	cout << "THere are this many annotation " << ids.size() << endl;

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
		
		cout << inputrecord.ref << " " << inputrecord.beginPos << " " << inputrecord.endPos << " " << ( std::hash<string>()(toCString(inputrecord.ref)) ^ std::hash<int>()(inputrecord.beginPos) ^ std::hash<int>()(inputrecord.endPos) ) << endl;
		ids[Name(toCString(inputrecord.ref), make_pair(inputrecord.beginPos, inputrecord.endPos))]++;
		sum = sum + inputrecord.score;
	}
	cout << "Sum " << sum << endl;

	for ( auto ii = ids.begin() ; ii != ids.end() ; ii++ )
		cout << ii->first.first << " " << ii->first.second.first << " " << ii->first.second.second << " " << ii->second << endl;

	return 0;
}
