#include "common.h"
#include <unordered_map>

struct ModifyStringOptions
{
        CharString inputFastFileName;
	CharString inputFilterFileName;
	CharString outputFileName;
	bool exclude;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("trim");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

	setShortDescription(parser, "A very specific trim for Shengbo");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i input.fastq -o output.fastq [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Trim the first 5bp from the first reads and then the last 21bp if the read is >50bp (after the first 5bp trim)");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFastFileName, parser, "input-file");
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

	//open the fasta/fastq file	
	SeqFileIn seqFileIn;
	if (!open(seqFileIn, toCString(options.inputFastFileName)))
	{
		cerr << "ERROR: Could not open file." << endl;
		return 1;
	}

	CharString id;
	Dna5String seq;
	CharString qual;

	SeqFileOut seqFileOut;
	if (!open(seqFileOut, toCString(options.outputFileName)))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

	while (!atEnd(seqFileIn))
	{
		readRecord(id, seq, qual, seqFileIn);

		if(length(seq) > 55)
		{	
			writeRecord(seqFileOut, id, infix(seq, 5, length(seq) - 21), infix(qual, 5, length(qual) - 21));
			
		} else {
			writeRecord(seqFileOut, id, infix(seq, 5, length(seq)), infix(qual, 5, length(qual)));
		}
	}

	return 0;
}
