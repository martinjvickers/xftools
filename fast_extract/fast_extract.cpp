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
	seqan::ArgumentParser parser("fast_extract");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("f", "input-filter-file", "Path to the input filter file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-filter-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");


	addOption(parser, seqan::ArgParseOption("e", "exclude", "Remove reads in fasta file that have IDs in the text file."));

	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "Jan 2017");
	addUsageLine(parser, "-i input.fastq -f filter.txt -o output.fastq [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Extract or exclude reads based on an input list text file of IDs.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFastFileName, parser, "input-file");
	getOptionValue(options.inputFilterFileName, parser, "input-filter-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	options.exclude = isSet(parser, "exclude");

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

	//map<CharString, int> map;
	unordered_map<string, int> map;

	//read in the list of IDs we wish to remove or keep and put into a hash table
	ifstream myfile;
	string line;
	myfile.open (toCString(options.inputFilterFileName));
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			//cout << line << '\n';
			map[line] = map[line] + 1;
		}
		myfile.close();
	}
	myfile.close();

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

		//this is annoying, the sam file we have in our text file cuts off anything after a space. So the IDs don't match.
		string myid = toCString(id);		
		string delimiter = " ";
		string token = myid.substr(0, myid.find(delimiter));

		//If the ID from the text file exists in fasta file, then print it
		if((map.count(token) > 0) && (options.exclude == false))
		{
			writeRecord(seqFileOut, id, seq, qual);
		}

		if((map.count(token) == 0) && (options.exclude == true))
                {
			writeRecord(seqFileOut, id, seq, qual);
                }
	}

	return 0;
}
