#include "common.h"

struct ModifyStringOptions
{
        CharString queryFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("methylation_tools");
	addOption(parser, seqan::ArgParseOption("q", "query-file", "Path to the query file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "query-file");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "July 2016");
	addUsageLine(parser, "-q query.gff3 [\\fIOPTIONS\\fP] ");
	addDescription(parser, "A tool to quickly deal with Bismark downstream analysis in a version controlled stable way");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.queryFileName, parser, "query-file");

	return seqan::ArgumentParser::PARSE_OK;

}

/*
Aim: calculate introns.

Current progress: It compiles!

*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	GffFileIn gffIn(toCString(options.queryFileName));

	GffRecord record;

	bool newgene = false;
	bool newmrna = false;
	bool newprotein = false;
	int currentgenestart = 0;
	int currentgeneend = 0;

	while (!atEnd(gffIn))
    	{	
		readRecord(record, gffIn);
	
		if(record.type == "gene"){
			newgene = true;
		} else if(record.type == "mRNA"){
			newmrna = true;
		} else if(record.type == "protein"){
			newprotein = true;
		}

		if(newprotein && newmrna && newgene){
			cout << "A new one!!" << endl;
			newgene = false;
			newmrna = false;
			newprotein = false;
		}

	}

	return 0;
}
