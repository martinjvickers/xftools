#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include <string>
#include <thread>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	int tss;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("tss_label");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	setShortDescription(parser, "XFTOOLS");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i reference.gff [\\fIOPTIONS\\fP] ");
	addOption(parser, seqan::ArgParseOption("tss", "tss-size", "Size of the TSS you wish to extract",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "tss-size", "200");

	addDescription(parser, "TSS Label. Labels the TSS region based upon the start of the gene and the TSS size parameter.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.tss, parser, "tss-size");

	return seqan::ArgumentParser::PARSE_OK;

}

int get_TSS(vector<GffRecord> gene, ModifyStringOptions options)
{
	int shortfall = options.tss;
	for(auto g : gene)
        {
		//if the gene is smaller than the TSS we want, print it out and we can go to the next one
		if((g.endPos - g.beginPos) < shortfall)
                {
			cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.beginPos+1 << "\t" << g.endPos << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=TSS" << endl;
			shortfall = shortfall - ((g.endPos - g.beginPos));
		}
		else if(shortfall > 0)
		{
			if(gene[0].strand == '+')
			{
				cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.beginPos+1 << "\t" << g.beginPos+1+shortfall << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=TSS" << endl;
				cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.beginPos+1+shortfall+1 << "\t" << g.endPos << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=Body" << endl;
				shortfall = 0;
			}
			else if(gene[0].strand == '-')
			{
				cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.endPos-shortfall << "\t" << g.endPos << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=TSS" << endl;
				cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.beginPos+1 << "\t" << g.endPos-shortfall-1 << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=Body" << endl;
				shortfall = 0;
			}
		}
		else 
		{
			cout << g.ref << "\t" << g.source << "\t"<< g.type << "\t" << g.beginPos+1 << "\t" << g.endPos << "\t.\t" << g.strand << "\t.\tParent=" << g.tagValues[0] << ";Extract=Body" << endl;
		}
	}
}

/*
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

	GffFileIn gffIn(toCString(options.inputFileName));

	GffFileOut gffOut(std::cout, Gff());
	GffRecord current_record, last_record;
	vector<GffRecord> gene;
	bool first = true;
	while (!atEnd(gffIn))
	{
		readRecord(current_record, gffIn);
		
		if(first != true)
		{
			if(current_record.tagValues[0] == last_record.tagValues[0])
			{
				gene.push_back(current_record);
			}
			else
			{
				get_TSS(gene, options);
				gene.clear();
				gene.push_back(current_record);
			}
		}
		else
		{
			gene.push_back(current_record);
			first = false;
		}

		last_record = current_record;
	}

	// do the last one
	get_TSS(gene, options);


	return 0;
}
