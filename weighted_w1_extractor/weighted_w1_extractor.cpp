#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>
#include <string>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        seqan::ArgumentParser parser("bam_2_w1_extractor");
        setShortDescription(parser, "XFTOOLS");
        setVersion(parser, "0.0.1");
        setDate(parser, "June 2017");
        addUsageLine(parser, "-i input.bam -o output.gff [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Create a w50 file from a w1 file.");

	// accept an input file
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");

        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        // If parsing was not successful then exit with code 1 if there were errors.
        // Otherwise, exit with code 0 (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
                return res;

	// parse the inputs into the options struct
        getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.outputFileName, parser, "output-file");

        return seqan::ArgumentParser::PARSE_OK;
}

typedef std::map<int, double> chrmap;
typedef std::map<int, chrmap> completemap;

void writeToFile(completemap &counter, BamFileIn &inFile, GffFileOut &gffOutFile)
{
	//completemap[rID][pos][count]

	for(auto i : counter)
	{
		for(auto j : i.second)
		{
			typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
			TBamContext const & bamContext = context(inFile);

			GffRecord record;
			record.ref = contigNames(bamContext)[i.first];
			record.source = "xftools";
			record.type = "label";
			record.beginPos = j.first;
		//	record.endPos = (j.first+1);
			record.endPos = j.first;
			record.strand = '.';
			record.score = (double)j.second;
			writeRecord(gffOutFile, record);

		}
	}
}

int main(int argc, char const ** argv)
{
	ModifyStringOptions options;
        seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
	//read in input bam file
	BamFileIn inFile;
	if (!open(inFile, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open " << options.inputFileName << " for reading.\n";
		return 1;
	}

	//create output file
	GffFileOut gffOutFile;
	if(!open(gffOutFile, toCString(options.outputFileName)))
	{
		std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
		return 1;
	}

	completemap megacounter;
	vector<BamAlignmentRecord> meh;

	BamHeader header;
	readHeader(header, inFile);
	BamAlignmentRecord record;
	int rID = -1;
	CharString qName;

	do
	{
		readRecord(record, inFile);

		// if first time
		if(rID < 0)
		{
			rID = record.rID;
			qName = record.qName;
		}

		// if the current record ID does not match the last one then we know it's a new read
		if(record.qName != qName)
		{
			//go through each record in vector
			for(auto rec : meh)
			{
				//go through each pos of the record
				for(int i = rec.beginPos; i < (rec.beginPos + length(record.seq)); i++)
				{
					double tmp = megacounter[rec.rID][i];
					tmp = (double)tmp + ((double)1.0/(double)meh.size());
					megacounter[rec.rID][i] = tmp;
	//				cout << "Adding " << rID << " " << i << " " << (double)megacounter[rID][i] << " " << ((double)1.0/(double)meh.size()) << endl;
				}
			}
			meh.clear(); // clear vector
			rID = record.rID;
			qName = record.qName;
		}

		//add currect record to current vector
		meh.push_back(record);
		
	} while(!atEnd(inFile));

	writeToFile(megacounter, inFile, gffOutFile);

	close(inFile);
	close(gffOutFile);

	return 0;
}
