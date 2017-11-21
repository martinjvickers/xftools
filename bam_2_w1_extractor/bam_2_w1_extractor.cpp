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
	bool cigar_aware = false;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        seqan::ArgumentParser parser("bam_2_w1_extractor");
        setShortDescription(parser, "XFTOOLS");
        setVersion(parser, "0.0.6");
        setDate(parser, "November 2017");
        addUsageLine(parser, "-i input.bam -o output.gff [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Create a w50 file from a w1 file.");

	// accept an input file
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	addOption(parser, seqan::ArgParseOption("c", "cigar-aware", "This will use the full CIGAR span length. NOTE: presently this does not account for a gap due to a reference insertion. So this may not be what you want. I may change this in the future"));

        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        // If parsing was not successful then exit with code 1 if there were errors.
        // Otherwise, exit with code 0 (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
                return res;

	// parse the inputs into the options struct
        getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	options.cigar_aware = isSet(parser, "cigar-aware");

        return seqan::ArgumentParser::PARSE_OK;
}

void writeToFile(map<int,int> &counter, BamFileIn &inFile, int &rID, GffFileOut &gffOutFile)
{
	for(auto i : counter)
	{
		typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
		TBamContext const & bamContext = context(inFile);

		GffRecord record;
		record.ref = contigNames(bamContext)[rID];
		record.source = "xftools";
		record.type = "label";
		record.beginPos = i.first+1;
		record.endPos = (i.first+1);
		record.strand = '.';
		record.score = GffRecord::INVALID_SCORE();
		appendValue(record.tagNames, "n");
		appendValue(record.tagValues, to_string(i.second));
		writeRecord(gffOutFile, record);
	}
}

int main(int argc, char const ** argv)
{
	ModifyStringOptions options;
        seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	cout << "about to read bam "<< endl;

	//read in input bam file
	BamFileIn inFile;
	if (!open(inFile, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open " << options.inputFileName << " for reading.\n";
		return 1;
	}

	cout << "finished reading bam"<< endl;

	//create output file
	GffFileOut gffOutFile;
	// out(std::cout, Gff());
	if(!open(gffOutFile, toCString(options.outputFileName)))
	{
		std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
		return 1;
	}

	//<position,count>
	map<int,int> counter;

	BamHeader header;
	readHeader(header, inFile);
	BamAlignmentRecord record;
	int rID = -1;

	cout << "reading out file" << endl;

	do
	{
		readRecord(record, inFile);
		if(rID < 0)
			rID = record.rID;

		if(record.rID != rID)
		{
			writeToFile(counter, inFile, rID, gffOutFile);
			rID = record.rID;
			counter.clear();
		}

		//this is the bit right? I need to work out the length from the CIGAR score
		if(options.cigar_aware == true)
		{
			
			int spanlength;
			_getLengthInRef(spanlength, record.cigar);
			int querylength = _getQueryLength(record.cigar);
			//cout << record.qName << " " << length(record.seq) << " " << spanlength << " " << querylength << endl;
			for(int i = record.beginPos; i < record.beginPos+spanlength; i++)
			{
				counter[i]++;
			}
		}
		else 
		{
			for(int i = record.beginPos; i < (record.beginPos + length(record.seq)); i++)
			{
				counter[i]++;
			}
		}
		
	} while(!atEnd(inFile));

	writeToFile(counter, inFile, rID, gffOutFile);
	close(inFile);
	close(gffOutFile);

	return 0;
}
