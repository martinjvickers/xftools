#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
#include <unordered_map>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
	CharString outputFileName;
	CharString referenceFileName;
	CharString insertFileName;
	CharString contig;
	int position;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("fasta_insert");
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the reference fasta file that your other file will be inserted into. This might be a reference genome.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "reference-file");
	addOption(parser, seqan::ArgParseOption("i", "insert-file", "Path to the fasta file you wish to insert into your reference file.", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "insert-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	
	addOption(parser, seqan::ArgParseOption("c", "contig", "The contig in the reference file you wish to insert your fasta into. This is case sensitive.", seqan::ArgParseArgument::STRING, "TEXT"));
	setRequired(parser, "contig");

	addOption(parser, seqan::ArgParseOption("p", "position", "This is the zero-based position in the reference file you wish to insert your fasta into. ", seqan::ArgParseArgument::INTEGER, "INT"));
	setRequired(parser, "position");

	setShortDescription(parser, "XFTOOLS");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-r reference.fasta -i insert.fasta -o result.fasta [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Allows you to insert a fasta file into another fasta file at a given position.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.referenceFileName, parser, "reference-file");
	getOptionValue(options.insertFileName, parser, "insert-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	getOptionValue(options.contig, parser, "contig");
	getOptionValue(options.position, parser, "position");

	return seqan::ArgumentParser::PARSE_OK;
}

IupacString insert_seq(IupacString current, IupacString to_insert, int position)
{
	IupacString modified = current;
	insert(modified, position, to_insert);
	return modified;
}

/*
Current progress: It compiles!
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//read in the reference file
	CharString id;
	IupacString seq;

	// reference file
	SeqFileIn seqFileIn;
	if (!open(seqFileIn, toCString(options.referenceFileName)))
	{
		std::cerr << "ERROR: Could not open the reference file: " << options.referenceFileName << "\n";
		return 1;
	}

	// output file
	SeqFileOut seqFileOut;
	if (!open(seqFileOut, toCString(options.outputFileName)))
	{
		std::cerr << "ERROR: Could not open the output file: " << options.outputFileName << "\n";
		return 1;
	}

	while (!atEnd(seqFileIn))
	{
		readRecord(id, seq, seqFileIn);
		
		if(id == options.contig)
		{
			//open up the insert!
			SeqFileIn insertSeqFileIn(toCString(options.insertFileName));
			StringSet<CharString> ids;
			StringSet<IupacString> seqs;
			readRecords(ids, seqs, insertSeqFileIn);

			if(length(seqs) > 1)
			{
				cerr << "ERROR: There is more than one contig in the insert-file: " << options.insertFileName << endl;
				cerr << "INFO : You can only insert a single contig at a time." << endl;
				return 1;
			}

			cout << id << " was length: " << length(seq) << endl;
			IupacString modified_seq = insert_seq(seq, seqs[0], options.position); // should always be seqs[0] as is length(seqs) > 1 it would have thrown an error.
			cout << id << " now length: " << length(modified_seq) << endl;
			writeRecord(seqFileOut, id, modified_seq);
		} else {
			writeRecord(seqFileOut, id, seq);
			cout << id << " length: " << length(seq) << endl;
		}
			//insert();
	}

	//only allow the options.insertFileName to have one contig



	//check that the specific chromosome and position exist
		//maybe have a lazy version? until the first space?

	//go through reference until contig is found
		//insert
	

	return 0;
}
