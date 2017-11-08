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
#include <boost/algorithm/string.hpp>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString annotationFileName;
	CharString outputFileName;
	int tss;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("pausing_index");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("a", "annotation-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "annotation-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	setShortDescription(parser, "XFTOOLS");
	setVersion(parser, "0.0.1");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i input.w1.gff -a reference.gff -tss 200 [\\fIOPTIONS\\fP] ");
	addOption(parser, seqan::ArgParseOption("tss", "tss-size", "Size of the TSS you wish to extract",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "tss-size", "200");

	addDescription(parser, "Calculates the pausing index of genes");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.annotationFileName, parser, "annotation-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	getOptionValue(options.tss, parser, "tss-size");

	return seqan::ArgumentParser::PARSE_OK;
}

struct geneElement
{
	CharString name;
	CharString chr;
	int start;
	int end;
	int TSS_end;
	int sum_exon;
	bool TSS_initial;
	int TSS_remaining;
	char strand;
	vector<GffRecord> exons;
	double score_for_tss;
	double score_for_gene;
};

// put the annotation into a datastructure
/*
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	5174	5326	.	+	.	Parent=AT1G01010.1
*/
void process_annotation(GffFileIn &gffAnnotationIn, map< CharString, map<CharString, geneElement>> &exons, ModifyStringOptions options)
{
	GffRecord record;
	while (!atEnd(gffAnnotationIn))
	{
		readRecord(record, gffAnnotationIn);
		CharString label = record.tagValues[0];

		// caps for the ref
		string curr_ref = toCString(record.ref);
		for (string::size_type i = 0; i < curr_ref.length(); ++i)
			boost::to_upper(curr_ref);

		if ( exons[curr_ref].find(label) == exons[curr_ref].end() )
		{
			exons[curr_ref][label].name = label;
			exons[curr_ref][label].start = record.beginPos;
			exons[curr_ref][label].end = record.endPos;

			// Set chromosome
			exons[curr_ref][label].chr = record.ref;

			// set flag to say this is the first time we've
			// encountered this label
			exons[curr_ref][label].TSS_initial = true;
			exons[curr_ref][label].TSS_remaining = options.tss;

			//set scores to zero
			exons[curr_ref][label].score_for_tss = 0.0;
			exons[curr_ref][label].score_for_gene = 0.0;
		} 
		else 
		{
			// Set start
			if(record.beginPos+1 < exons[curr_ref][label].start)
				exons[curr_ref][label].start = record.beginPos+1;

			// Set end
			if(record.endPos > exons[curr_ref][label].end)
				exons[curr_ref][label].end = record.endPos;
		}

		// Push back the record
		exons[curr_ref][label].exons.push_back(record);
		exons[curr_ref][label].strand = record.strand;
		exons[curr_ref][label].sum_exon = exons[curr_ref][label].sum_exon + (record.endPos - record.beginPos);

		if(exons[curr_ref][label].TSS_remaining > 0)
		{
			// if the bases remaining is less than the 
			if(exons[curr_ref][label].strand == '+')
			{
				// this decides if our current exon contains the end of our TSS
				if((record.beginPos + exons[curr_ref][label].TSS_remaining) < record.endPos)
				{
					exons[curr_ref][label].TSS_end = record.beginPos + exons[curr_ref][label].TSS_remaining;
					exons[curr_ref][label].TSS_remaining = 0;
				}
				else // if it's not within this one, simply update the remaining
				{
					exons[curr_ref][label].TSS_remaining = exons[curr_ref][label].TSS_remaining - (record.endPos-record.beginPos);
				}
			}
			else if (exons[curr_ref][label].strand == '-')
			{
				if((record.endPos - exons[curr_ref][label].TSS_remaining) > record.beginPos)
				{
					exons[curr_ref][label].TSS_end = record.endPos - exons[curr_ref][label].TSS_remaining;
					exons[curr_ref][label].TSS_remaining = 0;
				}
				else
				{
					exons[curr_ref][label].TSS_remaining = exons[curr_ref][label].TSS_remaining - (record.endPos-record.beginPos);
				}
			}
		}
	}
}

/*
        CharString name;
        CharString chr;
        int start;
        int end;
        int TSS_end;
        bool TSS_initial;
        int TSS_remaining;
        char strand;
        vector<GffRecord> exons;
        double score_for_tss;
        double score_for_gene;
*/
void print_out(map< CharString, map <CharString, geneElement>> &exons, double sum_score, ModifyStringOptions options)
{
	GffFileOut gffOutFile;
	if(!open(gffOutFile, toCString(options.outputFileName)))
		std::cerr << "ERROR: Could not open output.gff" " for reading.\n";

	for(auto chr : exons)
	{
		for(auto gene : chr.second)
		{
			double tss_rpkm = ((gene.second.score_for_tss*1000000000)/sum_score) / options.tss;
			double gene_rpkm = ((gene.second.score_for_gene*1000000000)/sum_score) / (gene.second.sum_exon - options.tss);

			GffRecord record;
			record.ref = gene.second.chr;
			record.source = "pausing_index";
			record.type = "gene";
			record.beginPos = gene.second.start;
			record.endPos = gene.second.end;
			record.strand = gene.second.strand;
			if(gene_rpkm == 0.0 || ((tss_rpkm / gene_rpkm) == 0))
				record.score = GffRecord::INVALID_SCORE();
			else
				record.score = tss_rpkm / gene_rpkm;

			appendValue(record.tagNames, "Parent");
			appendValue(record.tagValues, gene.first);
			appendValue(record.tagNames, "SumExonLenth");
                        appendValue(record.tagValues, to_string(gene.second.sum_exon));
			appendValue(record.tagNames, "SumTSS");
			appendValue(record.tagValues, to_string(gene.second.score_for_tss));
			appendValue(record.tagNames, "SumGene");
                        appendValue(record.tagValues, to_string(gene.second.score_for_gene));
			appendValue(record.tagNames, "rpkmTSS");
                        appendValue(record.tagValues, to_string(tss_rpkm));
                        appendValue(record.tagNames, "rpkmGene");
                        appendValue(record.tagValues, to_string(gene_rpkm));

			writeRecord(gffOutFile, record);
		}
	}

	close(gffOutFile);
}

void process_input(GffFileIn &gffFileIn, map< CharString, map<int, GffRecord>> &w1_file, double &sum_score)
{
	GffRecord record;
	while (!atEnd(gffFileIn))
        {
		readRecord(record, gffFileIn);

		// make upper case
		string ref = toCString(record.ref);
		for (string::size_type i = 0; i < ref.length(); ++i)
			boost::to_upper(ref);

		w1_file[ref][record.beginPos+1] = record;
		
		// sum the score for RPKM later
		sum_score += record.score;		
	}
}

void calculate_counts(map< CharString, map <CharString, geneElement>> &exons, map< CharString, map<int, GffRecord>> &w1_file)
{
	for(auto &chr : exons) // chr now is a map of elements
	{
		for(auto &gene : chr.second) // gene is now a particular geneElement
		{
			for(auto &exon : gene.second.exons)
			{
				// make the ref upper case
				string ref = toCString(exon.ref);
				for (string::size_type i = 0; i < ref.length(); ++i)
					boost::to_upper(ref);
				
				// let's find the finding
				for(int i = exon.beginPos+1; i <= exon.endPos; i++)
				{
					if ( w1_file[ref].find(i) != w1_file[ref].end() )
                			{
						if(gene.second.strand == '+')
						{
                                			if(i < gene.second.TSS_end)
							{
								gene.second.score_for_tss = gene.second.score_for_tss + w1_file[ref][i].score;
								
							}
							else
							{
								gene.second.score_for_gene = gene.second.score_for_gene + w1_file[ref][i].score;
							}
						}
						else if(gene.second.strand == '-')
						{
							if(i > gene.second.TSS_end)
							{
								gene.second.score_for_tss = gene.second.score_for_tss + w1_file[ref][i].score;
							}
							else
							{
								gene.second.score_for_gene = gene.second.score_for_gene + w1_file[ref][i].score;
							}
						}
					}
				}
			}
		}
	}
}

int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	GffFileIn gffFileIn, gffAnnotationIn;

	if (!open(gffFileIn, toCString(options.inputFileName)))
		return 1;
	if (!open(gffAnnotationIn, toCString(options.annotationFileName)))
		return 1;

	// Let's make a map of the input file
	double sum_score = 0.0;
	map< CharString, map<int, GffRecord>> w1_file;
	process_input(gffFileIn, w1_file, sum_score);
	close(gffFileIn);

	// now we put the exons into RAM and get our TSS/gene info
	map< CharString, map <CharString, geneElement>> exons;
	process_annotation(gffAnnotationIn, exons, options);
	close(gffAnnotationIn);

	// do the calculations
	calculate_counts(exons, w1_file);

	// write out
	print_out(exons, sum_score, options);

	cout.precision(17);
	cout << "The sum of the score (column 6) from " << options.inputFileName << " was " << sum_score << endl;

	return 0;
}
