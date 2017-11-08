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
	int tss;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("pausing_index");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("a", "annotation-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "annotation-file");
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
	bool TSS_initial;
	int TSS_remaining;
	char strand;
	vector<GffRecord> exons;
	double score_for_tss;
	double score_for_gene;
};

void calculate_values(GffRecord &record, auto &element)
{
	for(auto e : element.second.exons)
	{
		// count if it's within an exon
		if ( (record.beginPos+1 >= e.beginPos+1) && (record.endPos <= e.endPos) )
		{
			// here, we are within the exon, so we increment either
			// the tss or the gene depending on what side we're on
			if(element.second.strand == '+')
			{
				if(record.beginPos+1 < element.second.TSS_end)
				{
					element.second.score_for_tss = element.second.score_for_tss + record.score;
					cout << "Point:" << record.ref << " " << record.beginPos+1 << " with exon (" << element.second.name << ")[" << element.second.chr << "," << element.second.strand << " [" << e.beginPos+1 << "-" << e.endPos << "]" << " TSS_remain=" << element.second.TSS_remaining << " TSS score=" << element.second.score_for_tss << " " << element.second.TSS_end << " " << element.second.start << "-" << element.second.end << endl;
				}	
				else
				{
					element.second.score_for_gene = element.second.score_for_gene + record.score;
					cout << "Point:" << record.ref << " " << record.beginPos+1 << " with exon (" << element.second.name << ")[" << element.second.chr << "," << element.second.strand << " [" << e.beginPos+1 << "-" << e.endPos << "]" << " TSS_remain=" << element.second.TSS_remaining << " GENE score=" << element.second.score_for_gene << " " << element.second.TSS_end << " " << element.second.start << "-" << element.second.end << endl;
				}
			} 
			else if(element.second.strand == '-')
			{
				if(record.beginPos+1 > element.second.TSS_end)
				{
					element.second.score_for_tss = element.second.score_for_tss + record.score;
					cout << "Point:" << record.ref << " " << record.beginPos+1 << " with exon (" << element.second.name << ")[" << element.second.chr << "," << element.second.strand << " [" << e.beginPos+1 << "-" << e.endPos << "]" << " TSS_remain=" << element.second.TSS_remaining << " TSS score=" << element.second.score_for_tss << " " << element.second.TSS_end << " " << element.second.start << "-" << element.second.end << endl;
				}
				else
				{
					element.second.score_for_gene = element.second.score_for_gene + record.score;
                                        cout << "Point:" << record.ref << " " << record.beginPos+1 << " with exon (" << element.second.name << ")[" << element.second.chr << "," << element.second.strand << " [" << e.beginPos+1 << "-" << e.endPos << "]" << " TSS_remain=" << element.second.TSS_remaining << " GENE score=" << element.second.score_for_gene << " " << element.second.TSS_end << " " << element.second.start << "-" << element.second.end << endl;
				}
			}
		}
	}
}

// put the annotation into a datastructure
/*
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	exon	5174	5326	.	+	.	Parent=AT1G01010.1
*/
void process_annotation(GffFileIn &gffAnnotationIn, map< CharString, geneElement> &exons, map< CharString, map <CharString, geneElement>> &newExons)
{
	GffRecord record;
	while (!atEnd(gffAnnotationIn))
	{
		readRecord(record, gffAnnotationIn);
		CharString label = record.tagValues[0];

		newExons[record.ref][label].exons.push_back(record);

		if ( exons.find(label) == exons.end() )
		{
			exons[label].name = label;
			exons[label].start = record.beginPos;
			exons[label].end = record.endPos;

			// Set chromosome
			exons[label].chr = record.ref;

			// set flag to say this is the first time we've
			// encountered this label
			exons[label].TSS_initial = true;
			exons[label].TSS_remaining = 200;
		} 
		else 
		{
			// Set start
			if(record.beginPos+1 < exons[label].start)
				exons[label].start = record.beginPos+1;

			// Set end
			if(record.endPos > exons[label].end)
				exons[label].end = record.endPos;
		}

		// Push back the record
		exons[label].exons.push_back(record);
		exons[label].strand = record.strand;

		if(exons[label].TSS_remaining > 0)
		{
			// if the bases remaining is less than the 
			if(exons[label].strand == '+')
			{
				// this decides if our current exon contains the end of our TSS
				if((record.beginPos + exons[label].TSS_remaining) < record.endPos)
				{
					exons[label].TSS_end = record.beginPos + exons[label].TSS_remaining;
					exons[label].TSS_remaining = 0;
				}
				else // if it's not within this one, simply update the remaining
				{
					exons[label].TSS_remaining = exons[label].TSS_remaining - (record.endPos-record.beginPos);
				}
			}
			else if (exons[label].strand == '-')
			{
				if((record.endPos - exons[label].TSS_remaining) > record.beginPos)
				{
					exons[label].TSS_end = record.endPos - exons[label].TSS_remaining;
					exons[label].TSS_remaining = 0;
				}
				else
				{
					exons[label].TSS_remaining = exons[label].TSS_remaining - (record.endPos-record.beginPos);
				}
			}
		}

	}

	//// Now I've made this, I could transfer this into a better data structure 
	//// as we understand the extents of the data set

	// let's see our data structure
/*
	for(auto e : exons)
	{
		cout << e.first << "\t" << e.second.start << "\t" << e.second.end << "\t" << e.second.strand << "\t" << e.second.TSS_end << "\t" << e.second.TSS_remaining << endl;
		for(auto i : e.second.exons)
		{
			cout << "\t" << i.beginPos << "\t" << i.endPos << endl;
		}
	}
*/
}

/*
w1 file;

chr1	dzlab	window	5978	5978	1	.	.	n=1
chr1	dzlab	window	6000	6000	2	.	.	n=2
chr1	dzlab	window	6017	6017	1	.	.	n=1
chr1	dzlab	window	6020	6020	1	.	.	n=1
chr1	dzlab	window	6063	6063	3	.	.	n=3
chr1	dzlab	window	6064	6064	1	.	.	n=1
chr1	dzlab	window	6075	6075	2	.	.	n=2
chr1	dzlab	window	6192	6192	1	.	.	n=1
chr1	dzlab	window	6193	6193	2	.	.	n=2
chr1	dzlab	window	6194	6194	1	.	.	n=1


We are making a big assumption with this in that the input file and the 
annotation file are sorted correctly. I probably should do a check on 
this at the beginning.
*/
void process_data(GffFileIn &gffFileIn, map< CharString, geneElement> &exons)
{
	GffRecord record;
	int count = 0;
        while (!atEnd(gffFileIn))
        {
		readRecord(record, gffFileIn);

		//this is super lazy right now but i just want it
		//to work. I will consider a better alternative
		//when i properly rewrite WBA
		for(auto e : exons)
		{
			// get our contigs and check that they're the same
			// conver to uppercase because it's easier
			string data_ref = toCString(record.ref);
			string exon_ref = toCString(e.second.chr);
			for (string::size_type i = 0; i < data_ref.length(); ++i)
				boost::to_upper(data_ref);
			for (string::size_type i = 0; i < exon_ref.length(); ++i)
				boost::to_upper(exon_ref);

			
			if(data_ref == exon_ref)
			{
				// the data is within this gene. Now lets go through the exons
				if( (record.beginPos+1 >= e.second.start) && (record.endPos <= e.second.end))
				{
//					geneElement meh = e.second;
					calculate_values(record, e);
				}
			}
		}	
		count++;
	}
	cout << "Processed " << count << endl;
}

/*
	I'm not too happy with my implementation of this but time
	constraints mean I'm struggling to get this in a nice way
	so I'm just completing the task without considering 
	performance.
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

	GffFileIn gffFileIn, gffAnnotationIn;

	if (!open(gffFileIn, toCString(options.inputFileName)))
		return 1;
	if (!open(gffAnnotationIn, toCString(options.annotationFileName)))
		return 1;

	// get the exons into memory
	map< CharString, geneElement> exons;
	map< CharString, map <CharString, geneElement>> newExons;
	process_annotation(gffAnnotationIn, exons, newExons);

	// process through the w1 file
	process_data(gffFileIn, exons);

	return 0;
}
