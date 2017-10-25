#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <string>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString outputFileName;
	CharString label;
	int window_size;
	bool percentage = false;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        seqan::ArgumentParser parser("cytosine_abundance");
        setShortDescription(parser, "XFTOOLS");
        setVersion(parser, "0.0.1");
        setDate(parser, "October 2017");
        addUsageLine(parser, "-i input.fa -o output.gff [\\fIOPTIONS\\fP] ");
	addDescription(parser, "The purpose of this program is to calculate the cytosine abundance and cytosine context's in a reference genome.");

	// accept an input file
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	addOption(parser, seqan::ArgParseOption("l", "label", "Column 3 GFF output label. Useful if using SignalMap as GFFs with the same label will be merged.", seqan::ArgParseArgument::STRING, "TEXT"));
	addOption(parser, seqan::ArgParseOption("s", "window-size", "Size of window",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "window-size", "50");
	addOption(parser, seqan::ArgParseOption("p", "percentage", "Rather than calculating the number of C's, calculate the percentage of C's in the window."));
	
	setDefaultValue(parser, "label", "window");

        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        // If parsing was not successful then exit with code 1 if there were errors.
        // Otherwise, exit with code 0 (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
                return res;

	// parse the inputs into the options struct
        getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	getOptionValue(options.label, parser, "label");
	getOptionValue(options.window_size, parser, "window-size");
	options.percentage = isSet(parser, "percentage");

        return seqan::ArgumentParser::PARSE_OK;
}

unsigned int count_Cs (Dna5String &seq)
{
	unsigned int counter = 0;

	for(unsigned int i = 0; i < length(seq); i++)
		if(seq[i] == 'C')
			counter++;
	return counter;
}

void calculate_contexts(Dna5String &seq, map<string,int> &context_map)
{
	for(unsigned i = 2; i < length(seq)-2; i++)
	{
		Infix<String<Dna5> >::Type inf = infix(seq, i, i+3);
		CharString trin = inf;
		if(inf[0] == 'C')
		{
			if(inf[1] == 'G')
			{
				context_map["CG"]++;		//CG
				context_map[toCString(trin)]++;
				//cout << "CG\t" << inf << "\t" << trin<< endl;
			}
			else
			{
				if(inf[2] == 'G')
				{
					context_map["CHG"]++;	//CHG
					context_map[toCString(trin)]++;
					//cout << "CHG\t" << inf << "\t" << trin<< endl;
				}
				else
				{
					context_map["CHH"]++; 	//CHH
					context_map[toCString(trin)]++;
					//cout << "CHH\t" << inf << "\t" << trin<< endl;
				}
			}
		}
	}
}

void print_map(map<string,int> &context_map)
{
	for (auto i : context_map)
		cout << i.first << "=" << i.second << "\t";
	cout << endl;
}

void new_context_map(map<string,int> &context_map)
{
	context_map.insert( make_pair("CG", 0));
	context_map.insert( make_pair("CGA", 0));
	context_map.insert( make_pair("CGC", 0));
	context_map.insert( make_pair("CGG", 0));
	context_map.insert( make_pair("CGT", 0));
	context_map.insert( make_pair("CGN", 0));
	context_map.insert( make_pair("CHG", 0));
	context_map.insert( make_pair("CAG", 0));
	context_map.insert( make_pair("CCG", 0));
	context_map.insert( make_pair("CTG", 0));
	context_map.insert( make_pair("CNG", 0));
	context_map.insert( make_pair("CHH", 0));
	context_map.insert( make_pair("CAA", 0));
	context_map.insert( make_pair("CAC", 0));
	context_map.insert( make_pair("CAT", 0));
	context_map.insert( make_pair("CCA", 0));
	context_map.insert( make_pair("CCC", 0));
	context_map.insert( make_pair("CCT", 0));
	context_map.insert( make_pair("CTA", 0));
	context_map.insert( make_pair("CTC", 0));
	context_map.insert( make_pair("CTT", 0));
	context_map.insert( make_pair("CNN", 0));
}

void write_record(GffFileOut &gffOutFile, FaiIndex &faiIndex, unsigned int &contig, unsigned int &beginPos, int endPos, unsigned int &C_count, map<string,int> &context_map, ModifyStringOptions &options)
{
	GffRecord record;
	record.ref = sequenceName(faiIndex, contig);
	record.source = "xftools";
	record.type = options.label;
	record.beginPos = beginPos + 1; 
	record.endPos = endPos;
	record.score = C_count;
	record.strand = '.';
	record.phase = '.';
	StringSet<CharString> tagNames;
	StringSet<CharString> tagValues;
	//cout << beginPos << "\t" << endPos << endl;
	for (auto i : context_map)
	{
        	appendValue(record.tagNames, toCString(i.first));
		//CharString m = i.second;
		appendValue(record.tagValues, to_string(i.second));
		//appendValue(record.tagValues, toCString(i.second));
	}
	writeRecord(gffOutFile, record);
}

/*
	Method.
          Create index of reference
          Loop through the reference in windows
            count c's + contexts
            reverse compliment and repeat
            print to GFF
*/
int main(int argc, char const ** argv)
{
	ModifyStringOptions options;
        seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	// create output file
	GffFileOut gffOutFile;
	if(!open(gffOutFile, toCString(options.outputFileName)))
	{
		std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
		return 1;
	}

	//index the reference
	FaiIndex faiIndex;
	if (!build(faiIndex, toCString(options.inputFileName)))
    	{
        	std::cerr << "ERROR: Could not build FAI index for file " << options.inputFileName << ".\n";
        	return 1;
	}

	// I need to know the number of contigs
	Dna5String reference_seq;

	// for each contig
	for(unsigned int contig = 0; contig < numSeqs(faiIndex); contig++)
	{
		//cout << "Contig " << contig << "\t" << sequenceLength(faiIndex, contig) <<  endl;

		unsigned int beginPos = 0;
		while (beginPos+options.window_size < sequenceLength(faiIndex, contig))
		{
			map<string,int> context_map;
			new_context_map(context_map);

			readRegion(reference_seq, faiIndex, contig, beginPos, beginPos + options.window_size);

			// count the c's
			unsigned int total_C = count_Cs(reference_seq);
			reverseComplement(reference_seq);
			total_C = total_C + count_Cs(reference_seq);
			//cout << contig << "\t" << beginPos << "\t" << beginPos + options.window_size << "\t" << total_C << endl;

			// contexts
			if(beginPos > 0 && (beginPos+options.window_size+2 < sequenceLength(faiIndex, contig)) )
			{
				// get what we want with +2 on either side
				readRegion(reference_seq, faiIndex, contig, beginPos-2, beginPos + options.window_size+2);
				calculate_contexts(reference_seq, context_map);
				reverseComplement(reference_seq);
				calculate_contexts(reference_seq, context_map);
			//	print_map(context_map);
			}
			else 
			{
				if(beginPos < 2)	//	 we're at the start!
					readRegion(reference_seq, faiIndex, contig, beginPos, beginPos + options.window_size+2);
				Dna5String pad = "NN";
				pad += reference_seq;
				reference_seq = pad;
				calculate_contexts(reference_seq, context_map);
				reverseComplement(reference_seq);
				calculate_contexts(reference_seq, context_map);
			}

			write_record(gffOutFile, faiIndex, contig, beginPos, beginPos + options.window_size, total_C, context_map, options);
			/* For debugging
			cout << reference_seq << endl;
			reverseComplement(reference_seq);
			cout << reference_seq << endl;
			print_map(context_map);
			*/

			beginPos = beginPos + options.window_size;
		}
		// Do the last one too. Will work that out later but it's just this
                map<string,int> context_map;
                new_context_map(context_map);
		readRegion(reference_seq, faiIndex, contig, beginPos, sequenceLength(faiIndex, contig));
                unsigned int total_C = count_Cs(reference_seq);
                reverseComplement(reference_seq);
                total_C = total_C + count_Cs(reference_seq);
		readRegion(reference_seq, faiIndex, contig, beginPos-2, sequenceLength(faiIndex, contig));
		Dna5String pad = "NN";
		reference_seq += pad;
		calculate_contexts(reference_seq, context_map);
		reverseComplement(reference_seq);
		calculate_contexts(reference_seq, context_map);
		write_record(gffOutFile, faiIndex, contig, beginPos, (int)sequenceLength(faiIndex, contig), total_C, context_map, options);
	}

	close(gffOutFile);

	return 0;
}
