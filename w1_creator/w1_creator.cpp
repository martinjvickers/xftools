#include "common.h"

struct ModifyStringOptions
{
        CharString inputFileName;
	int window_length;
};

struct WindowValues
{
        int c;
        int t;
	int n;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("w1/(soon w50)_creator");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.6");
	setDate(parser, "November 2017");
	addUsageLine(parser, "-i CX_report.txt [\\fIOPTIONS\\fP] ");
	addOption(parser, seqan::ArgParseOption("l", "window-length", "Size of window",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "window-length", "50");

	addDescription(parser, "Create a w1 (and soon w50) file from a CX report.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.window_length, parser, "window-length");

	return seqan::ArgumentParser::PARSE_OK;

}

int roundUp(int numToRound, int multiple)
{
	if (multiple == 0)
		return numToRound;

	int remainder = numToRound % multiple;
	if (remainder == 0)
		return numToRound;

	return numToRound + multiple - remainder;
}



/*
Aim: Read in a bismark CX report and output a w1 file

Current progress: It compiles!

DATA FORMAT:
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
1	6	+	0	0	CHH	CCA
1	7	+	0	0	CHH	CAA
1	11	-	0	0	CHH	CTT
1	12	+	0	0	CHH	CCA
1	13	+	0	0	CHH	CAA
1	17	-	0	0	CHH	CTT
1	21	-	0	0	CHH	CAA
1	22	+	0	0	CHH	CAT
1	25	+	0	0	CHG	CAG
1	27	-	0	0	CHG	CTG
*/
int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	//open input stream
	ifstream data(toCString(options.inputFileName));
	string line;
    
	//open output files
	GffFileOut CGgffOut(toCString("CG.w1.gff"));
	GffFileOut CHGgffOut(toCString("CHG.w1.gff"));
	GffFileOut CHHgffOut(toCString("CHH.w1.gff"));

	while(getline(data, line))
	{
		stringstream lineStream(line);
		string cell;
		string row[7];
		int c = 0;

		GffRecord gffrecord;

		//put each cell of a row into the row array
		while(getline(lineStream, cell, '\t'))
		{
			row[c] = cell;
			c++;
		}

		//now process

		//only do something if it's methylated
		if((stoi(row[3]) == 0) && (stoi(row[4]) == 0)){

		} else {

			//calculate the score
			int count = stoi(row[3]) + stoi(row[4]);
			double score = stoi(row[3]) / count;

			//start making our GFF record
			gffrecord.ref = row[0];
			gffrecord.source = ".";
			gffrecord.type = row[5];
			//this is -1 because of the way SeqAn stores start positions from 0
			gffrecord.beginPos = stoi(row[1])-1;
			gffrecord.endPos = stoi(row[1]);

			//converting our strand to a char so that it can be added to the record
			char strand[row[2].length()];
			row[2].copy(strand,row[2].length());
			gffrecord.strand = strand[0];

			gffrecord.score = score;
			appendValue(gffrecord.tagNames, "c");
			appendValue(gffrecord.tagValues, row[3]);
			appendValue(gffrecord.tagNames, "t");
			appendValue(gffrecord.tagValues, row[4]);

			//write out record to correct file
			if(row[5] == "CG")
			{
				writeRecord(CGgffOut, gffrecord);
			} 
			else if(row[5] == "CHG")
			{
				writeRecord(CHGgffOut, gffrecord);
			}
			else if(row[5] == "CHH")
			{
				writeRecord(CHHgffOut, gffrecord);
			} else {
				cout << "Huh, what other contexts are there?" << endl;
			}

			//clear the gffrecord
			clear(gffrecord.tagNames);
			clear(gffrecord.tagValues);


		}	

	}


	return 0;
}
