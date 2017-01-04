#include "common.h"

//overload the SeqFileBuffer_ so that it uses Iupac String. In this way 
//the input file is checked against Iupac and any non-A/C/G/T is silently 
//converted into a N.
namespace seqan {
	template <typename TString, typename TSSetSpec, typename TSpec>
	struct SeqFileBuffer_<StringSet<TString, TSSetSpec>, TSpec>
	{
		typedef String<Iupac> Type;
	};
}

struct ModifyStringOptions
{
        CharString inputFileName;
	int window_length;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("w50_creator");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "July 2016");
	addUsageLine(parser, "-i sequence.fastq [\\fIOPTIONS\\fP] ");
	addOption(parser, seqan::ArgParseOption("l", "window-length", "Size of window",seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "window-length", "50");

	addDescription(parser, "Create a w50 file from a w1 file.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.window_length, parser, "window-length");

	return seqan::ArgumentParser::PARSE_OK;

}

struct ConvertCT :
    public std::unary_function<Iupac, Iupac>
{
    inline Iupac operator()(Iupac x) const
    {
        if (x == 'C') return 'T';

        return x;
    }

};

struct ConvertGA :
    public std::unary_function<Iupac, Iupac>
{
    inline Iupac operator()(Iupac x) const
    {
        if (x == 'G') return 'A';

        return x;
    }
};


/*
BBF-Mapper: build reference

Build reference for Bismark-But-Faster-Mapper	(Beta)

*/
int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	//read in all fastqfile
	CharString id;
	//Dna5String seq;
	IupacString seq;
	SeqFileIn seqFileIn(toCString(options.inputFileName));
	
	//open two fileout handlers for C2T and G2A
	SeqFileOut seqFileOutC2T;
        if (!open(seqFileOutC2T, toCString("ref.temp.C2T.fa")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

        SeqFileOut seqFileOutG2A;
        if (!open(seqFileOutG2A, toCString("ref.temp.G2A.fa")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }
	while(!atEnd(seqFileIn))
	{
        	try
        	{
        	        readRecord(id, seq, seqFileIn);
        	}
        	catch (Exception const & e)
        	{
        	        std::cout << "ERROR: " << e.what() << std::endl;
        	        return 1;
        	}
		//convert C's to T's
		typedef ModifiedString<IupacString, ModView<ConvertCT> > TModCT;
		TModCT modCT(seq);	
		writeRecord(seqFileOutC2T, id, modCT);

		//convert G's to A's
		typedef ModifiedString<IupacString, ModView<ConvertGA> > TModGA;
		TModGA modGA(seq);
		writeRecord(seqFileOutG2A, id, modGA);

	}

	close(seqFileOutG2A);
	close(seqFileOutC2T);
	//open each file and run bowtie-build on it;
	string const bowtie_cmd_G2A = string("bowtie2-build ref.temp.G2A.fa ") + toCString(options.inputFileName) + string("_G2A");
	system(bowtie_cmd_G2A.c_str());
	string const bowtie_cmd_C2T = string("bowtie2-build ref.temp.C2T.fa ") + toCString(options.inputFileName) + string("_C2T");
	system(bowtie_cmd_C2T.c_str());
	return 0;
}
