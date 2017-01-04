#include "common.h"

struct ModifyStringOptions
{
        CharString inputFileName;
	int window_length;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("w50_creator");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
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
    public std::unary_function<Dna, Dna>
{
    inline Dna operator()(Dna x) const
    {
        if (x == 'C') return 'T';

        return x;
    }

};

struct ConvertGA :
    public std::unary_function<Dna, Dna>
{
    inline Dna operator()(Dna x) const
    {
        if (x == 'G') return 'A';

        return x;
    }

};


/*
BBF-Mapper

Bismark-But-Faster-Mapper	(Beta)

Current progress:	It compiles and for a single-end fastq file it will create the C2T and G2A conversions.

Next steps, align and then compare.
*/
int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	cout << "hello?" << endl;

	//CT/GA output fastq files
        SeqFileOut seqFileOutG2A;
        if (!open(seqFileOutG2A, toCString(".temp.G2A.fastq")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

        SeqFileOut seqFileOutC2T;
        if (!open(seqFileOutC2T, toCString(".temp.C2T.fastq")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }


	//create the GA and CT conversions of our fasta files
	SeqFileIn seqFileIn;
	//if(!open(seqFileIn,toCString("XF-HGXF001A_S1.fastq.gz.temp.1_ambiguous_reads.fq")))
	if(!open(seqFileIn,toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
                return 1;
	}

	CharString id;
	Dna5String seq;
	CharString qual;

	cout << "jhello?" << endl;

	while(!atEnd(seqFileIn))
	{
        	try
        	{
        	        readRecord(id, seq, qual, seqFileIn);
        	}
        	catch (Exception const & e)
        	{
        	        std::cout << "ERROR: " << e.what() << std::endl;
        	        return 1;
        	}

		//convert C's to T's
		//cout << seq << endl;
		typedef ModifiedString<Dna5String, ModView<ConvertCT> > TModCT;
		TModCT modCT(seq);
		//cout << modCT << endl;
		writeRecord(seqFileOutC2T, id, modCT, qual);

		//convert G's to A's
		typedef ModifiedString<Dna5String, ModView<ConvertGA> > TModGA;
		TModGA modGA(seq);
		writeRecord(seqFileOutG2A, id, modGA, qual);

	}

	close(seqFileOutG2A);
	close(seqFileOutC2T);

	//run bowtie
	/*
	string const bowtie_refC2T_readG2A = string("bowtie2 -p 4 -x reference/Bisulfite_Genome/CT_conversion/BS_CT -U .temp.G2A.fastq | samtools view -bS - > .temp.refC2T_readG2A.bam");
	system(bowtie_refC2T_readG2A.c_str());
        string const bowtie_refC2T_readC2T = string("bowtie2 -p 4 -x reference/Bisulfite_Genome/CT_conversion/BS_CT -U .temp.C2T.fastq | samtools view -bS - > .temp.refC2T_readC2T.bam");
        system(bowtie_refC2T_readC2T.c_str());
	string const bowtie_refG2A_readG2A = string("bowtie2 -p 4 -x reference/Bisulfite_Genome/GA_conversion/BS_GA -U .temp.G2A.fastq | samtools view -bS - > .temp.refG2A_readG2A.bam");
        system(bowtie_refG2A_readG2A.c_str());
        string const bowtie_refG2A_readC2T = string("bowtie2 -p 4 -x reference/Bisulfite_Genome/GA_conversion/BS_GA -U .temp.C2T.fastq | samtools view -bS - > .temp.refG2A_readC2T.bam");
        system(bowtie_refG2A_readC2T.c_str());
	*/

	//samtools merge and sort
	/*
	string const samtools_merge = string("samtools merge -f .temp.hits.bam .temp.refC2T_readG2A.bam .temp.refC2T_readC2T.bam .temp.refG2A_readG2A.bam .temp.refG2A_readC2T.bam");
	system(samtools_merge.c_str());

	string const samtools_sort = string("samtools sort -n -o .temp.sorted.hits.bam .temp.hits.bam");
        system(samtools_sort.c_str());
	

	//now filter is all
	BamFileIn bamFile_hits;

	//files
	if (!open(bamFile_hits, toCString(".temp.sorted.hits.bam")))
    	{
        	std::cerr << "ERROR: Could not open " << std::endl;
        	return 1;
    	}

        BamHeader header;
        readHeader(header, bamFile_hits); //read the header

	vector<int> v;

	//begin reading the file
	BamAlignmentRecord current_record;
	BamAlignmentRecord working_record;
	readRecord(working_record, bamFile_hits);

	while (!atEnd(bamFile_hits))
	{

		readRecord(current_record, bamFile_hits);

		//if it's mapped
		if(!hasFlagUnmapped(current_record))
		{
			cout << current_record.qName << " mapped! Score = " << current_record.mapQ <<  endl;
		}

	}
	*/

	return 0;
}
