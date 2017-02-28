#include "common.h"

struct ModifyStringOptions
{
        CharString inputFileName;
	int cpu_cores;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("Mapper");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("p", "cpu-cores", "CPU Cores.", seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "cpu-cores", "1");
	setShortDescription(parser, "Methylation Tools");
	setVersion(parser, "0.0.1");
	setDate(parser, "Feb 2017");
	addUsageLine(parser, "-i sequence.fastq -p 4 [\\fIOPTIONS\\fP] ");
	addDescription(parser, "Maps BS reads to reference");
	
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.cpu_cores, parser, "cpu-cores");

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

int mapping(ModifyStringOptions options)
{
	string const bowtie_cmd1 = string("~/bin/bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.C2T.fa -U ref.temp.C2T.fastq.gz | samtools view -bS - > test_BS_C2T_C_TO_T.bam");
	cout << "Running: " << bowtie_cmd1 << endl;
	system(bowtie_cmd1.c_str());
	string const bowtie_cmd2 = string("~/bin/bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.G2A.fa -U ref.temp.C2T.fastq.gz | samtools view -bS - > test_BS_G2A_C_TO_T.bam");
	cout << "Running: " << bowtie_cmd2 << endl;
	system(bowtie_cmd2.c_str());
	string const bowtie_cmd3 = string("~/bin/bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.C2T.fa -U ref.temp.G2A.fastq.gz | samtools view -bS - > test_BS_C2T_G_TO_A.bam");
	cout << "Running: " << bowtie_cmd3 << endl;
	system(bowtie_cmd3.c_str());
	string const bowtie_cmd4 = string("~/bin/bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.G2A.fa -U ref.temp.G2A.fastq.gz | samtools view -bS - > test_BS_G2A_G_TO_A.bam");
	cout << "Running: " << bowtie_cmd4 << endl;
	system(bowtie_cmd4.c_str());

	return 0;
}

int create_reads(ModifyStringOptions options)
{

	//read in all fastqfile
	StringSet<CharString> id;
	StringSet<CharString> qual;
	StringSet<IupacString> seq;
	SeqFileIn seqFileIn;
	if(!open(seqFileIn, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}
	
	//open two fileout handlers for C2T and G2A
	SeqFileOut seqFileOutC2T;
        if (!open(seqFileOutC2T, toCString("ref.temp.C2T.fastq.gz")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

        SeqFileOut seqFileOutG2A;
        if (!open(seqFileOutG2A, toCString("ref.temp.G2A.fastq.gz")))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }
	while(!atEnd(seqFileIn))
	{
        	try
        	{
        	        readRecords(id, seq, qual, seqFileIn, 100000);
        	}
        	catch (Exception const & e)
        	{
        	        std::cout << "ERROR: " << e.what() << std::endl;
        	        return 1;
        	}
		
		for(int i = 0; i < length(id); i++){
			//convert C's to T's
			typedef ModifiedString<IupacString, ModView<ConvertCT> > TModCT;
			TModCT modCT(seq[i]);	
			writeRecord(seqFileOutC2T, id[i], modCT, qual[i]);
	
			//convert G's to A's
			typedef ModifiedString<IupacString, ModView<ConvertGA> > TModGA;
			TModGA modGA(seq[i]);
			writeRecord(seqFileOutG2A, id[i], modGA, qual[i]);
		}
		clear(qual);
		clear(id);
		clear(seq);
		cout << " Size " << length(id) << endl;
	}

	close(seqFileOutG2A);
	close(seqFileOutC2T);

	return 0;
}

/*
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//read in file
	StringSet<CharString> ids;
        StringSet<Dna5String> seqs;
        SeqFileIn seqFileIn;

	if (!open(seqFileIn, toCString(options.inputFileName)))
        {
                std::cerr << toCString(options.inputFileName) << std::endl;
                return 1;
        }

	readRecords(ids, seqs, seqFileIn);

	for(int i = 0; i < length(ids); i++)
	{
		//this is what bismark does
		//open ($fh->{fh},"$path_to_bowtie $bt2_options -x $fh->{bisulfiteIndex} -U $temp_dir$fh->{inputfile} |") or die "Can't open pipe to bowtie: $!";

		string const bowtie_cmd1 = string("~/bin/bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.C2T.fa -U ref.temp.C2T.fastq.gz | samtools view -bS - > test_BS_C2T_C_TO_T.bam");
	        cout << "Running: " << bowtie_cmd1 << endl;
        	system(bowtie_cmd1.c_str());
	}

	close(seqFileIn);

	return 0;
}
