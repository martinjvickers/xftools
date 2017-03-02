#include "common.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

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
	//make a bam/sam header
	typedef typename BamHeaderRecord::TTag	TTag;
	BamHeader header;
	
	//contig name store?
	StringSet<CharString> referenceNameStore;
	NameStoreCache<StringSet<CharString> > referenceNameStoreCache(referenceNameStore);
	BamIOContext<StringSet<CharString> > bamIOContext(referenceNameStore, referenceNameStoreCache);

	string const bowtie_cmd1 = string("bowtie2 -p ") + to_string(options.cpu_cores) + (" -k 2 -x ref.temp.C2T.fa -U ") + toCString(options.inputFileName);

	FILE *in;
    	char buff[1024];

	if(!(in = popen(bowtie_cmd1.c_str(), "r"))){
		return 1;
	} else {
		cout << "Running: " << bowtie_cmd1 << endl;
	}

	String<BamAlignmentRecord> alignments; // store our alignments

	//new thought, do i need to do all of this, just keep ploughing in the buffer
	while(fgets(buff, sizeof(buff), in)!=NULL){

		CharString input = buff; //convert buffer to CharString
		Iterator<CharString, Rooted>::Type iter = begin(input); //create iterator for the current buffer

	    	readHeader(header, bamIOContext, iter, Sam()); //read header from buffer (i assume this will just work if it's a header and if it's not it won't work?

    		while (!atEnd(iter)) //go over stuff in the iterator
		{
			resize(alignments, length(alignments) + 1);
			readRecord(back(alignments), bamIOContext, iter, Sam());
		}
	}

	//BamFileOut bamFileOut(bamIOContext, toCString("meh.bam"));
	//BamFileOut bamFileOut(bamIOContext);
	//BamFileOut bamFileOut(bamIOContext, std::cout, Bam());
	//writeHeader(bamFileOut, header);

	String<char> text;
	write(text, header, bamIOContext, Sam());

	for(auto& i : alignments)
	{
	//	writeRecord(bamFileOut, i);
		write(text, i, bamIOContext, Sam());
	}
	cout << text << endl;
/*
	for(auto& i : text)
	{
		cout << i << endl;
	}*/
	pclose(in);
	//close(bamFileOut);

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

/*
	if(create_reads(options) == 0)
	{
		cout << "Temp Files Created" << endl;
	} else {
		cerr << "ERROR: Creating temp files failed." << endl;
	}
*/
	if(mapping(options) == 0)
	{
		cout << "Mapping Files Created" << endl;
	} else {
		cerr << "ERROR: Mapping files failed." << endl;
	}

/*	//load in the bam files
	BamFileIn bamFileInC2T_C2T;
	BamFileIn bamFileInG2A_C2T;
	BamFileIn bamFileInC2T_G2A;
	BamFileIn bamFileInG2A_G2A;

	if (!open(bamFileInC2T_C2T, toCString("test_BS_C2T_C_TO_T.bam")))
	{
		std::cerr << "ERROR: Could not open " << toCString("test_BS_C2T_C_TO_T.bam") << std::endl;
		return 1;
	}

	if (!open(bamFileInG2A_C2T, toCString("test_BS_G2A_C_TO_T.bam")))
        {
                std::cerr << "ERROR: Could not open " << toCString("test_BS_G2A_C_TO_T.bam") << std::endl;
                return 1;
        }

	if (!open(bamFileInC2T_G2A, toCString("test_BS_C2T_G_TO_A.bam")))
        {
                std::cerr << "ERROR: Could not open " << toCString("test_BS_C2T_G_TO_A.bam") << std::endl;
                return 1;
        }

	if (!open(bamFileInG2A_G2A, toCString("test_BS_G2A_G_TO_A.bam")))
        {
                std::cerr << "ERROR: Could not open " << toCString("test_BS_G2A_G_TO_A.bam") << std::endl;
                return 1;
        }
*/
	/*lets assume (ass out of you and me)
	  that the fastq and the bam is in order by id. That way we can
	  read each record in the bam and then skip records in the fastq
	  until they match, thus skipping over reads that didn't map at all
	  and not taking a long time to run.
	*/
/*
	//read in our raw sequence data
	StringSet<CharString> ids;
	StringSet<Dna5String> seqs;
	SeqFileIn seqFileIn;
	if (!open(seqFileIn, toCString(options.inputFileName)))
	{
		std::cerr << toCString(options.inputFileName) << std::endl;
		return 1;
	}
	readRecords(ids, seqs, seqFileIn);

	vector<BamAlignmentRecord> read_bamrecords;

	//okay with everything loaded i need to cycle through them all together
	for(int i = 0; i < length(ids); i++)
	{
		//add all the alignments for a specific read to a vector and then process the vector
						
		

		//clear the vector

	}

	close(bamFileInC2T_C2T);
	close(bamFileInG2A_C2T);
	close(bamFileInC2T_G2A);
	close(bamFileInG2A_G2A);
	close(seqFileIn);
*/
	return 0;
}
