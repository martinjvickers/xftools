#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <string>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
        CharString inputFileName;
	CharString referenceFileName;
	int mismatches;
	double epsilon;
	bool debug = true;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        seqan::ArgumentParser parser("sRNA_mismatch_aligner");
        setShortDescription(parser, "XFTOOLS");
        setVersion(parser, "0.0.6");
        setDate(parser, "November 2017");
        addUsageLine(parser, "-i input.fa -o output.gff [\\fIOPTIONS\\fP] ");
	addDescription(parser, "The purpose of this program is to match sRNA's to the genome.");

	// accept an input file
	addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "reference-file");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the output file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("m", "num-mismatches", "Number of mismatches.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("e", "epsilon", ".", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));


        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        // If parsing was not successful then exit with code 1 if there were errors.
        // Otherwise, exit with code 0 (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
                return res;

	// parse the inputs into the options struct
        getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.referenceFileName, parser, "reference-file");
	setDefaultValue(parser, "num-mismatches", "0");
	getOptionValue(options.mismatches, parser, "num-mismatches");
	setDefaultValue(parser, "epsilon", "0.08");
	getOptionValue(options.epsilon, parser, "epsilon");

        return seqan::ArgumentParser::PARSE_OK;
}

int meh2(ModifyStringOptions options)
{
	// for the alignment not the frag store
        typedef String<char> TSequence;                 // sequence type
        typedef Align<TSequence, ArrayGaps> TAlign;     // align type

	typedef Row<TAlign>::Type TRow;
	typedef Iterator<TRow>::Type TRowIterator;


	// Create a custom fragment store which accepts Iupac characters
	struct NewStoreConfig: public FragmentStoreConfig<>
	{
		typedef String<Iupac> TReadSeq;
		typedef String<Iupac> TContigSeq;
		typedef Owner<>	TReadSeqStoreSpec;
		typedef Owner<> TReadNameStoreSpec;
	};

	BamFileOut bamFile("meh.bam");

	typedef FragmentStore<void, NewStoreConfig> TNewFragStore;
	typedef FragmentStore<void, NewStoreConfig>::TReadSeqStore TReadSeqStore;
	typedef Value<TReadSeqStore>::Type TReadSeq;
	typedef FragmentStore<void, NewStoreConfig>::TContigStore TContigStore;
	typedef Value<TContigStore>::Type TContigStoreElement;
	typedef TContigStoreElement::TContigSeq TContigSeq;
	typedef Index<TReadSeqStore, IndexQGram<UngappedShape<3>, OpenAddressing> > TIndex;
	typedef Pattern<TIndex, Swift<SwiftSemiGlobal> > TPattern;
	typedef Finder<TContigSeq, Swift<SwiftSemiGlobal> > TFinder;
	typedef FragmentStore<void, NewStoreConfig>::TAlignedReadStore TAlignedReadStore;
	typedef Value<TAlignedReadStore>::Type TAlignedRead;
	TNewFragStore fragStore;

	cout << "Loading contigs" << endl;

	if (!loadContigs(fragStore, toCString(options.referenceFileName)))
		return 1;

	cout << "Loading reads" << endl;

	if (!loadReads(fragStore, options.inputFileName))
		return 1;

	cout << "Creating index and pattern" << endl;
	TIndex index(fragStore.readSeqStore);
	TPattern pattern(index);
	
	for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
	{
		TFinder finder(fragStore.contigStore[i].seq);
		while (find(finder, pattern, options.epsilon))
		{
			//Finder<TContigSeq> verifyFinder(fragStore.contigStore[i].seq);
			//setPosition(verifyFinder, beginPosition(finder));
			//Pattern<TReadSeq, HammingSimple> verifyPattern(fragStore.readSeqStore[position(pattern).i1]);
			//unsigned readLength = length(fragStore.readSeqStore[position(pattern).i1]);

                        TAlign align;
                        resize(rows(align), 2);
                        assignSource(row(align, 0), fragStore.readSeqStore[position(pattern).i1]);
                        assignSource(row(align, 1), infix(finder));

			int gapExtend = -2;
			int gapOpen = -10;
			//int score = globalAlignment(align, Blosum62(gapExtend, gapOpen), AlignConfig<false, true, true, false>(), AffineGaps());
			int score = localAlignment(align, Score<int>(3, -3, -2, -2), AffineGaps());


//                        int score = globalAlignment(align, EditDistanceScore());

//			Score<int, Simple> scoringScheme(5, -3, -1, -5);
//			int score = globalAlignment(align, scoringScheme);

			TRow & row1 = row(align,0);
			TRow & row2 = row(align,1);
			typedef Iterator<TRow>::Type TRowIterator;
			TRowIterator itRow1 = begin(row1);
			TRowIterator itEndRow1 = end(row1);

	//		if(options.debug == true)
	//			cout << align << endl;

			int cutfront = countLeadingGaps(row1);
			setClippedBeginPosition(row1, cutfront);
			setClippedBeginPosition(row2, cutfront);

			int cutend = unclippedLength(row1) - countTrailingGaps(row1);
			setClippedEndPosition(row1, cutend);
			setClippedEndPosition(row2, cutend);

			int mismatches = 0;

			for(unsigned i = 0; i < length(row1); ++i)
			{
				if(row1[i] != row2[i])
					mismatches++;
			}

			if(mismatches <= options.mismatches)
			{
		//		if(options.debug == true)
	//				cout << "Mismatches " << mismatches << endl << "Position:" << i << "\t" << beginPosition(finder) << "\t" << endPosition(finder) << endl << fragStore.readSeqStore[position(pattern).i1] << endl << align<< endl;

	//			cout << "Position in readstore? " << position(pattern).i1 << "\t" << position(pattern).i2 << endl;

				TAlignedRead match(length(fragStore.alignedReadStore), position(pattern).i1, i, beginPosition(finder), endPosition(finder));
				appendValue(fragStore.alignedReadStore, match);
			}

		}
	}

	writeRecords(bamFile, fragStore);

	return 0;
}

/*
*/
int main(int argc, char const ** argv)
{
	ModifyStringOptions options;
        seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//index the reference
	FaiIndex faiIndex;
	if (!build(faiIndex, toCString(options.referenceFileName)))
    	{
        	std::cerr << "ERROR: Could not build FAI index for file " << options.referenceFileName << ".\n";
        	return 1;
	}

	SeqFileIn seqFileIn;
	CharString id, qual;
	Dna5String seq;
	if(!open(seqFileIn, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return 1;
	}

	meh2(options);

/*
	typedef StringSet<Dna5String> THaystacks;
	THaystacks haystacks;
	for(unsigned int contig = 0; contig < numSeqs(faiIndex); contig++)
	{
		Dna5String meh;
		readRegion(meh, faiIndex, contig, 0, sequenceLength(faiIndex, contig));
		appendValue(haystacks, meh);
	}
*/

/*
	typedef String<char> TSequence;                 // sequence type
	typedef Align<TSequence, ArrayGaps> TAlign;     // align type


	typedef Finder<Dna5String, Swift<SwiftSemiGlobal> > TSwiftFinder;
	//typedef Finder<StringSet<Dna5String>, Swift<SwiftSemiGlobal> > TSwiftFinder;
	typedef Dna5String TReadSeq_;
	typedef StringSet<TReadSeq_> TReadSet;
	typedef typename Value<TReadSet>::Type TReadSeq;
	typedef typename Value<TReadSeq>::Type TAlphabet;
	typedef Shape<TAlphabet, UngappedShape<10> > TShape;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> > TQGramIndex;
	typedef Pattern<TQGramIndex, Swift<SwiftSemiGlobal> > TSwiftPattern;

	Dna5String ref;
	readRegion(ref, faiIndex, 0, 0, sequenceLength(faiIndex, 0));

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

		clock_t begin = clock();
		int counter = 0;
		StringSet<Dna5String> readSeqs;
		appendValue(readSeqs, seq);
		TQGramIndex index(readSeqs);
		TSwiftPattern swiftPattern(index);
		//TSwiftFinder swiftFinder(ref);
		TSwiftFinder swiftFinder(haystacks);

		while (find(swiftFinder, swiftPattern, 0.1))
		{
			TAlign align;
			resize(rows(align), 2);
			assignSource(row(align, 0), seq);
			assignSource(row(align, 1), infix(swiftFinder));
			int score = globalAlignment(align, EditDistanceScore());
			//cout << "Score: " << score << endl << align << endl;
			counter++;

			if(abs(score) <= options.mismatches)
				cout << "Score: " << score << endl <<  align << endl;
		}
*/

/*
		while (find(finder, seq))
		{
			cout << beginPosition(finder) << "\t" << endPosition(finder) << endl;
		}
*/

/*
		for(unsigned int contig = 0; contig < numSeqs(faiIndex); contig++)
		{

			cout << "Chromosome " << contig << " / " << numSeqs(faiIndex) << endl;
			for(unsigned int beginPos = 0; beginPos + length(seq) < sequenceLength(faiIndex, contig); beginPos++)
			{

				IupacString refSeq;
				readRegion(refSeq, faiIndex, contig, beginPos, beginPos + length(seq));
				typedef String<Iupac> TSequence;
				typedef Align<TSequence, ArrayGaps> TAlign;
				TAlign align;
				resize(rows(align), 2);
				assignSource(row(align, 0), seq);
				assignSource(row(align, 1), refSeq);
				int score = globalAlignment(align, EditDistanceScore());

				if(abs(score) <= options.mismatches)
					cout << score << endl <<  align << endl;
			}
		}
*/
//		clock_t end = clock();
//		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//cout << "Time taken " << elapsed_secs << " " << counter << endl;
//	}

	return 0;
}
