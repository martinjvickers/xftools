#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/store.h> /* FragmentStore */
#include <string>
#include <seqan/bam_io.h>
#include <map>
#include <vector>
//#include "tbb/tbb.h"
//#include "tbb/atomic.h"
#include <stdlib.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString referenceFileName;
   CharString inputLeftFileName;
   CharString inputRightFileName;
   CharString outputBamFileName;
   CharString ambigOutputBamFileName = NULL;
   int cpu_cores = 1;
   int memory;
   bool stranded = true;
   bool debug = false;
   bool non_directional = false;
   bool singleEnd;
   bool compress = false;
   bool uncoveredGff = false;
   bool writeCXreport = false;
};

struct CX_report
{
   unsigned int count_methylated;
   unsigned int count_unmethylated;
};

GffFileOut gffOutCG, gffOutCHG, gffOutCHH;
ofstream cxOut;

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv)
{
   ArgumentParser parser("miniature-sniffle-mapper");
   addOption(parser, ArgParseOption("i", "input-file", "Input filename", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file", "Output filename", 
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("r", "reference-file", "Reference file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "reference-file");
   addOption(parser, ArgParseOption("CX", "CX-report", 
                                    "Write Out a CX report."));
   addOption(parser, ArgParseOption("gz", "compress", "Compress temporary files\
                                    , this will increase the run time"));
   addOption(parser, ArgParseOption("c", "num-cores", "Number of Bowtie CPU \
                                    cores to use. NOTE: This is per bowtie2 \
                                    instance (x2 for directional and x4 for \
                                    non-directional).", 
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "num-cores", "1");
   addOption(parser, ArgParseOption("m", "memory", "Indicate how much memory \
                                    (in Gb) can be used to store the per base \
                                    methylation before writing out to file \
                                    NOTE: This must be large enough to fit \
                                    your largest chromosome in memory. If it \
                                    is set too small then this will \
                                    automatically be changed. 2GB is sufficient\
                                    for most genomes.",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "memory", "2");
   setShortDescription(parser, "BS Mapper");
   setVersion(parser, "0.0.5");
   setDate(parser, "September 2018");
   addUsageLine(parser, "-i mapping.bam [\\fIOPTIONS\\fP] -r reference.fa");

   addDescription(parser, "Perform the mapping.");
   seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

   if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputBamFileName, parser, "output-file");
   getOptionValue(options.referenceFileName, parser, "reference-file");

   options.writeCXreport = isSet(parser, "CX-report");
   options.singleEnd = isSet(parser, "input-file"); 
   getOptionValue(options.cpu_cores, parser, "num-cores");
   getOptionValue(options.memory, parser, "memory");
   options.compress = isSet(parser, "compress");

   return seqan::ArgumentParser::PARSE_OK;
}

// Given the fasta index, the postition and contig, return the Trinucleotide 
// context.
// TODO: Maybe make this return 0/1 as a succeed/fail
Dna5String getTrinucleotideContext(FaiIndex &faiIndex, CharString ref, unsigned int pos)
{
   unsigned idx = 0;
   if (!getIdByName(idx, faiIndex, ref))
   {
      cerr << "Error: " << ref << " not found in Reference" << endl;
   }
   Dna5String sequenceInfix;
   readRegion(sequenceInfix, faiIndex, idx, pos, pos+3);
   return sequenceInfix;
}

CharString identifyContext(Dna5String getTrinucleotideContext)
{
   if(getTrinucleotideContext[0] == 'C')
   {
      if(getTrinucleotideContext[1] == 'G')
         return "CG";
      else if(getTrinucleotideContext[2] == 'G')
         return "CHG";
      else
         return "CHH";
   }
   else
   {
      return "NNN";
   }
}

int writeGFFRecord(map<CharString, map<unsigned int, CX_report>> &testMap,
                  unsigned pos, CharString contig, bool hasRC,
                  CharString context, Dna5String infix)
{
   GffRecord record;
   record.ref = contig;
   record.source = "";
   record.type = context;

   unsigned refPos = pos;
   record.strand = '+';
   if(hasRC == true)
   {
      record.strand = '-';
   }

   record.beginPos = refPos;
   record.endPos = refPos+1;

   if( testMap[contig].find(refPos) == testMap[contig].end() )
   {
      record.score = GffRecord::INVALID_SCORE();
      appendValue(record.tagNames, "c");
      appendValue(record.tagValues, to_string(0));
      appendValue(record.tagNames, "t");
      appendValue(record.tagValues, to_string(0));
      appendValue(record.tagNames, "n");
      appendValue(record.tagValues, to_string(0));
   }
   else
   {
      unsigned c, t, n;
      c = testMap[contig][refPos].count_methylated;
      t = testMap[contig][refPos].count_unmethylated;
      n = c + t;

      if(n == 0)
         record.score = 0;
      else
         record.score = (float)c / (float)n;
      appendValue(record.tagNames, "c");
      appendValue(record.tagValues, to_string(c));
      appendValue(record.tagNames, "t");
      appendValue(record.tagValues, to_string(t));
      appendValue(record.tagNames, "n");
      appendValue(record.tagValues, to_string(n));
   }

   if(context == "CG")
      writeRecord(gffOutCG, record);
   else if(context == "CHG")
      writeRecord(gffOutCHG, record);
   else if(context == "CHH")
      writeRecord(gffOutCHH, record);

   return 0;
}


int writeCXRecord(map<CharString, map<unsigned int, CX_report>> &testMap, 
                  unsigned pos, CharString contig, bool hasRC, 
                  CharString context, Dna5String infix)
{
   unsigned refPos = pos;
   CharString strand = '+';
   if(hasRC == true)
   {
      strand = '-';
      refPos = refPos + 2;
   }

   if( testMap[contig].find(refPos) == testMap[contig].end() )
   {
      cxOut << contig << "\t" << pos+1 << "\t" << strand;
      cxOut << "\t0\t0\t" << context << "\t" << infix << endl;
   } 
   else
   {
      cxOut << contig << "\t" << refPos+1 << "\t" << strand << "\t";
      cxOut << testMap[contig][refPos].count_methylated << "\t";
      cxOut << testMap[contig][refPos].count_unmethylated << "\t";
      cxOut << context << "\t" << infix << endl;
   }

   return 0;
}

// Writes out the CX report but includes the bases not covered
int writeOut(map<CharString, map<unsigned int, CX_report>> &testMap, 
             FaiIndex &faiIndex,ModifyStringOptions options)
{
   for(auto chromosome : testMap)
   {
      unsigned idx = 0;
      if (!getIdByName(idx, faiIndex, chromosome.first))
      {
         cerr << "Error: " << chromosome.first;
         cerr << " not found in Reference" << endl;
      }

      unsigned contigLen = sequenceLength(faiIndex, idx);

      for(unsigned i = 0; i < contigLen - 2; ++i)
      {
         Dna5String sequenceInfix;
         readRegion(sequenceInfix, faiIndex, idx, i, i+3);
         CharString context = identifyContext(sequenceInfix);

         if(context != "NNN")
         {
            writeCXRecord(testMap, i, chromosome.first, false,
                          context, sequenceInfix);
         }

         reverseComplement(sequenceInfix);
         context = identifyContext(sequenceInfix);

         if(context != "NNN")
         {
            writeCXRecord(testMap, i, chromosome.first, true,
                          context, sequenceInfix);
         }
      }

   }
   return 0;
}

int writeOutCoveredGFF(map<CharString, map<unsigned int, CX_report>> &testMap,
             FaiIndex &faiIndex, ModifyStringOptions options)
{
   for(auto chromosome : testMap)
   {
      for(auto pos : chromosome.second)
      {
         Dna5String sequenceInfix = getTrinucleotideContext(faiIndex, 
                                                            chromosome.first, 
                                                            pos.first);
         unsigned idx = 0;
         if (!getIdByName(idx, faiIndex, chromosome.first))
         {
            cerr << "Error: " << chromosome.first;
            cerr << " not found in Reference" << endl;
         }

         if(sequenceInfix[0] == 'C')
         {
            CharString context = identifyContext(sequenceInfix);
            writeGFFRecord(testMap, pos.first, chromosome.first, false, context,
                          sequenceInfix);
         }
         else if(sequenceInfix[0] == 'G')
         {
            sequenceInfix = getTrinucleotideContext(faiIndex, 
                                                    chromosome.first,
                                                    pos.first - 2);
            reverseComplement(sequenceInfix);
            CharString context = identifyContext(sequenceInfix);
            writeGFFRecord(testMap, pos.first, chromosome.first, true, context, 
                          sequenceInfix);
         } 
         else
         {
            CharString context = identifyContext(sequenceInfix);
            Dna5String bitmore = getTrinucleotideContext(faiIndex,
                                                         chromosome.first,
                                                         pos.first+1);
            Dna5String bitless = getTrinucleotideContext(faiIndex,
                                                         chromosome.first,
                                                         pos.first-1);
            Dna5String different;
            readRegion(different, faiIndex, idx, pos.first, pos.first+3);
            cerr << "WARNING ";
            cerr << sequenceInfix << " at pos " << chromosome.first << " ";
            cerr << pos.first << " " << context << "  " << bitmore << " " << bitless << " " << different << endl;
         }

      }
   }
}

int openOutputGff(ModifyStringOptions options)
{
   string gffFileNameCG = toCString(options.outputBamFileName) + string("CG.w1.gff.gz");
   string gffFileNameCHG = toCString(options.outputBamFileName) + string("CHG.w1.gff.gz");
   string gffFileNameCHH = toCString(options.outputBamFileName) + string("CHH.w1.gff.gz");

   if(!open(gffOutCG, toCString(gffFileNameCG)))
      return 1;

   if(!open(gffOutCHG, toCString(gffFileNameCHG)))
      return 1;

   if(!open(gffOutCHH, toCString(gffFileNameCHH)))
      return 1;

   return 0;
}

int openOutputCX(ModifyStringOptions options)
{
   string filename = toCString(options.outputBamFileName)
                     + string("_CX_report.txt");
   cxOut.open(filename);
}

int closeOutputCX(ModifyStringOptions options)
{
   cxOut.close();
}

int closeOutputGff(ModifyStringOptions options)
{
   close(gffOutCG);
   close(gffOutCHG);
   close(gffOutCHH);
}

bool doesBaseExistInRef(String<CigarElement<>> cigar, int seqPos)
{
   int runningTotal = 0;

   for(unsigned cigarOp = 0; cigarOp < length(cigar); cigarOp++)
   {  
      for(int i = 0; i < cigar[cigarOp].count; i++)
      {
         if(runningTotal == seqPos)
         {
            if(cigar[cigarOp].operation == 'I')
            {
               return false;
            }
            else
            {
               return true;
            }
         }
         if(cigar[cigarOp].operation == 'I' || cigar[cigarOp].operation == 'M')
            runningTotal++;
      }
   }

}

// D Deletion; the nucleotide is present in the reference but not in the read 
// I Insertion; the nucleotide is present in the read  but not in the reference.
int refPos(String<CigarElement<>> cigar, unsigned int seqPos)
{
   int seqCorrection = 0;
   unsigned runningTotal = 0;
   for(unsigned cigarOp = 0; cigarOp < length(cigar); cigarOp++)
   {
      if(cigar[cigarOp].operation == 'M')
      {
         runningTotal = runningTotal + cigar[cigarOp].count;
      }
      else if(cigar[cigarOp].operation == 'I')
      {
         seqCorrection = seqCorrection - cigar[cigarOp].count;
      }
      else if(cigar[cigarOp].operation == 'D')
      {
         seqCorrection = seqCorrection + cigar[cigarOp].count;
         runningTotal = runningTotal + cigar[cigarOp].count;
      }

      if((seqPos+seqCorrection) < runningTotal)
         return seqPos + seqCorrection;
   }
}

int process(ModifyStringOptions options)
{
   // Open BamFileIn for reading.
   BamFileIn inFile;
   if(!open(inFile, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open " << options.inputFileName;
      cerr << " for reading." << endl;
      return 1;
   }

   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.referenceFileName)))
   {
      cerr << "ERROR: Could not build reference index" << endl;
      return 1;
   }

   BamHeader header;
   readHeader(header, inFile);
   typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
   TBamContext const & bamContext = context(inFile);

   // initialise the mega vector
   map<CharString,vector<pair<unsigned int, unsigned int>>> test;

   for(unsigned i = 0; i < length(contigNames(bamContext)); ++i)
   {
      test[contigNames(bamContext)[i]].resize(contigLengths(bamContext)[i]);
   }

   // read the BAM file
   while(!atEnd(inFile))
   {
      // get the record
      BamAlignmentRecord record;
      readRecord(record, inFile);
      CharString ref = contigNames(bamContext)[record.rID];

      // get the tag
      BamTagsDict tagsDict(record.tags);
      unsigned tagIdx = 0;
      if(!findTagKey(tagIdx, tagsDict, "XM"))
      {
         cerr << "ERROR: Unknown key!\n";
      }
      CharString methylTag;
      if(!extractTagValue(methylTag, tagsDict, tagIdx))
      {
         cerr << "ERROR: " << record.qName;
         cerr << " There was an error extracting XM from tags!\n";
      }

      for(unsigned i = 0; i < length(methylTag); ++i)
      {
         if(doesBaseExistInRef(record.cigar, i) == true)
         {
            int refpos = record.beginPos + refPos(record.cigar, i);

            if(methylTag[i] != '.')
            {
               if(methylTag[i] == 'z' || methylTag[i] == 'x' || methylTag[i] == 'h')
               {
                  test[ref][refpos].second++;
               }
               else if(methylTag[i] == 'Z' || methylTag[i] == 'X' || methylTag[i] == 'H')
               {
                  test[ref][refpos].first++;
               }
            }
         }
      }
   }


   for(unsigned c = 0; c < length(contigNames(bamContext)); ++c)
   {
      vector<pair<unsigned int, unsigned int>> i = test[contigNames(bamContext)[c]];

      // let's get the whole of the chromosome into RAM
      IupacString chromosome;
      unsigned idx = 0;
      if(!getIdByName(idx, faiIndex, contigNames(bamContext)[c]))
      {
         cerr << "Error: " << contigNames(bamContext)[c] << " not found in Reference" << endl;
      }
      readSequence(chromosome, faiIndex, idx);

      unsigned pos = 0; // zero-based in RAM but +1 based when printed out
      for(auto j : i)
      {

         if(chromosome[pos] == 'C' && 
            ((pos+1 < contigLengths(bamContext)[c] && chromosome[pos+1]=='G') ||
             (pos+2 < contigLengths(bamContext)[c] && chromosome[pos+1]!='G')))
         {
            IupacString tri = chromosome[pos];
            tri += chromosome[pos+1];
            tri += chromosome[pos+2];
            cxOut << contigNames(bamContext)[c] << "\t" << (pos+1) << "\t+\t" << j.first << "\t" << j.second << "\t";
            cxOut << identifyContext(tri) << "\t" << tri << endl;
         }
         else if(chromosome[pos] == 'G' && 
                 ((pos-1 < contigLengths(bamContext)[c] && chromosome[pos-1]=='C') ||
                  (pos-2 < contigLengths(bamContext)[c] && chromosome[pos-1]!='C')))
         {
            IupacString tri = chromosome[pos-2];
            tri += chromosome[pos-1];
            tri += chromosome[pos];
            reverseComplement(tri);
            cxOut << contigNames(bamContext)[c] << "\t" << (pos+1) << "\t-\t" << j.first << "\t" << j.second << "\t";
            cxOut << identifyContext(tri) << "\t" << tri << endl;
         }
         pos++;
      }
   }
}

int main(int argc, char const ** argv)
{

   // Parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   if(options.writeCXreport == true)
      openOutputCX(options);
   else
      openOutputGff(options);

   process(options);

   if(options.writeCXreport == true)
      closeOutputCX(options);
   else
      closeOutputGff(options);

   return 0;
}
