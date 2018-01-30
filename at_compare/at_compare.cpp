#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName1;
   CharString inputAnnotationFileName2;
   CharString outputFileName1;
   CharString outputFileName2;
   bool exclude;
   bool lazyRef;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("at_compare");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("1", "input-annotation-file-1", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-annotation-file-1");
   addOption(parser, ArgParseOption("2", "input-annotation-file-2",
                                    "Path to the input filter file",
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-annotation-file-2");
   addOption(parser, ArgParseOption("o1", "output-file-1", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file-1");
   addOption(parser, ArgParseOption("o2", "output-file-2",
                                    "Path to the output file",
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file-2");
   addOption(parser, ArgParseOption("l", "lazy-ref", 
       "Internally it will capitalise both the input and annoation reference \
        names so that chr1, Chr1 and CHR1 will all match. The output GFF will \
        be of the same format as the annoation file."));
   setShortDescription(parser, "at_compare");
   setVersion(parser, "0.0.1");
   setDate(parser, "January 2018");
   addUsageLine(parser, "-i input.fa -1 annotation.gff -2 annotation2.gff -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and an input w1 \
                           file, give counts for each feature");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName1, 
                  parser, "input-annotation-file-1");
   getOptionValue(options.inputAnnotationFileName2,
                  parser, "input-annotation-file-2");
   getOptionValue(options.outputFileName1, parser, "output-file-1");
   getOptionValue(options.outputFileName2, parser, "output-file-2");
   options.lazyRef = isSet(parser, "lazy-ref");

   return ArgumentParser::PARSE_OK;
}

double calculate(Dna5String seq)
{
   unsigned count = 0;
   for(unsigned i = 0; i < length(seq) - 1; i++)
      if(seq[i] == 'A' || seq[i] == 'T')
         count++;

   return count / (double)length(seq);
}

double f1 = 0.0;
double f2 = 0.0;
unsigned f1c = 0;
unsigned f2c = 0;

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // create the fasta index
   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.inputFileName)))
      std::cout << "ERROR: Could not build the index!\n";

   // read in first file
   GffFileIn gffIn1;
   if(!open(gffIn1, toCString(options.inputAnnotationFileName1)))
   {
      std::cerr << "ERROR: Could not open example.gff" << std::endl;
      return 1;
   }

   // read in second file
   GffFileIn gffIn2;
   if(!open(gffIn2, toCString(options.inputAnnotationFileName2)))
   {
      std::cerr << "ERROR: Could not open example.gff" << std::endl;
      return 1;
   }

   // open output files
   GffFileOut gffOut1;
   if(!open(gffOut1, toCString(options.outputFileName1)))
   {
      std::cerr << "ERROR: Could not open example.gff" << std::endl;
      return 1;
   }

   GffFileOut gffOut2;
   if(!open(gffOut2, toCString(options.outputFileName2)))
   {
      std::cerr << "ERROR: Could not open example.gff" << std::endl;
      return 1;
   }

   GffRecord record;

   try
   {
      while(!atEnd(gffIn1))
      {
         readRecord(record, gffIn1);
         // translate sequence name into an index id
         unsigned idx = 0;
         if(!getIdByName(idx, faiIndex, record.ref))
         {
            cerr << "ERROR: Index does not know about sequence " << record.ref << "\n";
            return 1;
         }
         Dna5String sequenceInfix;
         readRegion(sequenceInfix, faiIndex, idx, record.beginPos, record.endPos);
         appendValue(record.tagNames, "AT");
         double at;
         if(record.strand == '-')
         {
            reverseComplement(sequenceInfix);
            at = calculate(sequenceInfix);
         }
         else
         {
            at = calculate(sequenceInfix);
         }
         f1 = f1 + at;
         f1c++;
         ostringstream strs;
         strs << at;
         appendValue(record.tagValues, strs.str());
         writeRecord(gffOut1, record);
      }
   }
   catch(Exception const & e)
   {
      cerr << "ERROR: " << e.what() << std::endl;
      return 1;
   }

   try
   {
      while(!atEnd(gffIn2))
      {
         readRecord(record, gffIn2);

         // translate sequence name into an index id
         unsigned idx = 0;
         if(!getIdByName(idx, faiIndex, record.ref))
         {
            cerr << "ERROR: Index does not know about sequence " << record.ref << "\n";
            return 1;
         }
         Dna5String sequenceInfix;
         readRegion(sequenceInfix, faiIndex, idx, record.beginPos, record.endPos);
         appendValue(record.tagNames, "AT");
         double at; 
         if(record.strand == '-')
         {
            reverseComplement(sequenceInfix);
            at = calculate(sequenceInfix);
         }
         else
         {
            at = calculate(sequenceInfix);
         }
         f2 = f2 + at;
         f2c++;
         ostringstream strs;
         strs << at;
         appendValue(record.tagValues, strs.str());
         writeRecord(gffOut2, record);
      }
   }
   catch(Exception const & e)
   {
      cerr << "ERROR: " << e.what() << std::endl;
      return 1;
   }

   cout << options.inputAnnotationFileName1 << " AT content " << f1/(float)f1c << endl;
   cout << options.inputAnnotationFileName2 << " AT content " << f2/(float)f2c << endl;

   return 0;
}
