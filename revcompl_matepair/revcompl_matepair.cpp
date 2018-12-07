#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <set>
#include <map>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
   bool exclude;
   bool lazyRef;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("revcompl_matepair");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   //setRequired(parser, "output-file");
   setShortDescription(parser, "Overlap");
   setVersion(parser, "0.0.1");
   setDate(parser, "December 2018");
   addUsageLine(parser, "-i input.fq.gz  -o revcompl.fq.gz \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, " \
                           ");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");

   return ArgumentParser::PARSE_OK;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   CharString id, qual;
   Dna5String seq;

   SeqFileIn inputFileIn;

      // Open fasta/fastq file
   if(!open(inputFileIn, (toCString(options.inputFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.inputFileName) << endl;
      return 1;
   }

   SeqFileOut fileOut;

   if(!open(fileOut, (toCString(options.outputFileName))))
   {
      cerr << "Error: could not open file ";
      cerr << toCString(options.outputFileName) << endl;
      return 1;
   }

   while(!atEnd(inputFileIn))
   {
      readRecord(id, seq, qual, inputFileIn);
      reverseComplement(seq);
      reverse(qual);
      writeRecord(fileOut, id, seq, qual);
   }
   close(inputFileIn);
   close(fileOut);

   return 0;
}
