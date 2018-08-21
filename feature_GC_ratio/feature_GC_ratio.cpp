#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputReferenceFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("feature_GC_ratio");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input feature file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("r", "input-reference-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-reference-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   setShortDescription(parser, "Calculate the GC content of a GFF of features");
   setVersion(parser, "0.0.1");
   setDate(parser, "August 2018");
   addUsageLine(parser, "-r reference.fa -i features.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and the reference \
                           fasta, the GC content is displayed");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputReferenceFileName, parser, "input-reference-file");
   getOptionValue(options.inputAnnotationFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");

   return ArgumentParser::PARSE_OK;
}

double GC_ratio(CharString &seq)
{
   unsigned int GC_count = 0;
   unsigned int tot_count = 0;

   for(unsigned int i = 0; i < length(seq); i++)
   {
      if(seq[i] == 'C' || seq[i] == 'G')
         GC_count++;
     
      tot_count++;
   }

   return (double)GC_count/(double)tot_count;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // load reference genome
   // create index
   // load annotation GFF

   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.inputReferenceFileName)))
   {
      cerr << "ERROR: Could not build reference index" << endl;
      return 1;
   }

   GffFileOut gffOut(toCString(options.outputFileName));

   GffFileIn gffIn(toCString(options.inputAnnotationFileName));
   GffRecord record;
   while (!atEnd(gffIn))
   {
      readRecord(record, gffIn);

      // get the id for the reference
      unsigned idx = 0;
      if(!getIdByName(idx, faiIndex, record.ref))
         std::cerr << "ERROR: FAI index has no entry for " << record.ref << "\n";

      CharString seq;
      readRegion(seq, faiIndex, idx, record.beginPos, record.endPos);
      appendValue(record.tagNames, "GC_ratio");
      appendValue(record.tagValues, to_string(GC_ratio(seq)));
      //cout << record.ref << "\t" << record.beginPos << "\t" << record.endPos << "\t" << GC_ratio(seq) << endl;
      writeRecord(gffOut, record);
   }

   return 0;
}
