#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>

using namespace seqan;
using namespace std;

#define LENGTH 0x1000

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
   CharString filename = "example_data/GSM1085191_mC_calls_Aa_0.tsv.gz";
   bool exclude;
   bool lazyRef;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("Overlap");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   //setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("a", "input-annotation-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   //setRequired(parser, "input-annotation-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   //setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("l", "lazy-ref", 
       "Internally it will capitalise both the input and annoation reference \
        names so that chr1, Chr1 and CHR1 will all match. The output GFF will \
        be of the same format as the annoation file."));
   setShortDescription(parser, "Overlap");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and an input w1 \
                           file, give counts for each feature");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   options.lazyRef = isSet(parser, "lazy-ref");

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

   //char lineText[100];
   //ifstream in("example_data/GSM1085191_mC_calls_Aa_0.tsv.gz");

   gzFile file = gzopen("example_data/GSM1085191_mC_calls_Aa_0.tsv.gz", "r");

   while(1)
   {
      int err;
      int bytes_read;
      unsigned char buffer[LENGTH];
      bytes_read = gzread (file, buffer, LENGTH - 1);
      buffer[bytes_read] = '\0';
    //  printf ("%s", buffer);

      if(bytes_read < LENGTH-1)
      {
         if(gzeof(file))
            break;
         else
            cerr << "ERROR: It's all gone wrong" << endl;
      }

   }

   return 0;
}
