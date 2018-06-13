#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("max_wba");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Input data file, e.g. w50s", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("a", "input-annotation-file", 
                                    "Input annotation file, e.g. DMRs",
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-annotation-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   setShortDescription(parser, "WBA, returning max score record");
   setVersion(parser, "0.0.1");
   setDate(parser, "June 2018");
   addUsageLine(parser, "-i input.w50.gff -a annotation.dmrs.gff -o output.gff\
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and an input GFF \
                           file, it will return the record in the input GFF \
                           with the maximum score that resides within the \
                           range of the annotation file. If more than one \
                           record is found with equal max score, it returns\
                           each with that max score. \
                           NOTE: this program is new and relatively untested");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
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

   GffFileIn gffIn(toCString(options.inputAnnotationFileName));
   GffFileOut gffOut(toCString(options.outputFileName));
   
   GffRecord record;
   while (!atEnd(gffIn))
   {
      readRecord(record, gffIn);

      GffFileIn gffDatIn(toCString(options.inputFileName));
      GffRecord qrecord;

      bool set = false;
      map<float, GffRecord> results;
      while (!atEnd(gffDatIn))
      {
         readRecord(qrecord, gffDatIn);
         if(record.ref == qrecord.ref && qrecord.beginPos >= record.beginPos && qrecord.endPos <= record.endPos)
         {
            results.insert( std::pair<float, GffRecord>(qrecord.score, qrecord) );
         }
      }
      close(gffDatIn);

      unsigned int count = 0;
      float last = NULL;

      for (auto it = results.rbegin(); it != results.rend(); ++it)
      {
         if(count == 0)
         {
            last = it->second.score;
            writeRecord(gffOut, it->second);
         }
         else if(last == it->second.score)
         {
            last = it->second.score;
            writeRecord(gffOut, it->second);
         }
         else
         {
            break;
         }

         count++;
      }
   }

   return 0;
}
