#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString outputFileName;
   CharString label;
   int window_size;
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
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("s", "window-size", "Size of window",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "window-size", "50");
   addOption(parser, ArgParseOption("l", "label", "Column 3 GFF output label.\
                                    Useful if using SignalMap", 
                                    ArgParseArgument::STRING, "TEXT"));
   setShortDescription(parser, "RPKM maker");
   setVersion(parser, "0.0.1");
   setDate(parser, "August 2019");
   addUsageLine(parser, "-i input.bam -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and an input w1 \
                           file, give counts for each feature");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.label, parser, "label");
   getOptionValue(options.window_size, parser, "window-size");
   return ArgumentParser::PARSE_OK;
}

int roundUp(int numToRound, int multiple)
{
   if(multiple == 0)
      return numToRound;

   int remainder = numToRound % multiple;
   
   if(remainder == 0)
      return numToRound;

   return numToRound + multiple - remainder;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // open 
   BamFileIn inFile;
   if (!open(inFile, toCString(options.inputFileName)))
   {
      std::cerr << "ERROR: Could not open " << options.inputFileName;
      std::cerr << " for reading." << std::endl;
      return 1;
   }

   // Read header.
   BamHeader header;
   readHeader(header, inFile);
   typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
   TBamContext const & bamContext = context(inFile);

   // check to see if it's sorted
   CharString tagValue;
   bool keyFound;

   for(int i = 0; i < length(header); i++)
   {
       keyFound = getTagValue(tagValue, "SO", header[i]);
       if(keyFound)
       {
          break;
       }
   }

   if(keyFound == true && tagValue == "coordinate")
   {
      // all is well!
   }
   else if(keyFound == false)
   {
      std::cerr << "ERROR: SO tag not found. We don't know if";
      std::cerr << " this BAM is sorted." << std::endl;
      return 1;
   }
   else if(tagValue != "coordinate")
   {
      std::cerr << "ERROR: BAM is not sorted." << std::endl;
      return 1;
   }

   BamAlignmentRecord record;
   map<int, map<int, int>> values;
   long long unsigned int num_reads = 0;

   while(!atEnd(inFile))
   {
      readRecord(record, inFile);

      if(hasFlagUnmapped(record))
      {
         // nothing to do as it's unmapped
      }
      else
      {
         int window = roundUp(record.beginPos+1, options.window_size);
         values[record.rID][window]++;
         num_reads++;
      }
   }

   close(inFile);

   cout << "Number of Mapped reads processed: " << num_reads << endl;

   GffFileOut gffOut;
   if(!open(gffOut, toCString(options.outputFileName)))
   {
      std::cerr << "ERROR: Could not open ";
      std::cerr << toCString(options.outputFileName)  << std::endl;
      return 1;
   }

   for(auto i : values)
   {
      for(auto j : i.second)
      {
         GffRecord record;
         record.ref = contigNames(bamContext)[i.first];
         record.source = "xftools";
         record.type = toCString(options.outputFileName);
         record.beginPos = j.first - options.window_size;
         record.endPos = j.first;
         record.strand = '.';
         record.score = (double)j.second * (double)1000000000 / (double)num_reads / (double)options.window_size;
         appendValue(record.tagNames, "n");
         appendValue(record.tagValues, to_string(j.second));
         writeRecord(gffOut, record);
      }
   }

   close(gffOut);

   return 0;
}
