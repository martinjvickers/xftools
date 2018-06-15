#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <vector>
#include <seqan/gff_io.h>

#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/misc/interval_tree.h>

using namespace seqan;
using namespace std;

typedef IntervalAndCargo<int, GffRecord> TInterval;

struct ModifyStringOptions 
{
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
   bool exclude;
   bool lazyRef;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) 
{

   ArgumentParser parser("Overlap");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
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

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != ArgumentParser::PARSE_OK)
     return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   options.lazyRef = isSet(parser, "lazy-ref");

   return ArgumentParser::PARSE_OK;
}

//int insertIntervals(String<TInterval> &intervals, GffFileIn &gffInFile) 
int insertIntervals(map<CharString, String<TInterval>> &intervals, GffFileIn &gffInFile)
{
   GffRecord record;

   try
   {
      while(!atEnd(gffInFile))
      {
         readRecord(record, gffInFile);
         appendValue(intervals[record.ref], TInterval(record.beginPos, record.endPos, record));
      }
   }
   catch (Exception const & e)
   {
      cerr << "ERROR: " << e.what() << std::endl;
      return 1;
   }

   close(gffInFile);

   return 0;
}

// Given two gff files, this program will find where file B 
// overlaps file A. I will start with a nieve approach.
int main(int argc, char const ** argv) 
{
   // Parse our input options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // open annotation file
   GffFileIn gffAnnotInFile;
   if(!open(gffAnnotInFile, toCString(options.inputAnnotationFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputFileName;
      cerr << " for reading.\n";
      return 1;
   }

   // create interval tree
   typedef IntervalAndCargo<int, GffRecord> TInterval;
   String<TInterval> intervals;

   map<CharString, String<TInterval>> full;

   if(insertIntervals(full, gffAnnotInFile) != 0)
   {
      cerr << "ERROR: Cannot put intervals into interval vector " << endl;
      return 1;
   }

   map<CharString, IntervalTree<int, GffRecord>> trees;
   for(auto i : full)
      trees[i.first] = i.second;

   // open input file
   GffFileIn gffInFile;
   if(!open(gffInFile, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputFileName;
      cerr << " for reading.\n";
      return 1;
   }

   // search for overlaps
   GffRecord record;
   while(!atEnd(gffInFile))
   {
      readRecord(record, gffInFile);
      String<GffRecord> results;
      findIntervals(results, trees[record.ref], record.beginPos, record.endPos);
      for(unsigned i = 0; i < length(results); ++i)
         cout << record.ref << "\t" << record.beginPos << "\t" << record.endPos << "\t" << results[i].ref << "\t" << results[i].beginPos << "\t" << results[i].endPos << endl;
   }

   close(gffInFile);
   close(gffAnnotInFile);

   return 0;
}
