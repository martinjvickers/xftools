#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include <string>
#include <thread>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString annotationFileName;
   CharString outputFileName;
   bool lazy_match = false;
   CharString type;
   CharString attribute;
};

struct AnnotationFeature
{
   unsigned int beginPos;
   unsigned int endPos;
   vector<GffRecord> records;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv)
{
   ArgumentParser parser("window_by_annotation");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input w1 file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("a", "annotation-file", 
                                    "Path to your annotation file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "annotation-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to your output file", 
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("t", "score-type",
                                    "Specifies the calculation shown in the \
                                    score column ($6) of the output GFF.", 
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "score-type", "sum methyl avg");
   setDefaultValue(parser, "score-type", "methyl");

   addOption(parser, ArgParseOption("attributes", "treat-attributes",
                                    "Addition type specifies the behaviour \
                                    when adding two or more rows together.",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "treat-attributes", "sum avg");
   setDefaultValue(parser, "treat-attributes", "sum");

   setShortDescription(parser, "XFTOOLS");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-i input.w1.gff -a reference.gff \
                        -o out.gff [\\fIOPTIONS\\fP] ");

   addDescription(parser, "Calculates the pausing index of genes");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.annotationFileName, parser, "annotation-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.type, parser, "score-type");
   getOptionValue(options.attribute, parser, "treat-attributes");
   
   return ArgumentParser::PARSE_OK;
}

// Read the entire single C file into RAM for fast look up.
// TODO: do this in chromosome chunks to save RAM in the future
int process_single_c(GffFileIn &gffFileIn, 
                     map<CharString, map<unsigned int, GffRecord>> &w1_file,
                     double &sum_score)
{
   while(!atEnd(gffFileIn))
   {
      GffRecord record;
      readRecord(record, gffFileIn);
      w1_file[record.ref][record.beginPos] = record;
   }

   return 0;
}

// Iterate through the annotation file and extract the relevent
// single-c records from the single-c map
int windowByAnnotation(GffFileIn &gffAnnotationIn,
                       map<CharString, map<unsigned int, GffRecord>> &w1_file,
                       vector<AnnotationFeature> &annotation,
                       ModifyStringOptions options)
{
   while(!atEnd(gffAnnotationIn))
   {
      GffRecord record;
      readRecord(record, gffAnnotationIn);
      AnnotationFeature currentFeature;

      for(unsigned int pos = record.beginPos; pos < record.endPos; ++pos)
      {
         currentFeature.records.push_back(w1_file[record.ref][pos]);
         currentFeature.beginPos = record.beginPos;
         currentFeature.endPos = record.endPos;
      }

      annotation.push_back(currentFeature);
   }

   return 0;

}

// collect all the c's and t's from the records associated with this feature
int calcMethyl(AnnotationFeature &feature)
{
   unsigned int c, t, n;
   for(GffRecord record : feature.records)
   {
      StringSet<CharString> tagNames = record.tagNames;
      StringSet<CharString> tagValues = record.tagValues;
      bool haveC, haveT, haveN = false;

      for(unsigned tagPos = 0; tagPos < length(tagNames); ++tagPos)
      {
         if(tagNames[tagPos] == 'c')
         {
         
         }
         else if(tagNames[tagPos] == 't')
         {

         }
         else if(tagNames[tagPos] == 'n')
         {

         }
      }

      if(haveC != true && haveT != true && haveN != true)
      {
         cerr << "ERROR: The single-c file has an issue with the C T N tags \
                  in column 9. You cannot calculate methylation score is \
                  this is not corrected. Please check your file." << endl;
         return 1;
      }
   }
}

void calcSum()
{}

void calcAvg()
{}

int processFeatureRecords(AnnotationFeature &feature, 
                          ModifyStringOptions options)
{
   if(options.type == "methyl")
      calcMethyl(feature);
   else if(options.type == "sum")
      calcSum();
   else if(options.type == "avg")
      calcAvg();
   else
      return 1;


   return 0;
}

// Iterates through the annotation features and performs calculation
int calculateFeatureValues(ModifyStringOptions options, 
                           vector<AnnotationFeature> &annotation)
{
   for(AnnotationFeature feature : annotation)
   {
      // perform some calculation
      processFeatureRecords(feature, options);

      // print to output file
   }
   return 0;
}

int main(int argc, char const ** argv)
{
   // Parse our options
   ModifyStringOptions options;
   seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;

   GffFileIn gffFileIn, gffAnnotationIn;
   if (!open(gffFileIn, toCString(options.inputFileName)))
      return 1;
   if (!open(gffAnnotationIn, toCString(options.annotationFileName)))
      return 1;

   // Let's make a map of the input file
   double sum_score = 0.0;
   map<CharString, map<unsigned int, GffRecord>> w1_file;
   process_single_c(gffFileIn, w1_file, sum_score);
   close(gffFileIn);

   // Now we make and process the annoation file
   vector<AnnotationFeature> annotation;
   windowByAnnotation(gffAnnotationIn, w1_file, annotation, options);
   close(gffAnnotationIn);
   w1_file.clear();

   calculateFeatureValues(options, annotation);

//	// now we put the exons into RAM and get our TSS/gene info
//	map< CharString, map <CharString, geneElement>> exons;
///	process_annotation(gffAnnotationIn, exons, options);
//	close(gffAnnotationIn);

	// do the calculations
//	calculate_counts(exons, w1_file, options);

	// write out
//	print_out(exons, sum_score, options);

//	cout.precision(17);
//	cout << "The sum of the score (column 6) from " << options.inputFileName << " was " << sum_score << endl;

   return 0;
}
