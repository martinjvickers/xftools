#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <vector>
#include <string>
#include <map>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString outputFileName;
   int window_size;
   CharString label;
   CharString program_name;
   CharString type;
};

struct WindowValues
{
   int c;
   int t;
   int n;
   float score;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
      int argc, 
      char const ** argv)
{
   ArgumentParser parser("w50_creator");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file",
                                    "Path to the output file",
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   setShortDescription(parser, "Methylation Tools");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-i input_w1.gff [\\fIOPTIONS\\fP] ");
   addOption(parser, ArgParseOption("s", "window-size", "Size of window",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "window-size", "50");
   addOption(parser, ArgParseOption("t", "addition-type", 
                                    "Addition type specifies the behaviour \
                                    when adding two or more rows together.", 
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "addition-type", "count methyl avg");
   setDefaultValue(parser, "addition-type", "methyl");
   addOption(parser, ArgParseOption("l", "label", "Column 3 GFF output label. \
                                    Useful if using SignalMap", 
                                    ArgParseArgument::STRING, "TEXT"));
   setDefaultValue(parser, "label", "window");
   addOption(parser, ArgParseOption("p", "program_name", "Column 2 GFF output \
                                    label. Useful if using SignalMap", 
                                    ArgParseArgument::STRING, "TEXT"));
   setDefaultValue(parser, "program_name", "methyl_tools");
	
   addDescription(parser, "Create a w50 file from a w1 file.");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.window_size, parser, "window-size");
   getOptionValue(options.label, parser, "label");
   getOptionValue(options.type, parser, "addition-type");
   getOptionValue(options.program_name, parser, "program_name");

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

// Search to see if we've already seen this chromosome 
int checkSorted(vector<CharString> &haveSeen, CharString &currentRef, 
                ModifyStringOptions options)
{
   for(auto const& value: haveSeen)
   {
      if(value == currentRef)
      {
         cerr << "ERROR: The input file ";
         cerr << toCString(options.inputFileName);
         cerr << " is not sorted by chromosome. Chromosome ";
         cerr << currentRef << " is out of order. " << endl;
         cerr << "To resolve this you need to sort your file by ";
         cerr << "Chromosome order, something like; sort -k1n ";
         cerr << toCString(options.inputFileName) << " > ";
         cerr << toCString(options.inputFileName);
         cerr << "_sorted.w1.gff" << endl;
         cerr << "And then rerun this program using the ";
         cerr << toCString(options.inputFileName);
         cerr << "_sorted.w1.gff file" << endl;

         return 1;
      }
   }

   return 0;
}

int writeToGFF(GffFileOut &gffFileOut, CharString ref, unsigned int binPos,
               unsigned int endPos, double score, ModifyStringOptions options, 
               StringSet<CharString> tagNames, StringSet<CharString> tagValues)
{  
   GffRecord record;
   record.ref = ref;
   record.source = options.program_name;
   record.type = options.label;
   record.beginPos = (binPos - (options.window_size));
   record.endPos = endPos;
   record.score = score;
   record.strand = '.';
   record.phase = '.';
   record.tagNames = tagNames;
   record.tagValues = tagValues;

   writeRecord(gffFileOut, record);

   return 0;  
}

int writeChromosomeBins(std::map<int, WindowValues> &bins, 
                        ModifyStringOptions options, GffFileOut &gffFileOut, 
                        int &largest, CharString &currentRef)
{
   for(std::map<int, WindowValues>::iterator iter = bins.begin(); 
       iter != bins.end(); ++iter)
   {
      float score;
      auto bin = *iter;
      int end = bin.first;

      // Fixes an odd situation where in the methylation input file you have
      // something like this. The DZLAB program just skips, so this does too
      // chr5 . CG 23576549 23576549 0.0000 . . c=0;t=0;n=1
      if((options.type == "methyl") && (bin.second.c == 0) && 
         (bin.second.t == 0))
         continue;

      if(iter == --bins.end())
         end = largest + 1;

      StringSet<CharString> tagNames;
      StringSet<CharString> tagValues;

      if(options.type == "methyl")
      {
         appendValue(tagNames, "c");
         appendValue(tagValues, to_string(bin.second.c));
         appendValue(tagNames, "t");
         appendValue(tagValues, to_string(bin.second.t));

         score = (float)bin.second.c / (float)(bin.second.c + bin.second.t);
      }
      else if(options.type == "count")
      {
         score = bin.second.n;
      }
      else if(options.type == "avg") 
      {
         score = bin.second.score / bin.second.n;
      }

      appendValue(tagNames, "n");
      appendValue(tagValues, to_string(bin.second.n));

      writeToGFF(gffFileOut, currentRef, bin.first, end, score,
                 options, tagNames, tagValues);

   }

   return 0;
}

int insertIntoMap(std::map<int, WindowValues> &bins,
                  ModifyStringOptions options, int &largest, 
                  CharString &currentRef, GffRecord &record)
{
   // SeqAn uses 0-based half-open intervals for internal representation
   // however all the previous dzlab stuff (at least w50 creation so far)
   // uses a 1-based half-open interval.

   // READ HERE
   // https://seqan.readthedocs.io/en/master/Tutorial/InputOutput/GffAndGtfIO.html
   int window = roundUp(record.beginPos + 1, options.window_size);

   // We expect the following format of tags and values  c=4;t=0;n=1
   // but we will not assume they are always in that order when reading 
   // them from the file. Since we have a struct to store the values,
   // we don't need to store the tag names themselves, just values.

   int c, t, n;
   WindowValues currWindowVal;
   bool haveN = false;

   // iterate through looking for c,t,n (will ignore anything else)
   for(int i = 0; i < length(record.tagNames); i++)
   {
      if(record.tagNames[i] == "c")
         currWindowVal.c = atoi(toCString(record.tagValues[i]));

      if(record.tagNames[i] == "t")
         currWindowVal.t = atoi(toCString(record.tagValues[i]));

      if(record.tagNames[i] == "n")
      {
         currWindowVal.n = atoi(toCString(record.tagValues[i]));
         haveN = true;
      }
   }

   currWindowVal.score = record.score;

   // means we have no n's in the input file, so we add one
   if(haveN == true)
      currWindowVal.n = 1;

   // get element in hash table, if not there, just put this one in it
   // if window exists, get values and add current Window to it
   WindowValues get = bins[window];
   currWindowVal.c = get.c + currWindowVal.c;
   currWindowVal.t = get.t + currWindowVal.t;
   currWindowVal.n = get.n + currWindowVal.n;
   currWindowVal.score = get.score + currWindowVal.score;

   // now insert back into map
   bins[window] = currWindowVal;

   // update currentRef
   currentRef = record.ref;
   largest = record.beginPos;

   return 0;
}

int main(int argc, char const ** argv)
{
   // Parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // Open out files
   GffFileIn gffFileIn;
   if (!open(gffFileIn, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }

   GffFileOut gffFileOut;
   if (!open(gffFileOut, toCString(options.outputFileName)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }

   CharString currentRef = NULL;
   GffRecord record;

   // Define our map to store bin
   std::map<int, WindowValues> bins;

   // Stores chromosomes we've seen
   vector<CharString> haveSeen;

   int gffLineNumber = 0;
   int largestBaseFound = 0;

   // Copy the file record by record.
   while (!atEnd(gffFileIn))
   {
      try
      {
         readRecord(record, gffFileIn);
         gffLineNumber++;
      }
      catch (Exception const & e)
      {
         cerr << "ERROR: " << e.what() << endl;
         cerr << "ERROR: " << record.ref << " ";
         cerr << currentRef << " " << gffLineNumber << endl;
         return 1;
      }

      // First Chromosome
      if(currentRef == NULL)
         currentRef = record.ref;		

      // At this point we will have dump the contents of the hash
      if(record.ref != currentRef)
      {
         if(checkSorted(haveSeen, currentRef, options) == 1)
            return 1;

         writeChromosomeBins(bins, options, gffFileOut, largestBaseFound, 
                             currentRef);

         // We should save the current reference as a check to see if we come 
         // across it again out of order
         haveSeen.push_back(currentRef);

         bins.clear();
         currentRef = record.ref;
         largestBaseFound = 0;
      }

      // here, we add to the hash
      if(record.ref == currentRef)
         insertIntoMap(bins, options, largestBaseFound, currentRef, record);

   }

   // read out very last chromosome
   // I hate this, there has to be a better way of doing this loop
   writeChromosomeBins(bins, options, gffFileOut, largestBaseFound, 
                       currentRef);

   return 0;
}
