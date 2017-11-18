#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
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
   setShortDescription(parser, "Methylation Tools");
   setVersion(parser, "0.0.3");
   setDate(parser, "Jan 2017");
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

int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   GffFileIn gffFileIn;
   if (!open(gffFileIn, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }

   bool meth = false;

   CharString currentRef = NULL;

   // Define our map to store windows
   map<int, WindowValues> map;

   cout << setprecision(4) << fixed;

   // Stores chromosomes we've seen
   vector<CharString> haveSeen;

   int count = 0;
   int largest = 0;

   // Copy the file record by record.
   GffRecord record;
   while (!atEnd(gffFileIn))
   {
      try
      {
         readRecord(record, gffFileIn);
         count++;
      }
      catch (Exception const & e)
      {
         cerr << "ERROR: " << e.what() << endl;
         cerr << "ERROR: " << record.ref << " ";
         cerr << currentRef << " "<< count << endl;
         return 1;
      }

      // First Chromosome
      if(currentRef == NULL)
      {
         currentRef = record.ref;		
      }

      // At this point we will have dump the contents of the hash
      if(record.ref != currentRef)
      {
         // Search to see if we've already seen this chromosome	
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
	
         for(auto iter = map.begin(); iter != map.end(); ++iter)
         {
            float score;
            auto p = *iter;

            // Write out GFF
            if(options.type == "methyl")
            {
               if((p.second.c > 0) || (p.second.t > 0))
               {
                  score = (float)p.second.c / (float)(p.second.c + p.second.t);
                  int end = p.first;
                  
                  if(iter == --map.end())
                  {
                     end = largest;
                  }

                  cout << currentRef << "\t";
                  cout << toCString(options.program_name) << "\t";
                  cout << toCString(options.label) << "\t";
                  cout << p.first - (options.window_size-1) << "\t";
                  cout << end << "\t"<< score << "\t.\t.\tc=" << p.second.c;
                  cout << ";t=" << p.second.t << ";n=" << p.second.n << endl;
               }
               else
               {
                  score = 0.0000;
               }
            } 
            else if(options.type == "count")
            {
               score = p.second.score;
               cout << currentRef << "\t" << toCString(options.program_name);
               cout << "\t" << toCString(options.label) << "\t";
               cout << p.first - (options.window_size - 1) << "\t";
               cout << p.first << "\t"<< score;
               cout << "\t.\t.\tn=" << p.second.n << endl;
            }
            else if(options.type == "avg") 
            {
               score = p.second.score / p.second.n;
               cout << currentRef << "\t" << toCString(options.program_name);
               cout << "\t" << toCString(options.label) << "\t";
               cout << p.first - (options.window_size - 1) << "\t";
               cout << p.first << "\t"<< score;
               cout << "\t.\t.\tn=" << p.second.n << endl;
            }
         }
	
         // We should save the current reference as a check to see if we come 
         // across it again out of order
         haveSeen.push_back(currentRef);

         map.clear();
         currentRef = record.ref;
         largest = 0;
      }

      // here, we add to the hash
      if(record.ref == currentRef)
      {
         //////ARGH!!!
         // SeqAn uses 0-based half-open intervals for internal representation
         // however all the previous dzlab stuff (at least w50 creation so far)
         // uses a 1-based half-open interval.

         // READ HERE
         // https://seqan.readthedocs.io/en/master/Tutorial/InputOutput/GffAndGtfIO.html

         // int window = roundUp(record.beginPos+1, options.window_size);
         // The above 0-based half open stuff is done by seqan when reading the record. 
         // This causes a bunch of problems when you have a base at position 0. I've
         // resolved this in SeqAn so now the representation should be exactly as the 
         // GFF is written. So no need for the above record.beginPos+1
         int window = roundUp(record.beginPos, options.window_size);

         // we expect the following format of tags and values  c=4;t=0;n=1
         // but we will not assume they are always in that order when reading 
         // them from the file. Since we have a struct to store the values,
         // we don't need to store the tag names themselves, just values.
			
         int c, t, n;
         WindowValues currWindowVal;
         bool haveN = 1;

         // iterate through looking for c,t,n (will ignore anything else)
         for(int i = 0; i < length(record.tagNames); i++)
         {	
            if(record.tagNames[i] == "c")
            {
               currWindowVal.c = atoi(toCString(record.tagValues[i]));
            }

            if(record.tagNames[i] == "t")
            {
               currWindowVal.t = atoi(toCString(record.tagValues[i]));
            }

            if(record.tagNames[i] == "n")
            {
               currWindowVal.n = atoi(toCString(record.tagValues[i]));
               haveN = 0;
            }
         }
         
         currWindowVal.score = record.score;

         // means we have no n's in the input file, so we add one
         if(haveN == 1)
         {
            currWindowVal.n = 1;
         }

         // get element in hash table, if not there, just put this one in it
         // if window exists, get values and add current Window to it
         WindowValues get = map[window];
         currWindowVal.c = get.c + currWindowVal.c;
         currWindowVal.t = get.t + currWindowVal.t;
         currWindowVal.n = get.n + currWindowVal.n;

         currWindowVal.score = get.score + currWindowVal.score;

         // now insert back into map
         map[window] = currWindowVal;

         // update currentRef
         currentRef = record.ref;
         largest = record.beginPos;
      }
   }

   // read out very last chromosome
   for(auto iter = map.begin(); iter != map.end(); ++iter)
   {
      // search to see if we've already seen this chromosome
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


      float score;
      auto p = *iter;

      // write out GFF 
      if(options.type == "methyl")
      {
         if((p.second.c > 0) || (p.second.t > 0))
         {
            score = (float)p.second.c / (float)(p.second.c + p.second.t);
            int end = p.first;

            if(iter == --map.end())
            {
               end = largest;
            }

            score = (float)p.second.c / (float)(p.second.c + p.second.t);
            cout << currentRef << "\t" << toCString(options.program_name);
            cout << "\t" << toCString(options.label) << "\t";
            cout << p.first - (options.window_size - 1) << "\t";
            cout << end << "\t"<< score << "\t.\t.\tc=";
            cout << p.second.c << ";t=" << p.second.t;
            cout << ";n=" << p.second.n << endl;
         }
         else
         {
            score = 0;
         }
      }
      else if(options.type == "count")
      {
         score = p.second.n;
         cout << currentRef << "\t" << toCString(options.program_name);
         cout << "\t" << toCString(options.label) << "\t";
         cout << p.first - (options.window_size - 1) << "\t";
         cout << p.first << "\t"<< score;
         cout << "\t.\t.\tn=" << p.second.n << endl;
      }
      else if(options.type == "avg")
      {
         score = p.second.score / p.second.n;
         cout << currentRef << "\t" << toCString(options.program_name);
         cout << "\t" << toCString(options.label) << "\t";
         cout << p.first - (options.window_size-1) << "\t";
         cout << p.first << "\t"<< score;
         cout << "\t.\t.\tn=" << p.second.n << endl;
      }
   }

   return 0;
}
