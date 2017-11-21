#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <vector>
#include <string>
#include <map>

#include "w50_functions.h"

using namespace seqan;
using namespace std;

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
         if(checkSorted(haveSeen, currentRef, options, true) == 1)
            return 1;

         writeChromosomeBins(bins, options, gffFileOut, largestBaseFound, 
                             currentRef);

         // We should save the current reference as a check to see if we come 
         // across it again out of order
         haveSeen.push_back(currentRef);

         // Now that we've written out this chromosome information we can reset
         bins.clear();
         currentRef = record.ref;
         largestBaseFound = 0;
      }

      // insert record into our map
      if(record.ref == currentRef)
         insertIntoMap(bins, options, largestBaseFound, currentRef, record);

   }

   // read out very last chromosome
   // I hate this, there has to be a better way of doing this loop
   writeChromosomeBins(bins, options, gffFileOut, largestBaseFound, 
                       currentRef);

   close(gffFileIn);
   close(gffFileOut);

   return 0;
}
