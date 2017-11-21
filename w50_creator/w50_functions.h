#ifndef W50_FUNCTIONS_H
#define W50_FUNCTIONS_H

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
      char const ** argv);

int roundUp(int numToRound, int multiple);

int checkSorted(vector<CharString> &haveSeen, CharString &currentRef,
                ModifyStringOptions options, bool printErrorMsg);

int writeToGFF(GffFileOut &gffFileOut, CharString ref, unsigned int binPos,
               unsigned int endPos, double score, ModifyStringOptions options,
               StringSet<CharString> tagNames, StringSet<CharString> tagValues);

int writeChromosomeBins(std::map<int, WindowValues> &bins,
                        ModifyStringOptions options, GffFileOut &gffFileOut,
                        int &largest, CharString &currentRef);

int insertIntoMap(std::map<int, WindowValues> &bins,
                  ModifyStringOptions options, int &largest,
                  CharString &currentRef, GffRecord &record);

#endif
