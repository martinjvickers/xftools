#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/misc/interval_tree.h>
#include <set>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputSLMFileName;
   CharString inputSpermFileName;
   CharString inputS1FileName;
   CharString inputS2FileName;
   CharString inputS3FileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("SLMextend");
   addOption(parser, ArgParseOption("slm", "slm-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "slm-file");

   addOption(parser, ArgParseOption("spm", "sperm-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "sperm-file");

   addOption(parser, ArgParseOption("s1", "soma1-file",                
                                    "Path to the input filter file",
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "soma1-file");

   addOption(parser, ArgParseOption("s2", "soma2-file",     
                                    "Path to the input filter file",
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "soma2-file");

   addOption(parser, ArgParseOption("s3", "soma3-file",
                                    "Path to the input filter file",
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "soma3-file");

   setShortDescription(parser, "SLMextend");
   setVersion(parser, "0.0.1");
   setDate(parser, "February 2019");
   addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given a list of SLMs along with several soma \
                           tissue samples, extend SLM given specific \
                           methylation criteria.");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputSLMFileName, parser, "slm-file");
   getOptionValue(options.inputSpermFileName, parser, "sperm-file");
   getOptionValue(options.inputS1FileName, parser, "soma1-file");
   getOptionValue(options.inputS2FileName, parser, "soma2-file");
   getOptionValue(options.inputS3FileName, parser, "soma3-file");

   return ArgumentParser::PARSE_OK;
}

typedef IntervalAndCargo<int, GffRecord> TInterval;

int insertIntervals(map<CharString, String<TInterval>> &intervals, 
                    GffFileIn &gffInFile)
{
   GffRecord record;

   try
   {
      while(!atEnd(gffInFile))
      {
         readRecord(record, gffInFile);
         appendValue(intervals[record.ref], 
                     TInterval(record.beginPos, record.endPos, record));
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

//create average!
int bedAverage(GffFileIn &gffS1InFile, GffFileIn &gffS2InFile,
               GffFileIn &gffS3InFile)
{
   
   return 0;
}

std::set<CharString> getList(GffFileIn &gffSLMInFile)
{
   std::set<CharString> chr;
   GffRecord record;

   while(!atEnd(gffSLMInFile))
   {
      readRecord(record, gffSLMInFile);
      chr.insert(record.ref);
   }

   cout << "Unique chrs " << chr.size() << endl;

   return chr;
}

/*
Brief from Jimmy

Take SLM.
Take the adjacent 100bp from either side of SLM

Calculate the % methylation for these 100bp regions for the sperm and three 
somatic tissues using c's and t's

Calculate the average % methylation for soma using the three individual %.

If the % methylation difference between sperm and the average soma is 
greater than x, then add the region onto the SLM. Repeat until the previous 
statement is untrue.
*/
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // read in the SLM file
   GffFileIn gffSLMInFile;
   if(!open(gffSLMInFile, toCString(options.inputSLMFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputSLMFileName;
      cerr << " for reading.\n";
      return 1;
   }

   // get a uniq list of chromosomes?
   std::set<CharString> chromos = getList(gffSLMInFile);
   close (gffSLMInFile);

   // reopen - because I'm not sure how else to do this;
   if(!open(gffSLMInFile, toCString(options.inputSLMFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputSLMFileName;
      cerr << " for reading.\n";
      return 1;
   }


   // load the sperm and soma files
   GffFileIn gffSpermInFile;
   if(!open(gffSpermInFile, toCString(options.inputSpermFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputSpermFileName;
      cerr << " for reading.\n";
      return 1;
   }

   GffFileIn gffS1InFile;
   if(!open(gffS1InFile, toCString(options.inputS1FileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputS1FileName;
      cerr << " for reading.\n";
      return 1;
   }

   GffFileIn gffS2InFile;
   if(!open(gffS2InFile, toCString(options.inputS2FileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputS2FileName;
      cerr << " for reading.\n";
      return 1;
   }

   GffFileIn gffS3InFile;
   if(!open(gffS3InFile, toCString(options.inputS3FileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputS3FileName;
      cerr << " for reading.\n";
      return 1;
   }
      
   // for each SLM

      // do

         // get <window and >window soma average methylation
         // get <window and >window sperm average methylation

         // if abs(sperm-soma) > x, extend SLM

      // if SLM extended, loop




   return 0;
}
