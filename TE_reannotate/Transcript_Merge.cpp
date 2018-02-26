#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
   bool exclude;
   bool lazyRef;
   int methyl_cutoff = 0.25;
   int merge_distance = 300;
   int cull_size = 200;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("TE_reannotate");
   addOption(parser, ArgParseOption("a", "input-annotation-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-annotation-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("l", "lazy-ref", 
       "Internally it will capitalise both the input and annoation reference \
        names so that chr1, Chr1 and CHR1 will all match. The output GFF will \
        be of the same format as the annoation file."));
   setShortDescription(parser, "TE_reannotate");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Given an annotation/feature file and an input w1 \
                           file, give counts for each feature");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   options.lazyRef = isSet(parser, "lazy-ref");

   return ArgumentParser::PARSE_OK;
}

map<CharString, vector<pair<unsigned int, unsigned int>>> blocks;
BamFileOut outBam;
GffFileOut outGFF;

int initial_merge(ModifyStringOptions options, GffFileIn &gffIn, vector<GffRecord> &merging)
{
   GffRecord record;
   GffRecord last;

   readRecord(last, gffIn);

   while (!atEnd(gffIn)) // loop through GFF file
   {
      readRecord(record, gffIn);
      unsigned int pos = 0;
      for(auto i : record.tagNames)
      {
         if(i == "gene_id")
         {
            break;
         }

         pos++;
      }
      
      if(record.tagValues[pos] == last.tagValues[pos] && record.strand == last.strand && record.ref == last.ref)
      {
         if(record.endPos > last.endPos)
            last.endPos = record.endPos;
         else
            record.endPos = record.endPos;
      }
      else
      {
         merging.push_back(last);
         last = record;
      }
   }

   merging.push_back(last);

   return 0;
}

int within_range(ModifyStringOptions options, vector<GffRecord> &merging)
{
   GffRecord last = merging[0];   
   vector<GffRecord> newMerging;

   for(unsigned int i = 1; i < merging.size(); i++)
   {

      GffRecord record = merging[i];
      if(record.strand == last.strand && record.ref == last.ref && (record.beginPos-last.endPos) <= 50)
      {
         last.endPos = record.endPos;
         if(record.endPos > last.endPos)
            last.endPos = record.endPos;
         else
            record.endPos = record.endPos;
      }
      else
      {
         newMerging.push_back(last);
         last = record;
      }

   }

   newMerging.push_back(last);

   merging = newMerging;

   return 0;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   // parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if(res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // Read GFF
   GffFileIn gffIn;
   if(!open(gffIn, toCString(options.inputAnnotationFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputAnnotationFileName << endl;
      return 1;
   }

   vector<GffRecord> merging;

   // initial merge which would combine them if they are the same gene_id
   initial_merge(options, gffIn, merging);

   for(auto i : merging)
      cout << i.ref << "\txftools\ttranscript\t" << i.beginPos << "\t" << i.endPos << "\t.\t" << i.strand << "\t.\t."<< endl;


   within_range(options, merging);

//   for(auto i : merging)
  //    cout << i.ref << "\txftools\ttranscript\t" << i.beginPos << "\t" << i.endPos << "\t.\t" << i.strand << "\t.\t."<< endl;

   close(outBam);
   return 0;

}
