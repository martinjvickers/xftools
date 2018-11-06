#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/seq_io.h>
#include <boost/algorithm/string.hpp>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
   CharString inputGenome;
   bool lazyRef;
   int extension;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("DMR_generator");
   addOption(parser, ArgParseOption("f", "feature-file", 
                                    "Path to the feature file to mimick", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "feature-file");
   addOption(parser, ArgParseOption("a", "input-gene-annotation", 
                                    "Path to the gene annotation file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   //setRequired(parser, "input-annotation-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("g", "genome-file", "Genome file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("l", "lazy-ref", 
       "Internally it will capitalise both the input and annoation reference \
        names so that chr1, Chr1 and CHR1 will all match. The output GFF will \
        be of the same format as the annoation file."));
   addOption(parser, seqan::ArgParseOption("e", "extension", "Extension around the feature file",
                                           seqan::ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "extension", "0");
   setShortDescription(parser, "DMR_generator");
   setVersion(parser, "0.0.1");
   setDate(parser, "November 2018");
   addUsageLine(parser, "-f slms.gff -a genes.gff -g genome.fa -o output.gff \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, " \
                           ");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "feature-file");
   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-gene-annotation");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.inputGenome, parser, "genome-file");
   options.lazyRef = isSet(parser, "lazy-ref");
   getOptionValue(options.extension, parser, "extension");

   return ArgumentParser::PARSE_OK;
}

typedef IntervalAndCargo<int, GffRecord> TInterval;

/*
   Inserts the records into the interval tree

 */ 
int insertIntervals(map<CharString, String<TInterval>> &intervals, 
                    GffFileIn &gffInFile, unsigned int extension, 
                    ModifyStringOptions options)
{
   GffRecord record;

   try
   {
      while(!atEnd(gffInFile))
      {
         readRecord(record, gffInFile);
         CharString ref = record.ref;
         if(options.lazyRef == true)
         {
            toUpper(ref);
         }
         appendValue(intervals[ref], TInterval(record.beginPos-extension, record.endPos+extension, record));
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

unsigned int devrand(void)
{
   int fn;
   unsigned int r;

   fn = open("/dev/urandom", O_RDONLY);
   if(fn == -1)
      exit(-1);

   if(read(fn, &r, 4) != 4)
      exit(-1);

   close(fn);
   return r;
}

unsigned int KISS()
{
   unsigned int x = devrand();
   unsigned int y;
   while(!(y = devrand()));
   unsigned int z = devrand();

   unsigned int c = devrand() % 698769069 + 1;

   unsigned long long t, a = 698769069ULL;
   x = 69069 * x + 12345;
   y ^= (y<<13); y ^= (y>>17); y ^= (y<<5);
   t = a*z+c; c = (t>>32);

   return x+y+(z=t);
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // I guess the best thing is to load the Gene and TE list into an interval tree for searching
   String<TInterval> intervals;
   map<CharString, String<TInterval>> full;

   // open annotation file
   GffFileIn gffAnnotInFile;
   if(!open(gffAnnotInFile, toCString(options.inputAnnotationFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputAnnotationFileName;
      cerr << " for reading.\n";
      return 1;
   }

   if(insertIntervals(full, gffAnnotInFile, options.extension, options) != 0)
   {
      cerr << "ERROR: Cannot put intervals into interval vector " << endl;
      return 1;
   }

   map<CharString, IntervalTree<int, GffRecord>> trees;
   for(auto i : full)
      trees[i.first] = i.second;

   // then iterate through the thing I need to mimick to find length, distance from gene/TE
      // find location that mimicks what we have requested in a) options b) loci to mimick
   GffFileIn gffFeatureInFile;
   if(!open(gffFeatureInFile, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << options.inputFileName;
      cerr << " for reading.\n";
      return 1;
   }

   // load genome
   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.inputGenome)))
   {
      cerr << "ERROR: Could not build reference index" << endl;
      return 1;
   }

   GffRecord record;
   while(!atEnd(gffFeatureInFile))
   {
      readRecord(record, gffFeatureInFile);
      String<GffRecord> results;

      // this gives us what we are trying to mimic
      findIntervals(results, trees[record.ref], record.beginPos, record.endPos);

      // so for each record, what i need to do is determine if it overlaps something
      /*
       I need to select a random chromosome/contig
       and then a position. So I probably need the genome?
      */
      String<GffRecord> currentOverlaps;
      unsigned int randChr, randPos;

      do
      {
         randChr = KISS() % numSeqs(faiIndex);
         randPos = KISS() % sequenceLength(faiIndex, randChr);
         CharString ref = sequenceName(faiIndex, randChr);
         if(options.lazyRef == true)
            toUpper(ref);
         findIntervals(currentOverlaps, trees[ref], randPos, randPos + (record.endPos - record.beginPos));
      }
      while(length(currentOverlaps) != length(results));

      cout << sequenceName(faiIndex, randChr) << "\t" << randPos << "\t";
      cout << randPos + (record.endPos - record.beginPos) << "\t";
      cout << length(currentOverlaps) << endl;

   }
   return 0;
}
