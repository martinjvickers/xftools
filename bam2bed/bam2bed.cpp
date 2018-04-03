#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <string>
using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString outputFileName;
   bool cigar_aware = false;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv)
{
   seqan::ArgumentParser parser("bam2bed");
   setShortDescription(parser, "XFTOOLS");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-i input.bam -o output.gff [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Creates a BED file from a BAM/SAM file.");

   // accept an input file
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("c", "cigar-aware", "This will use the full CIGAR span length. NOTE: presently this does not account for a gap due to a reference insertion. So this may not be what you want. I may change this in the future"));

   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

   // parse the inputs into the options struct
   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   options.cigar_aware = isSet(parser, "cigar-aware");

   return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
	
   //read in input bam file
   BamFileIn inFile;
   if (!open(inFile, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open " << options.inputFileName;
      cerr << " for reading.\n";
      return 1;
   }

   //get the bam context
   typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
   TBamContext const & bamContext = context(inFile);

   //create output file
   BedFileOut bedOut;
   // out(std::cout, Gff());
   if(!open(bedOut, toCString(options.outputFileName)))
   {
      std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
      return 1;
   }

   //<position,count>
   map<int,int> counter;

   BamHeader header;
   readHeader(header, inFile);
   BamAlignmentRecord record;
   int rID = -1;

   do
   {
      readRecord(record, inFile);

      if(options.cigar_aware == true)
      {
         int spanlength;
         _getLengthInRef(spanlength, record.cigar);
         int querylength = _getQueryLength(record.cigar);
         cout << contigNames(bamContext)[record.rID] << "\t";
         cout << record.beginPos << "\t" << record.beginPos+spanlength << endl;
      }
      else 
      {
      }
   } while(!atEnd(inFile));

   //writeToFile(counter, inFile, rID, gffOutFile);
   close(inFile);
   close(bedOut);

   return 0;
}
