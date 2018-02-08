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
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("TE_reannotate");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
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

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   options.lazyRef = isSet(parser, "lazy-ref");

   return ArgumentParser::PARSE_OK;
}

map<CharString, vector<pair<unsigned int, unsigned int>>> blocks;
BamFileOut outBam;
GffFileOut outGFF;


int getMethylationBlocks(ModifyStringOptions options, GffFileIn &gffIn)
{
   GffRecord record;
   GffRecord last;
   CharString ref = NULL;
   int beginBlock = -1;
   int endBlock = -1;
   bool inBlock = false;
   while (!atEnd(gffIn)) // loop through GFF file
   {
      readRecord(record, gffIn);
      if(inBlock == false && record.score > 0.0) // start of block
      {
         inBlock = true;
         ref = record.ref;
         beginBlock = record.beginPos;
         endBlock = record.endPos;
      }
      else if((inBlock == true && record.score == 0) || (record.ref != ref && inBlock == true)) // finished block
      {
         pair<unsigned int, unsigned int> block = make_pair(beginBlock, endBlock);
         blocks[ref].push_back(block);
         inBlock = false;
         beginBlock = -1;
         endBlock = -1;
      }
      else if(inBlock == true && record.score > 0.0) // within block
      {
         endBlock = record.endPos;
      }
   }

   if(inBlock == true) // finished block
   {
      pair<unsigned int, unsigned int> block = make_pair(beginBlock, endBlock);
      blocks[ref].push_back(block);
      inBlock = false;
      beginBlock = -1;
      endBlock = -1;
   }

   
   return 0;
}

bool overlap(BamAlignmentRecord record1, BamAlignmentRecord record2)
{
   if(max(record1.beginPos+length(record1.seq), record2.beginPos+length(record2.seq))
      - min(record1.beginPos, record2.beginPos)
      <
      (record1.beginPos+length(record1.seq) - record1.beginPos) +
      (record2.beginPos+length(record2.seq) - record2.beginPos)
     )
   {
      return true;
   }
   else
   {
      return false;
   }
}

bool overlap(BamAlignmentRecord record1, pair<unsigned int, unsigned int> record2)
{
   unsigned int record1_first = record1.beginPos;
   unsigned int record1_second = record1.beginPos+length(record1.seq);

   if(max(record1_second, record2.second)
      - min(record1_first, record2.first)
      <
      (record1_second - record1_first) +
      (record2.second - record2.first)
     )
   {
      return true;
   }
   else
   {
      return false;
   }
}


/*
   The way I'm going to tackle this is to try to make an unbroken
   overlap of reads from X bp of the 5' end until I make contact
   with the feature, wiping attempts until we have an unbroken one.

   I should do this twice, one for each strand.
*/
int search(BamFileIn &inFile, pair<unsigned int, unsigned int> block, CharString ref)
{
   BamAlignmentRecord record;
   readRecord(record, inFile);

   vector<BamAlignmentRecord> current_run;
   current_run.push_back(record); // this record shouldn't overlap with
                                  // the block

   bool blockOverlap = false;

   // do while the record
   do
   {
      if(atEnd(inFile))
         return 1;
      else
         readRecord(record, inFile);
      

      if(overlap(record, block) == true)
         blockOverlap = true;

      // if they overlap, then add record to the run
      if(overlap(current_run[current_run.size()-1], record))
      {
         current_run.push_back(record);
      }
      else
      {
         if(blockOverlap == true)
         {
            for(auto i : current_run)
            {
   /*            GffRecord record;
               record.ref = ref;
               record.source = "TE_reannoation";
               record.type = "CG_block";
               record.beginPos = block.first;
               record.endPos = block.second;
               record.strand = '.';
               record.score = GffRecord::INVALID_SCORE();
               writeRecord(outGFF, record);*/
               writeRecord(outBam, i);
            }
         }
         current_run.clear();

         if(record.beginPos < block.second)
            current_run.push_back(record);
      }

   }
   while(current_run.size() != 0 && record.beginPos < block.second);

}

int findReads(ModifyStringOptions options, BamIndex<Bai> &baiIndex, BamFileIn &inFile)
{
 //  BamHeader header;
   //readHeader(header, inFile);

   for(auto chromosome : blocks)
   {
      int rID = 0;
      if (!getIdByName(rID, contigNamesCache(context(inFile)), chromosome.first))
      {
         std::cerr << "ERROR: Reference sequence named " << chromosome.first << " not known.\n";
         return 1;
      }

      for(auto block : chromosome.second)
      {
         bool hasAlignments = false;
         if(!jumpToRegion(inFile, hasAlignments, rID, block.first, block.second, baiIndex))
         {
            std::cerr << "ERROR: Could not jump to " << block.first << ":" << block.second << "\n";
            return 1;
         }

         search(inFile, block, chromosome.first);

         /*
         cout << "Jumped to region " << chromosome.first << " " << block.first;
         cout << " " << block.second << endl;
         cout << "First read " << record.qName << " " << record.beginPos;
         cout << " " << record.beginPos+length(record.seq) << endl;
         */

      }
   }

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

   if(!open(outGFF, toCString("blocks.gff")))
   {
      cerr << "ERROR: Could not open ";
      cerr << "blocks.gff" << endl;
      return 1;
   }

   getMethylationBlocks(options, gffIn);

   for(auto i : blocks)
   {
      for(auto block : i.second)
      {
      GffRecord record;
      record.ref = i.first;
      record.source = "TE_reannoation";
      record.type = "CG_block";
      record.beginPos = block.first;
      record.endPos = block.second;
      record.strand = '.';
      record.score = GffRecord::INVALID_SCORE();
      writeRecord(outGFF, record);
      }
   }

   BamFileIn inFile;
   if (!open(inFile, toCString(options.inputFileName)))
   {
      std::cerr << "ERROR: Could not open " << options.inputFileName << " for reading.\n";
      return 1;
   }

   if (!open(outBam, toCString("reads.bam")))
   {
      return 1;
   }

   BamHeader header;
   readHeader(header, inFile);
   outBam.context = context(inFile);
   writeHeader(outBam, header);

   BamIndex<Bai> baiIndex;
   string baiString = toCString(options.inputFileName) + string(".bai");
   if (!open(baiIndex, toCString(baiString) ))
   {
      cerr << "ERROR: Could not read BAI index file " << "\n";
      cerr << "Ensure you've sorted and indexed your BAM file using samtools.";
      cerr << endl;
      return 1;
   }

   // find reads that overlap block
//   findReads(ModifyStringOptions options);
   findReads(options, baiIndex, inFile); 

   close(outBam);
   return 0;

}
