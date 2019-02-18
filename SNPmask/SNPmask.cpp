#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <sstream>

using namespace seqan;
using namespace std;

struct ModifyStringOptions {
   CharString inputFileName;
   CharString inputAnnotationFileName;
   CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("SNPmask");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("s", "input-snp-file", 
                                    "Path to the SNP input file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-snp-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");
   setShortDescription(parser, "SNPmask");
   setVersion(parser, "0.0.1");
   setDate(parser, "January 2019");
   addUsageLine(parser, "-i input.fa -s annotation.txt -o masked.fa \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.inputAnnotationFileName, parser, "input-snp-file");
   getOptionValue(options.outputFileName, parser, "output-file");

   return ArgumentParser::PARSE_OK;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // load fasta file
   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not build FAI index for file ";
      cerr << options.inputFileName << ".\n";
      return 1;
   }

   // Masked file to write out
   SeqFileOut seqFileOut(toCString(options.outputFileName));

   unsigned idx = 0;
   CharString seq;
   bool firstTime = true;

   fstream file;
   file.open(toCString(options.inputAnnotationFileName));

   string line;
   while(getline( file, line,'\n'))
   {
      string ref, snp, data;
      istringstream templine(line);
      int i = 0;
      int pos;
      while (getline( templine, data,'\t'))
      {
         if(i == 1)
            ref = data;
         else if(i == 2)
            pos = stoi(data);
         else if(i == 4)
            snp = data;
         i++;
      }

      unsigned currentIdx;
      if(!getIdByName(currentIdx, faiIndex, ref))
         std::cout << "ERROR: FAI index has no entry for " << ref << "\n";

      if(idx != currentIdx && firstTime == false)
      {
         writeRecord(seqFileOut, sequenceName(faiIndex, idx), seq);
      }

      if(idx != currentIdx || firstTime == true)
      {
         idx = currentIdx;
         readSequence(seq, faiIndex, idx);
         firstTime = false;
      }

      //if(snp[0] != seq[pos-1])
      //   cout << ref << "\t" << sequenceName(faiIndex, currentIdx) << "\t" << currentIdx << "\t" << snp << "\t" << pos << "\t" << seq[pos-1] << endl;
      seq[pos-1]  = 'N';
   }

   // write out the last one
   writeRecord(seqFileOut, sequenceName(faiIndex, idx), seq);

   close(seqFileOut);
   file.close();

   return 0;
}
