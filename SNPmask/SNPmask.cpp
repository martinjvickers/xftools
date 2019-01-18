#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
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

   fstream file;
   file.open(toCString(options.inputAnnotationFileName));

   string line;
   while(getline( file, line,'\n'))
   {
      istringstream templine(line);
      string data;
      int i = 0;
      string ref;
      int pos;
      string snp;
      while (getline( templine, data,'\t'))
      {
         if(i == 2)
            ref = data;
         else if(i == 3)
            pos = stoi(data);
         else if(i == 5)
            snp = data;
         i++;
      }
      cout << ref << "\t" << snp << "\t" << pos << endl;
   }
   file.close();

   return 0;
}
