#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>

using namespace boost::iostreams;

using namespace seqan;
using namespace std;

#define DEFAULT_BUF_LENGTH (16 * 16384)

struct ModifyStringOptions {
   CharString outputPrefix;
   CharString inputAnnotationFileName;
   CharString outputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("tsv_extract");
   addOption(parser, ArgParseOption("a", "input-annotation-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "input-annotation-file");
   addOption(parser, ArgParseOption("p", "output-prefix", 
             "Specify the prefix for your output files.", 
             ArgParseArgument::STRING, "TEXT"));
   setRequired(parser, "output-prefix");
   setShortDescription(parser, "Overlap");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2017");
   addUsageLine(parser, "-a annotation.gff -p outputPrefix \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Extract methylation data from TSV based on input GFF annotation.\
                           Ensure you pipe your tsv into tsv_extract. eg \
                           zcat input.tsv.gz | tsv_extract -a annotation.gff -p sampleName");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputAnnotationFileName, 
                  parser, "input-annotation-file");
   getOptionValue(options.outputPrefix, parser, "output-prefix");

   return ArgumentParser::PARSE_OK;
}

struct Methyl
{
   char strand;
   string context;
   int methylated_bases;
   int total_bases;
   bool methylation_call;
};

string calculate_context(string mc_class)
{
   if(mc_class[1] == 'G')
      return "CG";
   else if(mc_class[2] == 'G')
      return "CHG";
   else
      return "CHH";
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // chromosomes, position, details
   map<CharString, map<unsigned int, Methyl>> data;
   
   // read input file into memory
   std::string line;
   while (std::getline(std::cin, line))
   {
      istringstream iss( line );
      string word;
      int pos = 0;
      Methyl row;
      string chromosome;
      unsigned int base;
      /*
      0	chrom
      1	pos
      2 strand
      3 context
      4 methylated_bases
      5 total_bases
      6 methylation_call
      */
      while(getline(iss, word, '\t'))
      {
         stringstream convert(word);
         if(pos == 0)
         {
            chromosome = word;
         }
         else if(pos == 1)
         {
            convert >> base;
         }
         else if(pos == 2)
         {
            char strand;
            convert >> strand;
            row.strand = strand;
         }
         else if(pos == 3)
         {
            row.context = word;
         }
         else if(pos == 4)
         {
            int methylated_bases;
            convert >> methylated_bases;
            row.methylated_bases = methylated_bases;
         }
         else if(pos == 5)
         {
            int total_bases; 
            convert >> total_bases;
            row.total_bases = total_bases;
         }
         else if(pos == 6)
         {
            bool methylation_call;
            convert >> methylation_call;
            row.methylation_call = methylation_call;
         }

         pos++;
      }

      if(pos != 7)
      {
         cerr << "Error: line does not contain the correct number of columns";
         cerr << endl;
      }

      data[chromosome][base] = row;
   }

   // create three output GFF file handlers, for each context
   GffFileOut gffCGFileOut, gffCHGFileOut, gffCHHFileOut;
   string cg_out = toCString(options.outputPrefix) + string(".CG.gff");
   string chg_out = toCString(options.outputPrefix) + string(".CHG.gff");
   string chh_out = toCString(options.outputPrefix) + string(".CHH.gff");
   if (!open(gffCGFileOut, toCString(cg_out)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }
   if (!open(gffCHGFileOut, toCString(chg_out)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }
   if (!open(gffCHHFileOut, toCString(chh_out)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }

   string raw_out = toCString(options.outputPrefix) + string(".raw.txt");
   ofstream ofFile; 
   ofFile.open(raw_out);

   GffFileIn gffFileIn;
   if (!open(gffFileIn, toCString(options.inputAnnotationFileName)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }
   GffRecord record;

   while (!atEnd(gffFileIn))
   {
      try
      {
         readRecord(record, gffFileIn);
      }
      catch(Exception const & e)
      {
         cerr << "ERROR: " << e.what() << endl;
      }

      // zero counter for this GFF record
      int CG_tot = 0;
      int CHG_tot = 0;
      int CHH_tot = 0;
      int CG_meth = 0;
      int CHG_meth = 0;
      int CHH_meth = 0;

      for(int i = record.beginPos; i < record.endPos; i++)
      {
         if(data.find(record.ref) != data.end() && data[record.ref].find(i) != data[record.ref].end())
         {
            //work out the context and count total and methyl bases
            string context = calculate_context(data[record.ref][i].context);

            if(context == "CG")
            {
               CG_tot += data[record.ref][i].total_bases;
               CG_meth += data[record.ref][i].methylated_bases;
            }
            else if(context == "CHG")
            {
               CHG_tot += data[record.ref][i].total_bases;
               CHG_meth += data[record.ref][i].methylated_bases;
            }
            else if(context == "CHH")
            {
               CHH_tot += data[record.ref][i].total_bases;
               CHH_meth += data[record.ref][i].methylated_bases;
            }

            ofFile << record.ref << "\t" << i << "\t" << data[record.ref][i].strand;
            ofFile << "\t" << data[record.ref][i].context;
            ofFile << "\t" << data[record.ref][i].methylated_bases;
            ofFile << "\t" << data[record.ref][i].total_bases;
            ofFile << "\t" << data[record.ref][i].methylation_call;
            ofFile << endl;
            
         }
      }

      // calculate the record information
      GffRecord cgRecord = record;
      cgRecord.score = (float)CG_meth / (float)CG_tot;
      appendValue(cgRecord.tagNames, "c");
      appendValue(cgRecord.tagValues, to_string(CG_meth));
      appendValue(cgRecord.tagNames, "t");
      appendValue(cgRecord.tagValues, to_string(CG_tot-CG_meth));
      appendValue(cgRecord.tagNames, "n");
      appendValue(cgRecord.tagValues, to_string(CG_tot));
      cgRecord.type = toCString(options.outputPrefix) + string("_CG");
      writeRecord(gffCGFileOut, cgRecord);

      GffRecord chgRecord = record;
      chgRecord.score = (float)CHG_meth / (float)CHG_tot;
      appendValue(chgRecord.tagNames, "c");
      appendValue(chgRecord.tagValues, to_string(CHG_meth));
      appendValue(chgRecord.tagNames, "t");
      appendValue(chgRecord.tagValues, to_string(CHG_tot-CHG_meth));
      appendValue(chgRecord.tagNames, "n");
      appendValue(chgRecord.tagValues, to_string(CHG_tot));
      chgRecord.type = toCString(options.outputPrefix) + string("_CHG");
      writeRecord(gffCHGFileOut, chgRecord);

      GffRecord chhRecord = record;
      chhRecord.score = (float)CHH_meth / (float)CHH_tot;
      appendValue(chhRecord.tagNames, "c");
      appendValue(chhRecord.tagValues, to_string(CHH_meth));
      appendValue(chhRecord.tagNames, "t");
      appendValue(chhRecord.tagValues, to_string(CHH_tot-CHH_meth));
      appendValue(chhRecord.tagNames, "n");
      appendValue(chhRecord.tagValues, to_string(CHH_tot));
      chhRecord.type = toCString(options.outputPrefix) + string("_CHH");
      writeRecord(gffCHHFileOut, chhRecord);
   }

   return 0;
}
