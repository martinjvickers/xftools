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
   CharString chrom;
   double ratio;
   int window;
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

   addOption(parser, ArgParseOption("s", "window-size", "Size of window",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "window-size", "100");
   addOption(parser, ArgParseOption("r", "ratio", "Ratio",
                                    ArgParseArgument::DOUBLE, "DOUBLE"));
   setDefaultValue(parser, "ratio", "0.2");

   addOption(parser, ArgParseOption("c", "chr", "Which chromosome to work on",
                                    ArgParseArgument::STRING, "STR"));
   getOptionValue(options.chrom, parser, "chr");
   setDefaultValue(parser, "chr", "1");

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

   getOptionValue(options.window, parser, "window-size");
   getOptionValue(options.ratio, parser, "ratio");

   return ArgumentParser::PARSE_OK;
}

typedef IntervalAndCargo<int, GffRecord> TInterval;

std::set<CharString> getList(GffFileIn &gffSLMInFile, 
                             map<CharString, String<TInterval>> &intervals)
{
   std::set<CharString> chr;
   GffRecord record;

   while(!atEnd(gffSLMInFile))
   {
      readRecord(record, gffSLMInFile);
      chr.insert(record.ref);
      appendValue(intervals[record.ref], 
                  TInterval(record.beginPos, record.endPos, record));
   }

   cout << "Unique chrs " << chr.size() << endl;

   return chr;
}

int average(GffFileIn &gffSpermInFile, GffFileIn &s1, GffFileIn &s2, 
            GffFileIn &s3, CharString chr, 
            map<int, unsigned short int[4]> &average)
{
   bool s1Started = false, s2Started = false;
   bool s3Started = false, spmStarted = false;
   GffRecord record;

   cout << "[Starting S1]" << endl;

   while(!atEnd(s1))
   {
      readRecord(record, s1);

      if(record.ref == chr)
         s1Started = true;
      else if(record.ref != chr && s1Started == true)
         break;

      if(record.ref == chr)
      {
         int c, t;
         for(int i = 0; i < length(record.tagNames); i++)
         {
            if(record.tagNames[i] == "c")
               c = stoi(toCString(record.tagValues[i]));
            else if(record.tagNames[i] == "t")
               t = stoi(toCString(record.tagValues[i]));
         }
         average[record.beginPos][0] += c;
         average[record.beginPos][1] += t;
      } 
   }

   cout << "[Completed S1]" << endl;
   cout << "[Starting S2]" << endl;

   while(!atEnd(s2))
   {
      readRecord(record, s2);

      if(record.ref == chr)
         s2Started = true;
      else if(record.ref != chr && s2Started == true)
         break;

      if(record.ref == chr)
      {
         int c, t;
         for(int i = 0; i < length(record.tagNames); i++)
         {
            if(record.tagNames[i] == "c")
               c = stoi(toCString(record.tagValues[i]));
            else if(record.tagNames[i] == "t")
               t = stoi(toCString(record.tagValues[i]));
         }

         average[record.beginPos][0] += c;
         average[record.beginPos][1] += t;
      }
   }

   cout << "[Completed S2]" << endl;
   cout << "[Starting S3]" << endl;

   while(!atEnd(s3))
   {
      readRecord(record, s3);

      if(record.ref == chr)
         s3Started = true;
      else if(record.ref != chr && s3Started == true)
         break;

      if(record.ref == chr)
      {
         int c, t;
         for(int i = 0; i < length(record.tagNames); i++)
         {
            if(record.tagNames[i] == "c")
               c = stoi(toCString(record.tagValues[i]));
            else if(record.tagNames[i] == "t")
               t = stoi(toCString(record.tagValues[i]));
         }

         average[record.beginPos][0] += c;
         average[record.beginPos][1] += t;
      }
   }

   cout << "[Completed S3]" << endl;
   cout << "[Started spm]" << endl;

   while(!atEnd(gffSpermInFile))
   {
      readRecord(record, gffSpermInFile);

      if(record.ref == chr)
         spmStarted = true;
      else if(record.ref != chr && spmStarted == true)
         break;

      if(record.ref == chr)
      {
         int c, t;
         for(int i = 0; i < length(record.tagNames); i++)
         {
            if(record.tagNames[i] == "c")
               c = stoi(toCString(record.tagValues[i]));
            else if(record.tagNames[i] == "t")
               t = stoi(toCString(record.tagValues[i]));
         }

         average[record.beginPos][2] += c;
         average[record.beginPos][3] += t;
      }
   }

   cout << "[Completed spm]" << endl;

   return 0;
}

int perform_extend(GffFileIn &gffSLMInFile,
                   map<int, unsigned short int[4]> data,
                   double ratio, int window, CharString chrom)
{
   GffRecord record;
   while(!atEnd(gffSLMInFile))
   {
      readRecord(record, gffSLMInFile);
      if(record.ref == chrom)
      {
         cout << "[START]\t" << record.beginPos << "\t" << record.endPos << endl;
         
         bool front = false, back = false; // record if stopping
         int front_pos = record.beginPos;
         int back_pos = record.endPos;

         while(front == false)
         {

            int distance = window;
            int front_c_soma = 0, front_t_soma = 0;
            int front_c_spm = 0, front_t_spm = 0;
            for(int i = front_pos-distance; i < front_pos-distance+window; i++)
            {
               front_c_soma += data[i][0];
               front_t_soma += data[i][1];
               front_c_spm += data[i][2];
               front_t_spm += data[i][3];
            }
            double soma = (double)front_c_soma/((double)front_c_soma+(double)front_t_soma);
            double spm = (double)front_c_spm/((double)front_c_spm+(double)front_t_spm);
            if(abs(soma - spm) >= ratio)
            {
               front_pos -= distance;
               cout << "Extending" << endl;
            }
            else
            {
               front = true;
            }
         }

         while(back == false)
         {
            int back_c_soma = 0, back_t_soma = 0;
            int back_c_spm = 0, back_t_spm = 0;
            int distance = window;

            for(int i = back_pos; i < back_pos+distance; i++)
            {
               back_c_soma += data[i][0];
               back_t_soma += data[i][1];
               back_c_spm += data[i][2];
               back_t_spm += data[i][3];
            }
            double soma = (double)back_c_soma/((double)back_c_soma+(double)back_t_soma);
            double spm = (double)back_c_spm/((double)back_c_spm+(double)back_t_spm);
            if(abs(soma - spm) >= ratio)
            {
               back_pos += distance;
               cout << "Extending" << endl;
            }
            else
            {
               back = true;
            }
         }
         cout << "[END]\t"<< record.ref << "\t" << front_pos << "\t" << back_pos << endl;
      }
   }

   return 0;
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

   CharString meh = options.chrom;

   map<int, unsigned short int[4]> data;
   average(gffSpermInFile, gffS1InFile, gffS2InFile, gffS3InFile, 
           "1", data);

   perform_extend(gffSLMInFile, data, options.ratio, options.window, "1");

   return 0;
}
