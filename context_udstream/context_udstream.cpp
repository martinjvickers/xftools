#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <string>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString outputFileName;
   CharString refFileName;
   int upstream, downstream;
   bool suppress;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv)
{
   ArgumentParser parser("context_udstream");
   setShortDescription(parser, "XFTOOLS");
   setVersion(parser, "0.0.1");
   setDate(parser, "November 2019");
   addUsageLine(parser, "-r input.fa -i input_w1.gff -o output.gff");

   addDescription(parser, "Calculate up and downstream sequence of a context.");

   // accept an input file
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the w1 input file.",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");

   addOption(parser, ArgParseOption("r", "reference-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "reference-file");

   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");

   addOption(parser, ArgParseOption("u", "upstream", "Upstream of Context", 
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "upstream", "3");

   addOption(parser, ArgParseOption("d", "downstream", "Downstream of Context",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "downstream", "3");

   addOption(parser, ArgParseOption("s", "suppress-errors", 
             "Suppress error messages."));

   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

   // parse the inputs into the options struct
   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.refFileName, parser, "reference-file");
   getOptionValue(options.outputFileName, parser, "output-file");

   getOptionValue(options.upstream, parser, "upstream");
   getOptionValue(options.downstream, parser, "downstream");

   options.suppress = isSet(parser, "suppress-errors");

   return ArgumentParser::PARSE_OK;
}

void getInfo(ModifyStringOptions &options, GffRecord &r, 
             FaiIndex &faiIndex)
{
   Dna5String reference_seq, upstream_seq, downstream_seq, span_seq;
   
   unsigned idx = 0;
   if(!getIdByName(idx, faiIndex, r.ref))
   {
      if(options.suppress == true)
      {
         cerr << "ERROR: FAI index has no entry for ";
         cerr << r.ref << endl;
      }
   }
   else
   {
      int len = length(r.type);
      if(r.strand == '+')
      {
         readRegion(reference_seq, faiIndex, idx, r.beginPos,    
                    r.beginPos+len);
         readRegion(span_seq, faiIndex, idx, r.beginPos-options.upstream,
                    r.beginPos+len+options.downstream);

         upstream_seq = prefix(span_seq, options.upstream);
         downstream_seq   = suffix(span_seq, length(span_seq)-options.downstream);
      }
      else if(r.strand == '-')
      {
         readRegion(reference_seq, faiIndex, idx, r.beginPos-len+1,
                    r.beginPos+1);

         readRegion(span_seq, faiIndex, idx, r.beginPos-options.upstream-len+1,
                    r.beginPos+options.downstream+1);

         reverseComplement(reference_seq);
         reverseComplement(span_seq);

         downstream_seq = prefix(span_seq, options.downstream);
         upstream_seq   = suffix(span_seq, length(span_seq)-options.upstream);
         reverseComplement(downstream_seq);
         reverseComplement(upstream_seq);
         reverseComplement(span_seq);
      }
      else
      {
         cerr << "[ERROR] No strand information " << endl;
         // maybe I could add a function to work it out in the end
         // but not today.
      }

      cout << r.ref << "\t" << r.beginPos << "\t" << r.type << "\t";
      cout << r.score << "\t";
      cout << r.tagValues[0] << "\t" << r.tagValues[1] << "\t";
      cout << reference_seq << "\t";
      cout << r.strand << "\t" << upstream_seq << "\t";
      cout << downstream_seq << "\t" << span_seq << endl;
   }
}

int main(int argc, char const ** argv)
{
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
        
   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if(res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   //index the reference
   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.refFileName)))
   {
      cerr << "ERROR: Could not build FAI index for file ";
      cerr << options.refFileName << ".\n";
      return 1;
   }

   // create output file
   GffFileIn gffInFile;
   if(!open(gffInFile, toCString(options.inputFileName)))
   {
      cerr << "ERROR: Could not open ";
      cerr << toCString(options.inputFileName) << endl;
      return 1;
   }

   while(!atEnd(gffInFile))
   {
      GffRecord region;
      readRecord(region, gffInFile);

      getInfo(ref(options), ref(region), ref(faiIndex));
   }

   close(gffInFile);

   return 0;
}


/*
unsigned int count_Cs (Dna5String &seq)
{
   unsigned int counter = 0;

   for(unsigned int i = 0; i < length(seq); i++)
   if(seq[i] == 'C')
      counter++;
   return counter;
}

void calculate_contexts(Dna5String &seq, map<string,int> &context_map)
{
   for(unsigned i = 2; i < length(seq)-2; i++)
   {
      Infix<String<Dna5> >::Type inf = infix(seq, i, i+3);
      CharString trin = inf;
      if(inf[0] == 'C')
      {
         if(inf[1] == 'G')
	 {
            context_map["CG"]++; //CG
            context_map[toCString(trin)]++;
            //cout << "CG\t" << inf << "\t" << trin<< endl;
         }
         else
         {
            if(inf[2] == 'G')
            {
               context_map["CHG"]++; //CHG
               context_map[toCString(trin)]++;
               //cout << "CHG\t" << inf << "\t" << trin<< endl;
            }
            else
            {
               context_map["CHH"]++; //CHH
               context_map[toCString(trin)]++;
               //cout << "CHH\t" << inf << "\t" << trin<< endl;
            }
         }
      }
   }
}

void print_map(map<string,int> &context_map)
{
   for (auto i : context_map)
      cout << i.first << "=" << i.second << "\t";
   cout << endl;
}

void new_context_map(map<string,int> &context_map)
{
   context_map.insert( make_pair("CG", 0));
   context_map.insert( make_pair("CGA", 0));
   context_map.insert( make_pair("CGC", 0));
   context_map.insert( make_pair("CGG", 0));
   context_map.insert( make_pair("CGT", 0));
   context_map.insert( make_pair("CGN", 0));
   context_map.insert( make_pair("CHG", 0));
   context_map.insert( make_pair("CAG", 0));
   context_map.insert( make_pair("CCG", 0));
   context_map.insert( make_pair("CTG", 0));
   context_map.insert( make_pair("CNG", 0));
   context_map.insert( make_pair("CHH", 0));
   context_map.insert( make_pair("CAA", 0));
   context_map.insert( make_pair("CAC", 0));
   context_map.insert( make_pair("CAT", 0));
   context_map.insert( make_pair("CCA", 0));
   context_map.insert( make_pair("CCC", 0));
   context_map.insert( make_pair("CCT", 0));
   context_map.insert( make_pair("CTA", 0));
   context_map.insert( make_pair("CTC", 0));
   context_map.insert( make_pair("CTT", 0));
   context_map.insert( make_pair("CNN", 0));
}

void write_record(GffFileOut &gffOutFile, FaiIndex &faiIndex, 
                  unsigned int &contig, unsigned int &beginPos, 
                  int endPos, unsigned int &C_count, 
                  map<string,int> &context_map, 
                  ModifyStringOptions &options)
{
   GffRecord record;
   record.ref = sequenceName(faiIndex, contig);
   record.source = "xftools";
   record.type = options.label;
   record.beginPos = beginPos + 1; 
   record.endPos = endPos;
   if(options.percentage == true)
      record.score = (float)C_count / (float)(endPos - beginPos + 1);
   else
      record.score = C_count;
   record.strand = '.';
   record.phase = '.';
   StringSet<CharString> tagNames;
   StringSet<CharString> tagValues;
   //cout << beginPos << "\t" << endPos << endl;
   for (auto i : context_map)
   {
      appendValue(record.tagNames, toCString(i.first));
      //CharString m = i.second;
      appendValue(record.tagValues, to_string(i.second));
      //appendValue(record.tagValues, toCString(i.second));
   }
   writeRecord(gffOutFile, record);
}
*/


/*
	Method.
          Create index of reference
          Loop through the reference in windows
            count c's + contexts
            reverse compliment and repeat
            print to GFF
*/
/*
int main(int argc, char const ** argv)
{
   ModifyStringOptions options;
   seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if(res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;

   // create output file
   GffFileOut gffOutFile;
   if(!open(gffOutFile, toCString(options.outputFileName)))
   {
      std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
      return 1;
   }

   // create output file
   GffFileIn gffInFile;
   if(!open(gffInFile, toCString(options.inputRegionFileName)))
   {
      std::cerr << "ERROR: Could not open output.gff" " for reading.\n";
      return 1;
   }

   //index the reference
   FaiIndex faiIndex;
   if(!build(faiIndex, toCString(options.inputFileName)))
   {
      std::cerr << "ERROR: Could not build FAI index for file ";
      std::cerr << options.inputFileName << ".\n";
      return 1;
   }

   // I need to know the number of contigs
   Dna5String reference_seq;

   if(options.inputRegionFileName == NULL)
   {
      // for each contig
      for(unsigned int contig = 0; contig < numSeqs(faiIndex); contig++)
      {
         unsigned int beginPos = 0;
         while (beginPos+options.window_size < sequenceLength(faiIndex, contig))
         {
            map<string,int> context_map;
            new_context_map(context_map);
	
            readRegion(reference_seq, faiIndex, contig, beginPos, 
                       beginPos + options.window_size);

            // count the c's
            unsigned int total_C = count_Cs(reference_seq);
            reverseComplement(reference_seq);
            total_C = total_C + count_Cs(reference_seq);

            // contexts
            if(beginPos > 0 && 
               (beginPos+options.window_size+2 < sequenceLength(faiIndex, contig)))
            {
               // get what we want with +2 on either side
               readRegion(reference_seq, faiIndex, contig, beginPos-2, 
                          beginPos + options.window_size+2);
               calculate_contexts(reference_seq, context_map);
               reverseComplement(reference_seq);
               calculate_contexts(reference_seq, context_map);
            }
            else 
            {
               if(beginPos < 2)	// we're at the start!
                  readRegion(reference_seq, faiIndex, contig, beginPos, 
                             beginPos + options.window_size+2);

               Dna5String pad = "NN";
               pad += reference_seq;
               reference_seq = pad;
               calculate_contexts(reference_seq, context_map);
               reverseComplement(reference_seq);
               calculate_contexts(reference_seq, context_map);
            }

            write_record(gffOutFile, faiIndex, contig, beginPos, 
                         beginPos + options.window_size, total_C, 
                         context_map, options);

            beginPos = beginPos + options.window_size;
         }

         // Do the last one too. Will work that out later but it's just this
         map<string,int> context_map;
         new_context_map(context_map);

         readRegion(reference_seq, faiIndex, contig, beginPos, 
                    sequenceLength(faiIndex, contig));

         unsigned int total_C = count_Cs(reference_seq);
         reverseComplement(reference_seq);
         total_C = total_C + count_Cs(reference_seq);

         readRegion(reference_seq, faiIndex, contig, beginPos-2, 
                    sequenceLength(faiIndex, contig));

         Dna5String pad = "NN";
         reference_seq += pad;
         calculate_contexts(reference_seq, context_map);
         reverseComplement(reference_seq);
         calculate_contexts(reference_seq, context_map);

         write_record(gffOutFile, faiIndex, contig, beginPos, 
                      (int)sequenceLength(faiIndex, contig), 
                      total_C, context_map, options);
      }
   } 
   else 
   {
      while(!atEnd(gffInFile))
      {
         GffRecord region;
         readRecord(region, gffInFile);

         map<string,int> context_map;
         new_context_map(context_map);
         unsigned idx = 0;
         if(!getIdByName(idx, faiIndex, region.ref))
            std::cout << "ERROR: FAI index has no entry for " << endl;

         readRegion(reference_seq, faiIndex, idx, region.beginPos, region.endPos);

         unsigned int total_C = count_Cs(reference_seq);
         reverseComplement(reference_seq);
         total_C = total_C + count_Cs(reference_seq);
					
         if( region.beginPos > 0 && 
            (region.endPos+2 < sequenceLength(faiIndex, idx)) )
         {
            readRegion(reference_seq, faiIndex, idx, region.beginPos-2, 
                       region.endPos+2);
            calculate_contexts(reference_seq, context_map);
            reverseComplement(reference_seq);
            calculate_contexts(reference_seq, context_map);
         }
         else
         {
            if(region.beginPos < 2)
               readRegion(reference_seq, faiIndex, idx, region.beginPos, 
                          region.endPos+2);

            Dna5String pad = "NN";
            pad += reference_seq;
            reference_seq = pad;
            calculate_contexts(reference_seq, context_map);
            reverseComplement(reference_seq);
            calculate_contexts(reference_seq, context_map);
         }

         write_record(gffOutFile, faiIndex, idx, region.beginPos, region.endPos, 
                      total_C, context_map, options);
      }
   }

   close(gffOutFile);
   return 0;
}


*/
