#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>
#include <functional>
#include "boost/icl/interval.hpp"
#include "boost/icl/interval_map.hpp"
#include <boost/icl/split_interval_map.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include <set>

#include <chrono>
using namespace seqan;
using namespace std;

using namespace boost::icl;
using boost::bad_lexical_cast;

#include "feature.cpp"

/*
Map of map of split interval map. AKA 3D map/map/split-interval-map
The rationale here is that the complete map is able to separate out based
on strand and chromosome.
*/
typedef split_interval_map<int, Feature> featuremap;
typedef std::map<char, featuremap> strandedmap;
typedef std::map<string, strandedmap> completemap;

struct ModifyStringOptions
{
        CharString inputFileName;
        CharString inputAnnotationFileName;
        CharString outputFileName;
        bool exclude;
        bool lazyRef;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("Overlap");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	addOption(parser, seqan::ArgParseOption("a", "input-annotation-file", "Path to the input filter file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setRequired(parser, "input-annotation-file");
	addOption(parser, seqan::ArgParseOption("o", "output-file", "Path to the output file", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output-file");
	addOption(parser, seqan::ArgParseOption("l", "lazy-ref", "Internally it will capitalise both the input and annoation reference names so that chr1, Chr1 and CHR1 will all match. The output GFF will be of the same format as the annoation file."));

	setShortDescription(parser, "Overlap");
	setVersion(parser, "0.0.1");
	setDate(parser, "March 2017");
	addUsageLine(parser, "-i input.gff -a annotation.gff -o output.gff [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Given an annotation/feature file and an input w1 file, give counts for each feature");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");
	getOptionValue(options.inputAnnotationFileName, parser, "input-annotation-file");
	getOptionValue(options.outputFileName, parser, "output-file");
	options.lazyRef = isSet(parser, "lazy-ref");
	

	return seqan::ArgumentParser::PARSE_OK;
}

/*
Loads the 'base' GFF file and puts it into our completemap interval map
*/
completemap loadBase(ModifyStringOptions options)
{
	completemap results;

        //get our value and add value to map
        GffFileIn gffAnnotationIn;
        GffRecord annotationrecord;
        if (!open(gffAnnotationIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
        }

        //put all of our annotations into a vector and store the pairs in an interval map
        while (!atEnd(gffAnnotationIn))
        {
                readRecord(annotationrecord, gffAnnotationIn);
                Feature feat(0, annotationrecord.strand, annotationrecord.ref, annotationrecord.beginPos, annotationrecord.endPos);
                discrete_interval<int> inter_val = discrete_interval<int> (annotationrecord.beginPos, annotationrecord.endPos);
                string currRef;
                if(options.lazyRef == true){
                        //convert to upper
                        string tmp = toCString(annotationrecord.ref);
                        boost::to_upper(tmp);
                        currRef = boost::to_upper_copy(tmp);
                } else {
                        currRef = toCString(annotationrecord.ref);
                }
                results[currRef][annotationrecord.strand] += make_pair(inter_val, feat);
        }

        close(gffAnnotationIn);
	return results;
}

/*
With the base loaded into the completemap, it's time to look for an overlap
*/
int findOverlaps(ModifyStringOptions options, completemap results)
{
        //now load our new data
        GffFileIn gffRawIn;
        GffRecord to_bin_record;
        if (!open(gffRawIn, toCString(options.inputFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }

	bool differror = false;
	bool stranderror = false;
	int count = 0;

        //logic, this is where we bin or count etc with our overlaps.
        while (!atEnd(gffRawIn))
        {
                try
                {
			readRecord(to_bin_record, gffRawIn); //read record
			discrete_interval<int> key = discrete_interval<int> (to_bin_record.beginPos, to_bin_record.endPos); //create interval
			string currRef = toCString(to_bin_record.ref); //get the chromosome

			//convert chromosome to uppercase if treatment of reference/chromosome text is treated as lazy 
			if(options.lazyRef == true)
                                for (string::size_type i = 0; i < currRef.length(); ++i)
                                        boost::to_upper(currRef);

                        //checks to see if the chromosome/ref we are working with actually exists
                        if(results.find(currRef) == results.end())
                        {
                                if(differror == false){
                                        cerr << "Error: Your input file contains a reference [ " << currRef << " ] that does not exist in the annotation" << endl;
                                        cerr << "Error: References that exist in your annotation are :" << endl;
                                        cerr << "Error: ";
                                        for(auto& a: results)
                                                cerr << a.first << " ";
                                        cerr << endl;
                                        cerr << "Error: If this is simply that the annotation and input ref columns differ in the use of capital letters then you can rerun with the -l flag." << endl;
                                        cerr << "Error: The -l flag will capitalise both the input and the annotation internally." << endl;
                                        cerr << "Error: This will continue to run, but will not count these inputs." << endl;
                                        differror = true;
                                }
                        }

			//check if stranded
                        std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itresnew;

			//if we have more than one strand
                        if(results[currRef][to_bin_record.strand].size() > 0)
                        {
				//then we worry about strand data. 
				//Find the intervals in the results map that match out key range.
                                itresnew = results[currRef][to_bin_record.strand].equal_range(key);

				//now go through each of the intervals that match and then do something
                                for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
				{
					//cout << to_bin_record.beginPos << " " << to_bin_record.endPos << " matches " << (*itnew).second.startPos() << " " << (*itnew).second.endPos() << endl;
					cout << "Hello " << to_bin_record.beginPos << " " << to_bin_record.endPos << " matches " << (*itnew).second.startPos() << " " << (*itnew).second.endPos() << " " << (*itnew).second.strand() << " " << (*itnew).second.ref() << endl;
                                }
                                
                        } else {
                                if(stranderror == false){
                                        cout << "Error: Input data has no strand information but the reference does. The input data will be added to features in your annotation from both strands" << endl;
                                        stranderror = true;
                                }

                                for(auto& i : results[currRef])
                                {
                                        itresnew = i.second.equal_range(key);
                                        for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
                                        {
						cout << "Hello " <<  to_bin_record.beginPos << " " << to_bin_record.endPos << " matches " << (*itnew).second.startPos() << " " << (*itnew).second.endPos() << " " << (*itnew).second.strand() << " " << (*itnew).second.ref() << endl;
                                        }
                                }
                        }

			count++;

		} 
		catch (Exception const & e) 
		{
                        std::cerr << "Warning: " << e.what() << " and died on line " << count << std::endl;
                }

	}

	close(gffRawIn);

	return 0;
}

/*
Overlap. Given two gff files, this program will find where file B overlaps file A. 
*/
int main(int argc, char const ** argv)
{
	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	completemap results; // new way results[ref][strand].thingy

	results = loadBase(options); // put our 'Base' into the complete map

	findOverlaps(options, results); // using out complete map, let's search for overlaps

/*
	//now, let's reopen out feature file and get our feature counts
	GffFileIn gffFeatureIn;
	GffRecord featurerecord;
	if (!open(gffFeatureIn, toCString(options.inputAnnotationFileName)))
        {
                std::cerr << "ERROR: Could not open example.gff" << std::endl;
                return 1;
        }

	while (!atEnd(gffFeatureIn))
        {
		readRecord(featurerecord, gffFeatureIn);

                string currRef;
                if(options.lazyRef == true)
                {
			string tmp = toCString(featurerecord.ref);
                        boost::to_upper(tmp);
                        currRef = boost::to_upper_copy(tmp);
                } else {
                        currRef = toCString(featurerecord.ref);
                }

		discrete_interval<int> key = discrete_interval<int> (featurerecord.beginPos, featurerecord.endPos);

		//here is where we can take some user input to decide what to do with the 6th column. e.g. if we want to calculate methyl, average. count

		int score = 0;
		std::map<string, string> tagmap;
		std::pair<split_interval_map<int, Feature>::iterator, split_interval_map<int, Feature>::iterator> itresnew = results[currRef][featurerecord.strand].equal_range(key);
		for(auto itnew = itresnew.first; itnew != itresnew.second; ++itnew)
                {
			Feature mmm = (*itnew).second;
			score = score + mmm.score();
			tagmap = mmm.tagmap();
		}
		cout << featurerecord.ref << "\t" << featurerecord.source << "\t" << featurerecord.type << "\t" << featurerecord.beginPos << "\t" << featurerecord.endPos << "\t" << score << "\t" << featurerecord.strand << "\t" << featurerecord.phase << "\t";
	
		//now loop the tags at the end
		for(auto& m: tagmap)
		{
			cout << m.first << "=" << m.second <<";";
		}
		cout << endl;
	}

	close(gffFeatureIn);
*/
	return 0;
}
